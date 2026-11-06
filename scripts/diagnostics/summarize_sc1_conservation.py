#!/usr/bin/env python3
"""Read SC1 timestep audits; optionally report independent native-time endpoints."""
import argparse
import math
from pathlib import Path

PHASES=('liquid','bubbleGas','headspaceGas','gas','total')


def require_numbers(row, keys, context):
    for key in keys:
        if key not in row or not isinstance(row[key], (int, float)):
            raise ValueError(f'{context}: missing or nonnumeric {key}')


def records(text, prefix):
    rows=[]
    for line in text.splitlines():
        if not line.startswith(prefix+' '):
            continue
        row={}
        for item in line.split()[1:]:
            if '=' not in item:
                continue
            k,v=item.split('=',1)
            if k in row:
                raise ValueError(f'Duplicate {k} in {prefix}')
            try:
                n=float(v)
            except ValueError:
                row[k]=v
            else:
                if not math.isfinite(n):
                    raise ValueError(f'Nonfinite {k} in {prefix}')
                row[k]=n
        rows.append(row)
    return rows


def load_log(path):
    text=path.read_text()
    if not any(line.strip()=='End' for line in text.splitlines()):
        raise ValueError(f'{path}: no solver End; log incomplete')
    if 'FOAM FATAL' in text:
        raise ValueError(f'{path}: fatal solver error')
    method=records(text,'SC1_BALANCE_METHOD')
    if len(method)!=1 or method[0].get('version')!=1:
        raise ValueError(f'{path}: need one v1 timestep-audit session; legacy endpoint logs are insufficient')
    starts={}
    for row in records(text,'SC1_BALANCE_START'):
        phase=row.get('phase')
        if phase not in PHASES:
            raise ValueError(f'Unknown initial phase {phase}')
        require_numbers(row, ('t_s','mass_kg','oxygen_kg'), 'Initial inventory')
        if 'volume_m3' in row:
            require_numbers(row, ('volume_m3',), 'Initial inventory')
        if phase in starts:
            raise ValueError('Repeated restart/session in log')
        starts[phase]=row
    if set(starts)!=set(PHASES):
        raise ValueError('Missing initial phase inventories')
    if any(starts[p]['t_s']!=starts['total']['t_s'] for p in PHASES):
        raise ValueError('Initial phase start times differ within log')
    rows={p:[] for p in PHASES}
    for row in records(text,'SC1_BALANCE'):
        p=row.get('phase')
        if p not in rows:
            raise ValueError(f'Unknown phase {p}')
        required=('t_s','start_s','steps','mass_kg','oxygen_kg','massResidual_kg','oxygenResidual_kg',
                  'massExternalIntegral_kg','massTransferIntegral_kg','massInIntegral_kg','massOutIntegral_kg',
                  'oxygenTransferIntegral_kg','normalization_kg','equationResidualIntegral_kg')
        require_numbers(row, required, 'Balance record')
        if not row['steps'].is_integer():
            raise ValueError('Non-integer step count')
        previous=rows[p][-1] if rows[p] else {'t_s':starts[p]['t_s'],'steps':0}
        if row['t_s']<=previous['t_s'] or row['steps']<=previous['steps']:
            raise ValueError('Non-increasing time/steps')
        if row['start_s']!=starts[p]['t_s']:
            raise ValueError('Start time changed within log')
        # Independent mass closure recomputation checks parser, signs and output.
        expected=row['massExternalIntegral_kg']+row['massTransferIntegral_kg']+row['massInIntegral_kg']-row['massOutIntegral_kg']
        residual=row['mass_kg']-starts[p]['mass_kg']-expected
        if not math.isclose(residual,row['massResidual_kg'],rel_tol=1e-7,abs_tol=1e-9):
            raise ValueError(f'Inconsistent mass arithmetic: {p}')
        rows[p].append(row)
    if any(not rows[p] for p in PHASES):
        raise ValueError('No completed audit intervals')
    times=[r['t_s'] for r in rows['total']]
    if any([r['t_s'] for r in rows[p]]!=times for p in PHASES):
        raise ValueError('Missing phase records at report times')
    steps=[r['steps'] for r in rows['total']]
    if any([r['steps'] for r in rows[p]]!=steps for p in PHASES):
        raise ValueError('Phase step counts differ at report times')
    return starts,rows


def summarize(path, rows, native_times=False):
    print(f'\n{path}: cumulative every-timestep quadrature; NOT a fully discrete conservation certificate')
    end=rows['total'][-1]
    duration=end['t_s']-end['start_s']
    print(f'Own interval: start_s={end["start_s"]:.17g} endpoint_s={end["t_s"]:.17g}'
          f' duration_s={duration:.17g} steps={int(end["steps"])} reports={len(rows["total"])}')
    print('phase             mass residual (mg)   O2 residual (mg)   O2 transfer sum (mg)   normalization (mg)')
    for p in PHASES:
        r=rows[p][-1]
        print(f'{p:18} {1e6*r["massResidual_kg"]:18.7g} {1e6*r["oxygenResidual_kg"]:18.7g}'
              f' {1e6*r["oxygenTransferIntegral_kg"]:22.7g} {1e6*r["normalization_kg"]:20.7g}')
    print('Inspect total O2 transfer sum for implicit/explicit coupling lag; internal exchange should cancel at convergence.')
    if native_times:
        print('Native endpoint inventories; mean rates use only this run\'s own duration above.')
        print('phase                mass (kg)       O2 (kg)   mean mass source (mg/s)'
              '   mean mass transfer (mg/s)   mean O2 transfer (mg/s)   mean normalization (mg/s)')
        for p in PHASES:
            r=rows[p][-1]
            print(f'{p:18} {r["mass_kg"]:13.9g} {r["oxygen_kg"]:13.9g}'
                  f' {1e6*r["massExternalIntegral_kg"]/duration:25.9g}'
                  f' {1e6*r["massTransferIntegral_kg"]/duration:27.9g}'
                  f' {1e6*r["oxygenTransferIntegral_kg"]/duration:25.9g}'
                  f' {1e6*r["normalization_kg"]/duration:27.9g}')


def validate_initial(control, normal):
    cs,_=control; ns,_=normal
    for p in PHASES:
        if cs[p]['t_s']!=ns[p]['t_s']:
            raise ValueError(f'Initial t_s differs for {p}; not matched conditions')
        keys=['mass_kg','oxygen_kg']
        if 'volume_m3' in cs[p] or 'volume_m3' in ns[p]:
            if 'volume_m3' not in cs[p] or 'volume_m3' not in ns[p]:
                raise ValueError(f'Initial volume_m3 missing for {p}; cannot verify matched conditions')
            keys.append('volume_m3')
        for k in keys:
            if not math.isclose(cs[p][k],ns[p][k],rel_tol=1e-12,abs_tol=1e-14):
                raise ValueError(f'Initial {k} differs for {p}; not matched conditions')


def compare(control,normal,tolerance=1e-6):
    validate_initial(control,normal)
    _,cr=control; _,nr=normal
    if len(cr['total'])!=len(nr['total']):
        raise ValueError('Different number of saved comparisons')
    pairs=[]
    for c,n in zip(cr['bubbleGas'],nr['bubbleGas']):
        if abs(c['t_s']-n['t_s'])>tolerance:
            raise ValueError('Output times differ; do not silently interpolate')
        pairs.append((c['t_s'],1e6*(n['mass_kg']-c['mass_kg']),
                      1e6*(n['massResidual_kg']-c['massResidual_kg'])))
    return pairs


def native_endpoint_difference(control, normal):
    """Descriptive difference only: never pair or interpolate saved rows."""
    validate_initial(control, normal)
    c=control[1]['bubbleGas'][-1]; n=normal[1]['bubbleGas'][-1]
    return {
        'control_endpoint_s':c['t_s'], 'normal_endpoint_s':n['t_s'],
        'control_duration_s':c['t_s']-c['start_s'],
        'normal_duration_s':n['t_s']-n['start_s'],
        'endpoint_offset_s':n['t_s']-c['t_s'],
        'bubble_inventory_difference_mg':1e6*(n['mass_kg']-c['mass_kg']),
        'bubble_residual_difference_mg':1e6*(n['massResidual_kg']-c['massResidual_kg']),
    }


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('control',type=Path)
    p.add_argument('--normal',type=Path)
    p.add_argument('--native-times',action='store_true',help=
                   'Report independent native-time intervals and endpoint differences; '
                   'no interpolation or synchronized comparison')
    args=p.parse_args(argv)
    loaded={}; errors=[]
    # Report every valid completed run even if its peer or the comparison fails.
    for name,path in (('control',args.control),('normal',args.normal)):
        if path is None:
            continue
        try:
            loaded[name]=load_log(path)
            summarize(path,loaded[name][1],native_times=args.native_times)
        except (OSError,ValueError,KeyError) as e:
            errors.append(f'{name} ({path}): {e}')
    if errors:
        for error in errors:
            print(f'Audit interpretation stopped: {error}')
        return 1
    try:
        if args.normal:
            control=loaded['control']; normal=loaded['normal']
            if args.native_times:
                difference=native_endpoint_difference(control,normal)
                print('\nNative endpoint difference (normal minus control); NOT a synchronized comparison.')
                print('Independent timestamps and durations; no saved rows paired, no interpolation.')
                for key,value in difference.items():
                    print(f'{key}={value:.17g}')
            else:
                paired=compare(control,normal)
                print('\nTime(s)  normal-control bubble inventory(mg)  residual difference(mg)')
                for t,dm,dr in paired:
                    print(f'{t:.12g} {dm:.8g} {dr:.8g}')
        return 0
    except (OSError,ValueError,KeyError) as e:
        print(f'Audit interpretation stopped: {e}')
        if args.normal and not args.native_times:
            print('Use --native-times only for independently timed endpoint reporting; '
                  'it does not repair unmatched initial conditions.')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
