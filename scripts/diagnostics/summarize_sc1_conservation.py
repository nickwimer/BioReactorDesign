#!/usr/bin/env python3
"""Read new SC1 audit logs; never infer timestep integrals from legacy writes."""
import argparse
import math
from pathlib import Path

PHASES=('liquid','bubbleGas','headspaceGas','gas','total')


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
        if row['phase'] in starts:
            raise ValueError('Repeated restart/session in log')
        starts[row['phase']]=row
    if set(starts)!=set(PHASES):
        raise ValueError('Missing initial phase inventories')
    rows={p:[] for p in PHASES}
    for row in records(text,'SC1_BALANCE'):
        p=row['phase']
        if p not in rows:
            raise ValueError(f'Unknown phase {p}')
        required=('t_s','start_s','steps','mass_kg','oxygen_kg','massResidual_kg','oxygenResidual_kg',
                  'massExternalIntegral_kg','massTransferIntegral_kg','massInIntegral_kg','massOutIntegral_kg',
                  'oxygenTransferIntegral_kg','normalization_kg','equationResidualIntegral_kg')
        if any(k not in row for k in required):
            raise ValueError('Incomplete balance record')
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
    return starts,rows


def summarize(path, rows):
    print(f'\n{path}: cumulative every-timestep quadrature; NOT a fully discrete conservation certificate')
    print('phase             mass residual (mg)   O2 residual (mg)   O2 transfer sum (mg)   normalization (mg)')
    for p in PHASES:
        r=rows[p][-1]
        print(f'{p:18} {1e6*r["massResidual_kg"]:18.7g} {1e6*r["oxygenResidual_kg"]:18.7g}'
              f' {1e6*r["oxygenTransferIntegral_kg"]:22.7g} {1e6*r["normalization_kg"]:20.7g}')
    print('Inspect total O2 transfer sum for implicit/explicit coupling lag; internal exchange should cancel at convergence.')


def compare(control,normal,tolerance=1e-6):
    cs,cr=control; ns,nr=normal
    for p in PHASES:
        for k in ('t_s','mass_kg','oxygen_kg'):
            if not math.isclose(cs[p][k],ns[p][k],rel_tol=1e-12,abs_tol=1e-14):
                raise ValueError(f'Initial {k} differs for {p}; not matched conditions')
    if len(cr['total'])!=len(nr['total']):
        raise ValueError('Different number of saved comparisons')
    pairs=[]
    for c,n in zip(cr['bubbleGas'],nr['bubbleGas']):
        if abs(c['t_s']-n['t_s'])>tolerance:
            raise ValueError('Output times differ; do not silently interpolate')
        pairs.append((c['t_s'],1e6*(n['mass_kg']-c['mass_kg']),
                      1e6*(n['massResidual_kg']-c['massResidual_kg'])))
    return pairs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('control',type=Path)
    p.add_argument('--normal',type=Path)
    args=p.parse_args()
    try:
        control=load_log(args.control)
        normal=load_log(args.normal) if args.normal else None
        paired=compare(control,normal) if normal else None
        summarize(args.control,control[1])
        if normal:
            summarize(args.normal,normal[1])
            print('\nTime(s)  normal-control bubble inventory(mg)  residual difference(mg)')
            for t,dm,dr in paired:
                print(f'{t:.12g} {dm:.8g} {dr:.8g}')
        return 0
    except (OSError,ValueError,KeyError) as e:
        print(f'Audit interpretation stopped: {e}')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
