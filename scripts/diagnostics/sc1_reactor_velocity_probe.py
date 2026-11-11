#!/usr/bin/env python3
"""Reuse the verified v5 binary for three fresh steps, writing each for diagnosis."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import sc1_reactor_pressure_coupling as coupled
import sc1_reactor_velocity_fields as fields
from instrument_sc1_alpha_subcycles import _once

base=coupled.base
startup=coupled.startup
TAG='sc1_meml_velocity_location_probe_v1'
FAILED_LOG='4c469e3438576ad12139cad251de18b159dd35161cdcb7aede16124b193ecaac'
TIMES=('5e-05','0.0001','0.00015')
END=0.00015


def recipe():
    paths=(Path(__file__),Path(fields.__file__),base.HERE/'sc1_reactor_startup_case.py',
        base.HERE/'sc1_reactor_startup_flow.py',base.HERE/'sc1_wall_mesh.py')
    return {p.name:base.sha(p) for p in paths}


def control(root):
    root=root.resolve(strict=True);m,_=coupled.verify(root)
    if base.sha(root/'quiescent/log.solver')!=FAILED_LOG:
        raise ValueError('Requires the received unchanged v5 velocity failure')
    attempt=json.loads((root/'quiescent/attempt.json').read_text())
    if attempt!={k:m[k] for k in ('build','assets','ranks')}:raise ValueError('Control attempt/build mismatch')
    return m,dict(root=str(root),log_sha256=FAILED_LOG,
        manifest_sha256=base.sha(root/'startup_manifest.json'))


def output_control(text):
    text=_once(text,'stopAt endTime; endTime 0.01;','stopAt endTime; endTime 0.00015;')
    return _once(text,'writeControl timeStep; writeInterval 5;',
        'writeControl timeStep; writeInterval 1;')


def prepare(failed_root,out):
    failed_root=failed_root.resolve(strict=True);m,identity=control(failed_root)
    startup.require_allocation(m['ranks'])
    protected=[failed_root,base.REPO]+[Path(m[k]) for k in ('build','archive','regression_root')]
    protected += [Path(m[k]['root']) for k in ('baseline','failed_control','face_control','coupled_control','pressure_probe')]
    out=base.fresh(out,protected);out.mkdir(parents=True,exist_ok=False)
    source=failed_root/'quiescent';case=out/'quiescent'
    for folder in ('0','constant','system'):shutil.copytree(source/folder,case/folder)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):shutil.copy2(source/name,case/name)
    text=output_control((source/'system/controlDict').read_text())
    (case/'system/controlDict').write_text(text)
    expected=dict(m['inputs']);expected['system/controlDict']=hashlib.sha256(text.encode()).hexdigest()
    if base.case_input_hashes(case)!=expected or control(failed_root)!=(m,identity):
        raise ValueError('Control/inputs changed during preparation')
    manifest=dict(m,root=str(out),variant=TAG,inputs=expected,velocity_control=identity,
        velocity_recipe=recipe(),source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',manifest)
    print('VELOCITY_PROBE='+str(out)+' (fresh t=0; same compiled v5 binary; three steps with every-step output; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads((root/'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['velocity_recipe']!=recipe():raise ValueError('Wrong/changed spatial probe')
    parent,identity=control(Path(m['velocity_control']['root']))
    if identity!=m['velocity_control']:raise ValueError('Spatial probe control changed')
    expected=dict(parent['inputs'])
    text=output_control((Path(parent['root'])/'quiescent/system/controlDict').read_text())
    expected['system/controlDict']=hashlib.sha256(text.encode()).hexdigest()
    if m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Spatial probe inputs changed beyond output/duration')
    for key in parent:
        if key not in ('root','variant','inputs','source_commit') and m[key]!=parent[key]:
            raise ValueError('Spatial probe metadata changed: '+key)
    return m,parent


def records(text):
    """Select the measured three-step prefix; include no failure/shutdown records."""
    active=None;result=[]
    for line in text.splitlines():
        match=re.fullmatch(r'Time = (\S+)',line)
        if match:active=float(match[1])
        if active is None or active>END+1e-12:continue
        marker=line.split(' ',1)[0]
        if marker in ('SC1_TWO_STEP','SC1_QUIESCENT_STEP','SC1_COUPLED_FLUX') or (
                marker=='SC1_PRESSURE_FIELD' and 'term=pUnclipped ' in line):
            result.append(dict(marker=marker,**dict(re.findall(r'(\w+)=([^\s]+)',line))))
    return result


def compare_prefix(actual,reference):
    a=records(actual);b=records(reference)
    if len(a)!=33 or len(a)!=len(b):raise ValueError('Incomplete three-step audit/pressure/flux prefix')
    for left,right in zip(a,b):
        if left.keys()!=right.keys():raise ValueError('Prefix record shape changed')
        for key,value in left.items():
            try:x=float(value);y=float(right[key])
            except ValueError:
                if value!=right[key]:raise ValueError('Prefix record identity changed: '+key)
            else:
                if not math.isfinite(x+y) or not math.isclose(x,y,rel_tol=1e-9,abs_tol=1e-12):
                    raise ValueError('Spatial probe changed measured trajectory: '+key)
    return a


def analyze(case,m,parent):
    text=(case/'log.solver').read_text()
    times=re.findall(r'^Time = (\S+)$',text,re.M)
    if len(times)!=3 or any(abs(float(t)-float(expected))>1e-12 for t,expected in zip(times,TIMES)) or text.splitlines().count('End')!=1:
        raise ValueError('Spatial probe did not complete exactly three steps')
    from sc1_reactor_watch import Monitor
    monitor=Monitor()
    for line in text.splitlines():
        if monitor.feed(line):raise ValueError('Numerical stop in spatial probe: '+monitor.reason['reason'])
    selected=compare_prefix(text,(Path(parent['root'])/'quiescent/log.solver').read_text())
    measurements=fields.snapshots(case,TIMES,startup.CELLS,m['initial']['interface_midheight_m'])
    for snapshot in measurements:
        t=snapshot['time_s'];speed=snapshot['phases']['liquid']['present']['max_m_s']
        row=next(r for r in selected if r['marker']=='SC1_QUIESCENT_STEP' and abs(float(r['t_s'])-t)<1e-12)
        if speed>1e-3 or abs(speed-float(row['maxLiquidPresent_m_s']))>1e-10:
            raise ValueError('Reconstructed liquid speed disagrees with passed step guard')
        row=next(r for r in selected if r['marker']=='SC1_TWO_STEP' and r['phase']=='liquid' and abs(float(r['t_s'])-t)<1e-12)
        if abs(snapshot['liquidVolume_L']-float(row['volume_m3'])*1000)>1e-8:
            raise ValueError('Reconstructed liquid inventory disagrees with audit')
    hashes={time+'/'+name:base.sha(case/time/name) for time in TIMES for name in fields.FIELDS}
    return dict(suite=TAG,diagnostic_complete=True,source_control=m['velocity_control'],
        prefixAgreement=dict(records=len(selected),absoluteTolerance=1e-12,relativeTolerance=1e-9),
        snapshots=measurements,fieldHashes=hashes,log_sha256=base.sha(case/'log.solver'),
        meaning='Three-step spatial diagnostic only. The parent 10-ms startup failed at 0.2 ms. No startup, restart, reactor or biology acceptance.',
        patchMeaning='Adjacent means owner of a boundary face; it is not a distance-to-wall calculation.',
        pressureMeaning='Static pressure reconstructed as p_rgh + (alphaGas + 1000*alphaLiquid)*(-9.81*y), matching the unchanged initial geometry/gravity.')


def run(root):
    root=root.resolve(strict=True);m,parent=verify(root);startup.require_allocation(m['ranks'])
    case=root/'quiescent';allowed={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}
    if {p.name for p in case.iterdir()}!=allowed:raise ValueError('Refusing attempted spatial probe')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=base.environment(Path(json.loads((Path(m['build'])/'coupling_build_manifest.json').read_text())['base_build']))
    app=m['candidate_identity']['binary']
    try:
        if m['ranks']>1:
            base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
            command=['srun','--kill-on-bad-exit=1','--time=00:10:00','-n',str(m['ranks']),'--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated']
        else:command=[app,'-fileHandler','uncollated']
        from sc1_reactor_watch import execute
        execute(command,case,case/'log.solver',env,wall_seconds=600)
        # A failed solver is never reconstructed. Three normal written timesteps
        # are still diagnostic outputs; native field lengths/values are checked.
        if m['ranks']>1:
            for rank in range(m['ranks']):
                for time in TIMES:
                    for name in fields.FIELDS:
                        p=case/f'processor{rank}'/time/name
                        if p.is_symlink() or not p.is_file():raise ValueError('Missing/redirected processor output '+str(p))
            base.execute(['reconstructPar','-time','5e-05:0.00015','-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
        report=analyze(case,m,parent)
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Spatial probe initial inputs changed')
        base.write_json(case/'result.json',report)
        compact=dict(suite=TAG,diagnostic_complete=True,prefixAgreement=report['prefixAgreement'],
            snapshots=[dict(time_s=r['time_s'],liquidVolume_L=r['liquidVolume_L'],
                liquidPresent=r['phases']['liquid']['present'],liquidBulk=r['phases']['liquid']['bulk'],
                liquidWallAdjacent=r['phases']['liquid']['wallAdjacent'],volumeMixtureRms_m_s=r['volumeMixtureRms_m_s']) for r in report['snapshots']],
            meaning=report['meaning'])
        print(json.dumps(compact,indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VELOCITY_PROBE=COMPLETE (diagnostic only; parent startup remains failed)',flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='action',required=True)
    s=subs.add_parser('prepare');s.add_argument('--failed-root',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('run');s.add_argument('--root',type=Path,required=True);s.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='prepare':prepare(a.failed_root,a.out)
    else:run(a.root)
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 velocity probe stopped: '+str(error)+'; preserve all roots',file=sys.stderr);raise SystemExit(1)
