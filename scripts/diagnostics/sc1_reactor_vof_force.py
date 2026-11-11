#!/usr/bin/env python3
"""Build and run a passive two-step VOF force decomposition."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sc1_two_field as base
import sc1_reactor_startup as startup
import sc1_reactor_vof_nonorth as nonorth
import sc1_reactor_vof_force_build as builder
from instrument_sc1_alpha_subcycles import _once
import sc1_reactor_vof_saved as saved
import sc1_reactor_vof_force_report as reporting
from sc1_reactor_failure import regular

TAG='sc1_meml_vof_force_probe_v1'
CONTROL_LOG='d2630e82e732eb6d0ff333ca69a2ef47e90e267897e478a6da23efc5eb9dd8d1'


def output_control(text):
    text=_once(text,'application sc1InterFoamReference;','application '+builder.APP+';')
    return _once(text,'stopAt endTime; endTime 0.01;','stopAt endTime; endTime 0.0001;')


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(reporting.__file__))}


def control(root):
    root=root.resolve(strict=True);m=nonorth.verify(root)
    log=regular(root,'quiescent/log.solver')
    if base.sha(log)!=CONTROL_LOG:raise ValueError('Requires the received nonorthogonal control')
    if json.loads(regular(root,'quiescent/attempt.json').read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('Nonorthogonal control attempt changed')
    view,iterations=nonorth.reporting.pressure_coverage(log.read_text(),require_complete=False)
    if iterations['complete'] or iterations['observedStepRecords']!=168 or iterations['corrections']!=672 or iterations['convergedCorrections']!=672:
        raise ValueError('Wrong nonorthogonal failure/convergence evidence')
    saved.budget(view,m['initial'])
    return m,dict(root=str(root),manifest_sha256=base.sha(root/'startup_manifest.json'),log_sha256=CONTROL_LOG)


def expected_inputs(parent,source):
    text=output_control(regular(source,'system/controlDict').read_text())
    expected=dict(parent['inputs']);expected['system/controlDict']=hashlib.sha256(text.encode()).hexdigest()
    return expected,text


def prepare(control_root,build_root,out):
    control_root=control_root.resolve(strict=True)
    parent,identity=control(control_root);startup.require_allocation(parent['ranks'])
    build_root=build_root.resolve(strict=True);binary=builder.verify(build_root)
    bm=json.loads((build_root/'force_build_manifest.json').read_text())
    if bm['base_build']!=parent['build']:raise ValueError('Probe must derive from the control VOF build')
    protected=[control_root,build_root,Path(parent['gradient_control']['root']),Path(parent['gradient_control']['spatial_root']),base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','control')]
    out=base.fresh(out,protected);source=control_root/'quiescent';expected,text=expected_inputs(parent,source)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent'
    for rel in parent['inputs']:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(regular(source,rel),p)
    (case/'system/controlDict').write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(regular(source,name),case/name)
    if base.case_input_hashes(case)!=expected or control(control_root)!=(parent,identity):
        raise ValueError('Control or copied inputs changed during preparation')
    m=dict(parent,root=str(out),variant=TAG,inputs=expected,force_control=identity,
        build=str(build_root),assets=binary['assets'],force_identity=binary,
        force_build_manifest_sha256=base.sha(build_root/'force_build_manifest.json'),force_recipe=recipe(),source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_FORCE_PROBE='+str(out)+' (fresh t=0; two steps; passive face-force observer; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads(regular(root,'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['force_recipe']!=recipe():
        raise ValueError('Wrong/changed force probe')
    c=m['force_control'];parent,identity=control(Path(c['root']))
    binary=builder.verify(Path(m['build']))
    bm=json.loads((Path(m['build'])/'force_build_manifest.json').read_text())
    if binary!=m['force_identity'] or binary['assets']!=m['assets'] or bm['base_build']!=parent['build']:
        raise ValueError('Force build/control changed')
    if base.sha(Path(m['build'])/'force_build_manifest.json')!=m['force_build_manifest_sha256']:
        raise ValueError('Force build manifest changed')
    expected,_=expected_inputs(parent,Path(c['root'])/'quiescent')
    if c!=identity or m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Force probe/control identity changed')
    for key in parent:
        if key not in ('root','variant','inputs','source_commit','build','assets') and m[key]!=parent[key]:
            raise ValueError('Candidate metadata changed: '+key)
    return m


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted force probe')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=os.environ.copy();app=m['force_identity']['binary']
    try:
        base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
        from sc1_reactor_watch import execute
        execute(['srun','--kill-on-bad-exit=1','--time=00:10:00','-n',str(m['ranks']),
                 '--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated'],case,case/'log.solver',env,wall_seconds=600)
        control_log=regular(Path(m['force_control']['root']),'quiescent/log.solver')
        if base.sha(control_log)!=CONTROL_LOG:raise ValueError('Control log changed during force probe')
        report=reporting.analyze((case/'log.solver').read_text(),control_log.read_text())
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Force probe initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'),control=m['force_control'])
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k not in ('cells','faces','faceCoverage')},indent=2),flush=True)
        print(json.dumps(dict(firstAndLastPeaks=[r for r in report['cells'] if r['sequence'] in (1,8)]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_FORCE_PROBE=COMPLETE (diagnostic only; no startup acceptance)',flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('build');q.add_argument('--base-build',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('prepare');q.add_argument('--control-root',type=Path,required=True);q.add_argument('--build-root',type=Path,required=True)
    q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='build':builder.build(a.base_build,a.out)
    elif a.action=='prepare':prepare(a.control_root,a.build_root,a.out)
    else:run(a.root)
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('VOF force probe stopped: '+str(error)+'; preserve all roots',file=sys.stderr)
        raise SystemExit(1)
