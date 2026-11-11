#!/usr/bin/env python3
"""Fresh SC1 VOF nonorthogonal pressure convergence comparison."""
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
import sc1_reactor_vof_gradients as gradients
import sc1_reactor_vof_saved as saved
from sc1_reactor_vof_failure import scan
import sc1_reactor_vof_nonorth_report as reporting
from sc1_reactor_failure import regular

TAG='sc1_meml_vof_nonorthogonal_v3'
CONTROL_LOG='5c7601fc23d7e8261f8f02910104ed9ad3b0bc8179e026861ba7ba8b4a710fc8'


def solution(text):
    anchor='nNonOrthogonalCorrectors 0;'
    if text.count(anchor)!=1:raise ValueError('Unexpected control pressure settings')
    return text.replace(anchor,'nNonOrthogonalCorrectors 9;')


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(reporting.__file__))}


def control(root):
    root=root.resolve(strict=True);m=gradients.verify(root)
    log=regular(root,'quiescent/log.solver')
    if base.sha(log)!=CONTROL_LOG:raise ValueError('Requires the received failed gradient control')
    if json.loads(regular(root,'quiescent/attempt.json').read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('Gradient control attempt changed')
    evidence=scan(log)
    if evidence['malformedRecordCount'] or evidence['duplicateStepRecords'] or evidence['validStepRecordCount']!=168:
        raise ValueError('Unexpected gradient control evidence')
    if evidence['firstLimitEvents'].get('liquidMax_m_s',{}).get('t_s',0)<.008399999999:
        raise ValueError('Unexpected gradient control stop')
    saved.budget(log.read_text(),m['initial'])
    return m,dict(root=str(root),manifest_sha256=base.sha(root/'startup_manifest.json'),log_sha256=CONTROL_LOG)


def expected_inputs(parent,source):
    text=solution(regular(source,'system/fvSolution').read_text())
    expected=dict(parent['inputs']);expected['system/fvSolution']=hashlib.sha256(text.encode()).hexdigest()
    return expected,text


def prepare(control_root,out):
    control_root=control_root.resolve(strict=True)
    parent,identity=control(control_root);startup.require_allocation(parent['ranks'])
    protected=[control_root,Path(parent['gradient_control']['root']),Path(parent['gradient_control']['spatial_root']),base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','control')]
    out=base.fresh(out,protected);source=control_root/'quiescent';expected,text=expected_inputs(parent,source)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent'
    for rel in parent['inputs']:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(regular(source,rel),p)
    (case/'system/fvSolution').write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(regular(source,name),case/name)
    if base.case_input_hashes(case)!=expected or control(control_root)!=(parent,identity):
        raise ValueError('Control or copied inputs changed during preparation')
    m=dict(parent,root=str(out),variant=TAG,inputs=expected,nonorth_control=identity,
        nonorth_recipe=recipe(),source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_NONORTH='+str(out)+' (fresh t=0; ten pressure passes per correction; existing binary; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads(regular(root,'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['nonorth_recipe']!=recipe():
        raise ValueError('Wrong/changed nonorthogonal candidate')
    c=m['nonorth_control'];parent,identity=control(Path(c['root']))
    expected,_=expected_inputs(parent,Path(c['root'])/'quiescent')
    if c!=identity or m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Nonorthogonal candidate/control identity changed')
    for key in parent:
        if key not in ('root','variant','inputs','source_commit') and m[key]!=parent[key]:
            raise ValueError('Candidate metadata changed: '+key)
    return m


def finish_report(report,steps,control_report):
    # Reuse the verified snapshot comparison and gas-rest check; add convergence.
    report=gradients.finish_report(report,steps,control_report)
    report.pop('gradient_pass')
    report.update(nonorth_pass=report['reference_pass'] and report['gas_rest_pass']
        and report['pressureIteration']['pressure_iteration_pass'],
        meaning='10-ms nonorthogonal pressure comparison only. Original limits plus gas peak <=1 mm/s and converged final pressure passes. No operating-reactor, long-time equilibrium, aeration or biology acceptance.')
    return report


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted nonorthogonal candidate')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=os.environ.copy();app=m['reference_identity']['binary'];report=None
    try:
        base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
        from sc1_reactor_watch import execute
        execute(['srun','--kill-on-bad-exit=1','--time=00:30:00','-n',str(m['ranks']),
                 '--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated'],case,case/'log.solver',env)
        _,steps=reporting.budget((case/'log.solver').read_text(),m['initial'])
        for rank in range(m['ranks']):
            for i in range(5,startup.STEPS+1,5):
                for field in ('U','alpha.liquid','p_rgh'):
                    regular(case,f'processor{rank}/{i*startup.DT:.12g}/{field}')
        base.execute(['reconstructPar','-time','0.00025:0.01','-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
        report=reporting.analyze(case,m['initial'],startup.CELLS)
        control_path=Path(m['gradient_control']['spatial_root'])/'report.json'
        if base.sha(control_path)!=gradients.SPATIAL_SHA:raise ValueError('Spatial report changed during run')
        report=finish_report(report,steps,json.loads(control_path.read_text()))
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Candidate initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'),control=m['nonorth_control'])
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='snapshots'},indent=2),flush=True)
        print(json.dumps(dict(finalSnapshot=report['snapshots'][-1]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            try:
                _,iterations=reporting.pressure_coverage((case/'log.solver').read_text(),require_complete=False)
                iterations['meaning']='Pressure-iteration evidence only; a converged observed prefix is not a completed startup.'
                base.write_json(case/'nonorth_pressure_iterations.json',iterations)
                print(json.dumps(dict(nonorthPressureIterations=iterations),indent=2),flush=True)
            except (ValueError,KeyError) as error:
                print('Pressure-iteration summary unavailable: '+str(error),flush=True)
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_NONORTH='+('PASS' if report['nonorth_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('prepare');q.add_argument('--control-root',type=Path,required=True)
    q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='prepare':prepare(a.control_root,a.out)
    else:return 0 if run(a.root)['nonorth_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('VOF nonorthogonal comparison stopped: '+str(error)+'; preserve all roots',file=sys.stderr)
        raise SystemExit(1)
