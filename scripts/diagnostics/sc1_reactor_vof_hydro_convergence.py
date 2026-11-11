#!/usr/bin/env python3
"""Fresh hydrostatic startup with 24 pressure passes and unchanged acceptance limits."""
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
import sc1_reactor_vof_hydro as hydro
import sc1_reactor_vof_hydro_convergence_report as reporting
from instrument_sc1_alpha_subcycles import _once
from sc1_reactor_failure import regular

TAG='sc1_meml_vof_hydrostatic_convergence_v2'
CONTROL_LOG='156073fefe79f07866cbfd65c62e072cf7192145e5727148e9d3d9b4ee7cef4a'
CONTROL_MANIFEST='c8b27f64573dac4ed63d5936f71a2e9a89caa759135d600eb231cbfd13add9aa'
CONTROL_RESULT='6aca0527e3a02ffa9474d9840d9f43e8d075af4847f4575d1069cb87a557869a'


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(reporting.__file__))}


def control(root):
    root=root.resolve(strict=True)
    for name,digest in (('startup_manifest.json',CONTROL_MANIFEST),('quiescent/log.solver',CONTROL_LOG),('quiescent/result.json',CONTROL_RESULT)):
        if base.sha(regular(root,name))!=digest:raise ValueError('Requires the received hydrostatic control: '+name)
    m=hydro.verify(root)
    if json.loads(regular(root,'quiescent/attempt.json').read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('Hydrostatic control attempt changed')
    budget,steps=hydro.reporting.budget(regular(root,'quiescent/log.solver').read_text(),m['initial'])
    stored=json.loads(regular(root,'quiescent/result.json').read_text())
    if any(stored.get(k)!=v for k,v in budget.items() if k!='errors') or budget['errors'] or stored['errors']!=['Nonorthogonal pressure corrections did not converge']:
        raise ValueError('Control budget/evidence disagrees')
    if stored['hydrostatic_pass'] or budget['pressureIteration']['convergedCorrections']!=597 or budget['maxGas_m_s']>1e-3:
        raise ValueError('Unexpected hydrostatic control outcome')
    if len(stored['snapshots'])!=40:raise ValueError('Incomplete control snapshot report')
    return m,dict(root=str(root),manifest_sha256=CONTROL_MANIFEST,log_sha256=CONTROL_LOG,result_sha256=CONTROL_RESULT)


def expected_inputs(parent,source):
    text=_once(regular(source,'system/fvSolution').read_text(),
        'nNonOrthogonalCorrectors 9;',f'nNonOrthogonalCorrectors {reporting.PASSES-1};')
    expected=dict(parent['inputs']);expected['system/fvSolution']=hashlib.sha256(text.encode()).hexdigest()
    return expected,text


def prepare(control_root,out):
    control_root=control_root.resolve(strict=True);parent,identity=control(control_root)
    startup.require_allocation(parent['ranks'])
    protected=[control_root,base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','control','gradient_control')]
    protected += [Path(parent['gradient_control']['spatial_root']),Path(parent['hydro_control']['root']),
        Path(parent['hydro_control']['build']),Path(parent['hydro_control']['nonorth']['root'])]
    out=base.fresh(out,protected);source=control_root/'quiescent';expected,text=expected_inputs(parent,source)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent'
    for rel in parent['inputs']:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(regular(source,rel),p)
    (case/'system/fvSolution').write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(regular(source,name),case/name)
    if base.case_input_hashes(case)!=expected or control(control_root)!=(parent,identity):
        raise ValueError('Control or copied inputs changed during preparation')
    m=dict(parent,root=str(out),variant=TAG,inputs=expected,convergence_control=identity,
        convergence_recipe=recipe(),source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_HYDRO_CONVERGENCE='+str(out)+' (fresh t=0; 24 pressure passes; existing binary; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads(regular(root,'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['convergence_recipe']!=recipe():raise ValueError('Wrong/changed convergence candidate')
    parent,identity=control(Path(m['convergence_control']['root']))
    expected,_=expected_inputs(parent,Path(identity['root'])/'quiescent')
    if m['convergence_control']!=identity or m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Candidate/control input identity changed')
    for k in parent:
        if k not in ('root','variant','inputs','source_commit') and m[k]!=parent[k]:raise ValueError('Candidate metadata changed: '+k)
    return m


def compare(report,control_report):
    if len(report['snapshots'])!=40 or len(control_report['snapshots'])!=40:raise ValueError('Incomplete comparison snapshots')
    rows=[]
    for new,old in zip(report['snapshots'],control_report['snapshots']):
        if new['time_s']!=old['time_s']:raise ValueError('Comparison times differ')
        rows.append(dict(time_s=new['time_s'],liquidVolumeDifference_L=new['liquidVolume_L']-old['liquidVolume_L'],
            **{phase+'PeakDifference_m_s':new['phases'][phase]['present']['max_m_s']-old['phases'][phase]['present']['max_m_s'] for phase in ('liquid','gas')}))
    report['comparison']=rows
    report['comparisonMeaning']='Matched-time differences against the preserved ten-pass control; no tolerance is relaxed based on agreement.'
    return report


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted convergence candidate')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=os.environ.copy();report=None
    try:
        base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
        from sc1_reactor_watch import execute
        execute(['srun','--kill-on-bad-exit=1','--time=00:30:00','-n',str(m['ranks']),
                 '--cpus-per-task=1',m['hydro_identity']['binary'],'-parallel','-fileHandler','uncollated'],case,case/'log.solver',env)
        reporting.budget((case/'log.solver').read_text(),m['initial'])
        for rank in range(m['ranks']):
            for i in range(5,startup.STEPS+1,5):
                for field in ('U','alpha.liquid','p_rgh','p','pPert'):
                    regular(case,f'processor{rank}/{i*startup.DT:.12g}/{field}')
        base.execute(['reconstructPar','-time','0.00025:0.01','-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
        report=reporting.analyze(case,m['initial'],startup.CELLS)
        control_file=regular(Path(m['convergence_control']['root']),'quiescent/result.json')
        if base.sha(control_file)!=CONTROL_RESULT:raise ValueError('Control result changed during run')
        report=compare(report,json.loads(control_file.read_text()))
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Candidate initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'),control=m['convergence_control'])
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k not in ('snapshots','comparison')},indent=2),flush=True)
        print(json.dumps(dict(finalSnapshot=report['snapshots'][-1],finalComparison=report['comparison'][-1]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_HYDRO_CONVERGENCE='+('PASS' if report['hydrostatic_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('prepare');q.add_argument('--control-root',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='prepare':prepare(a.control_root,a.out)
    else:return 0 if run(a.root)['hydrostatic_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 hydrostatic convergence stopped: '+str(error)+'; preserve all roots',file=sys.stderr)
        raise SystemExit(1)
