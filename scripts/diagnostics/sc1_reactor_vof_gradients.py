#!/usr/bin/env python3
"""Fresh SC1 VOF scalar-gradient comparison using the existing verified binary."""
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
import sc1_reactor_vof as reference
import sc1_reactor_vof_saved as saved
import sc1_reactor_vof_report as reporting
from sc1_reactor_failure import regular

TAG='sc1_meml_vof_scalar_gradients_v2'
SPATIAL_SHA='f99b20ce7019e4d66a142532885dacc0761289d73174b20b2dc4eb332da0ddcf'
GAS_REST_LIMIT=1e-3


def schemes(text):
    anchor='gradSchemes { default Gauss linear; }'
    if text.count(anchor)!=1:raise ValueError('Unexpected control gradient dictionary')
    return text.replace(anchor,'''gradSchemes {
    default Gauss linear;
    grad(p_rgh) leastSquares;
    grad(rho) leastSquares;
    grad(alpha.liquid) leastSquares;
    nHat leastSquares;
}''')


def recipe():
    return {Path(__file__).name:base.sha(Path(__file__))}


def control(root,spatial_root):
    root=root.resolve(strict=True);spatial_root=spatial_root.resolve(strict=True)
    if base.sha(regular(root,'startup_manifest.json'))!=saved.MANIFEST_SHA or base.sha(regular(root,'quiescent/log.solver'))!=saved.LOG_SHA:
        raise ValueError('Requires the received VOF control')
    m=reference.verify(root)
    rp=regular(spatial_root,'report.json')
    if base.sha(rp)!=SPATIAL_SHA:raise ValueError('Requires the received saved-field report')
    report=json.loads(rp.read_text())
    if report['suite']!='sc1_vof_saved_spatial_v1' or report['sourceRoot']!=str(root):
        raise ValueError('Wrong spatial report source')
    for rel,digest in report['sourceHashes'].items():
        if base.sha(regular(root,rel))!=digest:raise ValueError('Spatial source identity changed: '+rel)
    for rel,digest in report['fieldHashes'].items():
        if base.sha(regular(root/'quiescent',rel))!=digest:raise ValueError('Saved control field changed: '+rel)
    if json.loads(regular(root,'quiescent/attempt.json').read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('VOF control attempt changed')
    return m,dict(root=str(root),manifest_sha256=saved.MANIFEST_SHA,log_sha256=saved.LOG_SHA,
        spatial_root=str(spatial_root),spatial_report_sha256=SPATIAL_SHA)


def expected_inputs(parent,source):
    text=schemes(regular(source,'system/fvSchemes').read_text())
    expected=dict(parent['inputs']);expected['system/fvSchemes']=hashlib.sha256(text.encode()).hexdigest()
    return expected,text


def prepare(control_root,spatial_root,out):
    control_root=control_root.resolve(strict=True);spatial_root=spatial_root.resolve(strict=True)
    parent,identity=control(control_root,spatial_root);startup.require_allocation(parent['ranks'])
    protected=[control_root,spatial_root,base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','control')]
    out=base.fresh(out,protected);source=control_root/'quiescent';expected,text=expected_inputs(parent,source)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent'
    for rel in parent['inputs']:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(regular(source,rel),p)
    (case/'system/fvSchemes').write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(regular(source,name),case/name)
    if base.case_input_hashes(case)!=expected or control(control_root,spatial_root)!=(parent,identity):
        raise ValueError('Control or copied inputs changed during preparation')
    m=dict(parent,root=str(out),variant=TAG,inputs=expected,gradient_control=identity,
        gradient_recipe=recipe(),source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_GRADIENTS='+str(out)+' (fresh t=0; four scalar-gradient entries; existing binary; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads(regular(root,'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['gradient_recipe']!=recipe():
        raise ValueError('Wrong/changed gradient candidate')
    c=m['gradient_control'];parent,identity=control(Path(c['root']),Path(c['spatial_root']))
    expected,_=expected_inputs(parent,Path(c['root'])/'quiescent')
    if c!=identity or m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Gradient candidate/control identity changed')
    for key in parent:
        if key not in ('root','variant','inputs','source_commit') and m[key]!=parent[key]:
            raise ValueError('Candidate metadata changed: '+key)
    return m


def finish_report(report,steps,control_report):
    by_time={round(s['time_s'],12):s for s in report['snapshots']}
    by_step={round(s['t_s'],12):s for s in steps}
    for t,snapshot in by_time.items():
        if abs(snapshot['phases']['gas']['present']['max_m_s']-by_step[t]['gasMax_m_s'])>1e-10:
            raise ValueError('Native gas speed disagrees with candidate audit')
    comparison=[]
    for old in control_report['snapshots']:
        new=by_time[round(old['time_s'],12)]
        comparison.append(dict(time_s=old['time_s'],phases={phase:dict(
            controlMax_m_s=old['phases'][phase]['present']['max_m_s'],
            candidateMax_m_s=new['phases'][phase]['present']['max_m_s'],
            controlRms_m_s=old['phases'][phase]['present']['rms_m_s'],
            candidateRms_m_s=new['phases'][phase]['present']['rms_m_s']) for phase in ('liquid','gas')}))
    gas_pass=report['maxGas_m_s']<=GAS_REST_LIMIT
    report.update(comparison=comparison,gasRestLimit_m_s=GAS_REST_LIMIT,gas_rest_pass=gas_pass,
        gradient_pass=report['reference_pass'] and gas_pass,
        meaning='10-ms scalar-gradient comparison only. Original startup limits retained; additionally require gas peak <=1 mm/s over all logged steps. No long-time equilibrium, aeration, oxygen-transfer or biology acceptance.')
    return report


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted gradient candidate')
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
        if base.sha(control_path)!=SPATIAL_SHA:raise ValueError('Spatial report changed during run')
        report=finish_report(report,steps,json.loads(control_path.read_text()))
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Candidate initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'),control=m['gradient_control'])
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='snapshots'},indent=2),flush=True)
        print(json.dumps(dict(finalSnapshot=report['snapshots'][-1]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_GRADIENTS='+('PASS' if report['gradient_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('prepare');q.add_argument('--control-root',type=Path,required=True)
    q.add_argument('--spatial-root',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='prepare':prepare(a.control_root,a.spatial_root,a.out)
    else:return 0 if run(a.root)['gradient_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('VOF gradient comparison stopped: '+str(error)+'; preserve all roots',file=sys.stderr)
        raise SystemExit(1)
