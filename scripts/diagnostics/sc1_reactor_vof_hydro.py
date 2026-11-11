#!/usr/bin/env python3
"""Fresh 10-ms VOF startup with a fixed hydrostatic pressure/density reference."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import sc1_two_field as base
import sc1_reactor_startup as startup
import sc1_reactor_startup_case as fields
import sc1_reactor_vof_force as force
import sc1_reactor_vof_hydro_build as builder
import sc1_reactor_vof_hydro_report as reporting
from instrument_sc1_alpha_subcycles import _once
from sc1_reactor_failure import regular
from sc1_two_field_cases import foam

TAG='sc1_meml_vof_hydrostatic_v1'
FORCE_LOG='91695d8e18a24fc6ee2e9e9d0d5d31f2c5db17f41b407bb3605c0c6e1f93c3ab'
FORCE_MANIFEST='9079a93cd5b1c338e2bdcafd66eef304e3fc8c77c29d154ed9fbd12fd1dc1e0c'


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(reporting.__file__))}


def control(root):
    root=root.resolve(strict=True)
    if base.sha(regular(root,'startup_manifest.json'))!=FORCE_MANIFEST or base.sha(regular(root,'quiescent/log.solver'))!=FORCE_LOG:
        raise ValueError('Requires the received completed force probe')
    fm=force.verify(root)
    if json.loads(regular(root,'quiescent/attempt.json').read_text())!={k:fm[k] for k in ('build','assets','ranks')}:
        raise ValueError('Force probe attempt changed')
    parent,identity=force.control(Path(fm['force_control']['root']))
    result=force.reporting.analyze(regular(root,'quiescent/log.solver').read_text(),
        regular(Path(identity['root']),'quiescent/log.solver').read_text())
    stored=json.loads(regular(root,'quiescent/result.json').read_text())
    if any(stored.get(k)!=v for k,v in result.items()):raise ValueError('Force result disagrees with logs')
    return parent,dict(root=str(root),build=fm['build'],manifest_sha256=FORCE_MANIFEST,log_sha256=FORCE_LOG,
        result_sha256=base.sha(root/'quiescent/result.json'),nonorth=identity)


def expected_inputs(parent,source):
    changes={}
    changes['system/controlDict']=_once(regular(source,'system/controlDict').read_text(),
        'application sc1InterFoamReference;','application '+builder.APP+';')
    changes['system/fvSchemes']=_once(regular(source,'system/fvSchemes').read_text(),
        'grad(p_rgh) leastSquares;','grad(p_rgh) leastSquares;\n    grad(pPert) leastSquares;')
    text=regular(source,'system/fvSolution').read_text()
    text=_once(text,'    p_rgh {','    pPert {')
    changes['system/fvSolution']=_once(text,'    p_rghFinal { $p_rgh; }','    pPertFinal { $pPert; }')
    changes['0/pPert']=fields.field('p_rgh','[1 -1 -2 0 0 0 0]','uniform 0').replace('object p_rgh;','object pPert;').replace('value uniform 100000;','value uniform 0;')
    initial=parent['initial']
    if initial['transition_width_m']!=fields.WIDTH:raise ValueError('Unexpected reference width')
    changes['constant/sc1HydrostaticReference']=foam('sc1HydrostaticReference',
        f"interfaceHeight {initial['interface_midheight_m']:.17g};\ntransitionWidth {fields.WIDTH:.17g};")
    expected=dict(parent['inputs']);expected.update({k:hashlib.sha256(v.encode()).hexdigest() for k,v in changes.items()})
    return expected,changes


def check_initial(source,initial):
    n=initial['cells'];c=fields.internal(source/'0/C',n,True);v=fields.internal(source/'0/V',n)
    alpha=fields.internal(source/'0/alpha.liquid',n);prgh=fields.internal(source/'0/p_rgh',n)
    y=c[:,1];h=initial['interface_midheight_m'];w=fields.WIDTH;x=y-h
    expected_alpha=np.clip((h+w/2-y)/w,0,1)
    integral=np.where(x<-w/2,x+w/8,np.where(x>w/2,w/8,x/2-x*x/(2*w)))
    ph=100000-fields.G*y+999*fields.G*(w/8-integral)
    if not np.allclose(alpha,expected_alpha,rtol=0,atol=1e-12) or not np.allclose(prgh-(1+999*alpha)*fields.G*y,ph,rtol=0,atol=1e-7):
        raise ValueError('Input is not the prescribed analytic hydrostatic fill')
    if np.any(v<=0) or abs(np.dot(alpha,v)-1)>1e-12 or np.max(np.abs(fields.internal(source/'0/U',n,True)))!=0:
        raise ValueError('Expected exactly 1000 L at rest')


def prepare(force_root,build_root,out):
    force_root=force_root.resolve(strict=True);parent,identity=control(force_root)
    startup.require_allocation(parent['ranks']);build_root=build_root.resolve(strict=True)
    binary=builder.verify(build_root);bm=json.loads((build_root/'hydro_build_manifest.json').read_text())
    if bm['base_build']!=parent['build']:raise ValueError('Candidate must derive from control VOF build')
    source=Path(identity['nonorth']['root'])/'quiescent';check_initial(source,parent['initial'])
    protected=[force_root,Path(identity['build']),source.parent,build_root,base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','control','gradient_control')]
    protected.append(Path(parent['gradient_control']['spatial_root']))
    out=base.fresh(out,protected);expected,changes=expected_inputs(parent,source)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent'
    for rel in parent['inputs']:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(regular(source,rel),p)
    for rel,text in changes.items():(case/rel).write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(regular(source,name),case/name)
    if base.case_input_hashes(case)!=expected or control(force_root)!=(parent,identity):
        raise ValueError('Control or copied inputs changed during preparation')
    m=dict(parent,root=str(out),variant=TAG,inputs=expected,hydro_control=identity,
        build=str(build_root),assets=binary['assets'],hydro_identity=binary,
        hydro_build_manifest_sha256=base.sha(build_root/'hydro_build_manifest.json'),hydro_recipe=recipe(),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_HYDROSTATIC='+str(out)+' (fresh 1000-L fill; fixed hydrostatic reference; maximum 10 ms; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads(regular(root,'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['hydro_recipe']!=recipe():raise ValueError('Wrong/changed hydrostatic candidate')
    parent,identity=control(Path(m['hydro_control']['root']));binary=builder.verify(Path(m['build']))
    bm=json.loads((Path(m['build'])/'hydro_build_manifest.json').read_text())
    if binary!=m['hydro_identity'] or binary['assets']!=m['assets'] or bm['base_build']!=parent['build'] or base.sha(Path(m['build'])/'hydro_build_manifest.json')!=m['hydro_build_manifest_sha256']:
        raise ValueError('Hydrostatic binary/build changed')
    expected,_=expected_inputs(parent,Path(identity['nonorth']['root'])/'quiescent')
    if identity!=m['hydro_control'] or m['inputs']!=expected or base.case_input_hashes(root/'quiescent')!=expected:
        raise ValueError('Hydrostatic/control inputs changed')
    for k in parent:
        if k not in ('root','variant','inputs','source_commit','build','assets') and m[k]!=parent[k]:
            raise ValueError('Candidate metadata changed: '+k)
    return m


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted hydrostatic candidate')
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
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Hydrostatic initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'),control=m['hydro_control'])
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='snapshots'},indent=2),flush=True)
        print(json.dumps(dict(finalSnapshot=report['snapshots'][-1]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_HYDROSTATIC='+('PASS' if report['hydrostatic_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('build');q.add_argument('--base-build',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('prepare');q.add_argument('--force-root',type=Path,required=True);q.add_argument('--build-root',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='build':builder.build(a.base_build,a.out)
    elif a.action=='prepare':prepare(a.force_root,a.build_root,a.out)
    else:return 0 if run(a.root)['hydrostatic_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 hydrostatic startup stopped: '+str(error)+'; preserve all roots',file=sys.stderr)
        raise SystemExit(1)
