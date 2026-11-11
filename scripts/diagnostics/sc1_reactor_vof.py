#!/usr/bin/env python3
"""Build and run a separate single-velocity SC1 free-surface reference."""
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
import sc1_reactor_velocity_probe as velocity
import sc1_reactor_vof_build as builder
import sc1_reactor_vof_case as cases
import sc1_reactor_vof_report as reporting

TAG='sc1_meml_vof_reference_v1'
CONTROL_LOG='c0cd5cfde42651b23262ad2cc545cf8f5b66336cc64d4e95e3bc1e11a0f967e8'
CONTROL_RESULT='269ebee5ab081b098caec86e6fcd7a56dcd1852f72ce99d12b5653508a49c21d'
COPIED=('0/alpha.liquid','0/p_rgh','0/p','0/C','0/V','constant/g','system/decomposeParDict')


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(cases.__file__),Path(reporting.__file__))}


def control(root):
    root=root.resolve(strict=True);m,parent=velocity.verify(root)
    if base.sha(root/'quiescent/log.solver')!=CONTROL_LOG or base.sha(root/'quiescent/result.json')!=CONTROL_RESULT:
        raise ValueError('Requires the received unchanged spatial diagnostic')
    report=json.loads((root/'quiescent/result.json').read_text())
    if not report['diagnostic_complete']:raise ValueError('Incomplete spatial control')
    if json.loads((root/'quiescent/attempt.json').read_text())!={k:m[k] for k in ('build','assets','ranks')}:
        raise ValueError('Spatial control attempt identity changed')
    velocity.compare_prefix((root/'quiescent/log.solver').read_text(),
        (Path(parent['root'])/'quiescent/log.solver').read_text())
    for name,digest in report['fieldHashes'].items():
        if base.sha(root/'quiescent'/name)!=digest:raise ValueError('Saved spatial field changed: '+name)
    return m,dict(root=str(root),manifest_sha256=base.sha(root/'startup_manifest.json'),
        log_sha256=CONTROL_LOG,result_sha256=CONTROL_RESULT)


def expected_inputs(parent):
    copied={k:h for k,h in parent['inputs'].items() if k in COPIED or k.startswith('constant/polyMesh/')}
    if not set(COPIED).issubset(copied):raise ValueError('Missing preserved reference input')
    return dict(copied,**{p:hashlib.sha256(t.encode()).hexdigest() for p,t in cases.case_files().items()})


def prepare(control_root,build_root,out):
    control_root=control_root.resolve(strict=True);build_root=build_root.resolve(strict=True)
    parent,identity=control(control_root);startup.require_allocation(parent['ranks'])
    binary=builder.verify(build_root)
    protected=[control_root,build_root,base.REPO]+[Path(parent[k]) for k in ('build','archive','regression_root')]
    protected += [Path(parent[k]['root']) for k in ('baseline','failed_control','face_control','coupled_control','pressure_probe','velocity_control')]
    out=base.fresh(out,protected);expected=expected_inputs(parent)
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent';source=control_root/'quiescent'
    shutil.copytree(source/'constant/polyMesh',case/'constant/polyMesh')
    for rel in COPIED:
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/rel,p)
    for rel,text in cases.case_files().items():
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):shutil.copy2(source/name,case/name)
    if base.case_input_hashes(case)!=expected or control(control_root)!=(parent,identity):
        raise ValueError('VOF preparation changed preserved inputs/control')
    m=dict(root=str(out),variant=TAG,build=str(build_root),assets=binary['assets'],ranks=parent['ranks'],
        archive=parent['archive'],regression_root=parent['regression_root'],baseline=parent['baseline'],
        initial=parent['initial'],mesh_hashes=parent['mesh_hashes'],inputs=expected,control=identity,
        reference_identity=binary,recipe=recipe(),build_manifest_sha256=base.sha(build_root/'vof_build_manifest.json'),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',m)
    print('VOF_REFERENCE='+str(out)+' (fresh unchanged mesh/fill; single velocity; maximum 10 ms; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads((root/'startup_manifest.json').read_text())
    if m['root']!=str(root) or m['variant']!=TAG or m['recipe']!=recipe():raise ValueError('Wrong/changed VOF case recipe')
    parent,identity=control(Path(m['control']['root']))
    binary=builder.verify(Path(m['build']))
    if identity!=m['control'] or binary!=m['reference_identity'] or binary['assets']!=m['assets']:
        raise ValueError('VOF reference/control identity changed')
    if base.sha(Path(m['build'])/'vof_build_manifest.json')!=m['build_manifest_sha256']:
        raise ValueError('VOF build manifest changed')
    for k in ('ranks','archive','regression_root','baseline','initial','mesh_hashes'):
        if m[k]!=parent[k]:raise ValueError('VOF reference changed control metadata: '+k)
    if m['inputs']!=expected_inputs(parent) or base.case_input_hashes(root/'quiescent')!=m['inputs']:
        raise ValueError('VOF reference case inputs changed')
    return m


def run(root):
    root=root.resolve(strict=True);m=verify(root);startup.require_allocation(m['ranks']);case=root/'quiescent'
    if {p.name for p in case.iterdir()}!={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}:
        raise ValueError('Refusing an attempted VOF reference')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=os.environ.copy();app=m['reference_identity']['binary'];report=None
    try:
        if m['ranks']>1:
            base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
            command=['srun','--kill-on-bad-exit=1','--time=00:30:00','-n',str(m['ranks']),'--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated']
        else:command=[app,'-fileHandler','uncollated']
        from sc1_reactor_watch import execute
        execute(command,case,case/'log.solver',env)
        reporting.budget((case/'log.solver').read_text(),m['initial'])
        if m['ranks']>1:
            for rank in range(m['ranks']):
                for i in range(5,startup.STEPS+1,5):
                    for field in ('U','alpha.liquid','p_rgh'):
                        p=case/f'processor{rank}'/f'{i*startup.DT:.12g}'/field
                        if p.is_symlink() or not p.is_file():raise ValueError('Missing/redirected VOF output '+str(p))
            base.execute(['reconstructPar','-time','0.00025:0.01','-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
        report=reporting.analyze(case,m['initial'],startup.CELLS)
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('VOF initial inputs changed')
        report.update(suite=TAG,log_sha256=base.sha(case/'log.solver'))
        base.write_json(case/'result.json',report)
        print(json.dumps({k:v for k,v in report.items() if k!='snapshots'},indent=2),flush=True)
        print(json.dumps(dict(finalSnapshot=report['snapshots'][-1]),indent=2),flush=True)
    finally:
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('VOF_REFERENCE='+('PASS' if report['reference_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='action',required=True)
    q=s.add_parser('build');q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('prepare');q.add_argument('--control-root',type=Path,required=True);q.add_argument('--build-root',type=Path,required=True);q.add_argument('--out',type=Path,required=True)
    q=s.add_parser('run');q.add_argument('--root',type=Path,required=True);q.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='build':builder.build(a.out)
    elif a.action=='prepare':prepare(a.control_root,a.build_root,a.out)
    else:return 0 if run(a.root)['reference_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 VOF reference stopped: '+str(error)+'; preserve all roots',file=sys.stderr);raise SystemExit(1)
