#!/usr/bin/env python3
"""Prepare/run one fresh audited, quiescent MEML-mesh startup; never restart archives."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import numpy as np
import sc1_two_field as base
import sc1_two_field_rates as rates
from sc1_two_field_rates_instrument import APP, RECIPE
from sc1_two_field_report import AuditContract, summarize_contract, parse
from sc1_two_field_rates_report import validate_rates
from sc1_wall_mesh import zone_adjacency
from sc1_reactor_startup_case import DT, STEPS, END, WRITE, PATCHES, case_files, internal, initialized_fields

TAG='sc1_meml_quiescent_startup_v1'
DRAG_TAG='sc1_meml_quiescent_segregated_drag_v2'
ARCHIVE=Path('/projects/microenviro/final_runs/reacting_Julia_08282026/meml/putida_sc1_final')
CELLS=281664
PINNED=dict(owner='b93ddd6e9322acd08b7c79cdf3fd346b4ba5bef8b92e803446c1fa91f5740b99',
    neighbour='57f913518890989d2bfabd90c71bb82c2e1713f05d0e9194bc250abd0a099149',
    cellZones='fede00ab610f79f958a06bd21519a8e47eebcac0b6259d0cf560955daa8d39b3',
    boundary='9e1312243892baaca3f5dff5553c17b6fc91df68e152181c9a969e1a1e08944e')


def tree_hashes(root):
    if root.is_symlink() or not root.is_dir():raise ValueError('Missing/redirected input tree')
    result={}
    for p in sorted(root.rglob('*')):
        if p.is_symlink():raise ValueError('Symlink in input tree: '+str(p))
        if p.is_file():result[p.relative_to(root).as_posix()]=base.sha(p)
    return result


def regression(root):
    m,original=rates.verify(root)
    result=rates.comparison(root,m,original)
    if not result['regression_pass']:raise ValueError('All eight rate regressions must pass first')
    return m,hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()


def prepare(regression_root, archive, out, ranks):
    regression_root=regression_root.resolve(strict=True); archive=archive.resolve(strict=True)
    generated=case_files(ranks)
    m,proof=regression(regression_root)
    mesh=archive/'constant/polyMesh'; original=tree_hashes(mesh)
    adjacency=zone_adjacency(mesh)
    if {k:v['sha256'] for k,v in adjacency['inputs'].items()}!=PINNED:
        raise ValueError('Mesh addressing differs from completed MEML wall evidence')
    if {p['patch']:p['type'] for p in adjacency['patches']}!=PATCHES:
        raise ValueError('Unexpected mesh boundaries')
    out=base.fresh(out,(archive,regression_root,Path(m['build']),base.REPO))
    out.mkdir(parents=True,exist_ok=False); case=out/'quiescent'
    for rel,text in generated.items():
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x') as f:f.write(text)
    shutil.copytree(mesh,case/'constant/polyMesh')
    if tree_hashes(mesh)!=original or tree_hashes(case/'constant/polyMesh')!=original:
        raise ValueError('Mesh changed while copied; preserve new root')
    env=base.environment(Path(m['build']))
    base.execute(['checkMesh','-fileHandler','uncollated'],case,case/'log.checkMesh',env)
    if 'Mesh OK.' not in (case/'log.checkMesh').read_text():raise ValueError('Mesh check failed')
    for function in ('writeCellCentres','writeCellVolumes'):
        base.execute(['postProcess','-func',function,'-time','0','-fileHandler','uncollated'],case,case/('log.'+function),env)
    centres=internal(case/'0/C',CELLS,True); volumes=internal(case/'0/V',CELLS)
    fields,initial=initialized_fields(centres,volumes)
    if not .99<initial['interface_midheight_m']<1.05:
        raise ValueError('Unexpected height for prescribed 1000-L fill')
    for rel,text in fields.items():(case/rel).write_text(text)
    manifest=dict(variant=TAG,root=str(out),archive=str(archive),mesh_hashes=original,
        build=m['build'],assets=m['assets'],regression_root=str(regression_root),regression_sha256=proof,
        ranks=ranks,initial=initial,inputs=base.case_input_hashes(case),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',manifest)
    print('REACTOR_STARTUP='+str(out)+' (native 1000-L fill; no CFD launched)',flush=True)
    print(json.dumps(initial,indent=2),flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True); m=json.loads((root/'startup_manifest.json').read_text())
    if m['variant'] not in (TAG,DRAG_TAG) or m['root']!=str(root):raise ValueError('Wrong startup manifest')
    if base.verify_build(Path(m['build']),recipe=RECIPE)['assets']!=m['assets']:
        raise ValueError('Startup build changed')
    _,proof=regression(Path(m['regression_root']))
    if proof!=m['regression_sha256']:raise ValueError('Regression evidence changed')
    if tree_hashes(Path(m['archive'])/'constant/polyMesh')!=m['mesh_hashes']:
        raise ValueError('Archived mesh changed')
    case=root/'quiescent'
    if case.is_symlink() or base.case_input_hashes(case)!=m['inputs']:
        raise ValueError('Startup inputs changed')
    generated=case_files(m['ranks'],segregated_drag=m['variant']==DRAG_TAG)
    centres=internal(case/'0/C',CELLS,True);volumes=internal(case/'0/V',CELLS)
    fields,initial=initialized_fields(centres,volumes);generated.update(fields)
    if initial!=m['initial']:raise ValueError('Startup initialization changed')
    for rel,text in generated.items():
        if base.sha(case/rel)!=hashlib.sha256(text.encode()).hexdigest():
            raise ValueError('Startup differs from its recipe: '+rel)
    if 'Mesh OK.' not in (case/'log.checkMesh').read_text():raise ValueError('Missing successful mesh check')
    if m['variant']==DRAG_TAG:
        parent,report=saved_startup(Path(m['baseline']['root']))
        if baseline_identity(Path(m['baseline']['root']),report)!=m['baseline']:
            raise ValueError('Saved startup control changed')
        expected=dict(parent['inputs'])
        expected['constant/phaseProperties']=hashlib.sha256(generated['constant/phaseProperties'].encode()).hexdigest()
        if m['inputs']!=expected or m['initial']!=parent['initial'] or any(
            m[k]!=parent[k] for k in ('ranks','build','assets','archive','mesh_hashes','regression_root','regression_sha256')):
            raise ValueError('Drag comparison changed more than its blending dictionary')
    return m


def baseline_identity(root,report):
    digest=base.sha(root/'quiescent/log.solver')
    if digest!=report['log_sha256']:raise ValueError('Control log changed after validation')
    return dict(root=str(root.resolve()),log_sha256=digest,
        result_sha256=base.sha(root/'quiescent/result.json'),manifest_sha256=base.sha(root/'startup_manifest.json'))


def saved_startup(root):
    """Validate the completed v1 control without rewriting any of its evidence."""
    root=root.resolve(strict=True)
    if root.is_symlink() or json.loads((root/'startup_manifest.json').read_text())['variant']!=TAG:
        raise ValueError('Comparison requires the original v1 startup control')
    m=verify(root);case=root/'quiescent'
    for name in ('attempt.json','log.solver','result.json'):
        if (case/name).is_symlink() or not (case/name).is_file():raise ValueError('Missing/redirected saved startup evidence')
    if json.loads((case/'attempt.json').read_text())!=dict(build=m['build'],assets=m['assets'],ranks=m['ranks']):
        raise ValueError('Startup attempt/build mismatch')
    report=analyze(case,m)
    if report!=json.loads((case/'result.json').read_text()):raise ValueError('Saved startup result/log/fields disagree')
    if not report['screen_pass']:raise ValueError('Startup control must have a valid mass screen')
    if base.sha(case/'log.solver')!=report['log_sha256'] or base.case_input_hashes(case)!=m['inputs']:
        raise ValueError('Startup control changed during validation')
    return m,report


def prepare_drag(baseline_root,out):
    baseline_root=baseline_root.resolve(strict=True);m,report=saved_startup(baseline_root)
    require_allocation(m['ranks'])
    identity=baseline_identity(baseline_root,report)
    out=base.fresh(out,(baseline_root,Path(m['archive']),Path(m['build']),Path(m['regression_root']),base.REPO))
    out.mkdir(parents=True,exist_ok=False);case=out/'quiescent';source=baseline_root/'quiescent'
    for folder in ('0','constant','system'):shutil.copytree(source/folder,case/folder)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):
        shutil.copy2(source/name,case/name)
    phase='constant/phaseProperties'
    text=case_files(m['ranks'],segregated_drag=True)[phase]
    (case/phase).write_text(text)
    expected=dict(m['inputs']);expected[phase]=hashlib.sha256(text.encode()).hexdigest()
    if base.case_input_hashes(case)!=expected or base.case_input_hashes(source)!=m['inputs'] or baseline_identity(baseline_root,report)!=identity:
        raise ValueError('Input/control changed while preparing comparison')
    manifest=dict(m,variant=DRAG_TAG,root=str(out),inputs=expected,baseline=identity,
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',manifest)
    print('REACTOR_STARTUP='+str(out)+' (fresh t=0 copy; only drag blending changed; no CFD launched)',flush=True)
    return out


def contract(initial):
    volumes=(initial['gas_volume_m3'],initial['liquid_volume_m3'])
    return AuditContract(dt=DT,steps=STEPS,volumes=volumes,
        volume_initial_tol=(1e-10,1e-10),mass_initial_tol=(1e-8,1e-7),
        volume_identity_tol=(1e-10,1e-10),mass_identity_tol=(1e-8,1e-7),
        frozen_mass_tol=(1e-8,1e-7),equation_l1_tol=(1e-10,1e-10),
        budget_arithmetic_tol=(1e-10,1e-9),
        budget_absolute=(1e-10+1e-8*volumes[0],1e-10+1e-8*1000*volumes[1]),
        budget_relative=(1e-3,1e-3),option_source_kg_s=(0,0))


def analyze(case, m):
    text=(case/'log.solver').read_text(); spec=contract(m['initial'])
    if m.get('variant',TAG)==DRAG_TAG and text.splitlines().count('Selecting drag blending method: none')!=1:
        raise ValueError('Missing unique runtime selection of segregated drag blending')
    report=summarize_contract(text,'meml_quiescent',spec,closed=False,reference=True)
    report=validate_rates(text,report,spec,source_expected=False)
    courant=[float(v) for v in re.findall(r'^Courant Number mean: \S+ max: (\S+)',text,re.M)]
    if not courant or not all(np.isfinite(courant)):raise ValueError('Missing/nonfinite Courant evidence')
    v=internal(case/'0/V',CELLS); a0=internal(case/'0/alpha.liquid',CELLS)
    measurements=[]
    for step in range(WRITE,STEPS+1,WRITE):
        time=f'{step*DT:.12g}'; folder=case/time
        a=internal(folder/'alpha.liquid',CELLS); b=internal(folder/'alpha.gas',CELLS)
        if np.min(a)<-1e-10 or np.max(a)>1+1e-10 or np.max(np.abs(a+b-1))>1e-10:
            raise ValueError('Invalid reconstructed phase fractions')
        velocity={}
        for phase,alpha in (('gas',b),('liquid',a)):
            speed=np.linalg.norm(internal(folder/('U.'+phase),CELLS,True),axis=1)
            velocity[phase]=dict(maxPresent_m_s=float(speed[alpha>.01].max(initial=0)),
                phaseWeightedRms_m_s=float(np.sqrt(np.dot(alpha*v,speed*speed)/np.dot(alpha,v))))
        measurements.append(dict(time_s=step*DT,liquidVolume_L=float(1000*np.dot(a,v)),
            alphaChangeL1_m3=float(np.dot(np.abs(a-a0),v)),velocity=velocity))
    accepted=parse(text)['SC1_TWO_STEP']
    for phase,density in (('gas',1),('liquid',1000)):
        final=internal(case/f'{END:.12g}'/('alpha.'+phase),CELLS)
        if abs(float(np.dot(final,v))*density-report['phases'][phase]['finalMass_kg'])>spec.mass_identity_tol[0 if phase=='gas' else 1]:
            raise ValueError('Reconstructed final inventory disagrees with accepted log')
        report['phases'][phase].update(alphaMinimum=min(r['alphaMin'] for r in accepted if r['phase']==phase),
            alphaMaximum=max(r['alphaMax'] for r in accepted if r['phase']==phase))
    quiet_errors=[]
    if max(courant)>.2:quiet_errors.append('Mixture Courant number exceeded 0.2')
    if max(r['alphaChangeL1_m3'] for r in measurements)>1e-5:
        quiet_errors.append('Quiescent alpha redistribution exceeded 0.01 L')
    if max(r['velocity']['liquid']['maxPresent_m_s'] for r in measurements)>1e-3:
        quiet_errors.append('Quiescent liquid speed exceeded 1 mm/s')
    report.update(suite=m.get('variant',TAG),initial=m['initial'],maxCourant=max(courant),snapshots=measurements,
        startup_pass=report['screen_pass'] and not quiet_errors,quiescent_errors=quiet_errors,
        log_sha256=base.sha(case/'log.solver'),
        meaning='10-ms actual-mesh startup screen, not equilibrium, stirred/aerated flow, oxygen or biology validation.')
    return report


def require_allocation(ranks):
    if not os.environ.get('SLURM_JOB_ID'):raise ValueError('Use an active Kestrel allocation')
    # srun also enforces the allocation; reject an obvious shortage before claiming a case.
    cpus=os.environ.get('SLURM_JOB_CPUS_PER_NODE','')
    if cpus and re.fullmatch(r'\d+(?:\(x\d+\))?(?:,\d+(?:\(x\d+\))?)*',cpus):
        available=sum(int(n)*int(k or 1) for n,k in re.findall(r'(\d+)(?:\(x(\d+)\))?',cpus))
        if available<ranks:raise ValueError('Requested ranks exceed allocated CPUs')


def run(root):
    root=root.resolve(strict=True);m=verify(root);require_allocation(m['ranks'])
    case=root/'quiescent';allowed={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}
    if {p.name for p in case.iterdir()}!=allowed:raise ValueError('Refusing an attempted/modified startup case')
    base.write_json(case/'attempt.json',dict(build=m['build'],assets=m['assets'],ranks=m['ranks']))
    env=base.environment(Path(m['build'])); app=str(Path(m['build'])/'bin'/APP)
    if m['ranks']>1:
        base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
        command=['srun','-n',str(m['ranks']),'--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated']
    else:command=[app,'-fileHandler','uncollated']
    base.execute(command,case,case/'log.solver',env)
    if m['ranks']>1:
        base.execute(['reconstructPar','-time',f'{WRITE*DT:.12g}:{END:.12g}',
            '-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
    if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Initial inputs changed during startup')
    report=analyze(case,m);base.write_json(case/'result.json',report)
    print(json.dumps(report,indent=2),flush=True)
    if m.get('variant')==DRAG_TAG:
        from sc1_reactor_startup_flow import compare_flow
        comparison=compare_flow(root,m,report)
        base.write_json(root/'drag_comparison.json',comparison)
        print(json.dumps(comparison,indent=2),flush=True)
    print('REACTOR_STARTUP='+('PASS' if report['startup_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='action',required=True)
    s=subs.add_parser('prepare');s.add_argument('--regression-root',type=Path,required=True)
    s.add_argument('--archive',type=Path,default=ARCHIVE);s.add_argument('--out',type=Path,required=True)
    s.add_argument('--ranks',type=int,default=8)
    s=subs.add_parser('prepare-drag');s.add_argument('--baseline-root',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('run');s.add_argument('--root',type=Path,required=True);s.add_argument('--run',action='store_true',required=True)
    s=subs.add_parser('analyze');s.add_argument('--root',type=Path,required=True)
    a=p.parse_args(argv)
    if a.action=='prepare':
        require_allocation(a.ranks);prepare(a.regression_root,a.archive,a.out,a.ranks)
    elif a.action=='prepare-drag':prepare_drag(a.baseline_root,a.out)
    elif a.action=='run':return 0 if run(a.root)['startup_pass'] else 2
    else:
        root=a.root.resolve(strict=True);m=verify(root)
        report=analyze(root/'quiescent',m)
        print(json.dumps(report,indent=2))
        if m['variant']==DRAG_TAG:
            from sc1_reactor_startup_flow import compare_flow
            print(json.dumps(compare_flow(root,m,report),indent=2))
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        print('SC1 reactor startup stopped: '+str(error)+'; preserve the startup root and logs',file=sys.stderr)
        raise SystemExit(1)
