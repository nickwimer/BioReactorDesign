#!/usr/bin/env python3
"""Build and run a restricted, private pressure-coupled 10-ms reactor startup."""
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
import sc1_reactor_startup as startup
import sc1_reactor_pressure_coupling_instrument as instrument
import sc1_reactor_pressure_probe as observation
from instrument_sc1_alpha_subcycles import _once
from sc1_two_field_rates_instrument import RECIPE, MARKER

TAG='sc1_meml_pressure_coupling_v5'
APP='birdPressureCoupled'
FAILED_LOG='d6f4e4b1345fdbe0029598c699e9be1f322479ce77152665ff3988c05052a9ea'
HEADERS=[base.HERE/'include'/name for name in ('sc1PressureProbe.H','sc1FaceDragCoupling.H','sc1FaceDragAlgebra.H')]
OBSERVED_LOG='d2c5d3b55ad096ee1fc8151a68353828fbfb05db9bc1e5c65af72fa0cc5ef656'


def recipe_hashes():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(instrument.__file__),Path(instrument.observation.__file__),*HEADERS)}


def build_sources(base_root,out):
    source=base_root/'src/multiphaseEulerFoam'
    changes={
        'multiphaseEulerFoam.C':instrument.solver((source/'multiphaseEulerFoam.C').read_text()),
        'pUf/pEqn.H':instrument.pressure((source/'pUf/pEqn.H').read_text()),
        'sc1PressureProbe.H':instrument.observer_header(HEADERS[0].read_text()),
        **{p.name:p.read_text() for p in HEADERS[1:]},
    }
    make=(source/'Make/files').read_text()
    make,count=re.subn(r'(?m)^EXE\s*=.*$',lambda _:'EXE = '+base.safe_make_path(out/'bin'/APP),make)
    if count!=1:raise ValueError('Unexpected candidate Make destination')
    changes['Make/files']=make
    options=(source/'Make/options').read_text()
    for directory in ('include','phaseSystems/lnInclude','interfacialModels/lnInclude','interfacialCompositionModels/lnInclude'):
        options=_once(options,'-I../'+directory,'-I'+base.safe_make_path(base_root/'src'/directory))
    changes['Make/options']=options
    return source,changes


def executable_identity(out,m):
    binary=out/'bin'/APP
    if binary.is_symlink() or not binary.is_file() or not os.access(binary,os.X_OK):
        raise ValueError('Missing/redirected/non-executable candidate')
    env=base.environment(Path(m['base_build']))
    result=subprocess.run(['ldd',str(binary)],env=env,text=True,capture_output=True,check=True,timeout=30)
    if 'not found' in result.stdout+result.stderr:raise ValueError('Candidate has unresolved libraries')
    mappings=dict(re.findall(r'^\s*(\S+)\s+=>\s+(\S+)\s+\(',result.stdout,re.M))
    expected={name:str(Path(m['base_build'])/'lib'/name) for _,name in base.COMPONENTS}
    if {k:v for k,v in mappings.items() if k.startswith(('libbird','libbirdeulerian'))}!=expected:
        raise ValueError('Candidate project libraries differ from verified rate build')
    symbols=subprocess.run(['nm','-D',str(binary)],text=True,capture_output=True,check=True,timeout=30).stdout
    if re.findall(r'^\s*([A-Za-z])\s+'+MARKER+r'\s*$',symbols,re.M)!=['U']:
        raise ValueError('Candidate lost mandatory rate-library dependency')
    return dict(binary=str(binary),sha256=base.sha(binary),dependencies=expected)


def build(base_root,out):
    startup.require_allocation(1)
    if os.environ.get('WM_PROJECT_VERSION')!='9' or not shutil.which('wmake'):
        raise ValueError('Source OpenFOAM 9 and make wmake available')
    base_root=base_root.resolve(strict=True)
    checked=base.verify_build(base_root,recipe=RECIPE)
    out=base.fresh(out,(base_root,base.REPO))
    source,changes=build_sources(base_root,out)
    before=base.source_hashes(source)
    out.mkdir(parents=True,exist_ok=False);(out/'bin').mkdir()
    shutil.copytree(source,out/'src',symlinks=False,ignore=base.ignore_generated)
    if base.source_hashes(out/'src')!=before:raise ValueError('Private application copy differs')
    for rel,text in changes.items():(out/'src'/rel).write_text(text)
    m=dict(variant=TAG,root=str(out),base_build=str(base_root),base_assets=checked['assets'],
        base_manifest_sha256=base.sha(base_root/'build_manifest.json'),recipe=recipe_hashes(),
        source_inputs=before,private_sources=base.source_hashes(out/'src'),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'coupling_build_manifest.json',m)
    print('PRESSURE_COUPLED_BUILD='+str(out),flush=True)
    base.execute(['wmake'],out/'src',out/'log.build.solver',base.environment(base_root))
    identity=verify_build(out,require_record=False)
    base.write_json(out/'coupling_build_verified.json',identity)
    print('Pressure-coupled application build/dependency verification: PASS (no CFD launched)',flush=True)
    print(json.dumps(identity,indent=2),flush=True)
    return out


def verify_build(out,require_record=True):
    out=out.resolve(strict=True);m=json.loads((out/'coupling_build_manifest.json').read_text())
    if m['variant']!=TAG or m['root']!=str(out) or m['recipe']!=recipe_hashes():
        raise ValueError('Wrong/changed pressure-coupling build recipe')
    base_root=Path(m['base_build'])
    if base.verify_build(base_root,recipe=RECIPE)['assets']!=m['base_assets'] or base.sha(base_root/'build_manifest.json')!=m['base_manifest_sha256']:
        raise ValueError('Verified rate build changed')
    source,changes=build_sources(base_root,out)
    expected=base.source_hashes(source)
    if expected!=m['source_inputs']:raise ValueError('Rate application sources changed')
    expected.update({p:hashlib.sha256(t.encode()).hexdigest() for p,t in changes.items()})
    if base.source_hashes(out/'src')!=expected or expected!=m['private_sources']:
        raise ValueError('Candidate sources differ from exact coupling recipe')
    identity=executable_identity(out,m)
    if require_record and json.loads((out/'coupling_build_verified.json').read_text())!=identity:
        raise ValueError('Verified candidate binary changed')
    return identity


def failed_control(root):
    return startup.failed_case(root,startup.COUPLED_TAG,FAILED_LOG,'received failed v4 coupled control')


def prepare(failed_root,candidate_build,out,probe_root):
    failed_root=failed_root.resolve(strict=True);candidate_build=candidate_build.resolve(strict=True)
    m,control=failed_control(failed_root);startup.require_allocation(m['ranks'])
    evidence=observed_probe(probe_root,m)
    identity=verify_build(candidate_build)
    pm=json.loads((candidate_build/'coupling_build_manifest.json').read_text())
    if pm['base_build']!=m['build'] or pm['base_assets']!=m['assets']:
        raise ValueError('Candidate does not use the failed control rate build')
    protected=[failed_root,candidate_build,probe_root.resolve(strict=True),base.REPO]+[Path(m[k]) for k in ('build','archive','regression_root')]
    protected += [Path(m[k]['root']) for k in ('baseline','failed_control','face_control')]
    out=base.fresh(out,protected);out.mkdir(parents=True,exist_ok=False)
    case=out/'quiescent';source=failed_root/'quiescent'
    for folder in ('0','constant','system'):shutil.copytree(source/folder,case/folder)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):shutil.copy2(source/name,case/name)
    if base.case_input_hashes(case)!=m['inputs'] or failed_control(failed_root)!=(m,control):
        raise ValueError('Inputs/control changed during candidate preparation')
    assets=dict(m['assets']);assets[identity['binary']]=identity['sha256']
    manifest=dict(m,variant=TAG,root=str(out),build=str(candidate_build),assets=assets,
        coupled_control=control,candidate_identity=identity,pressure_probe=evidence,
        coupling_build_manifest_sha256=base.sha(candidate_build/'coupling_build_manifest.json'),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',manifest)
    print('PRESSURE_COUPLED_STARTUP='+str(out)+' (fresh identical t=0 inputs; coupled pressure/drag; maximum 10 ms; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads((root/'startup_manifest.json').read_text())
    if m['variant']!=TAG or m['root']!=str(root):raise ValueError('Wrong pressure-coupled case')
    parent,control=failed_control(Path(m['coupled_control']['root']))
    identity=verify_build(Path(m['build']))
    if observed_probe(Path(m['pressure_probe']['root']),parent)!=m['pressure_probe']:
        raise ValueError('Observed pressure failure changed')
    if control!=m['coupled_control'] or identity!=m['candidate_identity'] or base.sha(Path(m['build'])/'coupling_build_manifest.json')!=m['coupling_build_manifest_sha256']:
        raise ValueError('Candidate control/build changed')
    pm=json.loads((Path(m['build'])/'coupling_build_manifest.json').read_text())
    assets=dict(parent['assets']);assets[identity['binary']]=identity['sha256']
    if pm['base_build']!=parent['build'] or pm['base_assets']!=parent['assets'] or m['assets']!=assets:
        raise ValueError('Candidate base build differs from control')
    for key in ('inputs','initial','ranks','archive','mesh_hashes','regression_root','regression_sha256','baseline','failed_control','face_control'):
        if m[key]!=parent[key]:raise ValueError('Candidate changed control metadata: '+key)
    if base.case_input_hashes(root/'quiescent')!=m['inputs']:raise ValueError('Candidate case inputs changed')
    return m,parent


def observed_probe(root,parent):
    root=root.resolve(strict=True)
    m,control_parent=observation.verify(root)
    if control_parent!=parent or base.sha(root/'quiescent/log.solver')!=OBSERVED_LOG:
        raise ValueError('Requires the received, unchanged pressure probe and same v4 control')
    return dict(root=str(root),log_sha256=OBSERVED_LOG,
        manifest_sha256=base.sha(root/'startup_manifest.json'))


def analyze(case,m):
    # Reuse unchanged conservation/quiescence limits and all 40 face snapshots.
    report=startup.analyze(case,dict(m,variant=startup.FACE_TAG))
    text=(case/'log.solver').read_text();lines=text.splitlines()
    if lines.count('SC1_PRESSURE_PROBE version=1 maxCompletedSteps=200')!=1 or lines.count(
            'SC1_PRESSURE_DONE completedSteps=200 meaning=requiresStartupAnalysis')!=1:
        raise ValueError('Missing pressure observer start/completion markers')
    if lines.count('SC1_FACE_COUPLING version=1 order=beforePressure')!=2*startup.STEPS:
        raise ValueError('Missing complete pre-pressure drag coupling coverage')
    if 'Inverting drag system:' in text or 'term=explicitDragPhi' in text:
        raise ValueError('Unexpected lagged/post-pressure drag path')
    fields=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('SC1_PRESSURE_FIELD ')]
    pressures=[float(r['min']) for r in fields if r.get('term')=='pUnclipped']
    if len(pressures)!=4*startup.STEPS or any(r.get('nonfinite')!='0' for r in fields):
        raise ValueError('Incomplete/nonfinite pressure evidence')
    if not all(np.isfinite(pressures)) or min(pressures)<1000:
        raise ValueError('Invalid unclipped pressure')
    fluxes=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('SC1_COUPLED_FLUX ')]
    if len(fluxes)!=4*startup.STEPS:raise ValueError('Incomplete coupled mixture flux checks')
    for row in fluxes:
        error=float(row['maxAbsError_m3_s']);scale=float(row['maxAbsMixture_m3_s'])
        if not np.isfinite(error+scale) or min(error,scale)<0 or error>1e-12+1e-10*scale:
            raise ValueError('Invalid coupled mixture flux identity')
    speeds=[dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('SC1_QUIESCENT_STEP ')]
    if len(speeds)!=startup.STEPS:raise ValueError('Incomplete accepted-step speed checks')
    for i,row in enumerate(speeds,1):
        t=float(row['t_s']);v=float(row['maxLiquidPresent_m_s'])
        if not np.isfinite(t+v) or abs(t-i*startup.DT)>1e-12 or not 0<=v<=1e-3:
            raise ValueError('Invalid accepted-step quiescent speed')
    report.update(suite=TAG,pressureCoupling=dict(calls=2*startup.STEPS,
        minimumUnclippedPressure_Pa=min(pressures),acceptedStepSpeedChecks=len(speeds),
        maxAcceptedLiquidPresent_m_s=max(float(r['maxLiquidPresent_m_s']) for r in speeds),
        mixtureFluxChecks=len(fluxes),maxMixtureFluxError_m3_s=max(float(r['maxAbsError_m3_s']) for r in fluxes),
        method='Exact symmetric two-phase face-drag inverse in both pressure coefficient and RHS'),
        pressure_probe=m['pressure_probe'])
    return report


def run(root):
    root=root.resolve(strict=True);m,parent=verify(root);startup.require_allocation(m['ranks'])
    case=root/'quiescent'
    allowed={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}
    if {p.name for p in case.iterdir()}!=allowed:raise ValueError('Refusing attempted/modified pressure-coupled startup')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=base.environment(Path(parent['build']));app=m['candidate_identity']['binary'];report=None
    try:
        if m['ranks']>1:
            base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
            command=['srun','--kill-on-bad-exit=1','--time=00:30:00','-n',str(m['ranks']),'--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated']
        else:command=[app,'-fileHandler','uncollated']
        from sc1_reactor_watch import execute
        execute(command,case,case/'log.solver',env)
        if m['ranks']>1:
            base.execute(['reconstructPar','-time',f'{startup.FACE_WRITE*startup.DT:.12g}:{startup.END:.12g}',
                '-fileHandler','uncollated'],case,case/'log.reconstructPar',env)
        if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Startup initial inputs changed')
        report=analyze(case,m);base.write_json(case/'result.json',report)
        print(json.dumps(report,indent=2),flush=True)
    finally:
        # Collect diagnostic evidence for success, numerical stop, or failed analysis.
        if (case/'log.solver').is_file():
            from sc1_reactor_failure import collect
            try:collect(root,root.with_name(root.name+'_evidence'))
            except Exception as error:print('Automatic evidence bundle unavailable: '+str(error),flush=True)
    print('PRESSURE_COUPLED_STARTUP='+('PASS' if report['startup_pass'] else 'FAIL'),flush=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='action',required=True)
    s=subs.add_parser('build');s.add_argument('--base-build',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('prepare');s.add_argument('--failed-root',type=Path,required=True);s.add_argument('--candidate-build',type=Path,required=True);s.add_argument('--out',type=Path,required=True);s.add_argument('--probe-root',type=Path,required=True)
    s=subs.add_parser('run');s.add_argument('--root',type=Path,required=True);s.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='build':build(a.base_build,a.out)
    elif a.action=='prepare':prepare(a.failed_root,a.candidate_build,a.out,a.probe_root)
    else:return 0 if run(a.root)['startup_pass'] else 2
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 pressure-coupled startup stopped: '+str(error),file=sys.stderr);raise SystemExit(1)
