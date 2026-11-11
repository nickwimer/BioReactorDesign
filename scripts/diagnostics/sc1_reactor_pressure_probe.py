#!/usr/bin/env python3
"""Build an application-only pressure probe and run a fresh, at-most-two-step diagnostic."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import sc1_two_field as base
import sc1_reactor_startup as startup
import sc1_reactor_pressure_instrument as instrument
from instrument_sc1_alpha_subcycles import _once
from sc1_two_field_rates_instrument import RECIPE, MARKER

TAG='sc1_meml_pressure_probe_v1'
APP='birdPressureProbe'
FAILED_LOG='d6f4e4b1345fdbe0029598c699e9be1f322479ce77152665ff3988c05052a9ea'
HEADER=base.HERE/'include/sc1PressureProbe.H'


def recipe_hashes():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(instrument.__file__),HEADER)}


def build_sources(base_root,out):
    source=base_root/'src/multiphaseEulerFoam'
    changes={
        'multiphaseEulerFoam.C':instrument.solver((source/'multiphaseEulerFoam.C').read_text()),
        'pUf/pEqn.H':instrument.pressure((source/'pUf/pEqn.H').read_text()),
        'sc1PressureProbe.H':HEADER.read_text(),
    }
    make=(source/'Make/files').read_text()
    make,count=re.subn(r'(?m)^EXE\s*=.*$',lambda _:'EXE = '+base.safe_make_path(out/'bin'/APP),make)
    if count!=1:raise ValueError('Unexpected probe Make destination')
    changes['Make/files']=make
    options=(source/'Make/options').read_text()
    for directory in ('include','phaseSystems/lnInclude','interfacialModels/lnInclude','interfacialCompositionModels/lnInclude'):
        options=_once(options,'-I../'+directory,'-I'+base.safe_make_path(base_root/'src'/directory))
    changes['Make/options']=options
    return source,changes


def executable_identity(out,m):
    binary=out/'bin'/APP
    if binary.is_symlink() or not binary.is_file() or not os.access(binary,os.X_OK):
        raise ValueError('Missing/redirected/non-executable probe')
    env=base.environment(Path(m['base_build']))
    result=subprocess.run(['ldd',str(binary)],env=env,text=True,capture_output=True,check=True,timeout=30)
    if 'not found' in result.stdout+result.stderr:raise ValueError('Probe has unresolved libraries')
    mappings=dict(re.findall(r'^\s*(\S+)\s+=>\s+(\S+)\s+\(',result.stdout,re.M))
    expected={name:str(Path(m['base_build'])/'lib'/name) for _,name in base.COMPONENTS}
    if {k:v for k,v in mappings.items() if k.startswith(('libbird','libbirdeulerian'))}!=expected:
        raise ValueError('Probe project libraries differ from verified rate build')
    symbols=subprocess.run(['nm','-D',str(binary)],text=True,capture_output=True,check=True,timeout=30).stdout
    if re.findall(r'^\s*([A-Za-z])\s+'+MARKER+r'\s*$',symbols,re.M)!=['U']:
        raise ValueError('Probe lost mandatory rate-library dependency')
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
    base.write_json(out/'probe_build_manifest.json',m)
    print('PRESSURE_PROBE_BUILD='+str(out),flush=True)
    base.execute(['wmake'],out/'src',out/'log.build.solver',base.environment(base_root))
    identity=verify_build(out,require_record=False)
    base.write_json(out/'probe_build_verified.json',identity)
    print('Pressure probe application build/dependency verification: PASS (no CFD launched)',flush=True)
    print(json.dumps(identity,indent=2),flush=True)
    return out


def verify_build(out,require_record=True):
    out=out.resolve(strict=True);m=json.loads((out/'probe_build_manifest.json').read_text())
    if m['variant']!=TAG or m['root']!=str(out) or m['recipe']!=recipe_hashes():
        raise ValueError('Wrong/changed pressure probe build recipe')
    base_root=Path(m['base_build'])
    if base.verify_build(base_root,recipe=RECIPE)['assets']!=m['base_assets'] or base.sha(base_root/'build_manifest.json')!=m['base_manifest_sha256']:
        raise ValueError('Verified rate build changed')
    source,changes=build_sources(base_root,out)
    expected=base.source_hashes(source)
    if expected!=m['source_inputs']:raise ValueError('Rate application sources changed')
    expected.update({p:hashlib.sha256(t.encode()).hexdigest() for p,t in changes.items()})
    if base.source_hashes(out/'src')!=expected or expected!=m['private_sources']:
        raise ValueError('Probe sources differ from exact observation recipe')
    identity=executable_identity(out,m)
    if require_record and json.loads((out/'probe_build_verified.json').read_text())!=identity:
        raise ValueError('Verified probe binary changed')
    return identity


def failed_control(root):
    return startup.failed_case(root,startup.COUPLED_TAG,FAILED_LOG,'received failed v4 coupled control')


def prepare(failed_root,probe_build,out):
    failed_root=failed_root.resolve(strict=True);probe_build=probe_build.resolve(strict=True)
    m,control=failed_control(failed_root);startup.require_allocation(m['ranks'])
    identity=verify_build(probe_build)
    pm=json.loads((probe_build/'probe_build_manifest.json').read_text())
    if pm['base_build']!=m['build'] or pm['base_assets']!=m['assets']:
        raise ValueError('Probe does not use the failed control rate build')
    protected=[failed_root,probe_build,base.REPO]+[Path(m[k]) for k in ('build','archive','regression_root')]
    protected += [Path(m[k]['root']) for k in ('baseline','failed_control','face_control')]
    out=base.fresh(out,protected);out.mkdir(parents=True,exist_ok=False)
    case=out/'quiescent';source=failed_root/'quiescent'
    for folder in ('0','constant','system'):shutil.copytree(source/folder,case/folder)
    for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):shutil.copy2(source/name,case/name)
    if base.case_input_hashes(case)!=m['inputs'] or failed_control(failed_root)!=(m,control):
        raise ValueError('Inputs/control changed during probe preparation')
    assets=dict(m['assets']);assets[identity['binary']]=identity['sha256']
    manifest=dict(m,variant=TAG,root=str(out),build=str(probe_build),assets=assets,
        probe_control=control,probe_identity=identity,
        probe_build_manifest_sha256=base.sha(probe_build/'probe_build_manifest.json'),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'startup_manifest.json',manifest)
    print('PRESSURE_PROBE='+str(out)+' (fresh identical t=0 inputs; at most two steps; no CFD launched)',flush=True)
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads((root/'startup_manifest.json').read_text())
    if m['variant']!=TAG or m['root']!=str(root):raise ValueError('Wrong pressure probe case')
    parent,control=failed_control(Path(m['probe_control']['root']))
    identity=verify_build(Path(m['build']))
    if control!=m['probe_control'] or identity!=m['probe_identity'] or base.sha(Path(m['build'])/'probe_build_manifest.json')!=m['probe_build_manifest_sha256']:
        raise ValueError('Probe control/build changed')
    pm=json.loads((Path(m['build'])/'probe_build_manifest.json').read_text())
    assets=dict(parent['assets']);assets[identity['binary']]=identity['sha256']
    if pm['base_build']!=parent['build'] or pm['base_assets']!=parent['assets'] or m['assets']!=assets:
        raise ValueError('Probe base build differs from control')
    for key in ('inputs','initial','ranks','archive','mesh_hashes','regression_root','regression_sha256','baseline','failed_control','face_control'):
        if m[key]!=parent[key]:raise ValueError('Probe changed control metadata: '+key)
    if base.case_input_hashes(root/'quiescent')!=m['inputs']:raise ValueError('Probe case inputs changed')
    return m,parent


def summary(text):
    stages=[];fields=[];stops=[];done=[]
    for line in text.splitlines():
        data=dict(re.findall(r'(\w+)=([^\s]+)',line))
        if line.startswith('SC1_PRESSURE_STAGE '):stages.append(data)
        if line.startswith('SC1_PRESSURE_FIELD ') and data.get('term') in ('pUnclipped','rAUfMixture','faceDiagonalExtra'):
            fields.append(data)
        if line.startswith('SC1_PRESSURE_STOP '):stops.append(data)
        if line.startswith('SC1_PRESSURE_DONE '):done.append(data)
    return dict(suite=TAG,probeMarkerCount=text.splitlines().count('SC1_PRESSURE_PROBE version=1 maxCompletedSteps=2'),
        stages=stages,selectedFields=fields,stops=stops,done=done,
        meaning='Diagnostic evidence only; no startup or conservation acceptance. Full field/location records are in the log.')


def run(root):
    root=root.resolve(strict=True);m,parent=verify(root);startup.require_allocation(m['ranks'])
    case=root/'quiescent'
    allowed={'0','constant','system','log.checkMesh','log.writeCellCentres','log.writeCellVolumes'}
    if {p.name for p in case.iterdir()}!=allowed:raise ValueError('Refusing attempted/modified probe')
    base.write_json(case/'attempt.json',{k:m[k] for k in ('build','assets','ranks')})
    env=base.environment(Path(parent['build']));app=m['probe_identity']['binary']
    error=None
    try:
        if m['ranks']>1:
            base.execute(['decomposePar','-fileHandler','uncollated'],case,case/'log.decomposePar',env)
            command=['srun','--kill-on-bad-exit=1','--time=00:30:00','-n',str(m['ranks']),'--cpus-per-task=1',app,'-parallel','-fileHandler','uncollated']
        else:command=[app,'-fileHandler','uncollated']
        from sc1_reactor_watch import execute
        execute(command,case,case/'log.solver',env)
    except (OSError,ValueError,subprocess.SubprocessError) as exc:error=exc
    if (case/'log.solver').is_file():
        report=summary((case/'log.solver').read_text());report['log_sha256']=base.sha(case/'log.solver')
        base.write_json(case/'result.json',report);print(json.dumps(report,indent=2),flush=True)
        from sc1_reactor_failure import collect
        collect(root,root.with_name(root.name+'_evidence'))
        if not error and (report['probeMarkerCount']!=1 or len(report['done'])!=1 or report['done'][0].get('completedSteps')!='2'):
            error=ValueError('Missing complete two-step probe coverage')
    elif error is None:error=ValueError('Missing probe solver log')
    if base.case_input_hashes(case)!=m['inputs']:raise ValueError('Probe inputs changed during execution')
    if error:raise ValueError(str(error)+'; preserve all probe/control roots')
    print('Pressure probe completed two diagnostic steps; this is not a startup pass.',flush=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='action',required=True)
    s=subs.add_parser('build');s.add_argument('--base-build',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('prepare');s.add_argument('--failed-root',type=Path,required=True);s.add_argument('--probe-build',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('run');s.add_argument('--root',type=Path,required=True);s.add_argument('--run',action='store_true',required=True)
    a=p.parse_args(argv)
    if a.action=='build':build(a.base_build,a.out)
    elif a.action=='prepare':prepare(a.failed_root,a.probe_build,a.out)
    else:run(a.root)
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print('SC1 pressure probe stopped: '+str(error),file=sys.stderr);raise SystemExit(1)
