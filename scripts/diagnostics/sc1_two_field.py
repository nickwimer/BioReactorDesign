#!/usr/bin/env python3
"""Prepare/build/run an isolated two-field numerical harness; never alter SC1 cases.

All project libraries are rebuilt privately. The build command launches no CFD.
Runs require a separate explicit --run and reject any previously attempted case.
continue-closed can validate/save the first result and launch only its fresh partner.
"""
import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from prepare_sc1_alpha_audit import ignore_generated, source_hashes, safe_make_path
from sc1_two_field_cases import CASES, case_files
from sc1_two_field_instrument import MARKER, phase_source, pressure_source, solver_source
from sc1_two_field_report import summarize

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO/'OFsolvers/birdmultiphaseEulerFoam'
APP = 'birdTwoFieldConservation'
COMPONENTS = (
    ('interfacialCompositionModels','libbirdeulerianInterfacialCompositionModels.so'),
    ('interfacialModels','libbirdeulerianInterfacialModels.so'),
    ('phaseSystems','libbirdphaseSystem.so'),
    ('multiphaseCompressibleMomentumTransportModels','libbirdmultiphaseMomentumTransportModels.so'),
    ('multiphaseReactions','libbirdmultiphaseReactions.so'),
    ('multiphaseEulerFoam/multiphaseSystems','libbirdmultiphaseSystems.so'),
)
COPY_DIRS = ('include','interfacialCompositionModels','interfacialModels','phaseSystems',
             'multiphaseCompressibleMomentumTransportModels','multiphaseReactions','multiphaseEulerFoam')


@dataclass(frozen=True)
class BuildRecipe:
    """Explicit private variant; default callers retain the original build contract."""
    tag: str
    app: str
    marker: str
    phase_transform: object
    solver_transform: object
    pressure_transform: object
    headers: tuple = ()


def build_recipe(recipe):
    return recipe or BuildRecipe('',APP,MARKER,phase_source,solver_source,pressure_source)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    with path.open('x') as f:
        json.dump(value,f,indent=2,allow_nan=False); f.write('\n')


def fresh(path, inputs):
    if path.exists() or path.is_symlink():
        raise ValueError('Destination already exists; preserve it: '+str(path))
    path=path.resolve()
    safe_make_path(path)
    for original in inputs:
        original=Path(original).resolve()
        if path==original or path in original.parents or original in path.parents:
            raise ValueError('Destination overlaps an input: '+str(original))
    return path


def selected_hashes(root):
    return {d:source_hashes(root/d) for d in COPY_DIRS}


def prepare_build(out, source=SOURCE, installed=None, *, recipe=None):
    recipe=build_recipe(recipe)
    source=source.resolve(strict=True)
    if installed is None:
        installed=Path(os.environ['FOAM_USER_LIBBIN']).resolve(strict=True)
    installed=installed.resolve(strict=True)
    out=fresh(out,(source,REPO,installed))
    before=selected_hashes(source)
    original_binaries={str(installed/name):sha(installed/name) for _,name in COMPONENTS if (installed/name).is_file()}
    # Validate source anchors and all Make destinations before creating anything.
    changes={
        'phaseSystems/phaseSystem/phaseSystemSolve.C':recipe.phase_transform((source/'phaseSystems/phaseSystem/phaseSystemSolve.C').read_text()),
        'multiphaseEulerFoam/multiphaseEulerFoam.C':recipe.solver_transform((source/'multiphaseEulerFoam/multiphaseEulerFoam.C').read_text()),
        'multiphaseEulerFoam/pEqnComps.H':recipe.pressure_transform((source/'multiphaseEulerFoam/pEqnComps.H').read_text()),
    }
    for directory, name in COMPONENTS:
        path=directory+'/Make/files'
        text=(source/path).read_text()
        target='LIB = '+safe_make_path(out/'lib'/name[:-3])
        text,count=re.subn(r'(?m)^LIB\s*=.*$',lambda _:target,text)
        if count!=1:
            raise ValueError('Unexpected Make library destination: '+path)
        changes[path]=text
    path='multiphaseEulerFoam/Make/files'
    text=(source/path).read_text()
    text,count=re.subn(r'(?m)^EXE\s*=.*$',lambda _:'EXE = '+safe_make_path(out/'bin'/recipe.app),text)
    if count!=1: raise ValueError('Unexpected Make executable destination')
    changes[path]=text
    out.mkdir(parents=True,exist_ok=False)
    (out/'src').mkdir(); (out/'lib').mkdir(); (out/'bin').mkdir()
    for directory in COPY_DIRS:
        shutil.copytree(source/directory,out/'src'/directory,symlinks=False,ignore=ignore_generated)
    if selected_hashes(out/'src')!=before or selected_hashes(source)!=before:
        raise ValueError('Source changed during copy; preserve partial tree')
    if any(p.is_symlink() for p in (out/'src').rglob('*')):
        raise ValueError('Unexpected symlink in private source copy')
    for path, text in changes.items(): (out/'src'/path).write_text(text)
    headers={
        'phaseSystems/phaseSystem/sc1PhaseFractionAudit.H':HERE/'include/sc1PhaseFractionAudit.H',
        'multiphaseEulerFoam/sc1TwoFieldStages.H':HERE/'include/sc1TwoFieldStages.H',
    }
    headers.update(dict(recipe.headers))
    for path, original in headers.items(): shutil.copyfile(original,out/'src'/path)
    manifest=dict(version=1,root=str(out),source=str(source),original_sources=before,
        private_sources=selected_hashes(out/'src'),installed_binaries=original_binaries,
        source_commit=subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip())
    if recipe.tag: manifest['variant']=recipe.tag
    write_json(out/'build_manifest.json',manifest)
    return out


def environment(root):
    env=os.environ.copy()
    env['FOAM_USER_LIBBIN']=str(root/'lib'); env['FOAM_USER_APPBIN']=str(root/'bin')
    env['LD_LIBRARY_PATH']=str(root/'lib')+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    env['PATH']=str(root/'bin')+':'+env.get('PATH','')
    return env


def verify_build(root, require_record=True, *, recipe=None):
    recipe=build_recipe(recipe)
    root=root.resolve(strict=True)
    m=json.loads((root/'build_manifest.json').read_text())
    if m['version']!=1 or m['root']!=str(root): raise ValueError('Wrong build manifest')
    if m.get('variant','')!=recipe.tag: raise ValueError('Wrong private build variant')
    if selected_hashes(Path(m['source']))!=m['original_sources']:
        raise ValueError('Original source tree changed since preparation')
    if selected_hashes(root/'src')!=m['private_sources']:
        raise ValueError('Private sources changed since preparation')
    for path, expected in m['installed_binaries'].items():
        if sha(path)!=expected: raise ValueError('Installed binary changed: '+path)
    env=environment(root)
    binaries=[root/'bin'/recipe.app]+[root/'lib'/name for _,name in COMPONENTS]
    assets={}
    for binary in binaries:
        if binary.is_symlink() or not binary.is_file(): raise ValueError('Missing/redirected private binary: '+str(binary))
        assets[str(binary)]=sha(binary)
        result=subprocess.run(['ldd',str(binary)],env=env,text=True,capture_output=True,check=True,timeout=30)
        if 'not found' in result.stdout+result.stderr:
            raise ValueError('Unresolved dependencies: '+str(binary)+'\n'+result.stdout)
        mappings=dict(re.findall(r'^\s*(\S+)\s+=>\s+(\S+)\s+\(',result.stdout,re.M))
        for name,path in mappings.items():
            if name.startswith(('libbird','libbirdeulerian')) and Path(path).resolve()!=root/'lib'/name:
                raise ValueError('Project library escaped private build: '+name+' => '+path)
        if binary.name in (recipe.app,'libbirdmultiphaseSystems.so') and mappings.get('libbirdphaseSystem.so')!=str(root/'lib/libbirdphaseSystem.so'):
            raise ValueError('Private phase-system dependency missing: '+str(binary))
    for binary, kind in ((root/'bin'/recipe.app,'U'),(root/'lib/libbirdphaseSystem.so','T')):
        result=subprocess.run(['nm','-D',str(binary)],text=True,capture_output=True,check=True,timeout=30)
        matches=re.findall(r'^\s*(?:[0-9a-fA-F]+\s+)?([A-Za-z])\s+'+recipe.marker+r'\s*$',result.stdout,re.M)
        if matches!=[kind]: raise ValueError('Mandatory observer symbol missing/wrong: '+str(binary))
    if not os.access(root/'bin'/recipe.app,os.X_OK): raise ValueError('Solver is not executable')
    if require_record and json.loads((root/'build_verified.json').read_text())['assets']!=assets:
        raise ValueError('Verified binaries changed')
    return dict(version=1,assets=assets,installed_unchanged=True,private_dependencies=True)


def execute(command, cwd, log, env):
    print('Executing: '+' '.join(map(str,command)),flush=True)
    with log.open('x') as output:
        result=subprocess.run(list(map(str,command)),cwd=cwd,env=env,stdout=output,stderr=subprocess.STDOUT)
    if result.returncode:
        lines=log.read_text(errors='replace').splitlines()
        errors=[line for line in lines if re.search(r'\b(?:fatal )?error:',line)]
        if errors:
            print('First compiler errors:\n'+'\n'.join(errors[:8]),flush=True)
        print('\n'.join(lines[-60:]),flush=True)
        raise ValueError(f'Command failed ({result.returncode}); preserve {log}')


def build(out, *, recipe=None):
    if os.environ.get('WM_PROJECT_VERSION')!='9': raise ValueError('Source the OpenFOAM 9 environment first')
    if not shutil.which('wmake'): raise ValueError('wmake unavailable')
    root=prepare_build(out,recipe=recipe)
    print('TWOFIELD_BUILD='+str(root),flush=True)
    env=environment(root)
    for directory, _ in COMPONENTS:
        execute(['wmake','libso'],root/'src'/directory,root/('log.build.'+directory.replace('/','_')),env)
    execute(['wmake'],root/'src/multiphaseEulerFoam',root/'log.build.solver',env)
    checked=verify_build(root,require_record=False,recipe=recipe)
    write_json(root/'build_verified.json',checked)
    print('SC1 two-field private build and runtime verification: PASS (no CFD launched)',flush=True)
    for path,digest in checked['assets'].items(): print(f'{path} sha256={digest}')


def prepare_cases(build_root, out):
    build_root=build_root.resolve(strict=True)
    checked=verify_build(build_root)
    out=fresh(out,(build_root,REPO))
    recipes={name:case_files(name) for name in CASES}
    # There must be precisely one dictionary-line difference within each pair.
    for pair in ('closed','source'):
        a=recipes[pair+'_no_reference']; b=recipes[pair+'_gas_reference']
        if set(a)!=set(b) or [k for k in a if a[k]!=b[k]]!=['constant/phaseProperties']:
            raise ValueError('Uncontrolled paired-case differences')
        if b['constant/phaseProperties'].replace('referencePhase gas;\n','')!=a['constant/phaseProperties']:
            raise ValueError('Unexpected paired phaseProperties difference')
    out.mkdir(parents=True,exist_ok=False)
    manifests={}
    for name, files in recipes.items():
        case=out/name; case.mkdir()
        for rel,text in files.items():
            path=case/rel; path.parent.mkdir(exist_ok=True)
            with path.open('x') as f: f.write(text)
        manifests[name]={rel:sha(case/rel) for rel in files}
    write_json(out/'cases_manifest.json',dict(version=1,root=str(out),build=str(build_root),assets=checked['assets'],cases=manifests))
    print('TWOFIELD_CASES='+str(out)+' (four fresh cases, no CFD launched)')


def case_input_hashes(case):
    result={}
    for folder in ('0','constant','system'):
        root=case/folder
        if root.is_symlink(): raise ValueError('Redirected case input directory')
        for p in sorted(root.rglob('*')):
            if p.is_symlink(): raise ValueError('Redirected case input')
            if p.is_file(): result[p.relative_to(case).as_posix()]=sha(p)
    return result


def verify_cases(root):
    root=root.resolve(strict=True)
    m=json.loads((root/'cases_manifest.json').read_text())
    if m['version']!=1 or m['root']!=str(root) or set(m['cases'])!=set(CASES):
        raise ValueError('Invalid cases manifest')
    if verify_build(Path(m['build']))['assets']!=m['assets']:
        raise ValueError('Cases belong to a different build')
    return m


def require_fresh_case(root, name, manifest):
    case=root/name
    if case.is_symlink() or any(p.name not in ('0','constant','system') for p in case.iterdir()):
        raise ValueError('Refusing attempted/modified case: '+str(case))
    if case_input_hashes(case)!=manifest['cases'][name]:
        raise ValueError('Case inputs changed: '+str(case))


def completed_report(root, name, manifest, *, summarize_log=summarize, require_pass=True):
    """Validate a saved run and return its report without writing or running CFD."""
    case=root/name
    if case.is_symlink() or not case.is_dir():
        raise ValueError('Missing/redirected completed case: '+str(case))
    for rel in ('attempt.json','initialized.json','log.blockMesh','log.checkMesh',
                'log.setFields','log.solver'):
        path=case/rel
        if path.is_symlink() or not path.is_file():
            raise ValueError('Missing/redirected completed-run evidence: '+str(path))
    expected_attempt=dict(case=name,build=manifest['build'],assets=manifest['assets'])
    if json.loads((case/'attempt.json').read_text())!=expected_attempt:
        raise ValueError('Completed attempt does not match the verified build')
    inputs=case_input_hashes(case)
    if inputs!=json.loads((case/'initialized.json').read_text()):
        raise ValueError('Completed case inputs differ from initialized snapshot')
    # blockMesh adds mesh files and setFields changes alpha; all other prepared
    # inputs must still match their original case-manifest hashes.
    for rel,digest in manifest['cases'][name].items():
        if rel not in ('0/alpha.gas','0/alpha.liquid') and inputs.get(rel)!=digest:
            raise ValueError('Completed case differs from prepared inputs: '+rel)
    if 'Mesh OK.' not in (case/'log.checkMesh').read_text():
        raise ValueError('Completed case lacks successful mesh check')
    raw=(case/'log.solver').read_bytes()
    report=summarize_log(raw.decode(),name)
    report['log_sha256']=hashlib.sha256(raw).hexdigest()
    if require_pass and not report['screen_pass']:
        raise ValueError('Completed case failed its numerical screen: '+str(report['errors']))
    if sha(case/'log.solver')!=report['log_sha256'] or case_input_hashes(case)!=inputs:
        raise ValueError('Completed evidence changed during validation')
    destination=case/'result.json'
    if destination.is_symlink() or (destination.exists() and
            (not destination.is_file() or json.loads(destination.read_text())!=report)):
        raise ValueError('Existing completed result differs; preserve it: '+str(destination))
    return report


def run_fresh_cases(root, manifest, names, *, summarize_log=summarize, app=APP):
    m=manifest
    env=environment(Path(m['build']))
    for name in names:
        # Recheck immediately before claiming the case with attempt.json.
        require_fresh_case(root,name,m)
        case=root/name
        write_json(case/'attempt.json',dict(case=name,build=m['build'],assets=m['assets']))
        execute(['blockMesh'],case,case/'log.blockMesh',env)
        execute(['checkMesh'],case,case/'log.checkMesh',env)
        if 'Mesh OK.' not in (case/'log.checkMesh').read_text(): raise ValueError('Mesh did not pass checkMesh')
        execute(['setFields'],case,case/'log.setFields',env)
        initialized=case_input_hashes(case)
        write_json(case/'initialized.json',initialized)
        execute([Path(m['build'])/'bin'/app],case,case/'log.solver',env)
        after=case_input_hashes(case)
        if after!=initialized: raise ValueError('Initial fields or dictionaries changed during solve')
        try:
            report=summarize_log((case/'log.solver').read_text(),name)
        except ValueError as error:
            raise ValueError(f'{error}; preserve {case / "log.solver"}; '
                             'reanalyze the saved log without rerunning CFD') from error
        report['log_sha256']=sha(case/'log.solver')
        write_json(case/'result.json',report)
        print(json.dumps(report,indent=2),flush=True)


def run_pair(root, pair):
    root=root.resolve(strict=True); m=verify_cases(root)
    names=[pair+'_no_reference',pair+'_gas_reference']
    if pair=='source':
        for name in ('closed_no_reference','closed_gas_reference'):
            result=json.loads((root/name/'result.json').read_text())
            if result['log_sha256']!=sha(root/name/'log.solver') or not summarize((root/name/'log.solver').read_text(),name)['screen_pass']:
                raise ValueError('Source pair requires both completed closed screens to pass')
    # Check BOTH destinations before attempting the first case.
    for name in names: require_fresh_case(root,name,m)
    run_fresh_cases(root,m,names)
    print('PAIR='+pair+' complete. Review both budgets before accepting a numerical formulation.')


def continue_closed(root):
    """Recover a passing first report and run only the unattempted closed partner."""
    root=root.resolve(strict=True); m=verify_cases(root)
    first, partner='closed_no_reference','closed_gas_reference'
    require_fresh_case(root,partner,m)
    report=completed_report(root,first,m)
    destination=root/first/'result.json'
    if not destination.exists(): write_json(destination,report)
    print('REUSED='+str(root/first)+' (validated saved run; no CFD rerun)',flush=True)
    print(json.dumps(report,indent=2),flush=True)
    run_fresh_cases(root,m,[partner])
    print('PAIR=closed complete. Review both budgets before the source pair.')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='action',required=True)
    b=subs.add_parser('build'); b.add_argument('--out',type=Path,required=True)
    v=subs.add_parser('verify'); v.add_argument('--root',type=Path,required=True)
    c=subs.add_parser('prepare-cases'); c.add_argument('--build-root',type=Path,required=True); c.add_argument('--out',type=Path,required=True)
    r=subs.add_parser('run'); r.add_argument('--root',type=Path,required=True); r.add_argument('--pair',choices=('closed','source'),required=True); r.add_argument('--run',action='store_true',required=True)
    r=subs.add_parser('continue-closed'); r.add_argument('--root',type=Path,required=True); r.add_argument('--run',action='store_true',required=True)
    a=subs.add_parser('analyze'); a.add_argument('--log',type=Path,required=True); a.add_argument('--case',choices=CASES,required=True)
    args=parser.parse_args(argv)
    if args.action in ('build','run','continue-closed') and not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Use an active Kestrel allocation for build/run')
    if args.action=='build': build(args.out)
    elif args.action=='verify': print(json.dumps(verify_build(args.root),indent=2))
    elif args.action=='prepare-cases': prepare_cases(args.build_root,args.out)
    elif args.action=='run': run_pair(args.root,args.pair)
    elif args.action=='continue-closed': continue_closed(args.root)
    else: print(json.dumps(summarize(args.log.read_text(),args.case),indent=2))


if __name__=='__main__':
    try: main()
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        print('SC1 two-field check stopped: '+str(error),file=sys.stderr)
        raise SystemExit(1)
