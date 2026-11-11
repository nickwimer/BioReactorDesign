"""Private Foundation interFoam application with observation-only changes."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sc1_two_field as base
import sc1_reactor_startup as startup
import sc1_reactor_vof_instrument as instrument

APP='sc1InterFoamReference'
FIXTURE=base.HERE/'fixtures/sc1_interfoam9_sources.json'
HEADER=base.HERE/'include/sc1VofAudit.H'


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(instrument.__file__),FIXTURE,HEADER)}


def native_sources(project):
    expected=json.loads(FIXTURE.read_text())['sources'];actual={}
    for rel,text in expected.items():
        p=project/rel
        if not p.is_file() or base.sha(p)!=hashlib.sha256(text.encode()).hexdigest():
            raise ValueError('Installed OF9 source differs from inspected upstream: '+str(p))
        actual[rel]=base.sha(p)
    return actual


def private_sources(project,out):
    native_sources(project);result={}
    for rel in json.loads(FIXTURE.read_text())['sources']:
        name=rel.split('interFoam/',1)[1] if '/interFoam/' in rel else Path(rel).name
        result[name]=(project/rel).read_text()
    result['interFoam.C']=instrument.solver(result['interFoam.C'])
    result['pEqn.H']=instrument.pressure(result['pEqn.H'])
    result['sc1VofAudit.H']=HEADER.read_text()
    result['Make/files']='interFoam.C\n\nEXE = '+base.safe_make_path(out/'bin'/APP)+'\n'
    return result


def identity(out):
    binary=out/'bin'/APP
    if binary.is_symlink() or not binary.is_file() or not os.access(binary,os.X_OK):
        raise ValueError('Missing/redirected VOF reference executable')
    r=subprocess.run(['ldd',str(binary)],capture_output=True,text=True,check=True,timeout=30)
    if 'not found' in r.stdout+r.stderr:raise ValueError('Unresolved VOF reference library')
    mapping=dict(re.findall(r'^\s*(\S+)\s+=>\s+(\S+)\s+\(',r.stdout,re.M))
    if any(name.startswith(('libbird','libbirdeulerian')) for name in mapping):
        raise ValueError('VOF reference unexpectedly links BiRD libraries')
    lib=Path(os.environ['FOAM_LIBBIN']).resolve(strict=True)
    for name in ('libOpenFOAM.so','libfiniteVolume.so','libimmiscibleIncompressibleTwoPhaseMixture.so','libtwoPhaseChangeModels.so'):
        if name not in mapping or Path(mapping[name]).resolve()!=(lib/name).resolve():
            raise ValueError('VOF reference not using active OF9 core library: '+name)
    assets={str(binary):base.sha(binary)}
    assets.update({str(Path(path).resolve(strict=True)):base.sha(path) for path in mapping.values()})
    return dict(binary=str(binary),sha256=assets[str(binary)],dependencies=mapping,assets=assets)


def build(out):
    startup.require_allocation(1)
    if os.environ.get('WM_PROJECT_VERSION')!='9' or not shutil.which('wmake'):
        raise ValueError('Source OpenFOAM 9 before building the reference')
    project=Path(os.environ['WM_PROJECT_DIR']).resolve(strict=True)
    out=base.fresh(out,(project,base.REPO));sources=private_sources(project,out)
    m=dict(root=str(out),project=str(project),recipe=recipe(),native_sources=native_sources(project),
        private_sources={p:hashlib.sha256(t.encode()).hexdigest() for p,t in sources.items()})
    out.mkdir(parents=True,exist_ok=False);(out/'bin').mkdir();(out/'src').mkdir()
    for rel,text in sources.items():
        p=out/'src'/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    base.write_json(out/'vof_build_manifest.json',m)
    print('VOF_BUILD='+str(out),flush=True)
    base.execute(['wmake'],out/'src',out/'log.build.solver',os.environ.copy())
    checked=verify(out,require_record=False);base.write_json(out/'vof_build_verified.json',checked)
    print('VOF reference application build/dependency verification: PASS (no CFD launched)',flush=True)
    print(json.dumps(dict(binary=checked['binary'],sha256=checked['sha256']),indent=2),flush=True)
    return out


def verify(out,require_record=True):
    out=out.resolve(strict=True);m=json.loads((out/'vof_build_manifest.json').read_text())
    if os.environ.get('WM_PROJECT_VERSION')!='9' or m['root']!=str(out) or m['recipe']!=recipe():
        raise ValueError('Wrong environment or changed VOF build recipe')
    project=Path(m['project'])
    if project!=Path(os.environ['WM_PROJECT_DIR']).resolve(strict=True) or native_sources(project)!=m['native_sources']:
        raise ValueError('Native OF9 source identity changed')
    expected={p:hashlib.sha256(t.encode()).hexdigest() for p,t in private_sources(project,out).items()}
    if base.source_hashes(out/'src')!=expected or m['private_sources']!=expected:
        raise ValueError('Private VOF source identity changed')
    checked=identity(out)
    if require_record and json.loads((out/'vof_build_verified.json').read_text())!=checked:
        raise ValueError('VOF binary or runtime dependencies changed')
    return checked
