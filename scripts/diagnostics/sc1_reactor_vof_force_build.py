"""Private two-step force observer derived from the verified VOF reference."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sc1_two_field as base
import sc1_reactor_startup as startup
import sc1_reactor_vof_force_instrument as instrument
import sc1_reactor_vof_build as reference

APP='sc1InterFoamForces'
HEADER=base.HERE/'include/sc1VofForceAudit.H'


def recipe():
    return {p.name:base.sha(p) for p in (Path(__file__),Path(instrument.__file__),HEADER)}


def private_sources(base_root,out):
    source=base_root/'src'
    result={rel:(source/rel).read_text() for rel in base.source_hashes(source)}
    result['interFoam.C']=instrument.solver(result['interFoam.C'])
    result['pEqn.H']=instrument.pressure(result['pEqn.H'])
    result['sc1VofForceAudit.H']=HEADER.read_text()
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


def build(base_root,out):
    startup.require_allocation(1)
    if os.environ.get('WM_PROJECT_VERSION')!='9' or not shutil.which('wmake'):
        raise ValueError('Source OpenFOAM 9 before building the force probe')
    base_root=base_root.resolve(strict=True);checked=reference.verify(base_root)
    out=base.fresh(out,(base_root,base.REPO,Path(os.environ['WM_PROJECT_DIR'])))
    sources=private_sources(base_root,out)
    m=dict(root=str(out),base_build=str(base_root),base_identity=checked,recipe=recipe(),
        private_sources={p:hashlib.sha256(t.encode()).hexdigest() for p,t in sources.items()})
    out.mkdir(parents=True,exist_ok=False);(out/'bin').mkdir();(out/'src').mkdir()
    for rel,text in sources.items():
        p=out/'src'/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    base.write_json(out/'force_build_manifest.json',m)
    print('VOF_FORCE_BUILD='+str(out),flush=True)
    base.execute(['wmake'],out/'src',out/'log.build.solver',os.environ.copy())
    verified=verify(out,require_record=False);base.write_json(out/'force_build_verified.json',verified)
    print('VOF force observer build/dependency verification: PASS (no CFD launched)',flush=True)
    print(json.dumps(dict(binary=verified['binary'],sha256=verified['sha256']),indent=2),flush=True)
    return out


def verify(out,require_record=True):
    out=out.resolve(strict=True);m=json.loads((out/'force_build_manifest.json').read_text())
    if os.environ.get('WM_PROJECT_VERSION')!='9' or m['root']!=str(out) or m['recipe']!=recipe():
        raise ValueError('Wrong environment or changed force-build recipe')
    base_root=Path(m['base_build'])
    if reference.verify(base_root)!=m['base_identity']:raise ValueError('Original VOF build changed')
    expected={p:hashlib.sha256(t.encode()).hexdigest() for p,t in private_sources(base_root,out).items()}
    if base.source_hashes(out/'src')!=expected or m['private_sources']!=expected:
        raise ValueError('Private force observer source identity changed')
    checked=identity(out)
    if require_record and json.loads((out/'force_build_verified.json').read_text())!=checked:
        raise ValueError('Force observer binary or runtime dependencies changed')
    # The observation-only application must retain the same shared libraries.
    if checked['dependencies']!=m['base_identity']['dependencies']:
        raise ValueError('Force observer runtime library mapping differs from control')
    for path,digest in m['base_identity']['assets'].items():
        if path!=m['base_identity']['binary'] and checked['assets'].get(path)!=digest:
            raise ValueError('Force observer runtime library content differs from control')
    return checked
