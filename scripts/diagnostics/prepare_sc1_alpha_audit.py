#!/usr/bin/env python3
"""Prepare a private SC1 solver and phase-system library; never compile or run CFD.

The existing conservation preparer retains the live solver's biology and other
local changes. Only the new phaseSystems symlink is replaced by copied sources.
An existing or partially prepared destination is never reused.
"""
import argparse
import fnmatch
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from prepare_sc1_conservation_solver import ROOT, replace_once, scratch

HERE=Path(__file__).resolve().parent
PHASE_HEADERS=ROOT/'OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseSystem'
HEADERS=('sc1AlphaSubcycleAudit.H','sc1AlphaSubcycleLedger.H')
DEPENDENCIES=('include','phaseSystems','interfacialModels',
              'interfacialCompositionModels','multiphaseCompressibleMomentumTransportModels')
GENERATED=('lnInclude','platforms','linux*','*.o','*.dep','*.so','*.a','build*.log','.git')
MARKER='sc1AlphaAuditLibraryVersion'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ignore_generated(directory, names):
    return {name for name in names if any(fnmatch.fnmatch(name,p) for p in GENERATED)}


def source_hashes(root):
    """Hash every copied input, following links but rejecting cycles/missing files."""
    result={}
    def visit(directory, relative, ancestors):
        resolved=directory.resolve(strict=True)
        if resolved in ancestors:
            raise ValueError(f'Source directory symlink cycle: {directory}')
        entries=sorted(directory.iterdir())
        ignored=ignore_generated(directory,[p.name for p in entries])
        for path in entries:
            if path.name in ignored:
                continue
            name=relative/path.name
            if path.is_dir():
                visit(path,name,ancestors|{resolved})
            elif path.is_file():
                result[name.as_posix()]=sha256(path)
            else:
                raise ValueError(f'Unreadable/non-file source entry: {path}')
    visit(root,Path(),set())
    return result


def safe_make_path(path):
    # Paths enter literal Make syntax; do not reinterpret whitespace/metacharacters.
    if not re.fullmatch(r'/[A-Za-z0-9_./+\-]+',str(path)):
        raise ValueError(f'Unsupported Make path: {path}')
    return str(path)


def run_base_preparer(source, out):
    completed=subprocess.run(
        [sys.executable,str(HERE/'prepare_sc1_conservation_solver.py'),
         '--source',str(source),'--out',str(out)],
        text=True,capture_output=True,check=False,
    )
    if completed.returncode:
        raise ValueError('Conservation preparation failed:\n'+completed.stdout+completed.stderr)


def bind_private_library(main, options, library_dir):
    if MARKER in main or 'SC1_ALPHA_PRIVATE_LIBRARY' in options:
        raise ValueError('Alpha library binding already/partially installed')
    main=replace_once(main,'#include "phaseSystem.H"',
                      '#include "phaseSystem.H"\nextern "C" int '+MARKER+'();')
    anchor='''        sc1MassStages.reset(new sc1MassStageAudit(fluid, LTS));
    }

    while (pimple.run(runTime))'''
    main=replace_once(main,anchor,'''        sc1MassStages.reset(new sc1MassStageAudit(fluid, LTS));
    }

    const int sc1AlphaLibraryVersion = sc1AlphaAuditLibraryVersion();
    if (sc1AlphaLibraryVersion != 1)
    {
        FatalErrorInFunction << "Wrong SC1 alpha-audit library version"
            << exit(FatalError);
    }
    Info<< "SC1_ALPHA_LIBRARY version=" << sc1AlphaLibraryVersion << endl;

    while (pimple.run(runTime))''')
    # Insert before the original assignment contents, never after other -L flags.
    pattern=re.compile(r'(?m)^([ \t]*EXE_LIBS[ \t]*=[ \t]*)')
    if len(pattern.findall(options))!=1:
        raise ValueError('Expected one EXE_LIBS assignment; inspect local options')
    options=pattern.sub(lambda m:m[1]+'-L'+safe_make_path(library_dir)+' ',options)
    options+='\n# SC1_ALPHA_PRIVATE_LIBRARY: private search path precedes installed libraries.\n'
    return main,options


def prepare(source, out, installed_library, base_prepare=run_base_preparer):
    """Prepare new paths; the CLI additionally restricts source/out to scratch."""
    if out.exists() or out.is_symlink():
        raise ValueError(f'Destination exists (including partial preparation): {out}')
    source=source.resolve(strict=True)
    out=out.resolve()
    if source==out or source in out.parents or out in source.parents:
        raise ValueError('Source and destination must not overlap')
    safe_make_path(out)
    if installed_library.name!='libbirdphaseSystem.so':
        raise ValueError('Installed library must be named libbirdphaseSystem.so')
    installed_library=installed_library.resolve(strict=True)
    if not installed_library.is_file():
        raise ValueError('Installed phase-system library is not a readable file')
    if out==installed_library or out in installed_library.parents or installed_library in out.parents:
        raise ValueError(f'Destination overlaps installed library: {installed_library}')
    installed_hash=sha256(installed_library)
    dependencies={}
    for name in DEPENDENCIES:
        target=(source.parent/name).resolve(strict=True)
        if not target.is_dir():
            raise ValueError(f'Missing source dependency: {target}')
        if out==target or target in out.parents or out in target.parents:
            raise ValueError(f'Destination overlaps source dependency: {target}')
        if name not in ('include','phaseSystems') and not (target/'lnInclude').is_dir():
            raise ValueError(f'Missing prebuilt dependency headers: {target}/lnInclude')
        dependencies[name]=target
    phase_source=dependencies['phaseSystems']
    original_hashes=source_hashes(phase_source)
    solve_name='phaseSystem/phaseSystemSolve.C'
    solve=(phase_source/solve_name).read_text()
    from instrument_sc1_alpha_subcycles import instrument_phase_system_solve
    modified_solve=instrument_phase_system_solve(solve)
    headers={name:(PHASE_HEADERS/name).read_bytes() for name in HEADERS}
    library_files=(phase_source/'Make/files').read_text()
    library_files,count=re.subn(r'(?m)^[ \t]*LIB[ \t]*=.*$',
                               'LIB = '+safe_make_path(out/'lib/libbirdphaseSystem'),library_files)
    if count!=1:
        raise ValueError('Expected one phaseSystems LIB target')
    # All available-input checks above precede creating the destination.
    base_prepare(source,out)
    private=out/'phaseSystems'
    if not private.is_symlink() or private.resolve()!=phase_source:
        raise ValueError('New conservation tree has an unexpected phaseSystems dependency')
    private.unlink()  # Unlink only the new tree's dependency, never its target.
    shutil.copytree(phase_source,private,symlinks=False,ignore=ignore_generated)
    copied_hashes=source_hashes(private)
    if copied_hashes!=original_hashes:
        raise ValueError('Source changed during phaseSystems copy; preserve this partial tree')
    if any(path.is_symlink() for path in private.rglob('*')):
        raise ValueError('Private phaseSystems copy unexpectedly contains symlinks')
    (out/'lib').mkdir(exist_ok=False)
    (private/solve_name).write_text(modified_solve)
    (private/'Make/files').write_text(library_files)
    for name,content in headers.items():
        (private/'phaseSystem'/name).write_bytes(content)
    solver=out/'solver'
    changed_files=('putidaFoam.C','Make/options','Make/files')
    solver_before={name:sha256(solver/name) for name in changed_files}
    main,options=bind_private_library((solver/'putidaFoam.C').read_text(),
                                     (solver/'Make/options').read_text(),out/'lib')
    executable=out/'bin/putidaFoamSC1AlphaAudit'
    files,count=re.subn(r'(?m)^[ \t]*EXE[ \t]*=.*$',
                       'EXE = '+safe_make_path(executable),(solver/'Make/files').read_text())
    if count!=1:
        raise ValueError('Expected one copied solver EXE target')
    (solver/'putidaFoam.C').write_text(main)
    (solver/'Make/options').write_text(options)
    (solver/'Make/files').write_text(files)
    if source_hashes(phase_source)!=original_hashes or sha256(installed_library)!=installed_hash:
        raise ValueError('Original phase source/library changed during preparation')
    manifest={
        'alpha_audit_version':1, 'prepared_only':True,
        'solver_source':str(source),'phase_source':str(phase_source),
        'executable':str(executable),'private_library':str(out/'lib/libbirdphaseSystem.so'),
        'installed_library':str(installed_library),'installed_library_sha256':installed_hash,
        'source_dependencies':{name:str(path) for name,path in dependencies.items()},
        'private_phase_source_sha256_before':copied_hashes,
        'private_phase_source_sha256_after':source_hashes(private),
        'instrumenter_sha256':sha256(HERE/'instrument_sc1_alpha_subcycles.py'),
        'solver_sha256_before_alpha':solver_before,
        'solver_sha256_after_alpha':{name:sha256(solver/name) for name in changed_files},
        'solver_source_sha256_after_alpha':source_hashes(solver),
        'runtime_requirement':'Prepend this root/lib to process-local LD_LIBRARY_PATH; '
            'verify executable and case-loaded libraries resolve only the private libbirdphaseSystem.so.',
    }
    base_manifest_path=out/'source_manifest.json'
    base_manifest=json.loads(base_manifest_path.read_text())
    base_manifest['dependencies']['phaseSystems']=str(private)
    base_manifest['alpha_audit_version']=1
    base_manifest['executable']=str(executable)
    base_manifest_path.write_text(json.dumps(base_manifest,indent=2)+'\n')
    # Completion marker written last; a failed partial destination is never reused.
    (out/'alpha_audit_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,
                        default=Path('/scratch/nwimer/sc1_bubble_audit_20261029/solver'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--installed-library',type=Path,
                        help='Default: $FOAM_USER_LIBBIN/libbirdphaseSystem.so')
    args=parser.parse_args(argv)
    try:
        if args.out.is_symlink():
            raise ValueError(f'Destination is a symlink: {args.out}')
        installed=args.installed_library
        if installed is None:
            libbin=os.environ.get('FOAM_USER_LIBBIN')
            if not libbin:
                raise ValueError('Set up OpenFOAM or provide --installed-library')
            installed=Path(libbin)/'libbirdphaseSystem.so'
        manifest=prepare(scratch(args.source),scratch(args.out),installed)
        print('Prepared private alpha audit: '+manifest['executable'])
        print('Private library target: '+manifest['private_library'])
        print('No compilation or CFD launched. Build phaseSystems with wmake libso, '
              'then solver with wmake; never use Allwmake.')
        return 0
    except (OSError,ValueError,ImportError) as error:
        print(f'SC1 alpha preparation stopped: {error}')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
