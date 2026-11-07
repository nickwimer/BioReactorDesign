#!/usr/bin/env python3
"""Read-only preflight for an isolated alpha-audit build; never launch the solver."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess

from prepare_sc1_alpha_audit import MARKER, sha256, source_hashes


def command_result(command, environment):
    result=subprocess.run(command,env=environment,text=True,capture_output=True,
                          timeout=30,check=False)
    return {'command':command,'returncode':result.returncode,
            'stdout':result.stdout,'stderr':result.stderr}


def marker_type(text):
    matches=re.findall(r'^\s*(?:[0-9a-fA-F]+\s+)?([A-Za-z?])\s+'
                       +MARKER+r'(?:@\S+)?\s*$',text,re.M)
    return matches


def private_mapping(text, private):
    lines=[line for line in text.splitlines() if 'libbirdphaseSystem' in line]
    if len(lines)!=1:
        raise ValueError('Expected exactly one libbirdphaseSystem.so dependency mapping')
    match=re.match(r'^\s*libbirdphaseSystem\.so\s+=>\s+(\S+)\s+\(',lines[0])
    if not match:
        raise ValueError('Unrecognized or absolute phase-system dependency: '+lines[0].strip())
    resolved=Path(match[1]).resolve(strict=True)
    if resolved!=private:
        raise ValueError(f'Phase-system fallback: {resolved}; expected {private}')
    return str(resolved)


def verify(root, loaded_libraries=(), extra_libraries=(), runner=command_result):
    root=root.resolve(strict=True)
    report={'root':str(root),'checks':[],'commands':[],'errors':[]}
    environment=os.environ.copy()
    environment['LD_LIBRARY_PATH']=str(root/'lib')+(
        ':'+environment['LD_LIBRARY_PATH'] if environment.get('LD_LIBRARY_PATH') else '')
    report['LD_LIBRARY_PATH']=environment['LD_LIBRARY_PATH']
    try:
        manifest=json.loads((root/'alpha_audit_manifest.json').read_text())
        if manifest.get('alpha_audit_version')!=1:
            raise ValueError('Unsupported alpha audit manifest version')
        executable=root/'bin/putidaFoamSC1AlphaAudit'
        private=root/'lib/libbirdphaseSystem.so'
        for name,path,key in (('executable',executable,'executable'),('private library',private,'private_library')):
            if Path(manifest[key])!=path or path.resolve(strict=True)!=path or not path.is_file():
                raise ValueError(f'Unexpected, missing or redirected {name}: {path}')
            report['checks'].append({'artifact':name,'path':str(path),'sha256':sha256(path)})
        if not os.access(executable,os.X_OK):
            raise ValueError(f'Audit executable is not executable: {executable}')
        installed=Path(manifest['installed_library'])
        installed_hash=sha256(installed)
        if installed_hash!=manifest['installed_library_sha256']:
            raise ValueError(f'Installed phase-system library changed: {installed}')
        report['checks'].append({'artifact':'unchanged installed library','path':str(installed),'sha256':installed_hash})
        for path,expected,label in (
            (Path(manifest['phase_source']),manifest['private_phase_source_sha256_before'],'original phase sources'),
            (root/'phaseSystems',manifest['private_phase_source_sha256_after'],'private phase sources'),
        ):
            if source_hashes(path)!=expected:
                raise ValueError(f'Changed {label}: {path}')
            report['checks'].append({'artifact':label,'path':str(path),'files':len(expected),'unchanged':True})
        solver_hashes=manifest['solver_source_sha256_after_alpha']
        if source_hashes(root/'solver')!=solver_hashes:
            raise ValueError('Changed prepared solver sources or diagnostic headers')
        report['checks'].append({'artifact':'prepared solver sources','path':str(root/'solver'),
                                 'files':len(solver_hashes),'unchanged':True})
        for binary,expected in ((private,['T']),(executable,['U'])):
            result=runner(['nm','-D',str(binary)],environment)
            report['commands'].append(result)
            if result['returncode'] or marker_type(result['stdout'])!=expected:
                raise ValueError(f'Missing/incorrect mandatory alpha marker in {binary}; expected {expected[0]}')
        binaries=[(executable,True)]
        binaries.extend((Path(path).resolve(strict=True),True) for path in loaded_libraries)
        binaries.extend((Path(path).resolve(strict=True),False) for path in extra_libraries)
        for binary,require_phase in binaries:
            if not binary.is_file():
                raise ValueError(f'Not a runtime binary: {binary}')
            result=runner(['ldd',str(binary)],environment)
            report['commands'].append(result)
            if result['returncode'] or 'not found' in result['stdout']+result['stderr']:
                raise ValueError(f'Unresolved runtime dependencies: {binary}')
            try:
                resolved=(private_mapping(result['stdout'],private)
                          if require_phase or 'libbirdphaseSystem' in result['stdout'] else None)
            except (ValueError,OSError) as error:
                raise ValueError(f'{binary}: {error}') from error
            report['checks'].append({'artifact':'runtime dependency','path':str(binary),
                                     'sha256':sha256(binary),'phase_library':resolved})
        report['coverage']=('Executable and explicitly supplied libraries only; '
                            'supply every applicable case-loaded extension.')
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        report['errors'].append(str(error))
    report['ok']=not report['errors']
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--loaded-library',type=Path,action='append',default=[],
                        help='Case-loaded phase-system-dependent library; repeat for each')
    parser.add_argument('--extra-library',type=Path,action='append',default=[],
                        help='Additional case library; any phase-system dependency must be private')
    parser.add_argument('--verbose-json',action='store_true',
                        help='Print the complete report, including ldd/nm output, as JSON')
    args=parser.parse_args(argv)
    try:
        report=verify(args.root,args.loaded_library,args.extra_library)
    except (OSError,ValueError) as error:
        report={'ok':False,'checks':[],'commands':[],'errors':[str(error)]}
    if args.verbose_json:
        print(json.dumps(report,indent=2))
    else:
        for check in report['checks']:
            detail=('sha256='+check['sha256']) if 'sha256' in check else f'unchanged files={check["files"]}'
            print(f'{check["artifact"]}: {check["path"]} {detail}')
            if check.get('phase_library'):
                print('  libbirdphaseSystem.so => '+check['phase_library'])
            elif 'phase_library' in check:
                print('  no phase-system dependency (extra-library check)')
        if report['errors']:
            for error in report['errors']:
                print('SC1 runtime verification stopped: '+error)
            if report['commands'] and report['commands'][-1]['command'][0]=='ldd':
                result=report['commands'][-1]
                print(result['stdout'],end='')
                print(result['stderr'],end='')
        else:
            print(report['coverage'])
            print('SC1 runtime verification: PASS (no solver launched)')
    return 0 if report['ok'] else 1


if __name__=='__main__':
    raise SystemExit(main())
