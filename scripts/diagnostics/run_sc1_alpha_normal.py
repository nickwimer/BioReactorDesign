#!/usr/bin/env python3
"""Verify the prepared SC1 normal case; launch only with explicit --run.

The preparation manifest remains the initial-state reference after control
finishes. Completed control Co files may differ ONLY in a proved outer-header
location. No file is normalized, copied, replaced or repaired by this tool.
"""
import argparse
from collections import deque
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from prepare_sc1_alpha_audit import sha256
from prepare_sc1_conservation_pair import CHECKPOINT, hashes, validate_settings
from prepare_sc1_conservation_solver import scratch
from sc1_checkpoint_metadata import COURANT_FILES, verify_copied_courant_location_metadata
import summarize_sc1_conservation as summary
from verify_sc1_alpha_runtime import verify

EXE_SHA='43aeeef4984dd016ba94664bfe64fb498458237a6ab925e15a2322357a49ca96'
LIB_SHA='bc6d9da6e63caf8392917fd318d41097f1a488826730b36be33f0e41e9b6dd15'
REQUIRED_FIELDS={'alpha.liquid','alpha.bubbleGas','alpha.headspaceGas',
                 'O2.liquid','O2.bubbleGas','O2.headspaceGas','f.bubbleGas'}


def differences(expected, actual):
    return {'changed':sorted(n for n in expected.keys() & actual.keys() if expected[n]!=actual[n]),
            'missing':sorted(expected.keys()-actual.keys()),
            'extra':sorted(actual.keys()-expected.keys())}


def preflight_pair(root):
    """Read-only, fail closed. Return complete metadata proofs when applicable."""
    root=Path(root).resolve(strict=True)
    control,normal=root/'control',root/'normal'
    for case in (control,normal):
        for path in (case,case/'processors52',case/'constant',case/'system'):
            if path.is_symlink() or not path.is_dir():
                raise ValueError(f'Expected an ordinary prepared directory: {path}')
    if any(normal.glob('log.*')):
        raise ValueError('Normal log already exists; preserve this case')
    manifest=json.loads((root/'comparison_manifest.json').read_text())
    if (manifest['checkpoint']!=CHECKPOINT or manifest['alpha_audit'] is not True
            or manifest['duration_s']!=0.1):
        raise ValueError('Expected the prepared 0.1 s alpha-audit pair')
    times=sorted(p.name for p in (normal/'processors52').iterdir()
                 if p.is_dir() and p.name[0].isdigit())
    if times!=[CHECKPOINT]:
        raise ValueError(f'Normal case has unexpected time directories: {times}')
    sections=manifest['initial_state_sha256']
    time_section=f'processors52/{CHECKPOINT}'
    if set(sections)!= {'processors52/constant',time_section}:
        raise ValueError('Unexpected initial-state manifest sections')
    if not (REQUIRED_FIELDS | COURANT_FILES)<=sections[time_section].keys():
        raise ValueError('Initial-state manifest lacks required fields')
    source=Path(manifest['initial_state_source'])
    if not source.is_absolute() or str(source)!=manifest['control_source']:
        raise ValueError('Unexpected common checkpoint source in manifest')
    for section,expected in sections.items():
        if not expected or any(not re.fullmatch(r'[0-9a-f]{64}',v) for v in expected.values()):
            raise ValueError(f'Invalid hash manifest: {section}')
        for case in (control,normal):
            if (case/section).is_symlink():
                raise ValueError(f'Redirected checkpoint/mesh: {case/section}')
        actual=hashes(normal/section)
        if actual!=expected:
            raise ValueError(f'Normal initial state changed: {section}: {differences(expected,actual)}')
    # Normal is proved pristine before any header exception is considered.
    proofs={}
    for section,expected in sections.items():
        actual=hashes(control/section)
        diff=differences(expected,actual)
        if diff['missing'] or diff['extra']:
            raise ValueError(f'Completed-control file set changed: {section}: {diff}')
        for name in diff['changed']:
            if section!=time_section or name not in COURANT_FILES:
                raise ValueError(f'Completed-control non-Courant input changed: {section}/{name}')
            try:
                proof=verify_copied_courant_location_metadata(
                    control/section/name,normal/section/name,source/section/name)
            except (OSError,ValueError) as error:
                raise ValueError(f'Completed-control {name} is NOT proved location-only: {error}') from error
            if proof['control_sha256']!=actual[name] or proof['normal_sha256']!=expected[name]:
                raise ValueError(f'Checkpoint changed during proof: {name}')
            proofs[name]=proof
    validate_settings(control,normal)
    if (control/'system/controlDict').read_bytes()!=(normal/'system/controlDict').read_bytes():
        raise ValueError('Prepared controlDict settings differ')
    starts,rows=summary.load_log(control/'log.conservation')
    if abs(starts['total']['t_s']-float(CHECKPOINT))>1e-12:
        raise ValueError('Completed control started at a different checkpoint')
    if abs(rows['total'][-1]['t_s']-float(CHECKPOINT)-0.1)>1e-4:
        raise ValueError('Completed control has an unexpected duration')
    with (control/'log.conservation').open() as stream:
        if not any(line.strip()=='SC1_ALPHA_LIBRARY version=1' for line in stream):
            raise ValueError('Completed control lacks the private alpha-library marker')
    return {'root':str(root),'normal_initial_files':sum(map(len,sections.values())),
            'control_location_proofs':proofs,'completed_control_steps':rows['total'][-1]['steps']}


def check_build(build):
    if sha256(build/'bin/putidaFoamSC1AlphaAudit')!=EXE_SHA:
        raise ValueError('Executable hash differs from the completed control')
    if sha256(build/'lib/libbirdphaseSystem.so')!=LIB_SHA:
        raise ValueError('Private library hash differs from the completed control')
    libbin=os.environ.get('FOAM_USER_LIBBIN')
    if not libbin:
        raise ValueError('Source OpenFOAM and load libfabric before checking/running')
    report=verify(build,(
        Path('/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so'),
        Path('/scratch/nwimer/sc1_headspace_entrainment_20261104/lib/libsc1HeadspaceEntrainmentModels.so'),
    ),(Path(libbin)/'libalphaLimitVelocity.so',))
    if not report['ok']:
        raise ValueError('Runtime verification failed: '+'; '.join(report['errors']))
    print('Verified unchanged binaries, original sources and private runtime dependencies.')


def launch(root, build, runner=subprocess.run):
    """Called only after preflight; exclusive log creation also prevents retries."""
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('No active SLURM allocation; solver not launched')
    normal=root/'normal'
    environment=os.environ.copy()
    environment['LD_LIBRARY_PATH']=str(build/'lib')+(
        ':'+environment['LD_LIBRARY_PATH'] if environment.get('LD_LIBRARY_PATH') else '')
    command=['srun','-n','52',str(build/'bin/putidaFoamSC1AlphaAudit'),
             '-parallel','-fileHandler','collated']
    print(f'Launching only {normal}; output: {normal}/log.conservation',flush=True)
    with (normal/'log.conservation').open('x') as log:
        result=runner(command,cwd=normal,env=environment,stdout=log,
                      stderr=subprocess.STDOUT,check=False)
    print(f'SC1 alpha-normal solver status: {result.returncode}',flush=True)
    return result.returncode


def report_run(root, status):
    log=root/'normal/log.conservation'
    summary_status=0
    if status==0:
        summary_status=summary.main([str(root/'control/log.conservation'),
                                    '--normal',str(log),'--native-times'])
    tail=deque(maxlen=30)
    with log.open(errors='replace') as stream:
        for line in stream:
            tail.append(line)
            if (line.startswith('SC1_') and
                (not line.startswith('SC1_BOUNDARY ') or 'patch=outlet ' in line)) or line.startswith(
                    ('BUBBLE_MASS_AUDIT ','BUBBLE_STAGE_AUDIT ','HEADSPACE_ENTRAINMENT_PB_AUDIT ')):
                print(line,end='')
    print(''.join(tail),end='')
    return status if status else summary_status


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pair-root',type=Path,
                        default=Path('/scratch/nwimer/sc1_alpha_pair_20261110T155635Z'))
    parser.add_argument('--audit-root',type=Path,
                        default=Path('/scratch/nwimer/sc1_alpha_audit_20261107T001119Z'))
    parser.add_argument('--run',action='store_true',help='Launch normal once after all checks; default is read-only')
    args=parser.parse_args(argv)
    try:
        root,build=scratch(args.pair_root),scratch(args.audit_root)
        check_build(build)
        report=preflight_pair(root)
        print(f'Normal initial state matches preparation manifest: {report["normal_initial_files"]} files.')
        for name,proof in report['control_location_proofs'].items():
            print(f'Proved completed-control outer-location-only change: {name}; '
                  f'unchanged payload SHA256={proof["payload_sha256"]}')
        print('Pair settings and completed control verified. No field files changed.')
        if not args.run:
            print('Read-only preflight PASS; no solver launched.')
            return 0
        return report_run(root,launch(root,build))
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print(f'SC1 normal stopped: {error}\nPreserve this pair; do not overwrite or restart it.')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
