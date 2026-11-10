#!/usr/bin/env python3
"""Prepare one fresh observer control; verify by default, launch only with --run.

Uses the pristine manifest checkpoint, not a completed case's rewritten Co files.
No equations, numerical controls, closures, biology or existing cases are edited.
"""
import argparse
from collections import Counter, deque
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess

from prepare_sc1_alpha_audit import sha256
from prepare_sc1_conservation_pair import CHECKPOINT, hashes, validate_settings
from prepare_sc1_conservation_solver import scratch
from prepare_sc1_headspace_entrainment_trial import enclosing_brace_end
from sc1_alpha_convergence import SECTIONS, assets, variants
from sc1_alpha_source_probe import verify_probe
import summarize_sc1_conservation as summary

BUILD = Path('/scratch/nwimer/sc1_alpha_source_probe_20261110T191459Z')
BASELINE = Path('/scratch/nwimer/sc1_alpha_pair_20261110T155635Z')
MANIFEST = 'source_probe_case_manifest.json'
EXE_SHA = '0f38c43c6c986810c6bf4578391cac6c71fefc22d73b5b4aef05432920a1c11b'
LIB_SHA = '672a40b44bda7e65e43f38e32a27b3a1c2217471c13ee768cb9564914573a20a'
EXTENSIONS = {
    '/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so':
        'f823c8a3a3c83067171882cae53dfce4f9988dbc34c07bae6f6d58668c77ab4d',
    '/scratch/nwimer/sc1_headspace_entrainment_20261104/lib/libsc1HeadspaceEntrainmentModels.so':
        '27c8ce353108690a041454d2a38957b9c88862e2215f0f75e56fb6a3b3f847fb',
}
LIMITER_SHA = 'e63f7be7cee3503d3d1c7b9bb5131c791e5645b07bff7e061a3cca9a5ce2d5ac'


def check_build(build):
    for relative, expected in (('bin/putidaFoamSC1AlphaAudit', EXE_SHA),
                               ('lib/libbirdphaseSystem.so', LIB_SHA)):
        if sha256(build / relative) != expected:
            raise ValueError(f'Binary differs from verified source-probe build: {relative}')
    libbin = os.environ.get('FOAM_USER_LIBBIN')
    if not libbin:
        raise ValueError('Source OpenFOAM and load libfabric first')
    limiter = Path(libbin) / 'libalphaLimitVelocity.so'
    for path, expected in [*[(Path(p), h) for p, h in EXTENSIONS.items()],
                           (limiter, LIMITER_SHA)]:
        if sha256(path) != expected:
            raise ValueError(f'Case extension differs from verified build: {path}')
    verify_probe(build, [Path(p) for p in EXTENSIONS], [limiter])
    print('Verified source-probe binary hashes, original sources and private runtime dependencies.')


def enable_probe(text):
    """Find the explicit alpha dictionary; preserve every existing byte."""
    masked = re.sub(r'/\*.*?\*/|//[^\n]*', lambda m: ' ' * len(m[0]), text, flags=re.S)
    token_re = re.compile(r'"(?:\\.|[^"\\])*"|[{};]|[^\s{};]+')

    def dictionaries(start, end):
        tokens = list(token_re.finditer(masked, start, end))
        found, i = [], 0
        while i < len(tokens):
            token = tokens[i]
            if token[0].startswith(('#', '$')):
                raise ValueError('Source probe requires explicit dictionaries')
            if token[0] == '{':
                if i == 0:
                    raise ValueError('Unnamed dictionary')
                close = enclosing_brace_end(text, token.start())
                if close >= end:
                    raise ValueError('Unbalanced dictionary')
                found.append((tokens[i-1][0].strip('"'), token.start(), close))
                while i < len(tokens) and tokens[i].start() <= close:
                    i += 1
                continue
            if token[0] == '}':
                raise ValueError('Unexpected closing brace')
            i += 1
        return found

    solvers = [d for d in dictionaries(0, len(text)) if d[0] == 'solvers']
    if len(solvers) != 1:
        raise ValueError('Expected one explicit solvers dictionary')
    _, start, end = solvers[0]
    candidates = [d for d in dictionaries(start+1, end) if 'alpha' in d[0]]
    if len(candidates) != 1 or candidates[0][0] not in ('alpha', 'alpha.*'):
        raise ValueError('Expected one explicit alpha solver dictionary')
    _, start, end = candidates[0]
    body = masked[start+1:end]
    if re.search(r'[#${}]|\bsc1AlphaSourceProbe\b', body):
        raise ValueError('Expected explicit, unprobed alpha controls without nested dictionaries')
    for key, expected in (('sc1AlphaAudit', 'true'), ('nAlphaCorr', '1'), ('nAlphaSubCycles', '6')):
        if re.findall(r'\b' + key + r'\s+(\S+)\s*;', body) != [expected]:
            raise ValueError(f'Unexpected alpha control: {key}')
    return text[:end] + '\n        sc1AlphaSourceProbe true;\n    ' + text[end:]


def prepare(baseline_pair, out, build):
    if out.exists() or out.is_symlink():
        raise ValueError(f'Destination exists; preserve it: {out}')
    baseline_pair, build = baseline_pair.resolve(strict=True), build.resolve(strict=True)
    baseline, out = baseline_pair / 'control', out.resolve()
    common = json.loads((baseline_pair / 'comparison_manifest.json').read_text())
    if common['checkpoint'] != CHECKPOINT or common['duration_s'] != 0.1 or common['alpha_audit'] is not True:
        raise ValueError('Expected the completed 0.1 s alpha pair')
    source = Path(common['initial_state_source']).resolve(strict=True)
    for original in (baseline_pair, source, build):
        if original == out or original in out.parents or out in original.parents:
            raise ValueError('Destination overlaps a preserved input/build')
    state = common['initial_state_sha256']
    if set(state) != set(SECTIONS[2:]) or not all(state.values()):
        raise ValueError('Incomplete initial-state manifest')
    for section, expected in state.items():
        if hashes(source / section) != expected:
            raise ValueError(f'Original common checkpoint/mesh changed: {section}')
    inputs = {section: hashes(baseline / section) for section in SECTIONS[:2]}
    validate_settings(baseline, baseline_pair / 'normal')
    if (baseline / 'system/controlDict').read_bytes() != (baseline_pair / 'normal/system/controlDict').read_bytes():
        raise ValueError('Completed pair controlDict settings differ')
    phase = re.sub(r'/\*.*?\*/|//[^\n]*', '', (baseline / 'constant/phaseProperties').read_text(), flags=re.S)
    foam = re.findall(r'\bfoamAreaMultiplier\s+(\S+)\s*;', phase)
    if re.findall(r'\bentrainmentEfficiency\s+(\S+)\s*;', phase) != ['0'] or not foam or any(float(v) != 0 for v in foam):
        raise ValueError('Expected zero entrainment and zero foam multiplier')
    control = (baseline / 'system/controlDict').read_text()
    solution = (baseline / 'system/fvSolution').read_text()
    variants(control, solution)  # Validate baseline; do not use either numerical variant.
    patched = enable_probe(solution)
    log = baseline / 'log.conservation'
    starts, rows = summary.load_log(log)
    if (abs(starts['total']['t_s'] - float(CHECKPOINT)) > 1e-12
            or abs(rows['total'][-1]['t_s'] - float(CHECKPOINT) - 0.1) > 1e-9
            or rows['total'][-1]['steps'] != 2000):
        raise ValueError('Completed baseline has a different interval or step count')
    with log.open() as stream:
        if not any(line.strip() == 'SC1_ALPHA_LIBRARY version=1' for line in stream):
            raise ValueError('Baseline lacks private alpha-library marker')
    check_build(build)
    model_assets = assets(baseline)
    manifest = {'version': 1, 'baseline_pair': str(baseline_pair),
                'baseline_log': str(log), 'baseline_log_sha256': sha256(log),
                'build_root': str(build), 'executable_sha256': EXE_SHA, 'library_sha256': LIB_SHA,
                'initial_state_source': str(source), 'initial_state_sha256': state,
                'checkpoint': CHECKPOINT, 'duration_s': 0.1, 'observer_only': True,
                'baseline_settings_sha256': inputs, 'assets_sha256': model_assets,
                'changed_entries': {'solvers/alpha/sc1AlphaSourceProbe': [None, True]}}
    out.mkdir(parents=True, exist_ok=False)
    for section in SECTIONS:
        origin = source if section.startswith('processors52/') else baseline
        shutil.copytree(origin / section, out / section)
    for filename in model_assets:
        shutil.copy2(baseline / filename, out / filename)
    (out / 'system/fvSolution').write_text(patched)
    current = {section: hashes(out / section) for section in SECTIONS}
    expected_system = dict(inputs['system'], fvSolution=sha256(out / 'system/fvSolution'))
    if (current['constant'] != inputs['constant'] or current['system'] != expected_system
            or any(current[s] != h for s, h in state.items()) or assets(out) != model_assets):
        raise ValueError('Copied inputs changed; preserve incomplete output')
    if ({section: hashes(baseline / section) for section in SECTIONS[:2]} != inputs
            or assets(baseline) != model_assets or sha256(log) != manifest['baseline_log_sha256']):
        raise ValueError('Baseline changed during preparation; preserve incomplete output')
    manifest['inputs_sha256'] = current
    (out / MANIFEST).write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def preflight(case):
    if case.is_symlink() or not case.is_dir() or any(case.glob('log.*')):
        raise ValueError('Case missing, redirected or already attempted; preserve it')
    if any(path.is_symlink() for path in case.rglob('*')):
        raise ValueError('Prepared case contains a redirected input')
    manifest = json.loads((case / MANIFEST).read_text())
    if (manifest['version'] != 1 or manifest['observer_only'] is not True
            or manifest['checkpoint'] != CHECKPOINT or manifest['duration_s'] != 0.1
            or manifest['executable_sha256'] != EXE_SHA or manifest['library_sha256'] != LIB_SHA):
        raise ValueError('Unexpected source-probe case manifest')
    times = sorted(p.name for p in (case / 'processors52').iterdir() if p.is_dir() and p.name[0].isdigit())
    if times != [CHECKPOINT]:
        raise ValueError(f'Unexpected case times: {times}')
    expected = manifest['inputs_sha256']
    if set(expected) != set(SECTIONS):
        raise ValueError('Incomplete case input manifest')
    for section, digests in expected.items():
        if hashes(case / section) != digests:
            raise ValueError(f'Prepared input changed: {section}')
    if assets(case) != manifest['assets_sha256']:
        raise ValueError('Prepared model assets changed')
    if sha256(Path(manifest['baseline_log'])) != manifest['baseline_log_sha256']:
        raise ValueError('Completed baseline log changed')
    check_build(Path(manifest['build_root']))
    return manifest


def legacy_records(log):
    with log.open() as stream:
        return Counter(line.strip() for line in stream
                       if (line.startswith('SC1_') and not line.startswith('SC1_SOURCE_PROBE'))
                       or line.startswith('BUBBLE_MASS_AUDIT '))


def regression(case, manifest):
    """Require complete new records and report exact printed legacy regression."""
    log, baseline = case / 'log.conservation', Path(manifest['baseline_log'])
    if sha256(baseline) != manifest['baseline_log_sha256']:
        raise ValueError('Completed baseline log changed')
    _, budgets = summary.load_log(log)  # Includes solver-End and finite-budget checks.
    endpoint = budgets['total'][-1]['t_s']
    keys, markers, methods = [], 0, 0
    with log.open() as stream:
        for line in stream:
            if line.strip() == 'SC1_SOURCE_PROBE_LIBRARY version=1':
                markers += 1
            if line.startswith('SC1_SOURCE_PROBE_METHOD '):
                methods += 1
            if not line.startswith('SC1_SOURCE_PROBE '):
                continue
            row = dict(token.split('=', 1) for token in line.split()[1:])
            keys.append((row['phase'], int(row['callInStep']), row['bin']))
            fields = ('phaseCellSubcycleSamples', 'rawBelowSamples', 'rawAboveSamples',
                      'pressureBelowSamples', 'pressureAboveSamples', 'rawSource_kg',
                      'pressureSource_kg', 'split_kg', 'implicit_kg', 'rawCorrection_kg',
                      'rawExcess_kg', 'pressureExcess_kg', 'identityRemainder_kg',
                      'identityRemainderL1_kg', 'maxDtAbsSp', 't_s')
            if any(not math.isfinite(float(row[k])) for k in fields) or abs(float(row['t_s']) - endpoint) > 1e-12:
                raise ValueError('Invalid source-probe row or reporting time')
    expected = {(p, c, b) for p in ('liquid', 'bubbleGas') for c in (1, 2)
                for b in ('lower', 'upper', 'interior', 'outside')}
    if markers != 1 or methods != 2 or len(keys) != 16 or set(keys) != expected:
        raise ValueError('Missing, duplicate or unexpected source-probe records')
    previous, observed = legacy_records(baseline), legacy_records(log)
    if previous != observed:
        missing, extra = previous - observed, observed - previous
        raise ValueError(f'Printed legacy regression differs: {sum(missing.values())} baseline / '
                         f'{sum(extra.values())} probe records unmatched; preserve logs for review, do not rerun')
    print('Observer regression PASS: all legacy SC1 and bubble-mass records match exactly at printed precision.')
    print('Complete source probe: 16 rows, two outer indices, four alpha bins per solved phase.')


def report(case, manifest, solver_status=0):
    status = solver_status
    if not status:
        status = summary.main([manifest['baseline_log'], '--normal', str(case / 'log.conservation'), '--native-times'])
        print('Reader label "normal" means the observer control; entrainment is ZERO.')
        try:
            regression(case, manifest)
        except (OSError, ValueError, KeyError) as error:
            print(f'Observer regression requires review: {error}')
            status = 1
    tail = deque(maxlen=30)
    with (case / 'log.conservation').open(errors='replace') as stream:
        for line in stream:
            tail.append(line)
            if ((line.startswith('SC1_') and (not line.startswith('SC1_BOUNDARY ') or 'patch=outlet ' in line))
                    or line.startswith('BUBBLE_MASS_AUDIT ')):
                print(line, end='')
    print(''.join(tail), end='')
    return status


def launch(case, manifest, runner=subprocess.run):
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('No active allocation; no solver launched')
    build = Path(manifest['build_root'])
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = str(build / 'lib') + (':' + env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    print(f'Launching only {case}; output: {case}/log.conservation', flush=True)
    with (case / 'log.conservation').open('x') as log:
        result = runner(['srun', '-n', '52', str(build / 'bin/putidaFoamSC1AlphaAudit'),
                         '-parallel', '-fileHandler', 'collated'], cwd=case, env=env,
                        stdout=log, stderr=subprocess.STDOUT, check=False)
    print(f'SC1 source-probe solver status: {result.returncode}', flush=True)
    return result.returncode


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare', help='Create one fresh case; never launch CFD')
    prep.add_argument('--baseline-pair', type=Path, default=BASELINE)
    prep.add_argument('--build-root', type=Path, default=BUILD)
    prep.add_argument('--out', type=Path, required=True)
    run = sub.add_parser('run', help='Read-only preflight unless --run is supplied')
    run.add_argument('--case', type=Path, required=True)
    run.add_argument('--run', action='store_true')
    read = sub.add_parser('report', help='Read a completed/attempted case; never launch CFD')
    read.add_argument('--case', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == 'prepare':
            if args.out.is_symlink():
                raise ValueError('Destination is a symlink')
            case = scratch(args.out)
            prepare(scratch(args.baseline_pair), case, scratch(args.build_root))
            print(f'Prepared PROBECASE={case}; no CFD launched.')
            print('Only alpha/sc1AlphaSourceProbe enabled; all physics, numerics and biology preserved.')
        else:
            if args.case.is_symlink():
                raise ValueError('Case is a symlink')
            case = scratch(args.case)
            if args.command == 'report':
                return report(case, json.loads((case / MANIFEST).read_text()))
            manifest = preflight(case)
            print('Pristine original checkpoint, case inputs and verified new build match.')
            if args.run:
                return report(case, manifest, launch(case, manifest))
            print('Read-only preflight PASS; no solver launched.')
        return 0
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        print(f'SC1 source-probe run stopped: {error}\nPreserve all existing/partial cases and builds.')
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
