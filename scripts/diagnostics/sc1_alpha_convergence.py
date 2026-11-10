#!/usr/bin/env python3
"""Prepare isolated numerical checks from the completed zero-entrainment control.

No solver rebuild or physics edits. Preparation creates two fresh cases;
run is read-only unless --run is supplied, and launches just one selected case.
"""
import argparse
from collections import deque
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

from prepare_sc1_conservation_pair import CHECKPOINT, hashes
from prepare_sc1_conservation_solver import scratch
from prepare_sc1_headspace_entrainment_trial import enclosing_brace_end
from prepare_sc1_headspace_o2_trial import set_control_entry, top_level_control_matches
from prepare_sc1_alpha_audit import sha256
from run_sc1_alpha_normal import check_build
import summarize_sc1_conservation as summary

VARIANTS=('four_outer','half_dt')
SECTIONS=('constant','system','processors52/constant',f'processors52/{CHECKPOINT}')


def value(text, key):
    matches=top_level_control_matches(text,key)
    if len(matches)!=1:
        raise ValueError(f'Expected one explicit top-level {key}')
    start,end=matches[0]
    return re.sub(r'/\*.*?\*/|//[^\n]*','',text[start:end],flags=re.S).strip()[len(key):].strip().rstrip(';').strip()


def pimple_bounds(text):
    masked=re.sub(r'/\*.*?\*/|//[^\n]*',lambda m:' '*len(m[0]),text,flags=re.S)
    matches=list(re.finditer(r'\bPIMPLE\s*\{',masked))
    if len(matches)!=1:
        raise ValueError('Expected one explicit PIMPLE dictionary')
    start=matches[0].end()-1
    return start,enclosing_brace_end(text,start)


def variants(control, solution):
    """Change only time-step controls or outer count, from the reviewed baseline."""
    start,end=pimple_bounds(solution)
    body=solution[start+1:end]
    if (int(value(body,'nOuterCorrectors'))!=2
            or value(body,'sc1ConservationAudit')!='true'):
        raise ValueError('Expected the audited two-outer-iteration baseline')
    if (float(value(control,'deltaT'))!=5e-5 or float(value(control,'maxDeltaT'))!=5e-5
            or value(control,'adjustTimeStep') not in ('yes','true','on')
            or value(control,'startFrom')!='latestTime'
            or abs(float(value(control,'endTime'))-float(CHECKPOINT)-0.1)>1e-9):
        raise ValueError('Expected the reviewed 0.1 s, 5e-5 s timestep baseline')
    if not re.search(r'\bsc1AlphaAudit\s+true\s*;',solution):
        raise ValueError('Baseline alpha audit is not enabled')
    four=solution[:start+1]+set_control_entry(body,'nOuterCorrectors','4')+solution[end:]
    half=set_control_entry(set_control_entry(control,'deltaT','2.5e-05'),'maxDeltaT','2.5e-05')
    return {'four_outer':(control,four),'half_dt':(half,solution)}


def assets(case):
    return {p.name:sha256(p) for p in sorted(case.iterdir()) if p.is_file() and
            (p.suffix in ('.pt','.txt','.onnx') or p.name.startswith('scaler_'))}


def prepare(baseline_pair, out):
    if out.exists() or out.is_symlink():
        raise ValueError(f'Destination exists; preserve it: {out}')
    baseline_pair=baseline_pair.resolve(strict=True)
    baseline=baseline_pair/'control'
    out=out.resolve()
    common=json.loads((baseline_pair/'comparison_manifest.json').read_text())
    if common['checkpoint']!=CHECKPOINT or common['duration_s']!=0.1 or common['alpha_audit'] is not True:
        raise ValueError('Expected the reviewed completed alpha pair')
    source=Path(common['initial_state_source']).resolve(strict=True)
    for original in (baseline_pair,source):
        if original==out or original in out.parents or out in original.parents:
            raise ValueError('Destination overlaps a preserved input')
    state=common['initial_state_sha256']
    if set(state)!={'processors52/constant',f'processors52/{CHECKPOINT}'}:
        raise ValueError('Unexpected common initial-state manifest')
    for section,expected in state.items():
        if hashes(source/section)!=expected:
            raise ValueError(f'Original common checkpoint/mesh changed: {source/section}')
    # Use the completed control's reviewed configuration, but never its rewritten
    # initial Co files. The pristine original source still matches the manifest.
    initial_inputs={section:hashes(baseline/section) for section in ('constant','system')}
    phase=(baseline/'constant/phaseProperties').read_text()
    phase=re.sub(r'/\*.*?\*/|//[^\n]*','',phase,flags=re.S)
    if re.findall(r'\bentrainmentEfficiency\s+(\S+)\s*;',phase)!=['0']:
        raise ValueError('Numerical checks require the zero-entrainment control')
    foam=re.findall(r'\bfoamAreaMultiplier\s+(\S+)\s*;',phase)
    if not foam or any(float(v)!=0 for v in foam):
        raise ValueError('Foam multiplier must remain zero')
    log=baseline/'log.conservation'
    starts,rows=summary.load_log(log)
    if (abs(starts['total']['t_s']-float(CHECKPOINT))>1e-12
            or abs(rows['total'][-1]['t_s']-float(CHECKPOINT)-0.1)>1e-4):
        raise ValueError('Completed baseline interval differs from the reviewed control')
    changes=variants((baseline/'system/controlDict').read_text(),
                     (baseline/'system/fvSolution').read_text())
    asset_hashes=assets(baseline)
    manifest={'version':1,'baseline_pair':str(baseline_pair),'baseline_log':str(log),
              'baseline_log_sha256':sha256(log),'initial_state_source':str(source),
              'checkpoint':CHECKPOINT,'duration_s':0.1,'initial_state_sha256':state,
              'baseline_settings_sha256':initial_inputs,'assets_sha256':asset_hashes,
              'physics':'same zero-entrainment control; biology schedule unchanged',
              'cases':{}}
    out.mkdir(parents=True,exist_ok=False)
    for name,(control,solution) in changes.items():
        case=out/name
        case.mkdir()
        for section in SECTIONS:
            origin=source if section.startswith('processors52/') else baseline
            shutil.copytree(origin/section,case/section)
        for filename in asset_hashes:
            shutil.copy2(baseline/filename,case/filename)
        (case/'system/controlDict').write_text(control)
        (case/'system/fvSolution').write_text(solution)
        current={section:hashes(case/section) for section in SECTIONS}
        for section,expected in state.items():
            if current[section]!=expected:
                raise ValueError(f'Copied state differs: {case/section}')
        if current['constant']!=initial_inputs['constant'] or assets(case)!=asset_hashes:
            raise ValueError('Copied physics/biology/assets changed')
        for key,expected in initial_inputs['system'].items():
            if key not in ('controlDict','fvSolution') and current['system'].get(key)!=expected:
                raise ValueError(f'Copied numerical setting changed unexpectedly: {key}')
        manifest['cases'][name]={'inputs_sha256':current,
            'changed_entries':({'PIMPLE/nOuterCorrectors':[2,4]} if name=='four_outer' else
                               {'deltaT':[5e-5,2.5e-5],'maxDeltaT':[5e-5,2.5e-5]})}
    if ({section:hashes(baseline/section) for section in ('constant','system')}!=initial_inputs
            or assets(baseline)!=asset_hashes or sha256(log)!=manifest['baseline_log_sha256']):
        raise ValueError('Baseline changed during preparation; preserve partial output')
    (out/'numerical_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return manifest


def preflight(root, name):
    manifest=json.loads((root/'numerical_manifest.json').read_text())
    if manifest['version']!=1 or set(manifest['cases'])!=set(VARIANTS):
        raise ValueError('Unexpected numerical experiment manifest')
    case=root/name
    if case.is_symlink() or not case.is_dir() or any(case.glob('log.*')):
        raise ValueError('Case missing, redirected or already attempted; preserve it')
    if any(path.is_symlink() for path in case.rglob('*')):
        raise ValueError('Prepared numerical case contains a redirected input')
    times=sorted(p.name for p in (case/'processors52').iterdir()
                 if p.is_dir() and p.name[0].isdigit())
    if times!=[CHECKPOINT]:
        raise ValueError(f'Unexpected case times: {times}')
    expected=manifest['cases'][name]['inputs_sha256']
    if set(expected)!=set(SECTIONS):
        raise ValueError('Incomplete numerical input manifest')
    for section,digests in expected.items():
        if (case/section).is_symlink() or hashes(case/section)!=digests:
            raise ValueError(f'Prepared input changed: {case/section}')
    if assets(case)!=manifest['assets_sha256']:
        raise ValueError('Prepared model assets changed')
    baseline_log=Path(manifest['baseline_log'])
    if sha256(baseline_log)!=manifest['baseline_log_sha256']:
        raise ValueError('Completed baseline log changed')
    summary.load_log(baseline_log)
    return manifest


def run_case(root, name, build, manifest, runner=subprocess.run):
    if not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('No active allocation; no solver launched')
    case=root/name
    env=os.environ.copy()
    env['LD_LIBRARY_PATH']=str(build/'lib')+(':'+env['LD_LIBRARY_PATH'] if env.get('LD_LIBRARY_PATH') else '')
    print(f'Launching {name} only: {case}/log.conservation',flush=True)
    with (case/'log.conservation').open('x') as log:
        result=runner(['srun','-n','52',str(build/'bin/putidaFoamSC1AlphaAudit'),
                       '-parallel','-fileHandler','collated'],cwd=case,env=env,
                      stdout=log,stderr=subprocess.STDOUT,check=False)
    print(f'SC1 {name} solver status: {result.returncode}',flush=True)
    report_status=0
    if result.returncode==0:
        report_status=summary.main([manifest['baseline_log'],'--normal',str(case/'log.conservation'),
                                    '--native-times'])
        print(f'Above reader label "normal" means numerical variant {name}; entrainment is ZERO.')
    tail=deque(maxlen=30)
    with (case/'log.conservation').open(errors='replace') as stream:
        for line in stream:
            tail.append(line)
            if (line.startswith('SC1_') and
                (not line.startswith('SC1_BOUNDARY ') or 'patch=outlet ' in line)) or line.startswith('BUBBLE_MASS_AUDIT '):
                print(line,end='')
    print(''.join(tail),end='')
    return result.returncode or report_status


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    prep=commands.add_parser('prepare',help='Create fresh cases; never run CFD')
    prep.add_argument('--baseline-pair',type=Path,default=Path('/scratch/nwimer/sc1_alpha_pair_20261110T155635Z'))
    prep.add_argument('--out',type=Path,required=True)
    run=commands.add_parser('run',help='Check a prepared case; launch only with --run')
    run.add_argument('--root',type=Path,required=True)
    run.add_argument('--case',choices=VARIANTS,required=True)
    run.add_argument('--audit-root',type=Path,default=Path('/scratch/nwimer/sc1_alpha_audit_20261107T001119Z'))
    run.add_argument('--run',action='store_true')
    args=parser.parse_args(argv)
    try:
        if args.command=='prepare':
            out=scratch(args.out)
            manifest=prepare(scratch(args.baseline_pair),out)
            print(f'Prepared {out}/four_outer and {out}/half_dt; no CFD launched.')
            print('Both retain the common original checkpoint, zero entrainment, all physics and biology settings.')
            print('Run four_outer first. Keep half_dt unrun until that result is reviewed.')
        else:
            root,build=scratch(args.root),scratch(args.audit_root)
            check_build(build)
            manifest=preflight(root,args.case)
            print('Numerical case, original checkpoint and baseline log verified.')
            print('Changed entries: '+json.dumps(manifest['cases'][args.case]['changed_entries']))
            if args.run:
                return run_case(root,args.case,build,manifest)
            print('Read-only preflight PASS; no solver launched.')
        return 0
    except (OSError,ValueError,KeyError,TypeError,subprocess.SubprocessError) as error:
        print(f'SC1 numerical check stopped: {error}\nPreserve all existing/partial cases; do not overwrite them.')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
