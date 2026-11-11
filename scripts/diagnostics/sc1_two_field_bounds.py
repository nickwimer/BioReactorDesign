#!/usr/bin/env python3
"""Prepare/run separate near-bound pairs using an already verified baseline build."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import sc1_two_field as base
from sc1_two_field_cases import CASES as BASELINE_CASES
from sc1_two_field_bounds_cases import CASES, LEVELS, VERSION, case_files
from sc1_two_field_bounds_report import summarize


def recipe_hashes(name):
    return {p:hashlib.sha256(text.encode()).hexdigest() for p,text in case_files(name).items()}


def baseline_evidence(root):
    """Read-only prerequisite: all four original tests and their provenance pass."""
    root=root.resolve(strict=True)
    manifest=base.verify_cases(root)
    logs={}
    for name in BASELINE_CASES:
        if not (root/name/'result.json').is_file():
            raise ValueError('Missing saved baseline result: '+name)
        logs[name]=base.completed_report(root,name,manifest)['log_sha256']
    return manifest,logs


def prepare(baseline_root, out):
    baseline_root=baseline_root.resolve(strict=True)
    baseline,logs=baseline_evidence(baseline_root)
    build=Path(baseline['build'])
    out=base.fresh(out,(baseline_root,build,base.REPO))
    recipes={name:case_files(name) for name in CASES}
    for level in LEVELS:
        a,b=(recipes[f'bound_{level}_{ref}'] for ref in ('no_reference','gas_reference'))
        if [p for p in a if a[p]!=b[p]]!=['constant/phaseProperties'] or \
                b['constant/phaseProperties'].replace('referencePhase gas;\n','')!=a['constant/phaseProperties']:
            raise ValueError('Uncontrolled near-bound reference comparison')
    out.mkdir(parents=True,exist_ok=False)
    for name,files in recipes.items():
        for rel,text in files.items():
            path=out/name/rel; path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('x') as f: f.write(text)
    manifest=dict(version=VERSION,suite='sc1_two_field_bounds',root=str(out),
        build=baseline['build'],assets=baseline['assets'],baseline_root=str(baseline_root),
        baseline_log_hashes=logs,cases={name:recipe_hashes(name) for name in CASES},
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip())
    base.write_json(out/'cases_manifest.json',manifest)
    print('NEARBOUND_CASES='+str(out)+' (four fresh cases; no CFD launched)')
    return out


def verify(root):
    root=root.resolve(strict=True)
    m=json.loads((root/'cases_manifest.json').read_text())
    if m['version']!=VERSION or m['suite']!='sc1_two_field_bounds' or m['root']!=str(root):
        raise ValueError('Wrong near-bound manifest')
    if m['cases']!={name:recipe_hashes(name) for name in CASES}:
        raise ValueError('Near-bound recipes changed since preparation')
    baseline,logs=baseline_evidence(Path(m['baseline_root']))
    if (m['build'],m['assets'],m['baseline_log_hashes'])!=(baseline['build'],baseline['assets'],logs):
        raise ValueError('Near-bound build or baseline evidence changed')
    return m


def saved_report(root, name, manifest):
    if not (root/name/'result.json').is_file():
        raise ValueError('A completed saved near-bound result is required: '+name)
    # A failed numerical screen is a valid input to timestep sensitivity.
    return base.completed_report(root,name,manifest,summarize_log=summarize,require_pass=False)


def run(root, level):
    root=root.resolve(strict=True); m=verify(root)
    if level not in LEVELS: raise ValueError('Unknown timestep level')
    if level=='fine':
        for ref in ('no_reference','gas_reference'):
            report=saved_report(root,'bound_coarse_'+ref,m)
            if not report['diagnostic_pass']:
                raise ValueError('Review coarse input/regime/accounting errors before fine runs: '+ref)
    names=['bound_'+level+'_'+ref for ref in ('no_reference','gas_reference')]
    for name in names: base.require_fresh_case(root,name,m)
    base.run_fresh_cases(root,m,names,summarize_log=summarize)
    print('NEAR_BOUND_LEVEL='+level+' complete. Review both numerical screens; failures are findings.')


def compare(root):
    root=root.resolve(strict=True); m=verify(root)
    reports={name:saved_report(root,name,m) for name in CASES}
    comparisons={}
    for ref in ('no_reference','gas_reference'):
        coarse=reports['bound_coarse_'+ref]; fine=reports['bound_fine_'+ref]
        comparisons[ref]={}
        for phase in ('gas','liquid'):
            a,b=coarse['phases'][phase],fine['phases'][phase]
            comparisons[ref][phase]=dict(
                coarseResidualL1_kg=a['residualL1_kg'],fineResidualL1_kg=b['residualL1_kg'],
                fineOverCoarseResidualL1=b['residualL1_kg']/a['residualL1_kg'] if a['residualL1_kg'] else None,
                fineMinusCoarseFinalMass_kg=b['finalMass_kg']-a['finalMass_kg'],
                fineMinusCoarseBoundaryOut_kg=b['boundaryOut_kg']-a['boundaryOut_kg'],
                coarseFinalOuterSplitVolume_m3=a['finalOuterSplitVolume_m3'],
                fineFinalOuterSplitVolume_m3=b['finalOuterSplitVolume_m3'],
                coarseFinalOuterRescaleVolume_m3=a['finalOuterRescaleVolume_m3'],
                fineFinalOuterRescaleVolume_m3=b['finalOuterRescaleVolume_m3'])
    return dict(suite='sc1_two_field_bounds',complete=True,comparisons=comparisons,
        cases={name:dict(screen_pass=r['screen_pass'],errors=r['errors'],log_sha256=r['log_sha256'])
               for name,r in reports.items()},
        meaning='Two timestep levels show sensitivity, not a demonstrated convergence order. '
                'Final-outer terms are iteration diagnostics, not extra physical mass.')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest='action',required=True)
    p=subs.add_parser('prepare'); p.add_argument('--baseline-root',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    p=subs.add_parser('run'); p.add_argument('--root',type=Path,required=True); p.add_argument('--level',choices=LEVELS,required=True); p.add_argument('--run',action='store_true',required=True)
    p=subs.add_parser('analyze'); p.add_argument('--log',type=Path,required=True); p.add_argument('--case',choices=CASES,required=True)
    p=subs.add_parser('compare'); p.add_argument('--root',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.action=='run' and not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Use an active Kestrel allocation for run')
    if args.action=='prepare': prepare(args.baseline_root,args.out)
    elif args.action=='run': run(args.root,args.level)
    elif args.action=='analyze': print(json.dumps(summarize(args.log.read_text(),args.case),indent=2))
    else: print(json.dumps(compare(args.root),indent=2))


if __name__=='__main__':
    try: main()
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        print('SC1 near-bound check stopped: '+str(error),file=sys.stderr)
        raise SystemExit(1)
