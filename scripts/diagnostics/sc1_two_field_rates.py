#!/usr/bin/env python3
"""Build and regress a private rate-split solver against completed controls."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import sc1_two_field as base
import sc1_two_field_bounds as bounds
from sc1_two_field_cases import CASES as BASELINE_CASES, case_files as baseline_files
from sc1_two_field_bounds_cases import CASES as BOUND_CASES, case_files as bound_files
from sc1_two_field_rates_instrument import APP, TAG, RECIPE
from sc1_two_field_rates_report import summarize

CASES=BASELINE_CASES+BOUND_CASES
PAIRS=(CASES[:2],CASES[2:4],CASES[4:6],CASES[6:])


def case_files(name):
    files=baseline_files(name) if name in BASELINE_CASES else bound_files(name)
    text=files['system/controlDict']
    anchor='application '+base.APP+';'
    if text.count(anchor)!=1:raise ValueError('Unexpected control application')
    files['system/controlDict']=text.replace(anchor,'application '+APP+';')
    return files


def recipes():
    return {n:{p:hashlib.sha256(t.encode()).hexdigest() for p,t in case_files(n).items()} for n in CASES}


def controls(root):
    root=root.resolve(strict=True); m=bounds.verify(root)
    baseline_root=Path(m['baseline_root']); baseline=base.verify_cases(baseline_root)
    reports={n:base.completed_report(baseline_root,n,baseline) for n in BASELINE_CASES}
    for n in BOUND_CASES:
        reports[n]=bounds.saved_report(root,n,m)
        if not reports[n]['diagnostic_pass']:
            raise ValueError('Existing near-bound controls have invalid diagnostics: '+n)
    return reports


def prepare(build_root, control_root, out):
    build_root=build_root.resolve(strict=True);control_root=control_root.resolve(strict=True)
    checked=base.verify_build(build_root,recipe=RECIPE)
    original=controls(control_root)
    out=base.fresh(out,(build_root,control_root,base.REPO))
    out.mkdir(parents=True,exist_ok=False)
    for n in CASES:
        for rel,text in case_files(n).items():
            p=out/n/rel;p.parent.mkdir(parents=True,exist_ok=True)
            with p.open('x') as f:f.write(text)
    base.write_json(out/'cases_manifest.json',dict(version=1,variant=TAG,root=str(out),
        build=str(build_root),assets=checked['assets'],control_root=str(control_root),
        control_logs={n:r['log_sha256'] for n,r in original.items()},cases=recipes(),
        source_commit=subprocess.check_output(['git','-C',str(base.REPO),'rev-parse','HEAD'],text=True).strip()))
    print('RATE_CASES='+str(out)+' (eight fresh cases; no CFD launched)')
    return out


def verify(root):
    root=root.resolve(strict=True);m=json.loads((root/'cases_manifest.json').read_text())
    if m['version']!=1 or m['variant']!=TAG or m['root']!=str(root) or m['cases']!=recipes():
        raise ValueError('Invalid/changed rate-split case manifest')
    if base.verify_build(Path(m['build']),recipe=RECIPE)['assets']!=m['assets']:
        raise ValueError('Rate cases belong to a different private build')
    original=controls(Path(m['control_root']))
    if {n:r['log_sha256'] for n,r in original.items()}!=m['control_logs']:
        raise ValueError('Completed control logs changed')
    return m,original


def saved(root,name,m):
    if not (root/name/'result.json').is_file():raise ValueError('Missing saved rate result: '+name)
    return base.completed_report(root,name,m,summarize_log=summarize,require_pass=False)


def comparison(root,m,original):
    result={}
    for n in CASES:
        r=saved(root,n,m)
        result[n]=dict(screen_pass=r['screen_pass'],diagnostic_pass=r['diagnostic_pass'],
            errors=r['errors'],log_sha256=r['log_sha256'],control_log_sha256=original[n]['log_sha256'],
            phases={p:dict(controlResidualL1_kg=original[n]['phases'][p]['residualL1_kg'],
                rateResidualL1_kg=r['phases'][p]['residualL1_kg'],
                screenTolerance_kg=r['phases'][p]['screenTolerance_kg'],
                finalMass_kg=r['phases'][p]['finalMass_kg'],boundaryOut_kg=r['phases'][p]['boundaryOut_kg'],
                allOuterSplitL1Volume_m3=r['phases'][p]['allOuterSplitL1Volume_m3'])
                for p in ('gas','liquid')})
    return dict(suite=TAG,complete=True,regression_pass=all(r['screen_pass'] for r in result.values()),
        cases=result,meaning='Numerical regression against unchanged thresholds. Reactor/interface validation remains required.')


def run_all(root):
    root=root.resolve(strict=True);m,original=verify(root)
    # Check every destination before claiming the first case.
    for name in CASES:base.require_fresh_case(root,name,m)
    for index,pair in enumerate(PAIRS):
        base.run_fresh_cases(root,m,pair,summarize_log=summarize,app=APP)
        reports=[saved(root,n,m) for n in pair]
        # Closed/mid-fraction regressions must pass. Near-bound mass failures
        # can be studied at half-dt only when the new diagnostics are valid.
        gate='screen_pass' if index<2 else 'diagnostic_pass'
        if not all(r[gate] for r in reports):
            raise ValueError('Regression stage failed '+gate+'; preserve '+str(root))
    result=comparison(root,m,original)
    print(json.dumps(result,indent=2),flush=True)
    print('RATE_REGRESSION='+('PASS' if result['regression_pass'] else 'FAIL'),flush=True)
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);subs=p.add_subparsers(dest='action',required=True)
    s=subs.add_parser('build');s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('prepare');s.add_argument('--build-root',type=Path,required=True);s.add_argument('--control-root',type=Path,required=True);s.add_argument('--out',type=Path,required=True)
    s=subs.add_parser('run-all');s.add_argument('--root',type=Path,required=True);s.add_argument('--run',action='store_true',required=True)
    s=subs.add_parser('compare');s.add_argument('--root',type=Path,required=True)
    s=subs.add_parser('analyze');s.add_argument('--case',choices=CASES,required=True);s.add_argument('--log',type=Path,required=True)
    a=p.parse_args(argv)
    if a.action in ('build','run-all') and not os.environ.get('SLURM_JOB_ID'):
        raise ValueError('Use an active Kestrel allocation')
    if a.action=='build':
        base.build(a.out,recipe=RECIPE)
        print('RATE_BUILD='+str(a.out.resolve()))
    elif a.action=='prepare':prepare(a.build_root,a.control_root,a.out)
    elif a.action=='run-all':return 0 if run_all(a.root)['regression_pass'] else 2
    elif a.action=='analyze':print(json.dumps(summarize(a.log.read_text(),a.case),indent=2))
    else:
        root=a.root.resolve(strict=True);m,original=verify(root)
        print(json.dumps(comparison(root,m,original),indent=2))
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,subprocess.SubprocessError) as error:
        print('SC1 rate-split experiment stopped: '+str(error),file=sys.stderr)
        raise SystemExit(1)
