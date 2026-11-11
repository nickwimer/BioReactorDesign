#!/usr/bin/env python3
"""Read a saved startup log and bundle evidence; never run CFD or modify the case.

This is deliberately independent of OpenFOAM, numpy and the success-only startup
analyzer. An incomplete/invalid trajectory is useful evidence, not a passed run.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import tarfile

VARIANTS=('sc1_meml_quiescent_startup_v1','sc1_meml_quiescent_segregated_drag_v2',
          'sc1_meml_quiescent_face_momentum_v3','sc1_meml_quiescent_coupled_drag_v4',
          'sc1_meml_pressure_probe_v1','sc1_meml_pressure_coupling_v5')
REPO=Path(__file__).resolve().parents[2]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def regular(root, relative, *, missing_ok=False):
    """Reject path escapes and symlinks, including intermediate components."""
    rel=Path(relative)
    if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe evidence path: '+str(rel))
    p=root
    for part in rel.parts:
        p=p/part
        if p.is_symlink():raise ValueError('Redirected evidence: '+str(p))
    if not p.is_file():
        if missing_ok and not p.exists():return None
        raise ValueError('Missing evidence: '+str(p))
    return p


def number(value):
    try:return float(value)
    except ValueError:return None


def scan(path):
    """Keep first event of each kind plus a bounded per-time chronology."""
    first={}; counts={}; time=None; outer=None; times=[]; paired={}; duplicates=0
    last_step={}; courant_max=None; largest_alpha=0.; lines=0; end=False
    malformed=0; marker_count=0; timeline=[]; current=None

    def event(kind,line,text,**extra):
        counts[kind]=counts.get(kind,0)+1
        first.setdefault(kind,dict(line=line,lastPrintedTime_s=time,outer=outer,text=text,**extra))

    with path.open(errors='replace') as f:
        for lines,raw in enumerate(f,1):
            text=raw.strip()
            match=re.fullmatch(r'Time = (\S+)',text)
            if match:
                time=number(match[1]);outer=None
                if time is None or not math.isfinite(time):
                    time=None;event('invalid_time',lines,text)
                times.append(time)
                current=dict(time_s=time,line=lines,alphaMin=None,alphaMax=None,
                             pressureIterationsMax=0,pressureLimitCalls=0)
                if len(timeline)<1000:timeline.append(current)
            match=re.fullmatch(r'PIMPLE: Iteration (\d+)',text)
            if match:outer=int(match[1])
            if text=='End':end=True
            if text=='Selecting drag blending method: none':marker_count+=1
            match=re.fullmatch(r'Courant Number mean: (\S+) max: (\S+)',text)
            if match:
                value=number(match[2])
                if value is None or not math.isfinite(value):event('nonfinite_courant',lines,text)
                else:
                    courant_max=value if courant_max is None else max(courant_max,value)
                    # OF prints Co before incrementing Time. Do not label it as the next step.
                    if value>.2:event('courant_above_0_2',lines,text,maximum=value)
            match=re.fullmatch(r'(gas|liquid) fraction, min, max = (\S+) (\S+) (\S+)',text)
            if match:
                values=[number(x) for x in match.groups()[1:]]
                if any(x is None or not math.isfinite(x) for x in values):
                    event('nonfinite_alpha',lines,text)
                else:
                    _,low,high=values;largest_alpha=max(largest_alpha,abs(low),abs(high))
                    if current is not None:
                        current['alphaMin']=low if current['alphaMin'] is None else min(low,current['alphaMin'])
                        current['alphaMax']=high if current['alphaMax'] is None else max(high,current['alphaMax'])
                    if low < -1e-10 or high > 1+1e-10:event('alpha_out_of_bounds',lines,text)
                    if max(abs(low),abs(high))>10:event('alpha_magnitude_above_10',lines,text)
            if 'Solving for p_rgh,' in text:
                match=re.search(r'No Iterations (\d+)',text)
                if match and current is not None:
                    current['pressureIterationsMax']=max(current['pressureIterationsMax'],int(match[1]))
                if re.search(r'\b(?:nan|inf)\b',text,re.I):event('nonfinite_pressure_solve',lines,text)
            if text.startswith('limitPressure:'):
                event('pressure_limit',lines,text)
                if current is not None:current['pressureLimitCalls']+=1
            if text.startswith('SC1_'):
                fields=dict(re.findall(r'(\w+)=([^\s]+)',text))
                if re.search(r'=(?:[+-]?(?:nan|inf))(?=\s|$)',text,re.I):
                    event('nonfinite_audit',lines,text)
                if text.startswith('SC1_PHASE_TERM ') and fields.get('term')=='clamp':
                    value=number(fields.get('volume_m3',''))
                    if value is not None and math.isfinite(value) and abs(value)>1e-10:
                        event('clamp_above_1e_10_m3',lines,text)
                if text.startswith('SC1_TWO_STEP '):
                    phase=fields.get('phase'); t=number(fields.get('t_s',''))
                    if phase not in ('gas','liquid') or t is None or not math.isfinite(t):malformed+=1
                    else:
                        # Printed SC1 times can differ from Time by ~1e-16 s.
                        key=round(t,12); seen=paired.setdefault(key,set())
                        if phase in seen:duplicates+=1
                        seen.add(phase);last_step[phase]=dict(line=lines,record=fields)
            lower=text.lower()
            trapping_banner=re.match(r'^sigfpe\s*:\s*enabling floating point exception trapping\b',lower)
            if not trapping_banner and ('floating point exception' in lower or 'sigfpe::sighandler' in lower):
                event('floating_point_exception',lines,text)
            if 'foam fatal' in lower:event('foam_fatal',lines,text)
            if 'due to time limit' in lower or 'exceeded its time limit' in lower:
                event('allocation_time_limit',lines,text)
    if not times and not first:raise ValueError('No recognizable startup runtime evidence')
    complete_pairs=sorted(t for t,p in paired.items() if p=={'gas','liquid'})
    # This is record coverage only: bounded accepted fields can hide catastrophic subcycles.
    return dict(lines=lines,endMarkerSeen=end,printedTimes=len(times),lastPrintedTime_s=time,
        pairedStepRecordCount=len(complete_pairs),lastPairedStepRecordTime_s=max(complete_pairs,default=None),
        duplicateStepRecords=duplicates,malformedStepRecords=malformed,lastStepRecords=last_step,
        dragNoneSelectionCount=marker_count,maxFiniteCourant=courant_max,
        largestFiniteAlphaMagnitude=largest_alpha,firstEvents=first,eventCounts=counts,
        timeline=timeline,timelineTruncated=len(times)>len(timeline),
        meaning='Diagnostic evidence only. End/paired records do not establish completion, conservation or startup acceptance.',
        timeContext='lastPrintedTime_s is the preceding Time line; Courant output occurs before the next Time increment.')


def context(path, events, out):
    selected=set()
    for e in events.values():selected.update(range(max(1,e['line']-8),e['line']+9))
    previous=0
    with path.open(errors='replace') as source,out.open('x') as dest:
        for n,line in enumerate(source,1):
            if n in selected:
                if n!=previous+1:dest.write('\n...\n')
                dest.write(f'{n}: {line}');previous=n


def collect(root,out):
    root=root.resolve(strict=True)
    manifest_path=regular(root,'startup_manifest.json')
    manifest_hash=sha(manifest_path)
    m=json.loads(manifest_path.read_text())
    if m.get('variant') not in VARIANTS or m.get('root')!=str(root):
        raise ValueError('Wrong startup manifest/root')
    log=regular(root,'quiescent/log.solver'); attempt=regular(root,'quiescent/attempt.json')
    attempt_hash=sha(attempt)
    expected_attempt={k:m[k] for k in ('build','assets','ranks')}
    if json.loads(attempt.read_text())!=expected_attempt:raise ValueError('Attempt/build identity mismatch')
    if out.exists() or out.is_symlink():raise ValueError('Destination exists; preserve it')
    out=out.resolve()
    protected=[root,REPO]+[Path(m[k]).resolve() for k in ('archive','build','regression_root') if k in m]
    if 'baseline' in m:protected.append(Path(m['baseline']['root']).resolve())
    for p in protected:
        if out==p or out in p.parents or p in out.parents:raise ValueError('Output overlaps an input')
    selected={'startup_manifest.json':manifest_path,'quiescent/attempt.json':attempt,'quiescent/log.solver':log}
    changed=[];missing=[];checked=0
    for rel,digest in m['inputs'].items():
        if not Path(rel).parts or Path(rel).parts[0] not in ('0','constant','system'):
            raise ValueError('Unexpected manifest input')
        p=regular(root,'quiescent/'+rel,missing_ok=True)
        if p is None:
            missing.append(rel);continue
        checked+=1
        if sha(p)!=digest:changed.append(rel)
        if rel.startswith('system/') or (rel.startswith('constant/') and not rel.startswith('constant/polyMesh/')):
            selected['quiescent/'+rel]=p
    extras=[]
    for folder in ('0','constant','system'):
        tree=root/'quiescent'/folder
        if tree.is_symlink():raise ValueError('Redirected input tree')
        for p in tree.rglob('*'):
            if p.is_symlink():raise ValueError('Redirected input: '+str(p))
            if p.is_file() and p.relative_to(root/'quiescent').as_posix() not in m['inputs']:
                extras.append(p.relative_to(root/'quiescent').as_posix())
    for name in ('log.checkMesh','log.decomposePar','result.json','guard_stop.json'):
        rel='quiescent/'+name
        if (root/rel).exists() or (root/rel).is_symlink():selected[rel]=regular(root,rel)
    hashes={rel:sha(p) for rel,p in selected.items()}
    if hashes['startup_manifest.json']!=manifest_hash or hashes['quiescent/attempt.json']!=attempt_hash:
        raise ValueError('Manifest/attempt changed during inspection')
    report=scan(log)
    folders={}
    case=root/'quiescent'
    for p in [case]+sorted(case.glob('processor[0-9]*')):
        if p.is_symlink() or not p.is_dir():raise ValueError('Redirected/missing processor directory')
        folders[p.relative_to(root).as_posix()]=sorted(q.name for q in p.iterdir()
            if q.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',q.name))
    report.update(root=str(root),variant=m['variant'],sourceCommit=m.get('source_commit'),
        solverLogSha256=hashes['quiescent/log.solver'],inputIdentity=dict(checked=checked,changed=changed,missing=missing,extra=extras),
        timeDirectories=folders,fieldAvailabilityMeaning='Directory names only; no assertion of complete written/reconstructable fields.',
        provenance='Local manifest/input/log inspection; build binaries and external regression/archive roots are not revalidated.')
    out.mkdir(parents=True,exist_ok=False)
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (out/'evidence_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    context(log,report['firstEvents'],out/'first_events.txt')
    bundle=out/'reactor_failure_evidence.tar.gz'
    with tarfile.open(bundle,'w:gz') as tar:
        for rel,p in selected.items():tar.add(p,arcname='case_evidence/'+rel,recursive=False)
        for name in ('report.json','evidence_sha256.json','first_events.txt'):
            tar.add(out/name,arcname=name,recursive=False)
    if {rel:sha(p) for rel,p in selected.items()}!=hashes:
        raise ValueError('Evidence changed during collection; bundle is not a stable snapshot; preserve output')
    compact={k:v for k,v in report.items() if k not in ('timeline','lastStepRecords','timeDirectories')}
    print(json.dumps(compact,indent=2,allow_nan=False))
    print('REPORT='+str(out/'report.json'))
    print('BUNDLE='+str(bundle))
    print('BUNDLE_SHA256='+sha(bundle))
    print('Evidence collected; no CFD launched and no case files changed.')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
    a=p.parse_args(argv);collect(a.root,a.out);return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,KeyError,TypeError,tarfile.TarError) as error:
        print('SC1 failure collection stopped: '+str(error),file=sys.stderr)
        raise SystemExit(1)
