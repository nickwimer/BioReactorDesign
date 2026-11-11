#!/usr/bin/env python3
"""Clone completed SC1 experiments at one explicit, identical checkpoint.

No CFD launch. Completed inputs are read-only. New copies retain the existing
physics and biological schedule. Only diagnostic output controls are changed.
Use --inspect-inputs without --out to report raw differences without copying.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

from prepare_sc1_conservation_solver import scratch
from prepare_sc1_headspace_entrainment_trial import enclosing_brace_end
from prepare_sc1_headspace_o2_trial import set_control_entry
from sc1_checkpoint_metadata import COURANT_FILES, verify_courant_location_metadata

CHECKPOINT = '31.199982372783058'


def hashes(root):
    """Stream hashes, including materialized symlink targets, without mutation."""
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink() and path.is_dir():
            raise ValueError(f'Directory symlink needs explicit provenance review: {path}')
        if path.is_file():
            before=path.stat()
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(4*1024*1024), b''):
                    digest.update(block)
            after=path.stat()
            if (before.st_size,before.st_mtime_ns,before.st_ino,before.st_dev) != (after.st_size,after.st_mtime_ns,after.st_ino,after.st_dev):
                raise ValueError(f'File changed while hashing: {path}')
            result[str(path.relative_to(root))] = digest.hexdigest()
    if not result:
        raise ValueError(f'Empty checkpoint/mesh: {root}')
    return result


def file_summary(path, digest):
    return {'size_bytes':path.stat().st_size, 'sha256':digest,
            'resolved_path':str(path.resolve()), 'symlink':path.is_symlink(),
            'compressed_suffix':path.suffix == '.gz'}


def first_difference(left, right):
    """Bounded raw-byte evidence, not a decoder or a header-only classifier."""
    offset=0
    with left.open('rb') as a, right.open('rb') as b:
        while True:
            x,y=a.read(1024*1024),b.read(1024*1024)
            if x != y:
                offset += next((i for i,(u,v) in enumerate(zip(x,y)) if u!=v), min(len(x),len(y)))
                break
            if not x:
                return None
            offset += len(x)
    result={'offset':offset}
    for label,path in (('control',left),('normal',right)):
        start=max(0,offset-32)
        with path.open('rb') as stream:
            stream.seek(start)
            context=stream.read(96)
        result[label]={'offset':start, 'hex':context.hex(),
                       'escaped_ascii':context.decode('ascii',errors='backslashreplace')}
    return result


def inspect_pair_inputs(control, active, checkpoint):
    """Read-only report; differences never authorize bypassing validation."""
    report={'inspection_only':True, 'checkpoint':checkpoint,
            'control_source':str(control), 'normal_source':str(active),
            'note':'Raw byte comparison only. An early header difference does not exclude later field differences. No files are normalized or copied.',
            'sections':{}, 'read_errors':0}
    for section in ('processors52/constant',f'processors52/{checkpoint}','constant','system'):
        try:
            left,right=control/section,active/section
            a,b=hashes(left),hashes(right)
            common=set(a)&set(b)
            changed=[]
            for name in sorted(n for n in common if a[n]!=b[n]):
                changed.append({'file':name,
                                'control':file_summary(left/name,a[name]),
                                'normal':file_summary(right/name,b[name]),
                                'first_difference':first_difference(left/name,right/name),
                                'classification':'undetermined'})
            report['sections'][section]={
                'files_control':len(a), 'files_normal':len(b),
                'identical_files':sum(a[n]==b[n] for n in common),
                'changed':changed,
                'only_control':[dict(file=n,**file_summary(left/n,a[n])) for n in sorted(set(a)-set(b))],
                'only_normal':[dict(file=n,**file_summary(right/n,b[n])) for n in sorted(set(b)-set(a))]}
        except (OSError,ValueError) as error:
            report['sections'][section]={'error':str(error)}
            report['read_errors']+=1
    return report


def enable_audit(text):
    # Strip comments before locating anchors; preserve character offsets.
    masked = re.sub(r'/\*.*?\*/|//[^\n]*', lambda m: ' '*len(m[0]), text, flags=re.S)
    matches = list(re.finditer(r'\bPIMPLE\s*\{', masked))
    if len(matches) != 1:
        raise ValueError('Expected one explicit PIMPLE dictionary')
    start=matches[0].end()-1
    end=enclosing_brace_end(text,start)
    body=text[start+1:end]
    if re.search(r'\bsc1ConservationAudit\b', masked[start+1:end]):
        raise ValueError('Conservation audit already configured')
    return text[:start+1]+body+'\n    sc1ConservationAudit true;\n'+text[end:]


def enable_alpha_audit(text):
    """Opt in only in the explicit alpha solver dictionary of a NEW case."""
    masked = re.sub(r'/\*.*?\*/|//[^\n]*', lambda m: ' '*len(m[0]), text, flags=re.S)
    token_re = re.compile(r'"(?:\\.|[^"\\])*"|[{};]|[^\s{};]+')

    def dictionaries(start, end):
        tokens = list(token_re.finditer(masked, start, end))
        found = []
        i = 0
        while i < len(tokens):
            token = tokens[i]
            if token[0].startswith(('#', '$')):
                raise ValueError('Alpha audit requires explicit dictionaries without includes or inheritance')
            if token[0] == '{':
                if i == 0:
                    raise ValueError('Unnamed dictionary in fvSolution')
                close = enclosing_brace_end(text, token.start())
                if close >= end:
                    raise ValueError('Unbalanced dictionary in fvSolution')
                found.append((tokens[i-1][0].strip('"'), token.start(), close))
                while i < len(tokens) and tokens[i].start() <= close:
                    i += 1
                continue
            if token[0] == '}':
                raise ValueError('Unexpected closing brace in fvSolution')
            i += 1
        return found

    solvers = [d for d in dictionaries(0, len(text)) if d[0] == 'solvers']
    if len(solvers) != 1:
        raise ValueError('Expected one explicit solvers dictionary')
    _, start, end = solvers[0]
    candidates = [d for d in dictionaries(start+1, end) if 'alpha' in d[0]]
    if len(candidates) != 1 or candidates[0][0] not in ('alpha', 'alpha.*'):
        raise ValueError('Expected one explicit alpha or "alpha.*" solver dictionary')
    _, start, end = candidates[0]
    if re.search(r'\bsc1AlphaAudit\b', masked[start+1:end]):
        raise ValueError('Alpha audit already configured')
    if re.search(r'[#\$]', masked[start+1:end]):
        raise ValueError('Alpha audit requires explicit alpha controls')
    return text[:end]+'\n        sc1AlphaAudit true;\n    '+text[end:]


def validate_pair(control, active, checkpoint, allow_courant_location_metadata=False,
                  metadata_proofs=None):
    evidence = {}
    for section in ('processors52/constant', f'processors52/{checkpoint}'):
        a,b=hashes(control/section),hashes(active/section)
        if a != b:
            if (not allow_courant_location_metadata or section!=f'processors52/{checkpoint}'
                    or set(a)!=set(b)):
                raise ValueError(f'Inputs are not identical: {section}; use --inspect-inputs (without --out) for a read-only file comparison')
            for name in sorted(n for n in a if a[n]!=b[n]):
                if name not in COURANT_FILES:
                    raise ValueError(f'Non-Courant checkpoint file differs: {name}')
                proof=verify_courant_location_metadata(control/section/name,active/section/name)
                if proof['control_sha256']!=a[name] or proof['normal_sha256']!=b[name]:
                    raise ValueError(f'Checkpoint changed between hashing and verification: {name}')
                if metadata_proofs is not None:
                    metadata_proofs[name]=proof
        evidence[section]=a
    for field in ('alpha.liquid','alpha.bubbleGas','alpha.headspaceGas',
                  'O2.liquid','O2.bubbleGas','O2.headspaceGas','f.bubbleGas'):
        if field not in evidence[f'processors52/{checkpoint}']:
            raise ValueError(f'Missing required checkpoint field: {field}')
    validate_settings(control, active)
    return evidence


def validate_settings(control, active):
    """Compare case settings without treating a completed checkpoint as pristine."""
    # Hold all constants, including microbeUpdateProperties, fixed except
    # the already established entrainment configuration in phaseProperties.
    a,b=hashes(control/'constant'),hashes(active/'constant')
    a.pop('phaseProperties',None); b.pop('phaseProperties',None)
    if a != b:
        raise ValueError('Constants/biology differ beyond phaseProperties')
    a,b=hashes(control/'system'),hashes(active/'system')
    for name in ('controlDict','fvSolution'):
        a.pop(name,None); b.pop(name,None)
    if a != b:
        raise ValueError('Numerical system files differ')
    # Ensure the controls relevant to the solve match exactly after enabling.
    if (control/'system/fvSolution').read_bytes() != (active/'system/fvSolution').read_bytes():
        raise ValueError('fvSolution differs; inspect before making a matched pair')
    normalized=[]
    for src,eff,criterion,downward in ((control,'0','legacyDownwardLiquid','true'),
                                       (active,'0.01','normalGasIntoLiquid','false')):
        p=(src/'constant/phaseProperties').read_text()
        p=re.sub(r'/\*.*?\*/|//[^\n]*','',p,flags=re.S)
        for key,value in (('entrainmentEfficiency',eff),('entrainmentVelocityCriterion',criterion),
                          ('requireDownwardLiquid',downward),('maxGlobalEntrainedKgS','0.002'),
                          ('nucleationDiameter','0.002')):
            matches=re.findall(r'\b'+key+r'\s+(\S+)\s*;',p)
            if key=='entrainmentVelocityCriterion' and not matches and criterion=='legacyDownwardLiquid':
                continue  # Original completed zero control predates this optional key.
            if matches != [value]:
                raise ValueError(f'Unexpected {key} in {src}: {matches}')
        foam=re.findall(r'\bfoamAreaMultiplier\s+(\S+)\s*;',p)
        if not foam or any(float(x)!=0 for x in foam):
            raise ValueError('Both cases must retain zero foam multiplier')
        for key in ('entrainmentEfficiency','entrainmentVelocityCriterion','requireDownwardLiquid'):
            p=re.sub(r'\b'+key+r'\s+\S+\s*;','',p)
        normalized.append(' '.join(p.split()))
    if normalized[0]!=normalized[1]:
        raise ValueError('phaseProperties differs beyond the three established entrainment switches')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    root=Path('/scratch/nwimer/sc1_headspace_entrainment_20261104')
    parser.add_argument('--control',type=Path,default=root/'run_control_05')
    parser.add_argument('--active',type=Path,default=root/'run_normalactive_05')
    parser.add_argument('--checkpoint',default=CHECKPOINT)
    parser.add_argument('--duration',type=float,choices=(0.1,0.5),default=0.5,
                        help='New diagnostic duration in CFD seconds; default preserves the completed 0.5 s comparison')
    parser.add_argument('--out',type=Path,help='Required for preparation; must not be used with --inspect-inputs')
    parser.add_argument('--inspect-inputs',action='store_true',help='Print raw input differences only; create no files or cases')
    parser.add_argument('--alpha-audit',action='store_true',
                        help='Enable private-library alpha subcycle diagnostics in new cases; requires putidaFoamSC1AlphaAudit')
    parser.add_argument('--allow-courant-location-metadata',action='store_true',
                        help='Require complete equality except proved Co outer location values; copy the full control checkpoint to both new cases')
    args=parser.parse_args()
    if args.inspect_inputs and args.out is not None:
        parser.error('--inspect-inputs must not be combined with --out')
    if args.inspect_inputs and args.allow_courant_location_metadata:
        parser.error('--inspect-inputs is a raw report; do not combine with --allow-courant-location-metadata')
    if args.inspect_inputs and args.duration != 0.5:
        parser.error('--duration only applies to new case preparation')
    if args.inspect_inputs and args.alpha_audit:
        parser.error('--alpha-audit only applies to new case preparation')
    if not args.inspect_inputs and args.out is None:
        parser.error('--out is required for case preparation')
    try:
        control,active=map(scratch,(args.control,args.active))
        # Explicit time names only; prohibit traversal and latestTime selection.
        if not re.fullmatch(r'\d+(?:\.\d+)?',args.checkpoint):
            raise ValueError('Checkpoint must be an explicit numeric directory')
        if args.inspect_inputs:
            report=inspect_pair_inputs(control,active,args.checkpoint)
            print(json.dumps(report,indent=2,ensure_ascii=True))
            return 2 if report['read_errors'] else 0
        out=scratch(args.out)
        if out.exists() or out.is_symlink():
            raise ValueError(f'Destination exists: {out}')
        for src in (control,active):
            if src==out or src in out.parents or out in src.parents:
                raise ValueError('Input/output paths overlap')
        metadata_proofs={}
        evidence=validate_pair(control,active,args.checkpoint,
                               args.allow_courant_location_metadata,metadata_proofs)
        contents=[]
        for src in (control,active):
            c=(src/'system/controlDict').read_text()
            for key,value in (('startFrom','latestTime'),('endTime',f'{float(args.checkpoint)+args.duration:.14g}'),
                              ('writeControl','runTime'),('writeInterval','0.1')):
                c=set_control_entry(c,key,value,insert_if_missing=(key=='writeControl'))
            f=enable_audit((src/'system/fvSolution').read_text())
            if args.alpha_audit:
                f=enable_alpha_audit(f)
            contents.append((c,f))
        # Don't equate trials with different hidden runtime control settings.
        if contents[0] != contents[1]:
            raise ValueError('Patched controlDict/fvSolution differ; review before running')
        out.mkdir(parents=True,exist_ok=False)
        for name,src,(c,f) in zip(('control','normal'),(control,active),contents):
            dest=out/name
            dest.mkdir()
            for section in ('constant','system','processors52/constant',f'processors52/{args.checkpoint}'):
                # The complete control state is copied unchanged into BOTH
                # new cases, only after the other source has been verified.
                origin=control if section.startswith('processors52/') else src
                shutil.copytree(origin/section,dest/section)
            # Copy runtime model assets; omit prior logs and derived outputs.
            for path in src.iterdir():
                if path.is_file() and (path.suffix in ('.pt','.txt','.onnx') or path.name.startswith('scaler_')):
                    shutil.copy2(path,dest/path.name)
            (dest/'system/controlDict').write_text(c)
            (dest/'system/fvSolution').write_text(f)
            for section,expected in evidence.items():
                if hashes(dest/section)!=expected:
                    raise ValueError(f'Copied state failed verification: {dest/section}')
        manifest={'checkpoint':args.checkpoint,'control_source':str(control),'normal_source':str(active),
                  'initial_state_source':str(control),
                  'courant_location_metadata_proofs':metadata_proofs,
                  'initial_state_sha256':evidence,'duration_s':args.duration,
                  'alpha_audit':args.alpha_audit,
                  'physics':'unchanged from completed experiments','run_order':['control','normal']}
        (out/'comparison_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(f'Prepared matched, scratch-only cases: {out}/control and {out}/normal')
        if metadata_proofs:
            print('Verified complete Co files differ only in outer location values: '+', '.join(sorted(metadata_proofs)))
        print(f'Both complete initial states copied unchanged from: {control}')
        print('Initial collated state and mesh hashes match. No CFD launched.')
        if args.alpha_audit:
            print('Alpha subcycle audit enabled: requires putidaFoamSC1AlphaAudit and its private libbirdphaseSystem.')
        print('Run control first; review conservation and regression before normal.')
        return 0
    except (OSError,ValueError) as e:
        print(f'SC1 pair preparation stopped: {e}')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
