#!/usr/bin/env python3
"""Clone completed SC1 experiments at one explicit, identical checkpoint.

No CFD launch. Completed inputs are read-only. New copies retain the existing
physics and biological schedule. Only diagnostic output controls are changed.
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

CHECKPOINT = '31.199982372783058'


def hashes(root):
    """Stream hashes, including materialized symlink targets, without mutation."""
    result = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink() and path.is_dir():
            raise ValueError(f'Directory symlink needs explicit provenance review: {path}')
        if path.is_file():
            digest = hashlib.sha256()
            with path.open('rb') as stream:
                for block in iter(lambda: stream.read(4*1024*1024), b''):
                    digest.update(block)
            result[str(path.relative_to(root))] = digest.hexdigest()
    if not result:
        raise ValueError(f'Empty checkpoint/mesh: {root}')
    return result


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


def validate_pair(control, active, checkpoint):
    evidence = {}
    for section in ('processors52/constant', f'processors52/{checkpoint}'):
        a,b=hashes(control/section),hashes(active/section)
        if a != b:
            raise ValueError(f'Inputs are not identical: {section}')
        evidence[section]=a
    for field in ('alpha.liquid','alpha.bubbleGas','alpha.headspaceGas',
                  'O2.liquid','O2.bubbleGas','O2.headspaceGas','f.bubbleGas'):
        if field not in evidence[f'processors52/{checkpoint}']:
            raise ValueError(f'Missing required checkpoint field: {field}')
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
    return evidence


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    root=Path('/scratch/nwimer/sc1_headspace_entrainment_20261104')
    parser.add_argument('--control',type=Path,default=root/'run_control_05')
    parser.add_argument('--active',type=Path,default=root/'run_normalactive_05')
    parser.add_argument('--checkpoint',default=CHECKPOINT)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    try:
        control,active,out=map(scratch,(args.control,args.active,args.out))
        if out.exists() or out.is_symlink():
            raise ValueError(f'Destination exists: {out}')
        for src in (control,active):
            if src==out or src in out.parents or out in src.parents:
                raise ValueError('Input/output paths overlap')
        # Explicit time names only; prohibit traversal and latestTime selection.
        if not re.fullmatch(r'\d+(?:\.\d+)?',args.checkpoint):
            raise ValueError('Checkpoint must be an explicit numeric directory')
        evidence=validate_pair(control,active,args.checkpoint)
        contents=[]
        for src in (control,active):
            c=(src/'system/controlDict').read_text()
            for key,value in (('startFrom','latestTime'),('endTime',f'{float(args.checkpoint)+0.5:.14g}'),
                              ('writeControl','runTime'),('writeInterval','0.1')):
                c=set_control_entry(c,key,value,insert_if_missing=(key=='writeControl'))
            contents.append((c,enable_audit((src/'system/fvSolution').read_text())))
        # Don't equate trials with different hidden runtime control settings.
        if contents[0] != contents[1]:
            raise ValueError('Patched controlDict/fvSolution differ; review before running')
        out.mkdir(parents=True,exist_ok=False)
        for name,src,(c,f) in zip(('control','normal'),(control,active),contents):
            dest=out/name
            dest.mkdir()
            for section in ('constant','system','processors52/constant',f'processors52/{args.checkpoint}'):
                shutil.copytree(src/section,dest/section)
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
                  'initial_state_sha256':evidence,'duration_s':0.5,
                  'physics':'unchanged from completed experiments','run_order':['control','normal']}
        (out/'comparison_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(f'Prepared matched, scratch-only cases: {out}/control and {out}/normal')
        print('Initial collated state and mesh hashes match. No CFD launched.')
        print('Run control first; review conservation and regression before normal.')
        return 0
    except (OSError,ValueError) as e:
        print(f'SC1 pair preparation stopped: {e}')
        return 1


if __name__=='__main__':
    raise SystemExit(main())
