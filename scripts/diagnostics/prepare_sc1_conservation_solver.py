#!/usr/bin/env python3
"""Clone and instrument the CURRENT scratch solver; never rebuild production.

Only prepares a new source tree and unique executable destination. Does not
compile, launch CFD, change a case or overwrite an existing destination.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / "OFsolvers/birdmultiphaseEulerFoam/putidaFoam"


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Expected one source anchor: {old!r}")
    return text.replace(old, new)


def instrument(main, species):
    if "sc1Conservation" in main or "sc1Conservation" in species:
        raise ValueError("Audit already/partially installed; use an unmodified source")
    main = replace_once(main, '#include "phaseSystem.H"',
                        '#include "phaseSystem.H"\n#include "sc1ConservationAudit.H"')
    main = replace_once(main, '    while (pimple.run(runTime))', '''    autoPtr<sc1ConservationAudit> sc1Conservation;
    if (pimple.dict().lookupOrDefault<Switch>("sc1ConservationAudit", false))
    {
        sc1Conservation.reset(new sc1ConservationAudit(fluid, LTS));
    }

    while (pimple.run(runTime))''')
    main = replace_once(main, '        runTime++;', '''        if (sc1Conservation.valid()) sc1Conservation->beforeAdvance();
        runTime++;''')
    pattern = re.compile(r"(?m)^([ \t]*)fluid\.solve\(rAUs, rAUfs\);\s*\n[ \t]*fluid\.correct\(\);")

    def pb_hooks(match):
        indent = match[1]
        guard = 'sc1Conservation.valid() && pimple.finalPimpleIter()'
        return '\n'.join(indent + line for line in (
            f'if ({guard}) sc1Conservation->beforePB();',
            'fluid.solve(rAUs, rAUfs);',
            f'if ({guard}) sc1Conservation->afterPBSolve();',
            'fluid.correct();',
            f'if ({guard}) sc1Conservation->afterPBCorrect();'))
    main, count = pattern.subn(pb_hooks, main)
    if count != 2:
        raise ValueError(f"Expected two adjacent solve/correct calls, found {count}")
    main = replace_once(main, '        runTime.write();', '''        if (sc1Conservation.valid())
        {
            sc1Conservation->sample
            (
                runTime.writeTime()
             || runTime.value() + SMALL >= runTime.endTime().value()
            );
        }
        runTime.write();''')

    # Retain the exact source matrix used in the existing equation. Never
    # substitute the repository YEqns.H for a differing cluster version.
    source = re.compile(r'fvModels\.source\(alpha,\s*rho,\s*Y\[i\]\)')
    if len(source.findall(species)) != 1:
        raise ValueError("Unrecognized YEqns fvModels source; inspect local source")
    species = source.sub('sc1ExternalSource()', species)
    species = replace_once(species, '            fvScalarMatrix YiEqn', '''            const tmp<fvScalarMatrix> sc1ExternalSource
            (
                fvModels.source(alpha, rho, Y[i])
            );
            fvScalarMatrix YiEqn''')
    species = replace_once(species, '            YiEqn.relax();', '''            if
            (
                sc1Conservation.valid() && pimple.finalPimpleIter()
             && Y[i].member() == "O2"
            )
            {
                sc1Conservation->equation
                (
                    phase, Y[i], YiEqn,
                    *specieTransfer[Y[i].name()], sc1ExternalSource()
                );
            }
            YiEqn.relax();''')
    species = replace_once(species, '    fluid.correctSpecies();', '''    if (sc1Conservation.valid() && pimple.finalPimpleIter())
        sc1Conservation->beforeSpeciesCorrection();
    fluid.correctSpecies();
    if (sc1Conservation.valid() && pimple.finalPimpleIter())
        sc1Conservation->afterSpeciesCorrection();''')
    return main, species


def scratch(path):
    path = path.resolve()
    allowed = (Path('/scratch/nwimer').resolve(), Path('/kfs3/scratch/nwimer').resolve())
    if not any(root in path.parents for root in allowed):
        raise ValueError(f"Must be below nwimer scratch: {path}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/scratch/nwimer/sc1_bubble_audit_20261029/solver'))
    parser.add_argument('--out', type=Path, required=True, help='New scratch root, not solver subdirectory')
    args = parser.parse_args()
    try:
        src, out = scratch(args.source), scratch(args.out)
        if out.exists() or out.is_symlink():
            raise ValueError(f"Destination exists: {out}")
        if src == out or src in out.parents or out in src.parents:
            raise ValueError('Source and destination must not overlap')
        originals = {name: (src/name).read_text() for name in ('putidaFoam.C', 'YEqns.H', 'Make/files', 'Make/options')}
        modified = instrument(originals['putidaFoam.C'], originals['YEqns.H'])
        make, count = re.subn(r'(?m)^\s*EXE\s*=.*$',
                             f'EXE = {out}/bin/putidaFoamSC1ConservationAudit', originals['Make/files'])
        if count != 1:
            raise ValueError('Unrecognized Make/files EXE entry')
        # Make/options still refers to ../include, ../phaseSystems, etc.
        dependencies = {}
        for name in ('include', 'phaseSystems', 'interfacialModels',
                     'interfacialCompositionModels', 'multiphaseCompressibleMomentumTransportModels'):
            target = (src.parent/name).resolve()
            if target.is_dir():
                dependencies[name] = target
            elif name != 'multiphaseCompressibleMomentumTransportModels':
                raise FileNotFoundError(f'Missing dependency: {target}')
        headers = {n: (AUDIT/n).read_text() for n in ('sc1BudgetLedger.H', 'sc1ConservationAudit.H')}
        out.mkdir(parents=True, exist_ok=False)
        shutil.copytree(src, out/'solver', ignore=shutil.ignore_patterns('linux*', '*.o', '*.dep', 'build.*.log'))
        (out/'bin').mkdir()
        for name, target in dependencies.items():
            (out/name).symlink_to(target, target_is_directory=True)
        (out/'solver/putidaFoam.C').write_text(modified[0])
        (out/'solver/YEqns.H').write_text(modified[1])
        (out/'solver/Make/files').write_text(make)
        for name, content in headers.items():
            (out/'solver'/name).write_text(content)
        manifest = {'source': str(src), 'dependencies': {k: str(v) for k,v in dependencies.items()},
                    'source_sha256': {k: hashlib.sha256(v.encode()).hexdigest() for k,v in originals.items()}}
        (out/'source_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
        print(f'Prepared {out}/solver; unique executable: {out}/bin/putidaFoamSC1ConservationAudit')
        print('No compilation or CFD launched. Review local source and compile before preparing cases.')
        return 0
    except (OSError, ValueError) as e:
        print(f'SC1 solver preparation stopped: {e}')
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
