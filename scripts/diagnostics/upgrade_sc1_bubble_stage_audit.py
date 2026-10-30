#!/usr/bin/env python3
"""Upgrade ONLY the scratch-copy RHEL9 putidaFoam bubble audit.

The script is a child process of the interactive shell: errors print warnings;
it never changes the production solver or CFD case, and does not invoke MPI.
Usage:
  python3 scripts/diagnostics/upgrade_sc1_bubble_stage_audit.py
"""
from pathlib import Path
import re
import shutil


def main():
    scratch = Path("/scratch/nwimer").resolve()
    solver = Path("/scratch/nwimer/sc1_bubble_audit_20261029/solver").resolve()
    root = Path(__file__).resolve().parents[2]
    audit_dir = root / "OFsolvers/birdmultiphaseEulerFoam/putidaFoam"
    src_path = solver / "putidaFoam.C"
    header = audit_dir / "bubbleMassAudit.H"
    stage_header = audit_dir / "bubbleMassTransferStageAudit.H"

    if scratch not in solver.parents:
        print("CHECK FAILED: target solver is outside canonical scratch.")
        return

    required = (src_path, header, stage_header)
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        print("CHECK FAILED: missing required file(s):", *missing, sep="\n  ")
        return

    src = src_path.read_text()
    if src.count('#include "bubbleMassAudit.H"') != 1:
        print("CHECK FAILED: existing scratch solver audit hook is unexpected.")
        return

    # Only the two explicit calls already present in local putidaFoam.C.
    pattern = re.compile(r"(?m)^([ \t]*)fluid\.correct\(\);$")
    count = len(pattern.findall(src))
    if count != 2:
        print(f"CHECK FAILED: expected two fluid.correct() calls; found {count}.")
        return

    guard = 'bubbleMassTransferStageAudit.H'
    if guard in src:
        if src.count(f'#include "{guard}"') != 2:
            print("CHECK FAILED: stage include already partially installed.")
            return
        new_src = src
        print("Stage hooks already present; refreshing audit headers.")
    else:
        def add_stage(match):
            indent = match.group(1)
            lines = [
                match.group(0),
                indent + "// Measure transfer before the pressure-velocity correction.",
                indent + "if",
                indent + "(",
                indent + "    pimple.dict().lookupOrDefault<Switch>",
                indent + "    (",
                indent + '        "bubbleMassAudit",',
                indent + "        false",
                indent + "    )",
                indent + " && runTime.writeTime()",
                indent + ")",
                indent + "{",
                indent + '    #include "bubbleMassTransferStageAudit.H"',
                indent + "}",
            ]
            return "\n".join(lines)

        new_src, patched = pattern.subn(add_stage, src)
        if patched != 2:
            print("CHECK FAILED: unable to safely instrument both calls.")
            return

    # All operations below are restricted to the separately copied solver.
    shutil.copy2(header, solver / "bubbleMassAudit.H")
    shutil.copy2(stage_header, solver / "bubbleMassTransferStageAudit.H")
    if new_src != src:
        src_path.write_text(new_src)

    print("Scratch source updated; production solver and case untouched.")
    print("Stage include count:", new_src.count(f'#include "{guard}"'))
    print("Active deposition lookup:")
    for line in header.read_text().splitlines():
        if 'fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas"))' in line:
            print(" ", line.strip())
    print("Rebuild from:", solver)
    print("Next: (cd", solver, "&& wmake > build.bubbleAudit.stage.log 2>&1)")


if __name__ == "__main__":
    main()
