# SC1 VOF saved-field inspection

## Purpose and verified evidence

The [VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped at 8.4 ms
because its liquid peak reached 1.003187 mm/s. Its phase mass budgets reconcile,
and all 672 pressure records remain finite. Gas peak speed nevertheless rises
to 0.972733 m/s. Raising the liquid threshold would hide a failed rest test.
The log cannot locate the gas peak or identify its driving force.

The existing processor directory listings extend to 8.25 ms, before the rejected
8.4-ms state. Inspect these existing native fields first. The new command reads
four snapshots at 0.25, 1, 4 and 8.25 ms. It does not run OpenFOAM, reconstruct
fields, compile a utility or launch CFD. No allocation is required by the script.

## What the command verifies and measures

`sc1_reactor_vof_saved.py` requires the exact received source manifest and log
hashes. It checks all manifest-listed initial input hashes and the attempt
identity, then reads each processor's `cellProcAddressing`. These maps must form
one complete permutation of all 281664 global cells. Mapped processor time-zero
C, alpha.liquid, U and p_rgh must agree with the preserved global initial fields.
This checks the map against native initial data rather than trusting directory names.

At each selected time, all eight processors must contain complete, finite ASCII
alpha.liquid, U and p_rgh internal fields with the correct object/class headers.
The gathered liquid inventory, alpha redistribution and both phase-present
velocity maxima must agree with the saved log. Fractions and reconstructed
static pressure must remain within the original bounds. Missing, redirected,
changed, malformed or inconsistent inputs stop analysis; no replacement fields
are synthesized. No file in the source case is written.

For liquid and gas separately, the report gives phase-weighted RMS velocity,
mean vertical velocity and peak location/vector in phase-present, bulk,
interface, wall-adjacent, propeller-adjacent and vent-adjacent cells. Peaks
include radius, height relative to the initial interface, initial/current alpha,
pressure and pressure change. Their immediate internal neighbours retain raw
native values for inspecting spatial patterns. Adjacency means owning a boundary
face, not distance from a wall. These stencils are not an OpenFOAM force budget.

The observed 168-step prefix is independently reconciled from phase-volume
changes and actual transport fluxes, including the final guard-rejected record.
This does not accept the failed startup. No external executable, dependency or
prior control-chain revalidation is claimed for this read-only analysis.

All selected input hashes are checked again after reading. A fresh output root
contains a JSON report and compact archive with the report, source manifest,
attempt identity and full source log. The report records hashes of the native
fields it read; the large native arrays remain in the original case.

## Next Kestrel command

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
VOF_CASE=/scratch/nwimer/sc1_vof_reference_20261111T191346Z
VOF_SPATIAL="/scratch/nwimer/sc1_vof_spatial_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -x "$PY" && -r "$VOF_CASE/quiescent/log.solver" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_saved.py" \
        --root "$VOF_CASE" --out "$VOF_SPATIAL"
    echo "Saved VOF inspection status: $?"
else
    echo "Inspection setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed `BUNDLE`. Preserve every existing root.
There is no new CFD result until a later, justified candidate is run and checked.

## Local verification and next decision

The 29 focused tests cover the VOF reference, incomplete-log collector and saved
reader. New checks exercise known nonzero flux budgets, partition-order mapping,
duplicate/missing/out-of-range addressing, native/log disagreement, wrong mapping,
missing and redirected fields, source mutation during reading, output overlap,
preserved input hashes, and compact archive contents. The actual received log
also reconciles with this code. No local OpenFOAM or Kestrel field access exists.

Use the measured spatial pattern to choose the next specific force-balance or
discretization correction. The reference removed independent phase slip but
did not establish rest. The current evidence does not select a new phase count,
prove a capillary or gravity cause, justify artificial damping, or validate
rotation, aeration, PBM, oxygen transfer or biology.
