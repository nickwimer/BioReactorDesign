# SC1 VOF saved-field inspection

## Completed Kestrel inspection; spatial archive pending

The command completed with status 0 at `8307261`, writing
`/scratch/nwimer/sc1_vof_spatial_20261111T195107Z`. It reports successful
processor/global addressing, mapped initial fields, finite native fields and
agreement with the logged inventory/velocity at all four selected times.
No CFD, reconstruction or case changes occurred. The [console result and
calculated cylindrical components](sc1_reactor_vof_spatial_console_2026-11-11.json)
are saved; the spatial archive has not yet been received or independently verified.

| At 8.25 ms | Gas peak | Liquid peak |
|---|---:|---:|
| Global cell | 199440 | 204379 |
| Speed, m/s | 0.955493898 | 0.000985278 |
| Radius, m | 0.270273556 | 0.270564788 |
| Height from initial midplane, mm | +12.617592 | -12.491097 |
| Liquid fraction | 2.01606e-8 | 0.9999999999999986 |
| Tangential velocity, m/s | +0.953835836 | -0.000977908 |
| Radial velocity, m/s | -0.056265181 | +0.000120287 |
| Vertical velocity, m/s | -0.000122349 | +7.14034e-7 |
| Adjacent boundary patches | None | None |

Here radial velocity is `(x*Ux+z*Uz)/r` and tangential velocity is
`(-z*Ux+x*Uz)/r`, with `r=sqrt(x*x+z*z)`. The reconstructed component norm
matches each reported speed. Tangential speed accounts for 99.8265% of the gas
peak and 99.2520% of the liquid peak. Both peaks are predominantly circumferential,
in nearly pure-phase cells on opposite sides of the initial transition layer,
and neither cell owns a boundary face. They are different cells/azimuths;
these values do not establish colocated counterflow, a global vortex or a mesh
join. The shared radius is a spatial clue, not a force attribution.

Gas peak speed rises from 0.029174 m/s at 0.25 ms to 0.116583 at 1 ms,
0.465039 at 4 ms and 0.955494 at 8.25 ms. Liquid stays near 1000 L.
The failed quiescent startup remains failed; the saved-field command's status 0
means successful inspection only.

**Next action: attach the existing archive; do not repeat the command below.**

`/scratch/nwimer/sc1_vof_spatial_20261111T195107Z/vof_spatial_evidence.tar.gz`

Reported SHA256:
`957cecdafd5b34a3220364da9709079bed7231fcbe71d1a87e312bf3311ff4ab`.

The full report contains phase/region statistics, native-field hashes and exact
peak-neighbour pressure/alpha/velocity stencils omitted from the compact console.
Verify the archive and those source identities, then inspect the stencil before
choosing a force-balance or discretization correction. No rebuild, new run,
relaxed threshold or new collection is needed to obtain this existing evidence.

## Purpose and verified evidence

The [VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped at 8.4 ms
because its liquid peak reached 1.003187 mm/s. Its phase mass budgets reconcile,
and all 672 pressure records remain finite. Gas peak speed nevertheless rises
to 0.972733 m/s. Raising the liquid threshold would hide a failed rest test.
The log alone cannot locate the gas peak or identify its driving force; the
completed saved-field inspection supplies the locations above.

The existing processor directory listings extend to 8.25 ms, before the rejected
8.4-ms state. The completed command reads
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

## Historical Kestrel command — completed; do not repeat

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

The console is received. Attach the existing `BUNDLE` listed above and preserve
every root; no new inspection is required.
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
