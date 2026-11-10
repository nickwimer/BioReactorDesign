# SC1 wall/MRF preflight: implementation and completed result

**Completed on Kestrel. Do not rerun this preflight.** The uploaded transcript
`Pasted text(20261110-212403).txt` records execution of commit
`26c2aed226694920f50248b9eca47635c290af9a` and successful collection in
`/scratch/nwimer/sc1_wall_preflight_20261110T212232Z` (canonical root
`/kfs3/scratch/nwimer`). Both meshes and all four existing VTKs were analyzed.
No mesh-reader, boundary-query or VTK errors are listed. The review items are
unexpanded directives in each case's fvModels, phaseProperties, liquid thermo
properties and controlDict. These require inspection, not a repeat run.

## Observed result and next action

| Generated `rotor` adjacency | Historical EMM | Original MEML |
|---|---:|---:|
| Zone cells | 43,680 | 35,629 |
| `walls` faces adjacent / total | 0 / 14,392 | 0 / 17,200 |
| `outlet` faces adjacent / total | 0 / 1,920 | 0 / 1,920 |
| `propeller` faces adjacent / total | 4,952 / 7,928 | 4,500 / 9,064 |

Both literal MRF dictionaries specify empty `nonRotatingPatches`, but there
are no adjacent `walls` faces to receive the MRF rotating-boundary assignment.
This rules out that specific direct-wall-forcing mechanism for the audited
constant meshes. Do not add `walls` to the exclusion list as a supposed fix
for these results. It does not establish correct wall velocities from every
other boundary/source treatment or validate the entire MRF configuration.

Only 62.46% (EMM) and 49.65% (MEML) of faces named `propeller` are zone-adjacent.
These are face-count fractions, not area, torque or blade-coverage fractions.
The [earlier geometry audit](sc1_three_phase_headspace_audit.md) already shows
that this patch combines shaft/hub/blade surfaces and has different extents:
EMM y=0.208–1.397 m versus MEML y=0.218–1.741 m. MEML restricts its rotating
zone to two axial windows. Therefore partial adjacency is not proof that half
the blades are stationary. Inspect actual U boundary conditions and the
location/purpose of nonadjacent faces before deciding whether shaft/rotor
treatment is incorrect. Mesh/geometry and initialization differences still
prevent attributing the EMM/MEML behavior solely to phase count.

Sector measurements below use the highest downward alphaL=0.5 crossing:

| Case/time | Liquid (L) | Core r=0–0.05 m crossing range (m) | Wall-side r=0.50–0.53 m crossing range (m) |
|---|---:|---:|---:|
| EMM 9.9 s | 1151.076505 | 1.2026–1.2070 | 1.2825–1.3216 |
| EMM 20 s | 1142.971378 | 1.1911–1.1960 | 1.2313–1.3242 |
| MEML 9.9 s | 999.812137 | 0.6491–0.7399 | 1.1577–1.2328 |
| MEML 20 s | 999.824467 | 0.6483–0.7470 | 1.1413–1.2249 |

EMM loses 8.105127 L between the sampled times, whereas MEML gains 0.012330 L.
Neither difference is a closed mass budget without boundary/source accounting.
EMM has an elevated bulk level plus additional wall-side rise; MEML retains
both elevated wall-side liquid and a deep central depression. No physical
surface validation or wall-climbing remedy follows from either observation.
All reported alpha values are bounded, but boundedness alone is not conservation
or hydrodynamic validation.

All sectors have a detected 0.5 crossing. MEML's core has multiple downward
crossings in 6/12 sectors at 9.9 s and 8/12 at 20 s; EMM has multiple crossings
in 3/12 sectors of r=0.45–0.50 m at 9.9 s. These are binned profiles, not proof
of connected films, droplets or a unique sharp surface. EMM occupies 44/50
height bins per sector, versus MEML's 50/50. Inspect the full profiles and mesh
coordinates before interpreting missing sampling bins as a geometric limit.

**Immediate next step:** retrieve the existing `evidence.tar.gz` from this
completed report directory. The pasted terminal summary does not contain
the saved U/alpha/pressure boundary entries, full phase/numerical settings,
initialization dictionaries or individual sector profiles. The bundle has not
yet been provided to the assistant. Read it before requesting another Kestrel
diagnostic or implementing a supposed fix; resolve referenced includes if
their contents are not captured. The t=0 inventory and historical integrated
flux/source budgets also remain unestablished. No new simulation is needed
to retrieve this evidence, and no production settings are changed.

This implements the read-only checks in the
[model-selection decision](SC1_MODEL_SELECTION_AND_PROBE_RESULTS_2026-11-10.md#wall-climbing-is-a-separate-acceptance-test).
No solver compilation, mesh conversion, CFD launch, entrainment change or
biological-schedule change is needed for this step. Do not rerun completed
alpha/source probes or overwrite any previous report or case.

## What the collector measures

`scripts/diagnostics/collect_sc1_wall_preflight.py` reads the historical EMM
and original MEML archives and existing scratch VTKs at 9.9/20 s. It writes
only a newly created scratch directory. `sc1_wall_mesh.py` counts boundary
faces whose owner cells belong to the generated `rotor` cell zone, using
actual owner/neighbour/cellZones/boundary addressing. ASCII, binary labels
of declared width/endianness, and gzip are supported; unsupported addressing
is reported as unavailable rather than silently approximated.

Counts are **zone adjacency, not proof of a wall defect**. Compare patch
purposes and actual MRF exclusions. Rotating propeller faces are expected;
stationary tank/baffle faces need scrutiny. Empty/coupled patches have different
MRF treatment. Literal MRF dictionaries are printed for review, and complete
configuration evidence is saved. Includes/substitutions are not expanded.

Boundary fields from available initial directories and the 9.9 s root
checkpoint are queried with Foundation OF9 `foamDictionary -entry boundaryField`.
The [OF9 source](https://cpp.openfoam.org/v9/foamDictionary_8C_source.html)
disables function entries unless `-expand` is selected. The collector never
passes expansion or mutation options, hashes source fields before/after each
query, and runs the utility from the new report directory. Queries are disabled
unless the environment identifies OpenFOAM version 9. Python bytecode writing
is disabled in the command below.

The surface analysis uses 12 azimuthal sectors about the y axis, radial bands
0–0.05, 0.35–0.45, 0.45–0.50 and 0.50–0.53 m, and the earlier audit's 0.025 m
height bins. It retains every up/down crossing at alphaL=0.1/0.5/0.9, cell
counts, sampled volumes, missing bins, alpha bounds, total liquid and liquid
above 1.02 m. There is no interpolation across missing bins and no clipping
of input alpha. Summary ranges use highest downward alphaL=0.5 crossings;
the full JSON retains the other evidence.

These remain cell-center-binned volume averages, not a reconstructed connected
interface or thin-film resolution claim. Missing crossings are not evidence
of no liquid. VTK provenance relies on the existing exports already used in
the historical audit; input hashes are recorded, but the collector does not
prove an export equals its source field. It flags time-local meshes because
constant-mesh adjacency would then be insufficient. This step does not recover
t=0 liquid inventory or integrated historical flux/source budgets.

## Historical invocation (completed; do not repeat)

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_WALLROOT="/scratch/nwimer/sc1_wall_preflight_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -z "${SLURM_JOB_ID:-}" || ! -x "$PY" ]]; then
    echo "Stopped: active allocation or Python missing."
elif [[ "$(git -C "$REPO" branch --show-current)" != "codex/sc1-headspace-bubble-audit-20261029" ]]; then
    echo "Stopped: repository is on a different branch."
elif git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only origin codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && git -C "$REPO" rev-parse HEAD \
    && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/collect_sc1_wall_preflight.py" \
        --out "$SC1_WALLROOT"; then
    echo "REPORT=$SC1_WALLROOT"
    echo "Collection finished. Review all reported gaps; no CFD was launched."
else
    echo "Collection stopped. Preserve any created report and paste the output."
fi
```

The Python process can terminate without terminating the allocation shell.
The command contains no allocation-release command or top-level `set -e`.
The first useful return is the printed MRF dictionaries and summary. The
`evidence.tar.gz` bundle contains the summary, full JSON and both configuration
reports if additional inspection is needed. Do not rerun to replace a partial
report; retain it and inspect the reported gap.

## Validation and next decision

Fourteen local tests cover known patch-adjacency counts, label widths and
endianness, whitespace-valued first binary bytes, compressed inputs, malformed
addressing, volume weighting without clipping, sector variation, missing bins,
multiple crossings, volume/patch VTK selection, input/output preservation and
retention of an inspectable evidence bundle when VTKs are unavailable.
The pure-Python tests do not run OpenFOAM or PyVista. The subsequent Kestrel
execution is recorded above; full configuration/bundle interpretation remains
pending, and successful collection is not a physical-validation PASS.

If a stationary wall is incorrectly included in MRF rotation, prepare and audit
an isolated correction first. Otherwise use the configuration and sector data
to choose the smallest targeted pressure/interface or rotation benchmark.
Either route must still resolve conservation and species clipping before a
nonreacting checkpoint is accepted. Three fields remain optional; no model
change is accepted solely because it suppresses a visible symptom.
