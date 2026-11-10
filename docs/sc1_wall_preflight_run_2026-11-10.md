# SC1 wall/MRF preflight: next Kestrel action

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

## Run once in the active allocation

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
The pure-Python tests do not run OpenFOAM or PyVista; the actual Kestrel mesh,
utility and VTK reader execution remain to be checked by this preflight.

If a stationary wall is incorrectly included in MRF rotation, prepare and audit
an isolated correction first. Otherwise use the configuration and sector data
to choose the smallest targeted pressure/interface or rotation benchmark.
Either route must still resolve conservation and species clipping before a
nonreacting checkpoint is accepted. Three fields remain optional; no model
change is accepted solely because it suppresses a visible symptom.
