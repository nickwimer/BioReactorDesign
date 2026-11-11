# SC1 VOF passive face-force probe

## Received result and next action

The probe completed on Kestrel at `f951b53` in
`/scratch/nwimer/sc1_vof_force_20261111T205814Z`. The attached archive and all
16 included file hashes verify, and all force vectors reconstruct independently
from their face records. The [full result](sc1_reactor_vof_force_result_2026-11-11.json)
identifies a predominantly horizontal gravity/pressure residual at the selected
peaks; capillarity and matrix-versus-fresh pressure differences are much smaller.
The 83-record prefix comparison passed on Kestrel; its external parent log is
not in the uploaded archive. This is diagnostic completion, not startup success.

Next use the [hydrostatic pressure perturbation candidate](SC1_REACTOR_VOF_HYDROSTATIC_2026-11-11.md).
The command below is historical; do not repeat the completed probe.

## Original motivation

The [nonorthogonal comparison](sc1_reactor_vof_nonorth_result_2026-11-11.json)
reports all 672 observed final pressure passes converged, yet the same 8.4-ms
liquid-speed stop and nearly identical gas speed remain. This strongly disfavors
explicit pressure-correction lag as the dominant error in that comparison.
The console does not establish a discrete force balance or identify the cause.
Its archive has not been received locally; preparation checks the exact reported
log hash, revalidates the full control/build chain and reconciles its budget prefix.

The next action is a two-step observation of the actual velocity correction,
not another parameter sweep. Build a separate application from the verified
`sc1InterFoamReference` source tree. Insert an observer immediately after the
native final-pressure velocity assignment. The original equation statements,
phase transport, physical forces and startup guards remain unchanged.

## Measurements and interpretation

The native update is `U = HbyA + rAU*reconstruct((phig + pressureFlux)/rAUf)`.
At each of eight pressure corrections across two timesteps, record:

- Gravity, capillary, matrix-pressure and net reconstructed contributions at
  the global gas-present and liquid-present peak-speed cells (fraction >0.01).
- Predictor velocity, actual velocity, alpha, density, rAU, pressure and location.
- Every nonempty adjacent face's area vector, centre, local identity, patch/coupled
  status, gravity potential, density gradient and signed normal force terms.
- Pressure from the matrix flux versus a freshly evaluated pressure snGrad.
- Global maximum closure errors for the component sum and native velocity update.

Cell contributions include the native rAU and have velocity units (m/s).
Face contributions are pressure gradients (Pa/m). They are not multiple
independently evolved phase velocities. The predictor includes effects already
assembled in the momentum equation, including velocity history and transport.

Also log the diagnostic difference between face-centre gravity potential and
interpolated cell potential. This produces a **frozen-pressure diagnostic delta**
only: it is never applied to the solver or pressure equation, and cannot by
itself establish that changing gravity discretization would restore balance.

The [OF9 reconstruction operator](https://cpp.openfoam.org/v9/fvcReconstruct_8C_source.html)
uses the area-normal tensor and face contributions. The Python report independently
reconstructs each printed peak-cell contribution from all its face records and
checks agreement with the native result. Large opposing face terms are compared
with an explicit floating-point roundoff allowance. These arithmetic checks do
not alter any physical startup limit.

The observer elects a single rank for each global peak and reduces its data
before master-only printing. Cell/face labels are processor-local. Processor
copies are not a global face integral. Empty faces are skipped; the probe requires
the existing three-dimensional reactor mesh. Internal-face dC is owner-to-neighbour;
boundary-face dC is zero and must not be interpreted as a neighbour location.

## Scope and evidence checks

Fresh t=0 copy of `/scratch/nwimer/sc1_vof_nonorth_20261111T203602Z`.
Only application name and end time (0.0001 s) change in `system/controlDict`.
Retain dt=5e-5 s, eight ranks, ten pressure passes, all initial fields, mesh,
gradient schemes, boundaries, forces and material properties. The private build
has the same runtime libraries as the original VOF application. Neither build
nor any completed/failed control is overwritten. Solver wall time is capped at
ten minutes and the application has an additional two-step bound.

Completion requires exactly two times, a normal End marker, eight force samples,
16 peak-cell records and complete unique face coverage. All 83 original audit
records (one initial, 80 pressure, two steps) must agree with the failed control's
two-step prefix within relative 1e-9 / absolute 1e-12 tolerance. This checks the
logged trajectory; it is not a bitwise comparison of every field value.
The final peak speeds at each timestep must match the original step audit.

No reconstruction, restart or new physical correction is performed. Existing
five-step output cadence is retained; the observation is delivered in the log
and result JSON. Automatic evidence collection includes the full force result
and log on normal completion and the available log on failure. Console output
shows the first and last peak records; the bundle contains every observation.

## Historical Kestrel command

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
VOF_BUILD=/scratch/nwimer/sc1_vof_build_20261111T191346Z
NONORTH_CONTROL=/scratch/nwimer/sc1_vof_nonorth_20261111T203602Z
SC1_FORCE_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
FORCE_BUILD="/scratch/nwimer/sc1_vof_force_build_${SC1_FORCE_STAMP}"
FORCE_CASE="/scratch/nwimer/sc1_vof_force_${SC1_FORCE_STAMP}"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_force.py" \
        build --base-build "$VOF_BUILD" --out "$FORCE_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_force.py" \
        prepare --control-root "$NONORTH_CONTROL" --build-root "$FORCE_BUILD" --out "$FORCE_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_force.py" \
            run --root "$FORCE_CASE" --run
        echo "VOF force-probe status: $?"
    else
        echo "Build/preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed BUNDLE. Preserve the nonorthogonal and
gradient failure archives as well. No repeat of either failed startup is requested.

## Validation and remaining work

All 51 focused Python tests pass across the force observer workflow, nonorthogonal
and gradient comparisons, VOF reference, saved-field inspector and failure collector.
Tests include independent reconstruction on rotated faces, cancellation roundoff,
missing/corrupt observations, passive-prefix disagreement, preserved native
equation statements, private build destination/source validation, fresh-case
copying, attempt reuse rejection and automatic failure evidence.
OpenFOAM is unavailable locally. Compilation and the two-step diagnostic have now completed on Kestrel; see the received result above.

`VOF_FORCE_PROBE=COMPLETE` means the diagnostic completed and its checks passed.
It does not mean the rest failure was corrected. Use the measured term vectors
and face geometry to identify the next force-discretization change, then repeat
the unchanged physical startup checks. No operating-reactor, aeration,
oxygen-transfer or biology acceptance is implied; the four-hour biological
update schedule remains unchanged.
