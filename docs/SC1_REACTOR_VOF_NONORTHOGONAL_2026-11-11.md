# SC1 VOF nonorthogonal pressure comparison

## Kestrel result: converged pressure corrections, failed rest state

At `ef01ab8`, `/scratch/nwimer/sc1_vof_nonorth_20261111T203602Z` stopped
at 0.0084 s. All 672 observed correction groups converged on their final pass;
the largest final-pass initial residual is 9.98761133388736e-13. The summary's
`pressure_iteration_pass: false` reflects an incomplete 200-step run, not a
failed observed correction. Gas peak remains 0.972526435841 m/s, changing only
-2.60328e-8 m/s from the gradient control. Liquid peak remains 0.001003115908 m/s
and again exceeds 1 mm/s. There are 168 step and 6720 pressure records, no End,
and unchanged initial-input hashes. See the [received console](sc1_reactor_vof_nonorth_result_2026-11-11.json).
The new failure archive remains on Kestrel and has not been locally verified.

This weighs strongly against explicit pressure-correction lag as the dominant
cause in this comparison. It does not establish facewise force balance.
**Next: the [passive two-step force probe](SC1_REACTOR_VOF_FORCE_PROBE_2026-11-11.md).**
Do not repeat this failed startup or increase the pressure iteration count again.

## Result motivating this test

The [gradient candidate console](sc1_reactor_vof_gradients_result_2026-11-11.json)
reports the same 8.4-ms failure as the original VOF reference. Gas speed changes
from 0.9727328354 to 0.9725264619 m/s and liquid speed from 0.001003187005 to
0.001003115918 m/s. The gradient changes alone have not remedied the rest state.
The candidate archive has not been received locally. Preparation on Kestrel
checks its exact reported log hash and independently reconciles the 168-step
prefix before allowing the new run. That prefix includes the rejected last step.

The existing `nNonOrthogonalCorrectors 0` gives one pressure solve per correction.
In the pinned OF9 source, the pressure matrix's explicit nonorthogonal contribution
is rebuilt within `correctNonOrthogonal()`. A small final linear residual does
not establish convergence after that explicit contribution is updated.
The verified original VOF log contains 672 pressure solves, all with initial
residual above 1e-12 and at least one linear iteration. Its first four initial
residuals are 4.53561e-4, 1.59681e-7, 3.53284e-9 and 4.06759e-11.
These are separate PIMPLE pressure corrections; their changing predictors mean
the numbers alone do not measure the nonorthogonal defect or prove its cause.

[Liu, Tolle and Marić](https://doi.org/10.1140/epjs/s11734-025-01758-x)
show that explicit nonorthogonal correction can prevent force balance despite
consistent gradient schemes. Their stationary-water tests support checking
correction convergence, rather than treating a final linear residual as sufficient.
This motivates the comparison; their results do not establish SC1's root cause.

## Exact change and acceptance

Start from a fresh t=0 copy of the failed gradient candidate. Change only
`system/fvSolution`, replacing:

```text
nNonOrthogonalCorrectors 0;
```

with:

```text
nNonOrthogonalCorrectors 9;
```

The [OF9 loop](https://cpp.openfoam.org/v9/nonOrthogonalSolutionControl_8C_source.html)
and [final-pass condition](https://cpp.openfoam.org/v9/nonOrthogonalSolutionControlI_8H_source.html)
make this ten passes, not nine. This fixed choice is a bounded convergence test,
not the paper's adaptive algorithm or a guaranteed sufficient iteration count.

The existing binary is reused without compilation. Retain the four least-squares
scalar gradients, all physical properties and forces, boundaries, initial fields,
mesh, pressure tolerance, two outer/two pressure corrections, timestep 5e-5 s,
end time 0.01 s, eight ranks and 30-minute solver cap. There is no added damping,
MRF or aeration. Both the original and gradient controls remain untouched.

After normal completion, validate all 8000 pressure records and matching linear
solver records, their time/order, finite pressure and unchanged pressure limits.
Each final pass of the 800 correction groups must enter with residual <=1e-12
and perform zero linear iterations. Failure of that test remains a failed
convergence comparison even if the physical screening numbers happen to pass.

Only after checking every pass, an in-memory view containing the final pressure
observation from each group is passed to the unchanged mass-budget parser.
No CFD log is rewritten, and no intermediate pressure state is exempted from
validation. Require 200 steps, a normal End marker and all requested processor
fields before reconstruction. Native inventory and both phase-present velocity
maxima must agree with the audit at all 40 written snapshots.

Original phase-budget, pressure, alpha/complement, liquid speed <=1 mm/s,
Co <=0.2 and redistribution <=0.01-L criteria remain. Maximum logged gas-present
speed must also be <=1 mm/s. `nonorth_pass` requires these criteria AND pressure
iteration convergence. The unchanged binary retains its original early guards;
gas/convergence checks are post-run requirements, not newly added solver guards.

On failure, print and bundle any valid observed pressure-iteration prefix as
`nonorth_pressure_iterations.json`. A converged prefix is not a completed run.
The normal snapshot comparison uses the original VOF spatial report, whose
fields are available; it is not a new spatial comparison with the gradient case.

## Historical Kestrel command (already failed)

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
GRADIENT_CONTROL=/scratch/nwimer/sc1_vof_gradients_20261111T201422Z
NONORTH_CASE="/scratch/nwimer/sc1_vof_nonorth_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_nonorth.py" \
        prepare --control-root "$GRADIENT_CONTROL" --out "$NONORTH_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_nonorth.py" \
            run --root "$NONORTH_CASE" --run
        echo "VOF nonorthogonal comparison status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed BUNDLE. Preserve the gradient evidence
archive as well; no repeat of either failed control is requested.

## Validation and next decision

All 45 focused Python tests pass across this candidate, gradient and VOF controls,
saved-field inspection and failure collection. Coverage includes strict complete
pressure passes, rejected intermediate states, explicit versus linear convergence,
incomplete-prefix reporting, unchanged budgets, native phase-speed agreement,
single-dictionary copying, attempted-root rejection and evidence on failure.
Original recipe/build files are unchanged. No local OpenFOAM runtime is available.

A passing 10-ms result would permit longer-rest and timestep-sensitivity checks
before controlled stirring. If the final pressure passes converge but the speeds
do not improve, that weighs against explicit pressure-correction lag as the
dominant problem and calls for a discrete face-force analysis. It does not justify
relaxing the rest limits or disabling physical forces in an operating case.
No operating-reactor, aeration, oxygen-transfer or biology acceptance is implied.
The four-hour biological update schedule is unchanged.
