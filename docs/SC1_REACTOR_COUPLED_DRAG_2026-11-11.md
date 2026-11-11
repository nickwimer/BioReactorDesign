# SC1 coupled-drag startup comparison

## Received face-run evidence and selected next change

The v3 face-momentum bundle is now received as
`reactor_failure_evidence.tar(1).gz`. Its SHA256 matches
`8eac9f4f97692121eea0c7ff35efdc8a94e0fd0d8890134dbb6ad35a6abd7bfb`.
All 19 bundled evidence hashes match; the 13 included initial dictionaries also
match the startup manifest. Native initial fields and mesh remain on Kestrel;
the full 59-input match is the collector's recorded check.

The [verified result](sc1_reactor_face_failure_result_2026-11-11.json) now
includes all four pressure solves and the two alpha-call records. The log
confirms two face-momentum assemblies and no coupled-drag inversion.

| t=0.00005 s | Initial pressure residual | Final residual | Iterations | Outcome |
|---|---:|---:|---:|---|
| Outer 1, correction 1 | 0.00138420 | 9.1420e-13 | 244 | No pressure limiter |
| Outer 1, correction 2 | 6.8244e-7 | 9.8066e-13 | 185 | No pressure limiter |
| Outer 2, correction 1 | 0.00417225 | 9.6382e-13 | 345 | Pressure minimum −2499398.2044661008 Pa; guard triggers |
| Outer 2, correction 2 | 2.6109e-5 | 9.6738e-13 | 301 | Printed while the step was being terminated |

All three printed rate-state records have zero physical rates. Both alpha calls
report zero raw-dgdt, pressure-dilatation, source-split, clamp and rescale terms.
The second liquid alpha call changes the field by only 1.14518505493816e-11 m3
in L1, with an equation-residual L1 of 1.80398486558592e-18 m3. Accepted phase
fractions remain in [0,1], their sum error is zero, and the single-step mass
residuals are tiny. None of this validates the pressure state or the full run.

The failure appears after the second outer iteration assembles momentum/drag,
despite converged pressure linear solves and essentially unchanged alpha. This
supports testing the existing coupled-drag solve. It does not isolate the cause:
the log contains no drag coefficients, phase velocities or phase Courant numbers.
No claim is made that tighter tolerances, a smaller timestep, or this new option
has already resolved it.

## Numerical mechanism

The current face path includes self drag in each phase's diagonal and the
other phase's flux explicitly. The local
[partialEliminationf implementation](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.C)
removes that explicit contribution, forms a coupled matrix, and solves it.
For exactly two phases its normalized matrix is

```text
diagonal entries: 1, 1
off-diagonal entries: -rAUf_gas*Kf, -rAUf_liquid*Kf
determinant: 1 - (rAUf_gas*Kf)*(rAUf_liquid*Kf)
```

After inversion, `setMixturePhi` applies the same flux correction to both
phases to restore the pressure equation's mixture flux for complementary
fractions. It retains relative phase flux. This is an existing OpenFOAM/BiRD
solver path, already compiled into the verified private binary.

The next variant, `sc1_meml_quiescent_coupled_drag_v4`, changes **only
`PIMPLE/partialElimination true`** relative to the failed v3 input. Face momentum
remains enabled. No physical coefficient, phase guard, timestep, pressure
tolerance, corrector count, snapshot interval or acceptance limit changes.
This tests simultaneous drag coupling rather than adding a new physical closure.

## Preparation, evidence and stopping behavior

`prepare-coupled` verifies the exact received v3 log hash
`083bbd2d39499fff0c19210e2280f3c13d9b0ec55508752953a8d3fe8151654c`, its
attempt/input identity, and the preserved v2/v1/build/regression chain. It copies
only original `0`, `constant`, `system` and mesh-preparation logs to a fresh
root. Only `system/fvSolution` differs. No failed fields, processor output,
solver logs or attempt files are copied. No CFD runs during preparation.

The run retains eight ranks, dt=5e-5 s, 200 steps through 0.01 s, snapshots
every 0.00025 s, complete phase/source audits, and all prior quiescent criteria.
The same monitor stops pressure limiting, mixture Co >0.2, substantial alpha
excursions, nonfinite values, fatal errors or a 30-minute timeout. It now also
stops a malformed, nonfinite or nonpositive logged face-drag determinant.

On a complete run, reporting requires 400 face-momentum assembly records and
800 coupled-drag determinant records (two outer iterations and two pressure
corrections per timestep). Every determinant must be finite and positive;
the minimum is reported. No arbitrary small-positive cutoff is imposed.
Positivity alone does not establish matrix conditioning or physical validity.

The usual 40 reconstructed snapshots, mass budgets and present-phase speed
limits must also pass. Final flow statistics and the immediate v3 control
identity are saved in `coupled_comparison.json`. A stopped solver produces
the existing automatic failure bundle in `<new-root>_failure`; all earlier
roots remain unchanged. No automatic extension, rotation or aeration follows.

## Next Kestrel command

This uses the existing verified build and the current allocation. No rebuild
or rerun of the small numerical suites is needed.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
FAILED_FACE=/scratch/nwimer/sc1_reactor_face_20261111T163440Z
COUPLED_STARTUP="/scratch/nwimer/sc1_reactor_coupled_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
        prepare-coupled --failed-root "$FAILED_FACE" --out "$COUPLED_STARTUP"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
            run --root "$COUPLED_STARTUP" --run
        echo "Coupled-drag startup status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach a `BUNDLE` file if one is printed. A startup pass
remains a 10-ms numerical/quiet-flow result, not equilibrium, operating-reactor,
oxygen or biology validation. The four-hour biological schedule is unchanged.

## Local verification

All 107 focused baseline, near-bound, rate-split, startup, drag, failure, face
and coupled tests pass. Six new tests enforce the one-dictionary change,
preserve all controls, reject changed control hashes/extra parameter changes,
require complete positive determinant records, and retain the monitored MPI
runner and early reconstruction interval.

A scalar C++ test compiles the actual decomposition/back-substitution loops
extracted from `partialEliminationf`, checks them against an independent
two-by-two analytical solution over weak and strong positive couplings, and
checks the actual `setMixturePhi` loops for mixture-flux restoration and retained
relative flux. This verifies isolated algebra, not the full field/pressure solve.
OpenFOAM/MPI execution is unavailable locally; the v4 CFD run remains pending.
