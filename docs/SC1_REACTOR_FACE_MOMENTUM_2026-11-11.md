# SC1 face-momentum startup comparison

## Evidence and decision

The failed drag bundle has been received and verified: SHA256
`cbe15f3cc08564aed19324ee0c397ae6296ee0475f99dea0ba9ffb290b60adf4`, all 18
included source-file hashes match, and solver-log SHA256 is
`4c55d6b46ec8738a4c3f67179843044b2b45b6101bc5a32e51a0c1eb4bba214f`.
The corrected [parsed result](sc1_reactor_drag_failure_result_2026-11-11.json)
retains the chronology. The 59-input identity result comes from the Kestrel
collector; those large initial fields and mesh are not in the uploaded bundle.

The 23 pressure solves through the first limiter event all finish below
1e-12 residual (maximum 9.985295722848731e-13). At 0.0003 s, the first outer
iteration starts pressure corrections at residuals 0.0098453 and 0.0253179.
The second outer iteration then reaches `p min -73364.513746858589` Pa.
Its preceding logged mixture Co is only 0.0144981. Subsequent Co maxima grow
to 0.0557276, 0.153297 and 0.373199. All 36 physical-rate records are zero.
The first actual FPE indication is line 2273 at 0.0009 s, not the startup
trapping banner at line 34. The source-floor correction is inactive here.

This evidence points to a coupled momentum/pressure instability rather than
failure to reach the requested linear-solver tolerance. It does not prove a
unique cause: phase Courant numbers, force imbalance and velocity locations
were not logged, and the failure precedes the first scheduled field write.
No change to tolerances or physical coefficients is justified by this log.

**The next comparison enables only `PIMPLE/faceMomentum true` in the numerical
equations.** It also writes snapshots more frequently and monitors failures.
The old startup inherited cell momentum from the small uniform test fixture.
The repository's [nonreacting reactor setup](../scripts/putida_stir/nonreact/system/fvSolution)
already selects face momentum. This supplies a concrete alternative to test,
not proof that it will fix this quiescent vessel.

## What changes in the solver

Both paths are compiled into the existing verified `birdTwoFieldRates` binary.
The [entry point](../OFsolvers/birdmultiphaseEulerFoam/multiphaseEulerFoam/multiphaseEulerFoam.C)
selects `pU` or `pUf` from `faceMomentum`. The default is false.

In the [cell path](../OFsolvers/birdmultiphaseEulerFoam/multiphaseEulerFoam/pU/pEqn.H),
the reciprocal momentum diagonal is formed in cells, then interpolated for
pressure/flux coupling; cell velocity and face flux are corrected separately.
In the [face path](../OFsolvers/birdmultiphaseEulerFoam/multiphaseEulerFoam/pUf/pEqn.H),
the reciprocal diagonal is formed on faces from old-time phase mass,
interpolated momentum diagonal and face drag. Corrected phase flux supplies
the reconstructed velocity. The face momentum assembly refreshes the same
drag coefficients. Both paths include the same instrumented `pEqnComps.H`,
so the physical-source checks and accepted-step audit remain active.

`partialElimination` remains at its existing false default; no simultaneous
coupled-drag inversion is introduced. The official Foundation OF9
[damBreak4phase fvSolution](https://github.com/OpenFOAM/OpenFOAM-9/blob/89839ae3b8cd7b7c8d044f6898d4991177b1db41/tutorials/multiphase/multiphaseEulerFoam/laminar/damBreak4phase/system/fvSolution)
omits both switches, so face momentum is not claimed to be a universal
requirement for segregated drag.

## Case contract and execution

The variant is `sc1_meml_quiescent_face_momentum_v3`. `prepare-face` validates
the failed v2 case's original inputs, exact received log hash and attempt/build
identity, plus its completed v1 control and the existing build/regression
provenance. It copies only the v2 `0`, `constant`, `system` and mesh-preparation
logs into a fresh destination. It copies no failed fields or processor output.

Only two dictionary files differ from v2:

| File | Change |
|---|---|
| `system/fvSolution` | Add `faceMomentum true` to PIMPLE |
| `system/controlDict` | Write every 5 steps instead of 40: 0.00025 s instead of 0.002 s |

Mesh, exact 1000-L fill, hydrostatic initial fields, stationary walls/impellers,
gravity, vent, densities, viscosities, restored segregated drag, surface tension,
compression, source-free physics, dt=5e-5 s, 200 steps, two outer/two pressure
correctors, six alpha subcycles, pressure tolerances and eight ranks all remain
unchanged. No source code or private binary is rebuilt. The earlier small
regressions are revalidated from saved evidence, not rerun; they are not claimed
to have tested this newly selected face path with nonzero sources.

The runner saves/reconstructs 40 snapshots through 0.01 s and applies the same
inventory, velocity and interface-redistribution thresholds at every snapshot.
It requires exactly 400 `Constructing face momentum equations` runtime records
and the unique drag-blending marker. Accepted steps and all alpha-call budget
identities retain their existing complete-coverage checks. The final report
also records present/bulk/interface velocities, peak locations and phase slip
in `face_comparison.json`.

The new run uses a live Python monitor and `srun --kill-on-bad-exit=1` with a
30-minute step limit. It stops its own child step on the first detected:

- Pressure-limiter activation in the quiescent vessel.
- Logged mixture Co >0.2, matching the existing startup criterion.
- Alpha subcycle outside [0,1] by more than 1e-8, or nonfinite alpha/Co/audit/pressure residual.
- Fatal solver error or FPE (excluding the trapping-enabled banner).
- 30-minute wall-clock limit.

The subcycle threshold tolerates insignificant printed roundoff; it does not
relax the existing accepted-field or mass thresholds. Detection is polled at
0.2-second intervals, not an atomic check inside a solver iteration. The monitor
signals only the new child session; it never cancels the allocation or another
job. Srun's [signal and failed-rank handling](https://slurm.schedmd.com/srun.html)
provides the corresponding MPI-step termination. Logs and already written
fields remain intact; an interrupted write is not a valid checkpoint.

On monitored solver failure the runner attempts a small automatic evidence
bundle in a fresh sibling `<startup-root>_failure` directory, including
`guard_stop.json`. It does not reconstruct or analyze the incomplete trajectory
as a successful run. The bundle omits large fields, which remain on Kestrel.
If collection itself fails, the original solver failure and log remain reported.

## Next Kestrel command

Use the existing allocation. This prepares a new t=0 case and runs the single
10-ms comparison; the earlier roots remain unchanged.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
FAILED=/scratch/nwimer/sc1_reactor_drag_20261111T144728Z
FACE_STARTUP="/scratch/nwimer/sc1_reactor_face_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
        prepare-face --failed-root "$FAILED" --out "$FACE_STARTUP"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
            run --root "$FACE_STARTUP" --run
        echo "Face-momentum startup status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output. If the monitor stops the solver, attach its printed `BUNDLE`
file. Status 0 is a completed startup pass; 2 is a completed failed screen;
1 is a stopped/incomplete attempt. No automatic extension, rotation or aeration
follows. A pass is a short startup result, not equilibrium or reactor validation.

## Verification

All 101 focused numerical, startup, drag, failure-collection and face tests
pass locally. The new tests enforce the two-file change, preserve both controls,
reject altered inputs/logs, require runtime branch coverage and every early
snapshot, and verify MPI arguments and reconstruction order. Live child-process
tests exercise pressure-triggered termination, partial final log lines, quick
success/failure, timeout and no overwrite. A stopped run bundles evidence and
does not enter reconstruction/reporting. Replaying the received failed log through
the monitor selects line 832 at 0.0003 s, before the large alpha excursions.
Legacy report/recipe tests remain green.
OpenFOAM and MPI CFD execution are not available locally; the new case has not
yet run on Kestrel.
