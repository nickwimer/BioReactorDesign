# SC1 restored-drag startup: numerical failure and saved-log inspection

## Observed result

The matched drag comparison at `0876843` prepared successfully in
`/kfs3/scratch/nwimer/sc1_reactor_drag_20261111T144728Z`, decomposed to eight
ranks, and launched the existing verified `birdTwoFieldRates` binary. The
solver command returned 143; the startup wrapper returned 1. No completed
startup result, reconstructed five-snapshot series or drag comparison was
reported. Preserve the entire attempted root. Do not resume it as a valid
checkpoint or repeat its unchanged configuration.

The uploaded console transcript is `Pasted text(20261111-155434).txt`, SHA256
`bd2414190730f64171db029b7d07e1e668f02492975d16f591f2c971aee65c25`.
This is a console excerpt, **not the complete solver log**. It contains both
numerical divergence and allocation cancellation. It cannot establish the
earliest failure or the exact ordering of the FPE and Slurm signals.

| Evidence in the excerpt | Consequence |
|---|---|
| At t=0.0009 s, liquid alpha-call 1 divergence volume is 9.56981804819008e123 m3; clamp volume has the same printed value | An inner transport calculation is catastrophically invalid even if final clipping returns bounded fractions |
| At the same time, alpha-call 2 prints liquid min/max −1.488591507618836e142 / 9.952198638737628e141 | This is numerical divergence, not merely an interrupted valid trajectory |
| The displayed `SC1_RATE_STATE` has all four physical rates zero; displayed raw-dgdt and source-split terms are zero | The tail does not implicate the previously corrected nonnegative physical source split; earlier records still need inspection |
| Pressure reaches the 1000-Pa limiter and PCG reports zero iterations with residuals around 4.34e−31 | A small normalized pressure residual does not establish a valid coupled state after transport has diverged |
| Rank 2 reports an FPE in `phaseSystem::solve`; Slurm reports cancellation due to time limit | Both events must remain in the failure record; a fresh allocation alone is not a numerical remedy |

The first scheduled field output is at 0.002 s, later than the displayed
failure. Do not assume evolving velocity fields exist at 0.0009 s. The user
has already obtained allocation 19829576 on `x1005c3s4b1n1` for 24 hours; the
next command needs neither an allocation nor the OpenFOAM environment.

## What the source establishes and what remains unknown

The [drag-weight correction](SC1_REACTOR_STARTUP_2026-11-10.md) fixes a real
dictionary cancellation. It has **not** produced an acceptable reactor startup.
The preserved v1 control still has valid mass budgets and excessive velocities.
The eight small rate-split regressions remain evidence only for their tested
numerical contract.

The current startup recipe omits `faceMomentum` and `partialElimination`.
The local [solver entry point](../OFsolvers/birdmultiphaseEulerFoam/multiphaseEulerFoam/multiphaseEulerFoam.C)
defaults both to false, selecting cell momentum. In
[MomentumTransferPhaseSystem.C](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.C),
`momentumTransfer()` puts self drag in the matrix diagonal while `KdUByAs()`
and `phiKdPhis()` supply the other phase's contribution explicitly.
`partialElimination()` can invert the coupled drag system; it is not selected
by this recipe.

Drag coupling, phase/mixture Courant growth and discrete pressure/interface
balance need to be distinguished using the start of the saved log. No new
coefficient, timestep, tolerance or coupling option is changed by this update.
The next CFD candidate should follow that evidence and stop on numerical
divergence instead of running until the allocation expires.

## Next Kestrel command: collect existing evidence only

The collector reads the saved root, verifies its attempt against the startup
manifest, reports changed/missing/extra initial inputs, and finds the first
logged Courant breach, alpha excursion, appreciable clamp, pressure limiting,
nonfinite audit, FPE and time-limit cancellation. Every event carries a log
line number and preceding printed time/outer iteration. Courant lines are
explicitly labelled as preceding the next `Time` increment.

`report.json` includes a bounded per-time alpha/pressure chronology, last
gas/liquid step records and available time-directory names. These are evidence
coverage, **not acceptance or proof of complete saved fields**. The collector
does not call the success-only startup analyzer, `ldd`, a build, mesh utilities,
reconstruction or CFD. It does not revalidate external binaries/regressions.

The new bundle contains the complete solver log, manifest, attempt, available
mesh/decomposition logs, small constant/system dictionaries, report, first-event
contexts and SHA256 hashes. Initial fields, mesh, processors and later field
outputs remain in place and are not copied. Changed inputs remain findings;
path redirection, attempt mismatch and evidence changing during collection stop
the collector. A pre-existing output destination is never overwritten.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
FAILED=/scratch/nwimer/sc1_reactor_drag_20261111T144728Z
EVIDENCE="/scratch/nwimer/sc1_reactor_failure_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -x "$PY" && -r "$FAILED/quiescent/log.solver" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_failure.py" \
        --root "$FAILED" --out "$EVIDENCE"
    echo "Evidence collection status: $?"
else
    echo "Collection setup stopped. Preserve existing files and paste the output."
fi
```

Paste the console report and attach the printed `BUNDLE` file. A collection
status 0 means evidence was collected, never that the failed simulation passed.
The completed v1 control, failed v2 run, numerical suites, production cases,
solver source and biological schedule remain unchanged.

## Local verification

All 92 focused baseline, near-bound, rate-split, startup, drag and failure-report
tests pass. The new tests cover first-event chronology rather than just the
giant failure tail, Courant time labelling, pre-clamp failures hidden by bounded
step records, nonfinite values, duplicate records, exact full-log inclusion,
input change reporting, path/overwrite protection, concurrent log-change
rejection and no mutation of the saved case. Runtime errors now say `First
command errors:`; actual `wmake` failures retain `First compiler errors:`.
Kestrel log collection and the next CFD remedy remain unexecuted locally.
