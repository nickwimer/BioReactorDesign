# SC1 restored-drag startup: numerical failure and saved-log inspection

## Bundle received and verified

The bundle below is now received. Its SHA256 and all 18 bundled source-file
hashes match. The corrected parser places the actual FPE at line 2273,
t=0.0009 s, outer 2. All 23 pressure solves through the first limiter event
reach residuals below 1e-12; all 36 physical-rate records are zero. See the
[parsed result](sc1_reactor_drag_failure_result_2026-11-11.json) and the
[subsequent face-momentum comparison](SC1_REACTOR_FACE_MOMENTUM_2026-11-11.md),
which failed on its first timestep. That new bundle is received and inspected;
the subsequent [coupled-drag comparison](SC1_REACTOR_COUPLED_DRAG_2026-11-11.md)
also failed at the same first-step outer iteration. Its bundle is now verified. The
[short pressure probe](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md) has now reproduced
the pressure collapse and stopped before limiting; its bundle is now verified.
The [private pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
stopped at 0.2 ms on the liquid-speed guard; its bundle is verified.
The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md) completed; its verified bundle also reveals horizontal gas motion near the interface. The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped on the liquid-speed guard at 8.4 ms; its bundle is verified. Converged nonorthogonal pressure corrections also failed to restore rest; next use the [two-step force probe](SC1_REACTOR_VOF_FORCE_PROBE_2026-11-11.md).
No recollection or repetition of the failed run is needed.

The collector at `264d209` completed with status 0 in
`/kfs3/scratch/nwimer/sc1_reactor_failure_20261111T161441Z`.
All 59 recorded initial inputs match, with none missing or extra. The log has
2,284 lines, 18 printed times and 17 paired step records through 0.00085 s;
there is no `End` marker. These records do not imply valid accepted states.

| First reported event | Time context | Evidence |
|---|---|---|
| Alpha excursion beyond 1e-10 | 0.0002 s, outer 2, line 545 | Min −1.232699447738956e−10, max 1.0000000005446585; initially a very small excursion |
| Pressure limiter | 0.0003 s, outer 2, line 832 | `p min -73364.513746858589` Pa |
| Mixture Co >0.2 | After printed 0.0004 s, line 1089 | Max 0.37319948362935979; Courant is printed before the next Time increment |
| Clamp volume >1e-10 m3 | 0.00065 s, outer 1, line 1617 | −1.36104293969645e−9 m3 |
| Alpha magnitude >10 | 0.00065 s, outer 2, line 1651 | Min −16.379226889599494, max 14.623227456604468 |

The full-log summary reports maximum finite Co 2.7252723295830327e29 and
largest finite alpha magnitude 7.433392801508894e144. Pressure limiting
precedes the first logged mixture-Co breach; this warrants inspecting the
early pressure/momentum sequence, but does not by itself establish causality
or rule out a phase-specific Courant problem.

The collector incorrectly called the startup banner `sigFpe : Enabling floating
point exception trapping (FOAM_SIGFPE).` an exception at line 34. That banner
only enables trapping. The detector now excludes it, with a regression test
that still detects the later stack trace and actual crash. The reported first
FPE line/count must not be used as the failure onset. Other reported events
are unaffected by this string-matching correction.

The user first pasted the collection summary, then attached this existing file:

`/scratch/nwimer/sc1_reactor_failure_20261111T161441Z/reactor_failure_evidence.tar.gz`

- Reported bundle SHA256: `cbe15f3cc08564aed19324ee0c397ae6296ee0475f99dea0ba9ffb290b60adf4`.
- Reported solver-log SHA256: `4c55d6b46ec8738a4c3f67179843044b2b45b6101bc5a32e51a0c1eb4bba214f`.

Those hashes have been verified and the bundled log analyzed with the corrected
parser locally. The pressure sequence motivates the linked face-momentum
comparison with the timestep and tolerances held fixed. The small
negative-alpha excursion alone is not evidence of the dominant mechanism.

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

## Historical collection command (completed; do not repeat)

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
Kestrel log collection is complete and its bundle has been inspected locally;
the next CFD comparison remains unexecuted locally.
After the Kestrel collection exposed the banner misclassification, all eight
failure-collector tests passed with the corrected detector. The full 92-test
run above predates this one-line parser correction and added regression test;
no solver code or startup settings changed.
