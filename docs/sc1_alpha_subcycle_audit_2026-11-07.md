# SC1 stage-control result and private alpha-subcycle audit

## Completed control and regression

The user-supplied `Pasted text(20261106-235025).txt` records the successful
52-rank control in `/scratch/nwimer/sc1_stage_pair_20261106T233034Z/control`.
It ran from 31.1999823727831 to 31.2999823727864 s in 2,000 steps, reported
solver status 0, `End`, and `Finalising parallel run`. ExecutionTime was
712.81 s. The normal case in this pair remains unrun; preserve both cases.

The existing revision-2 executable is
`/scratch/nwimer/sc1_stage_audit_20261106T231907Z/bin/putidaFoamSC1ConservationAudit`,
SHA256 `32262b1b49d57bc423a00a0cd7c393881a336d3c9471657e506401900db2fc42`.
Initial-state equality and all four Courant location metadata proofs passed.
Every printed v1 `SC1_BALANCE` field for all five phase/aggregate rows, and the
legacy bubble audit, match the original control's first report exactly at
displayed precision. This supports observational regression of revision 2;
it does not prove bitwise equality of every field.

## Oxygen residual attribution

All amounts below are cumulative over the actual approximately 0.1 s interval,
in mg. These are final-PIMPLE species ledgers, not sums of every outer solve.

| Phase | Previous O2 residual | Post-species inventory change | Stabilization contribution |
|---|---:|---:|---:|
| Liquid | +0.0004269842 | approximately zero | +0.0004269195 |
| Bubble | −0.4962267893 | −0.4974157010 | +0.0011889096 |
| Headspace | +0.5437204864 | +0.5437684760 | −0.0000480001 |
| Combined gas | +0.0474936971 | +0.0463527748 | +0.0011409096 |
| Total | +0.0479206813 | +0.0463527760 | +0.0015678291 |

The total remaining unexplained term is approximately −1.40e-9 mg; matrix
reconstruction mismatch L1 is approximately 1.22e-11 mg. The old-time and DDT
inventory defects are at rounding scale. Gas residuals are primarily inventory
changes after the species solve, associated with later pressure/thermodynamic
density updates while species mass fractions are held fixed. They are not
evidence of unexplained species leakage. Coupled convergence and physical
validity remain separate questions. Do not compare these final-iteration O2
changes directly with all-outer mass-stage totals, or multiply mass changes by
a phase-average mass fraction to infer an exact species change.

Liquid requested removal was **49.7847586202 mg**, with **9.67866235089 mg**
of interphase O2 transfer. The measured correction adds **40.1135511624525 mg**;
the independent integral of negative liquid O2 predicts
**40.1135511624520 mg**. Thus clipping cancels **80.57396%** of requested removal.
The requested-minus-correction amount, 9.67120745775 mg, is a numerical net
removal, not experimentally validated biological uptake. The final species
solve has 130,850 negative liquid cells, with minimum Y approximately
−2.50243e-11. The cumulative count 261,822,202 counts repeated phase-cell
samples, not unique cells. No gas clipping was measured.

These cases retain frozen OUR and the intended biological schedules. They are
not nonreacting checkpoints. No biological timestep, trigger, or sink physics
is changed by this result or the alpha diagnostic.

## Mass discrepancy localized to alpha allocation

Each of the five mass stages was called 4,000 times: two outer calls per
timestep. Unassigned mass changes are zero for every phase and aggregate.
The following cumulative allocation includes all outer iterations, in mg:

| Phase | Alpha/PB block | Energy block | Pressure block |
|---|---:|---:|---:|
| Liquid | +123.513743802 | −0.000002160 | 0 |
| Bubble | +201.453634226 | −0.049824432 | +7.757003407 |
| Headspace | −183.336301604 | −0.185520803 | −19.110761407 |
| Combined gas | +18.117332626 | −0.235345232 | −11.353758007 |
| Total | +141.631071415 | −0.235342100 | −11.353758282 |

Model-correction and species blocks contribute zero mass change. Alpha/PB
changes are entirely changes in alpha, with zero density contribution. The
liquid inventory rises **123.513741642 mg**. The sampled final-stage source
budget predicts **304.7386569 mg** (295.060 mg base source plus 9.6786569 mg
transfer), leaving **−181.2249153 mg**. Density cannot explain that difference.
The corresponding legacy bubble, headspace, combined-gas, and total mass
residuals are approximately +77.635748, −79.350487, −1.714739, and −182.939654 mg.

This localizes the remaining discrepancy to alpha/source/flux treatment. It
does not yet establish a mass leak: the previous ledger samples physical
sources after the solve, whereas alpha uses assembled Su/Sp, stored pressure
dilatation, limited fluxes, subcycles, and old-time resets. The repository's
`phaseSystemSolve.C` matches upstream OpenFOAM 9. No alpha-equation physics has
been changed. Increasing entrainment or foam coefficients is not justified.

## Diagnostic implementation and accounting

`prepare_sc1_alpha_audit.py` first uses the existing conservation solver
preparer, then replaces only the new tree's `phaseSystems` symlink with a full
materialized source copy. Generated build products are excluded. It injects
observations into the copied `phaseSystemSolve.C`; the tracked original and
installed sources remain unchanged. The private library retains the basename
`libbirdphaseSystem.so`, with an absolute scratch output target and the
private search path first in the new executable's link options.

The unique executable is `putidaFoamSC1AlphaAudit`. It requires an exported
`sc1AlphaAuditLibraryVersion()` marker, so falling back to the original
installed library cannot silently run an uninstrumented audit. No class
layout or public phase-system ABI is changed. The diagnostic is disabled by
default and enabled only with `sc1AlphaAudit true` in the alpha solver controls.
The case preparer's explicit `--alpha-audit` inserts that switch in the actual
`solvers/"alpha.*"` dictionary of fresh cases. All constants and the complete
initial state remain unchanged. Preparation never launches CFD.

For fixed mesh, global Euler stepping and geometricOneField MULES, the exact
subcycle identity is

\[
\alpha_a-\alpha_b=(\alpha_{old}-\alpha_b)
+\Delta t_{sub}(Su+Sp\alpha_a-\nabla\cdot\alpha\Phi).
\]

The audit measures current-before, actual oldTime, current-after, assembled
Su/Sp, and the limited flux used by each explicit solve. The oldTime reset
term is essential because later outer iterations reuse the same physical old
state. Repeated outer source actions must not be interpreted as additional
physical-time injection without their reset terms.

With alpha at assembly `a`, original mixture-divergence source `q`, raw
cross-phase dilatation `G`, and stored pressure input `D = phase.divU()`, it
also separates

\[
Su+Sp\alpha_a = D + (q+G-D) + (Su+Sp a-q-G) + Sp(\alpha_a-a).
\]

These identify stored pressure input, mixture coupling, source splitting at
assembly, and implicit evolution with frozen coefficients. Actual subcycle
durations are used. The ledger also captures the optional implicit alpha
correction, reference-headspace reconstruction, final bounds clamp, and final
averaged-flux versus subcycle-flux difference. It checks a full-solve inventory
telescope and equation residuals, including L1 norms. Source rows are emitted
only for explicitly solved phases; the reference phase is a complement, not
an independently solved MULES equation.

Each term is reported as volume and as mass weighted with the cell density
frozen at entry to that alpha solve. This exactly isolates alpha allocation
from density changes. Density-weighted divergence includes internal density
gradients; it is not mislabeled as a physical boundary mass flux. The unweighted
volume divergence is independently compared with physical boundary volume
flux. Primitive accumulators survive between calls; no registered field is
kept beyond the helper's mesh lifetime. Reductions occur once per complete
alpha solve rather than once per cell/subcycle term.

The audit rejects unsupported configurations: moving mesh, local time stepping,
non-Euler or per-field DDT schemes, more than one alpha correction, stationary
phases, a different reference phase, or a different phase count. The supplied
SC1 setup satisfies these assumptions. It does not reevaluate physical source
models. If stored `divU` already disagrees with the physical source budget,
the next attribution point is the existing assembled pressure source terms and
their timing, not an inferred defect in `massSource` or an arbitrary alpha fix.

## Build gate and subsequent run gate

All **105 local SC1 diagnostic tests pass**, including 32 new alpha
source-transformation, isolation/provenance, linkage and runtime-verifier
checks, plus the extended case-copy suite. The standalone arithmetic helper is
compiled and executed locally. This environment has no matching OpenFOAM 9
build; **private-library compilation and MPI validation are pending on Kestrel**.

Use a new `/scratch/nwimer/sc1_alpha_audit_<UTC>` destination with
`prepare_sc1_alpha_audit.py --out ...`. Build the copied `phaseSystems` using
`wmake -j 4 libso`, then the copied solver using `wmake`. **Do not use Allwmake**:
the inherited script can build shared installed targets outside this audit.
Capture each build to a new log and stop on failure while preserving the tree.

Run `verify_sc1_alpha_runtime.py --root ...` after the build, including the
film and entrainment libraries with repeatable `--loaded-library`, and the
velocity limiter with `--extra-library`. It checks provenance, original
installed-library integrity, exported/required marker, and dependency resolution
under a process-local private `LD_LIBRARY_PATH`. A phase-system dependency must
resolve to this exact private library. Missing dependencies are a failure.
Do not globally export the private library path. A standalone `ldd -r` on the
phase library is not a suitable gate: some of its symbols are supplied by the
full executable dependency graph.

After reviewing successful build/provenance output, prepare a fresh pair with
`--allow-courant-location-metadata --duration 0.1 --alpha-audit`. Launch only
its control with the private path scoped to the solver process, e.g. the
process portion `srun -n 52 env LD_LIBRARY_PATH="$ALPHAROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "$ALPHAROOT/bin/putidaFoamSC1AlphaAudit" -parallel -fileHandler collated`.
This is a later gate, not authorization to run the currently prepared stage-2
normal case. Check the library marker and actual method row, all legacy control
fields, alpha equation residuals, inventory telescopes, source decomposition,
and subcycle/outer counts. Only after the control explains the discrepancy
should the matched normal case run. Free-surface/foam validation against
experimental evidence and the nonreacting 9.9 s checkpoint remain subsequent
work; neither closure is approved for production.
