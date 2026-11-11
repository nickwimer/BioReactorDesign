# SC1 MEML oxygen-transfer investigation — new-session handoff

**Updated:** 2026-11-11

**Research system:** *Pseudomonas putida* SC1 stirred-tank simulation, OpenFOAM 9 / BiRD, NLR Kestrel  
**GitHub branch:** `codex/sc1-headspace-bubble-audit-20261029`  
**Pull request:** https://github.com/nickwimer/BioReactorDesign/pull/4  
**Status:** Experimental closures compile and run in isolated scratch cases. The oxygen deficit remains unresolved; no experimental closure has been approved for production.

**User direction (2026-11-10): three fields are optional.** The assistant had
introduced the gas split as a proposed conservation remedy. The user explicitly
requires retaining it only if needed. The objective is the simplest physically
validated, conservative nonreacting 9.9 s checkpoint, with a common hydrodynamic
basis for MM/EMM/MEML biology. Do not treat a three-field implementation as the
deliverable or expand its closures before evaluating a simpler candidate.

**Wall-climbing defect remains an explicit requirement.** The user reiterated
that the old liquid interface climbed artificially up the tank sides. A near-
1000 L inventory does not establish that the three-field model resolved it.
Check actual MRF-zone intersections with stationary walls/baffles and their
boundary treatment, then local wall films and azimuthal surface variation using
existing outputs. Annular means can hide the reported artifact. See the
[wall-climbing assessment](SC1_MODEL_SELECTION_AND_PROBE_RESULTS_2026-11-10.md#wall-climbing-is-a-separate-acceptance-test).
Neither a phase-count change nor coefficient tuning is an established remedy.

**Latest result: wall/MRF preflight completed.** The [completed report](sc1_wall_preflight_run_2026-11-10.md)
is `/scratch/nwimer/sc1_wall_preflight_20261110T212232Z`. Both cases have zero
`walls` faces adjacent to `rotor`, ruling out the specific direct MRF
rotating-wall assignment hypothesis for these meshes. Their partial
`propeller` adjacency needs component/velocity interpretation, not an automatic
zone expansion. EMM retains an elevated bulk level plus wall-side rise;
MEML retains elevated wall-side liquid and a deep central depression despite
near-1000 L inventory. Phase count remains optional and neither surface is
physically validated. Do not rerun the preflight. **The full evidence bundle
has now been received and inspected.** See [the configuration findings](SC1_WALL_BUNDLE_FINDINGS_2026-11-10.md).
EMM lacks a reference phase (selecting the inspected solver's alpha-rescaling
branch) and its empty compression table also omits the resolved capillary-force
path. MEML additionally changes the liquid outlet to slip, turbulence, source
layout and fill initialization. These are concrete mechanisms/confounders,
not a quantitative attribution of the old 151 L anomaly or proof of a fix.
**The prescribed-fill check is now complete, status 0**, at
`/kfs3/scratch/nwimer/sc1_prescribed_fill_20261110T215018Z`; do not rerun it.
EMM's archived prescription gives 999.182377 L, with cell-selection bounds
999.182377–1010.106301 L, versus 1151.076505 L in the 9.9 s field. Thus the
observed inventory exceeds even that upper bound by 140.970204 L. MEML's
estimate is 999.741417 L, bounds 993.247714–1004.405140 L, versus 999.812137 L
at 9.9 s. These are not recovered time-zero values or closed conservation
budgets; VTK integration error is outside the selection bounds. Prescribed
overfilling does not explain EMM under these assumptions. The rescaling,
boundary/source and pressure/density contributions still need measurement.
**The isolated two-field test harness builds and passes runtime-library
verification on Kestrel. All four baseline cases now pass all 100 steps of
the numerical screen: the closed pair at `56a84d7` and the source pair at
`ef74cf6`. Preserve the completed suite; do not rerun it. The separate
near-bound suite is also complete: coarse at `78ebbb8`, fine and comparison
at `75dc652`. All four diagnostics pass but all phase mass screens fail;
halving dt leaves the approximately 0.900-microgram liquid deficit unchanged
to within 0.016%.** See the [completed comparison and source-algebra review](SC1_TWO_FIELD_NEAR_BOUND_TEST_2026-11-10.md).

**Latest result:** the [private rate-split candidate](SC1_TWO_FIELD_RATE_SPLIT_TEST_2026-11-10.md)
built and passed all eight Kestrel CFD regressions at `d83bc83`, status 0.
Build `/scratch/nwimer/sc1_two_field_rates_build_20261111T021516Z`; cases
`/scratch/nwimer/sc1_two_field_rates_cases_20261111T021516Z`. Closed/mid-fraction
results are retained; near-bound liquid residual L1 drops by about 1,860–116,567×.
Original thresholds are unchanged. All diagnostics pass. Preserve every
completed control/candidate case and both builds; do not rerun these suites.

**Completed vessel control:** the [first MEML-mesh quiescent startup](SC1_REACTOR_STARTUP_2026-11-10.md)
completed at `e586c47` in `/scratch/nwimer/sc1_reactor_startup_20261111T045638Z`.
Its mass/diagnostic screens pass but quiescence fails (status 2): maximum
liquid speed 0.0963573 m/s, gas speed 2.320774 m/s, maximum mixture Co 0.623991.
Liquid stays at 1000 L to printed roundoff; its physical residual L1 is
3.344917609101449e-12 kg. Preserve this completed control and do not extend it.

**Earlier vessel attempt:** the restored-drag comparison at `0876843` in
`/scratch/nwimer/sc1_reactor_drag_20261111T144728Z` diverged by t=0.0009 s;
liquid alpha extrema reached about 1e142 before an FPE/time-limit cancellation.
The command returned 143, wrapper status 1. This is not simply a timeout and
not a valid checkpoint. The user already has a new 24-hour allocation.
Preserve both vessel roots; do not repeat either launch.

**Collection completed at `264d209`, status 0:** all 59 initial inputs match.
Pressure limiting begins at 0.0003 s, before the first logged mixture-Co >0.2
after printed 0.0004 s. The collector falsely classified the initial FPE-trapping
banner as a crash; the parser is corrected. See [the chronology and hashes](SC1_REACTOR_DRAG_FAILURE_2026-11-11.md).
**The bundle is now received and verified**, including all 18 bundled hashes.
All 23 pressure solves through the first limiter event reached residuals below
1e-12; all 36 physical-rate records are zero.

**Earlier vessel attempt:** the [face-momentum comparison](SC1_REACTOR_FACE_MOMENTUM_2026-11-11.md)
at `7604521` in `/scratch/nwimer/sc1_reactor_face_20261111T163440Z` failed on
its first timestep: t=0.00005 s, outer 2, `p min -2499398.2044661008` Pa.
The monitor stopped the step (137; wrapper status 1) and collected evidence.
All 59 inputs match; alpha remains bounded and max logged mixture Co is only
0.005974812253287854. There is one paired step record, no End and no reached
0.00025-s snapshot. Face momentum alone is not the remedy; preserve all roots.
**The face-run bundle is now received and verified:** all 19 included hashes
match; the full log confirms converged pressure solves, tiny alpha changes and
zero source/clamp terms before the second-outer pressure failure. See the
[verified result](sc1_reactor_face_failure_result_2026-11-11.json).
**Earlier vessel attempt:** the [coupled-drag startup](SC1_REACTOR_COUPLED_DRAG_2026-11-11.md)
at `8b81ede` in `/scratch/nwimer/sc1_reactor_coupled_20261111T165652Z` also
failed at t=0.00005 s, outer 2: `p min -2499398.1871999367` Pa. This differs
from v3 by only 0.0173 Pa. The monitor stopped the step (137; wrapper status 1).
The collector reports all 59 inputs unchanged, two paired records through
0.0001 s, no End, and max logged mixture Co 0.005975836748875787. No scheduled
snapshot was reached; this is not a valid restart or a passed startup.
**The v4 bundle is now received and verified:** all 19 hashes and 13 included
initial dictionaries match. The [verified result](sc1_reactor_coupled_failure_result_2026-11-11.json)
confirms active coupled drag, first-step minimum determinant falling from almost
1 to 0.0135944981275532, eight converged pressure solves and five zero-rate records.
**Latest vessel diagnostic:** the [short pressure probe](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md)
at `1dbe432` in `/scratch/nwimer/sc1_pressure_probe_20261111T172942Z` reaches
the intended pre-limit stop at the identical v4 pressure minimum, first timestep
outer 2. The face drag diagonal maximum rises from 0.112828 to 1.45118e9;
minimum mixture rAUf falls by a factor of 73.55. These global extrema need not
be colocated. All 59 inputs match; no completed-step records or restart exist.
The [bundle and full log](sc1_reactor_pressure_probe_result_2026-11-11.json)
are now received and verified: all 20 file hashes and 13 included dictionaries
match. The pressure guard precedes all eight fatal errors; return 139 arises
with MPI_Abort/sigSegv shutdown records afterward. At the minimum-pressure
cell alphaGas=0 and alphaLiquid=1, with pre-correction gas speed 0.387579 m/s
and liquid speed 1.30908e-7 m/s. Outer 1 produces absent-gas speeds up to
0.545059 m/s; outer-2 explicit liquid drag flux dominates its predictor range.

**Latest vessel attempt:** the [pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
built successfully at `408028e`, with binary SHA256
`e5414f2afb8192d2fe88cd4f881ac554a7605326b1938ccbfe5fc62291937dd4`
and the same six verified rate libraries. Build
`/scratch/nwimer/sc1_pressure_coupled_build_20261111T180325Z`; case
`/scratch/nwimer/sc1_pressure_coupled_20261111T180325Z`.
The run stopped at 0.0002 s, outer 2 (guard, then MPI_Abort/139, wrapper status 1).
The [bundle is now verified](sc1_reactor_pressure_coupling_result_2026-11-11.json):
archive, all 19 file hashes and 13 included dictionaries match. The exact stop
is `Quiescent liquid speed exceeded 1 mm/s`; peak liquid speed grows through
0.245, 0.529, 0.800 and 1.060 mm/s across four steps. All 16 pressure solves
converge and remain between 99982.986886 and 111847.796701 Pa; all 16 mixture
flux checks agree to <5e-18 m3/s. No nonfinite values in 624 field reductions.
The previous pressure collapse does not recur in this interval, but quiescence
still fails. No 0.25-ms native field snapshot or valid restart was produced.

The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md)
completed at `fa0636a` in `/scratch/nwimer/sc1_velocity_probe_20261111T183832Z`,
status 0. It reused the exact compiled v5 build, changing only endTime and
writeInterval in fresh t=0 inputs. All 33 compared records match the saved
v5 prefix. Liquid peaks occur at y=1.018662536 m, alphaLiquid approximately
0.49485, with downward liquid and upward gas motion. At 0.15 ms, present-liquid
peak is 0.800102 mm/s and bulk-liquid peak is 0.070844 mm/s. The last two global
peaks are not direct boundary-face owners; wall-adjacent peaks in the same
interface-height layer are almost equal. This points toward interface-normal
counterflow, without establishing a discrete force budget or accepted remedy.

The spatial-diagnostic bundle at
`/scratch/nwimer/sc1_velocity_probe_20261111T183832Z_evidence/reactor_failure_evidence.tar.gz`,
with SHA256 `5ba88937370743795665315d0f266522071bd592d813072a773875ada326191d`,
is now verified: archive, all 19 bundled hashes, 13 initial dictionaries and
33 prefix records match. The [full result is saved](sc1_reactor_velocity_location_result_2026-11-11.json).
It reveals predominantly horizontal gas peaks of 40.0167 mm/s in the
interface and 34.8828 mm/s in the bulk-gas layer above it at 0.15 ms. Neither
peak owns a boundary face. Normal counterflow at the liquid maximum is only
part of the disturbance. Whole-domain volume-mixture RMS is 1.482643 mm/s.

**Latest result: the [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md)
build passed; CFD stopped at printed time 0.0084 s.** The logging API error
in the first build (`Info.precision(17)`) was corrected to `Info().precision(17)`.
The fresh build at `a8349f1`, `/scratch/nwimer/sc1_vof_build_20261111T191346Z`,
passed compilation/dependency verification. Preserve it and the earlier failed
build `/scratch/nwimer/sc1_vof_build_20261111T190909Z`; do not rebuild.
The reference case is `/scratch/nwimer/sc1_vof_reference_20261111T191346Z`.
There are 168 printed time labels, no End, solver return 139 and wrapper status 1.
The fatal banner appears at line 4962, outer 2; its explanation is absent from
the console. All 49 initial inputs are reported unchanged. The maximum printed
Co is 0.001778043328599279, which does not establish pressure or speed acceptance.
The old collector's zero pair count/alpha magnitude cannot be interpreted for
VOF. It did not parse VOF audit/fraction lines; the read-only parser now does.

**Next action: receive and inspect the already-created evidence archive.**
`/scratch/nwimer/sc1_vof_reference_20261111T191346Z_evidence/reactor_failure_evidence.tar.gz`
has reported SHA256 `ab4e882c60e137ca99f4a646cf288c3324c2dfb6dd9348f3186f9fde3e785cb8`.
The [console result](sc1_reactor_vof_result_2026-11-11.json) is saved; the archive
has not yet been received or independently verified. Do not infer the cause,
rerun, rebuild, restart or recollect. Inspect the fatal context and reconcile the
VOF budgets from the existing full log first. Thirty-nine focused tests pass
for the reporting update. No CFD equations, guards or build recipes changed.
No rotation, aeration, PBM, oxygen-transfer or biology validation is implied.

The startup dictionary's symmetric 0.3/0.7 blend made the unordered segregated
drag weight `1-f1-f2` nearly zero. The attempted fresh t=0 copy added only
`blending/drag { type none; continuousPhase none; }`, restoring weight one.
It kept the same verified binary, initial fields, mesh, timestep, ranks and
acceptance limits. Sixteen startup tests pass, including actual C++ blending
functions and backward compatibility of the saved v1 control. Restored drag
alone has not produced a working startup.
No rotation, aeration, PBM, oxygen transfer or biology is enabled.

It compares no-reference and gas-reference settings
on a small closed periodic transport case, then on a known liquid-source case
with pressure-compatible outflow. Each pair differs only in the reference
setting. A new dynamically sized alpha observer measures rescaling as well as
the other existing alpha stages; a separate timestep observer captures actual
phase masses, boundary mass fluxes and assembled pressure sources. All six
required project libraries are rebuilt privately with the nonreacting solver.
**The first Kestrel build failed** at `a0f1929` in
`/kfs3/scratch/nwimer/sc1_two_field_build_20261110T222049Z`, after the first two
library commands completed. The third failed in the observer's chained field
initialization: OpenFOAM `Field::operator=` returns void. The assignments are
now separate; the two new headers were checked for similar chains. A targeted
C++ regression reproduces the old failure and passes with the actual corrected
initialization. All 208 local SC1 tests pass, including 26 harness tests.
**The corrected build at `a5989c1` subsequently passed all six library builds,
the solver build and runtime verification** in
`/kfs3/scratch/nwimer/sc1_two_field_build_20261110T224950Z`. The executable
SHA-256 is `1483c42b8bb3aeb37be743eafbd70f1b3e59d3ed053388ccd54f7c63a044f1f5`;
all seven binary hashes are recorded in the test document. No CFD launched
during that build. Preserve both builds and **do not rebuild**.

**The first closed-pair attempt** at
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T231637Z` completed mesh
generation, mesh checks and initialization for `closed_no_reference`, then
stopped at solver startup: OF9 rejects `pMin` in `fvSolution/PIMPLE`. The
generator now retains the pressure reference in PIMPLE and puts the same
1000 Pa floor in `system/fvConstraints` using `limitPressure`. No solver source
or binary changed. All 27 harness tests passed at that correction.

**The second closed-pair attempt**, at `c78199e`, created
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T232902Z`. The pressure
constraint loaded successfully, and initial gas/liquid masses were 0.0005/0.5 kg.
At time label 0.001 s, first PIMPLE iteration, `closed_no_reference` stopped
because `fvSolution/solvers` had only literal `alpha`, not a match for
`alpha.gas`. No timestep completed. The momentum-transfer implicit-pressure
check and MULES both look up phase-specific alpha fields. The generator now
uses `"alpha.*"`, matching `alpha`, `alpha.gas` and `alpha.liquid` with the same
controls. Only that dictionary key changed; all 29 harness tests passed, including
lookup coverage and reproduction of the old missing-phase-key failure.

**The third attempt**, at `d7b1b2e`, created
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T233431Z`.
`closed_no_reference` returned from the solver with status 0, then reporting
rejected `Incomplete or unexpected initial/step/alpha-call records`.
The harness reached this check after verifying unchanged initial inputs and
finding its library marker and solver `End`. No `result.json` was written.
The other closed case and both source cases were not launched.

An exact-float timestamp comparison reproduces this rejection on an otherwise
complete synthetic stream. The reporter now matches times only within 1e-12 s
of the prescribed grid, retains raw values, rejects matched-key duplicates,
and prints coverage counts plus missing/unexpected examples.

**Saved-log analysis at `e7825bd` passes:** the user supplied `complete: true`,
`screen_pass: true`, no errors and 100 steps. Gas cumulative absolute mass
residual is `5.769354784863862e-16` kg (limit `1.5e-12`); liquid is
`5.770176179219254e-13` kg (limit `5.01e-10`). Initial/final liquid mass is
0.5 / 0.499999999999939 kg. Both phases have zero reported external flux and
sources, with negligible final-outer rescaling. Of 3,360 adjusted time keys,
the maximum correction was `1.1102230246251565e-16` s, confirming timestamp
roundoff as the reporting failure. This validates only the small closed
constant-density transport fixture; it does not explain archived EMM excess
inventory or select the reactor model. The raw log remains on Kestrel.

**Do not rebuild or rerun either closed case.** The `continue-closed` command
revalidates build, attempt provenance, initialized/prepared inputs, mesh-check
evidence and the completed passing log, saves the missing `result.json` with
its log hash, then launches only the fresh gas-reference partner. Existing
inconsistent results, altered inputs or an attempted partner are refused.
All 39 local harness tests passed before this latest runtime result, including
mocked execution and preservation checks.

**Closed pair completed at `56a84d7`:** the user supplied the `REUSED=` record,
both passing 100-step reports, both log hashes, and `PAIR=closed complete`.
Gas-reference cumulative absolute residuals are `4.894737643911179e-16` kg gas
and `4.895160184386637e-13` kg liquid, below `1.5e-12` / `5.01e-10` kg limits.
Final gas/liquid masses are 0.000500000000000093 / 0.499999999999907 kg;
source, boundary and rescale terms are zero in that variant. Both closed
variants conserve within screening tolerances, so this fixture does not
select between them or explain archived EMM excess liquid.

The saved solver-log hashes are:

- no-reference: `9793a67274cf9260750c7b822f24bc0e9cf66a00d04f4b0cd707700042c212df`
- gas-reference: `1e563cb5179927771ffea8ccf691a89488f71cb871c500dd785117fd7e3e17ac`

**Source pair completed at `ef74cf6`:** both variants report complete, passing
100-step screens with no errors, followed by `PAIR=source complete`. They ran
in `/scratch/nwimer/sc1_two_field_cases_20261110T233431Z` with the same build.
Each received 0.0001 kg liquid over 0.1 s, retained about 9.04191e-05 kg and
exported 9.58092e-06 kg liquid. Liquid cumulative absolute residuals are
`6.991205047567109e-11` kg (no reference) and `3.8064394580038755e-11` kg
(gas reference), below the `1.001e-07` kg limit. Gas residuals are
`6.122027062821979e-15` / `3.79419149386864e-14` kg, below approximately
`1.90419e-10` kg. Gas inventory loss is accounted for by outlet export.
No interphase transfer was enabled. Source-split diagnostics are effectively
zero; no-reference final-outer rescaling is only about 3.1e-14 m3 per phase.

Source log hashes:

- no-reference: `b3124354624ce8ea8ff56e5e7663e89d3e36c385f01e0396e723314576fa6d96`
- gas-reference: `487bec0cc44e837d2de8ae23c075511a0433344afd27736b6d8e14a3ecff7b7d`

All four `result.json` files now exist with their log hashes. The user supplied
the runtime reports; raw evidence remains on Kestrel. Both two-field reference
settings passed that baseline; the near-bound failures below remain unresolved.
The uniform 50/50 source fixture does not
exercise the near-bound source-splitting mechanism diagnosed in the earlier
three-field probe or explain archived EMM excess liquid. See the [completed
source budgets](SC1_TWO_FIELD_CONSERVATION_TEST_2026-11-10.md#next-implementation-test-source-splitting-near-phase-bounds).
Preserve all completed and failed roots. **No rebuild, repeat of these four
cases, or reactor launch is the next action.** The [separate near-bound
harness](SC1_TWO_FIELD_NEAR_BOUND_TEST_2026-11-10.md) now creates four fresh
cases with alphaL=0.99999 / alphaG=0.00001: both reference settings at dt=0.001
and 0.0005 s, each to 0.1 s. It derives the existing verified build from the
completed baseline and revalidates every saved baseline result before launch.
All 53 focused tests passed locally, including unchanged baseline report output
and mocked preservation checks. No compiled source or baseline recipe changed.
Both pairs and the comparison have now completed. `diagnostic_pass` requires valid setup,
regime and accounting; `screen_pass` also requires the phase mass thresholds.
A complete mass-screen failure can proceed to timestep sensitivity only if
its diagnostics pass. No tolerances are automatically loosened. The existing
observer checks accepted-step bounds and split terms, not source-assembly
alpha bins.

**Near-bound coarse pair completed at `78ebbb8`**, in
`/kfs3/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z`, using the unchanged
verified build. Both 100-step reports have `complete: true`,
`diagnostic_pass: true`, `sourceSplitObserved: true` and `screen_pass: false`.
The only errors are the two phase-mass thresholds. Liquid signed residuals
are -9.001408870364855e-10 / -8.996987667297552e-10 kg (no reference / gas
reference), approximately 0.900 micrograms of budget shortfall and 8.91 times
the 1.01e-10 kg limit. Gas residuals are +8.996994139032305e-13 /
+8.99699284789516e-13 kg, about 156.75 times their limits. Liquid inventory
still increases; its deficit is relative to measured feed and export.
All steps stay in the target alpha regime, above the gas residual-alpha floor.

Liquid final-outer source splitting is -8.996993607492838e-10 /
-8.996993257360911e-10 kg in frozen-density units, closely matching the
accepted-step deficit. No-reference liquid rescaling is only
-4.41516614935189e-13 kg; gas reference has zero rescaling but retains the
mass-budget failure. This supports source splitting as the dominant measured
contribution here, not a quantitative attribution of the reactor anomaly.
Gas-reference reconstruction does not ensure a closed gas mass budget.
The raw logs and results remain on Kestrel; the user supplied the reports.
Their hashes are recorded in the linked test document. Preserve both coarse
cases. Both fine cases subsequently ran to the same endpoint as recorded below.

**Near-bound fine pair and comparison completed at `75dc652`.** The user
uploaded `Pasted text(20261111-011306).txt`; both 200-step reports have valid
diagnostics and the same two mass-screen errors. Fine liquid signed residuals
are -9.000017554931859e-10 / -8.997481055502226e-10 kg (no reference / gas
reference); fine/coarse L1 ratios are 0.9998454335923372 /
1.0000548392665323. Gas ratios are 1.00006636640375 /
1.0000545349602175. Source splitting still closely matches the liquid deficit.
The read-only comparison succeeded and reproduced all four log hashes.
The full parsed reports, comparison and transcript provenance are saved in
[the evidence JSON](sc1_two_field_bounds_fine_results_2026-11-10.json).
Raw solver logs and initialized inputs remain on Kestrel.

All eight baseline/near-bound cases are now complete; preserve them and the
verified build. No further timestep/reference sweep or reactor launch is the
next action. The source review shows that the 1e-4 denominator retains only
about 10% of the coupled dilatation at the tested 1e-5 distance from a bound.
The suppressed component predicts approximately 9e-10 kg liquid deficit in
this uniform-feed fixture. This is not a 90% loss of total feed and does not
quantify the old reactor anomaly.

The [implemented candidate](SC1_TWO_FIELD_RATE_SPLIT_TEST_2026-11-10.md)
is a guarded private two-field source-split variant: for complementary
fractions and nonnegative rates qi,qj, use increments Su=qi and
Sp=-(qi+qj), preserving `(1-a)*qi-a*qj` without a small-fraction denominator.
Source-only algebra/bounds and guards were checked locally, and all eight
pressure/flux-coupled CFD regressions now pass at `d83bc83`. The implementation preserves the
original divergence source and all diagnostics, uses the distinct executable
`birdTwoFieldRates` and library marker, and rejects unsupported states. Only
the private copied split changes. The linked historical command completed
the fresh build and all eight regression cases with unchanged thresholds.
Do not repeat it. The actual-mesh quiescent startup is the next action. No production sources,
original build or completed cases are overwritten. Do not reduce the floor
or add the measured deficit back after solving.
The existing three-field executable is not used for these cases.
See the [completed fill result and test requirements](SC1_WALL_BUNDLE_FINDINGS_2026-11-10.md#5-prescribed-fill-check-completed-the-emm-prescription-is-near-1000-l).

**Git transport correction:** the earlier hostname-resolution failure was
followed by a successful pull on the same allocated node; do not assert a
compute-node network restriction. That successful transcript still used an
HTTPS fetch URL. The user explicitly requires SSH for Git operations; future
Kestrel instructions must honor this and need not move the user to a login
node on the basis of the one failed lookup.

**Preceding source-probe result:** The probe in
`/scratch/nwimer/sc1_source_probe_control_20261110T194214Z` completed successfully:
status 0, 2,000 steps over 0.1 s, ExecutionTime 749.47 s. All legacy audit records
match baseline exactly at printed precision. The 16 new rows pass completeness
checks and their bin sums independently match the same-outer alpha ledger.
Final-outer liquid splitting is −96.438397 mg (chiefly near alpha=1), and
implicit evolution is −82.095139 mg (chiefly interior alpha). The corresponding
raw-target bound excess is only 0.000110285 mg, but coupled phase-sum/pressure
feasibility remains unproved. No source addback or physical correction is approved.

Read [the completed probe and reopened model-selection decision](SC1_MODEL_SELECTION_AND_PROBE_RESULTS_2026-11-10.md)
first. The archive/wall review and numerical source-split investigation above
are complete. The corrected private candidate passes all eight small-case
regressions; do not repeat them. The first quiescent MEML startup completed
with passing mass budgets but excessive velocities. Its drag-blending
cancellation is now diagnosed. Run the linked fresh `prepare-drag` comparison
before any longer quiescent interval or stirred flow. The
aerated 9.9-s reactor/oxygen/biology checkpoint remains unvalidated.
Do not repeat the completed probe. The historical two-field excess inventory
does not establish that phase count caused it; t=0 and boundary/source budgets
were not established. Retain three fields only for a demonstrated physical
need that the simpler formulation cannot satisfy.

**Preceding numerical checks:** Both cases in
`/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z` are complete. The half-dt
case completed 4,000 steps/8,000 alpha calls with status 0 over the identical
0.1 s interval. Halving dt worsened the liquid mass-budget shortfall from
181.224915 to **319.010443 mg** (+76.03%); actual liquid mass decreases
14.272415 mg despite positive net source. Final-outer source splitting
contributes −179.906786 mg and implicit evolution −135.674459 mg. Alpha
clamping is negligible. The prior four-outer case left a 154.384700 mg
shortfall and substantial iteration-to-iteration variation. Neither numerical
setting is an accepted fix. Mean OTR remains about 96.8 mg/s and oxygen
clipping cancels about 80.57% of requested removal.

The [half-timestep result and source-feasibility implementation](sc1_alpha_half_dt_results_2026-11-10.md)
documents that stage. The alpha/species accounting closes; the unresolved issue is the
source/flux/pressure formulation's physical consistency and coupling, not an
unattributed ledger error. No experimental closure is approved.

**Preserve the completed observer build.** The source-probe
build `/scratch/nwimer/sc1_alpha_source_probe_20261110T191459Z` compiled and
passed runtime verification, including unchanged installed library/original
sources and private linkage of the solver and trial extensions. The
[build/run record](sc1_alpha_source_probe_run_2026-11-10.md) is now historical;
do not repeat its command. The launcher pins
the successful binaries, copies the original manifest checkpoint into one
fresh case, retains all baseline physics/numerics/biology, enables only the
observer and checks legacy-record regression. All **162 local SC1 tests pass**.
The observer classifies source deficits and whether restoring the raw source
at fixed flux would violate phase bounds. Individual feasibility is necessary
but not proof of a conservative coupled correction.

**External-reference follow-up:** The [OpenFOAM benchmark and literature review](SC1_EXTERNAL_BENCHMARKS_2026-11-10.md)
identifies pinned upstream BiRD OF13 MRF/rotating-mesh oxygen/PBM tutorials,
newer HZDR code/case archives, stirred-tank PBM validation and free-surface
experiments. Newer alpha source code retains the capped Su/Sp structure, so a
version upgrade is not an established fix. The closest hybrid paper explicitly
left gas morphology transfer and froth to future work. Use the references to
guide isolated numerical and hydrodynamic benchmarks; none validates SC1 or
requires three fields. No code or simulation change accompanies
this research follow-up.

Do not rerun the completed `four_outer`, `half_dt`, or either case in
`/scratch/nwimer/sc1_alpha_pair_20261110T155635Z`. Preserve the existing verified
build `/scratch/nwimer/sc1_alpha_audit_20261107T001119Z` and the new verified
source-probe build and completed probe control. Its new measurements and
observational regression are recorded above. No further blind dt/outer/coefficient sweep
is requested. Older launch blocks are historical; allocation availability must
be checked rather than inferred from transcripts.

The earlier [full conservation pair](sc1_conservation_pair_results_2026-11-06.md)
remains authoritative for the 0.5 s normal/control comparison: own-duration
mean OTR improves only 0.353%, and normal's native final bubble inventory is
125.23 mg lower. The diagnostic continuations retain frozen OUR; they are not
the intended nonreacting 9.9 s checkpoint. Neither experimental closure is
approved for production. Historical follow-ups below record earlier gates;
the latest state above supersedes their pending-action instructions.

> **Instruction to the next assistant:** Read this handoff and the linked audit documents before changing anything. Review the actual BiRD source and trial outputs before forming new conclusions. Make **every code or documentation change as a Git commit pushed to the GitHub feature branch**. Do not leave work as an uncommitted patch, and do not modify the production CFD runs. Do not rerun already completed experiments without a reason.

**Conservation instrumentation follow-up (2026-11-06):** An opt-in timestep
ledger, final-PIMPLE O2 equation/source capture, PB pre/post-correction moments,
fresh solver/case preparers and log comparison tool are now implemented. Nine
local regression tests pass, including compilation/execution of the independent
C++ accumulator. **The corrected audit at `6ac7e36` now compiles on Kestrel;
all nine diagnostic and ten analytical checks pass there, and executable
dependencies resolve. There are no new control/active balance logs or
conservation conclusions.**
Read [the timestep audit design and build gate](sc1_timestep_conservation_audit.md)
next. This work explicitly distinguishes timestep quadrature from fully
discrete closure; full alpha/pressure/PB equation-level closure remains open.
All existing completed experiments and the biological schedule are preserved.

**Reliability investigation follow-up (2026-11-06):** Read the
[source review and recommended validation sequence](SC1_RELIABILITY_INVESTIGATION_2026-11-06.md).
It identifies omitted simultaneous gas counterexchange in the net-only model,
the need to validate deposition as free-surface disengagement, time-level and
clipping contributions to conservation residuals, and exact configuration
checks for capillarity and diffusive-transfer blending. Ten additional
analytical counterexample tests pass; they do not validate the OpenFOAM
implementation. The diagnostic stream API and terminal-interval predicate
were corrected. No physical closure or production case was changed, and no
new Kestrel simulation or experimental comparison has been performed.

**First build follow-up (2026-11-06):** The uploaded transcript confirms the
seven then-existing diagnostic checks and ten analytical checks passed on
Kestrel, but compilation failed at audit diffusion access. The corrected
lookup uses the registered phase transport model, with its required include
and link flags appended to the preserved local build options. Two additional
local checks pass. Preserve the failed scratch tree
`/scratch/nwimer/sc1_conservation_audit_20261106T182403Z`; rebuild in a fresh
tree. See the [build record](sc1_timestep_conservation_audit.md#first-kestrel-build-and-diffusion-access-correction-2026-11-06).

**Successful build follow-up (2026-11-06):** The next uploaded transcript
confirms `wmake` success at commit `6ac7e36`. Use the existing executable
`/scratch/nwimer/sc1_conservation_audit_20261106T183120Z/bin/putidaFoamSC1ConservationAudit`;
no further rebuild is currently needed. `ldd` reports no missing executable
dependencies, including the added transport library. This does not inspect
experimental libraries loaded later by `controlDict`. Next prepare a fresh
matched pair and review case settings, runtime assets and loaded libraries
before running the zero control. The matched audit retains the original
31.199982372783058 state and biology settings; it is not yet the nonreacting
9.9 s release checkpoint.

**Checkpoint comparison follow-up (2026-11-06):** The first matched-case
preparation stopped on a raw input mismatch in
`processors52/31.199982372783058`, after the collated mesh comparison passed.
No destination/cases were created and no CFD launched. The mismatch has not
been classified as metadata, file-set or physical-state differences. Use
`prepare_sc1_conservation_pair.py --inspect-inputs` (without `--out`) to obtain
file-level hashes, sizes and bounded byte evidence. The equality requirement
remains unchanged; do not strip headers, discard fields or silently replace
either checkpoint. See the [inspection instructions](sc1_timestep_conservation_audit.md#checkpoint-mismatch-before-case-preparation-2026-11-06).

**Mismatch narrowed (2026-11-06):** The uploaded inspection reports 36/36
collated mesh files identical, 55/59 checkpoint files identical, all 12 system
files identical, and only the expected `phaseProperties` difference among
51 constant files. The four changed checkpoint files are `Co.liquid`,
`Co.bubbleGas`, `Co.headspaceGas` and `Co.mixture`; each normal file is five
bytes longer and first differs in its outer-header run-directory path.
Their remaining bytes have not yet been independently verified on Kestrel.
The new opt-in `--allow-courant-location-metadata` proves equality of the
complete files except that one validated outer-header value, or refuses
preparation. Only after proof, it copies the entire control checkpoint to
both NEW cases and verifies every raw file hash. No original is altered and
no field is omitted. Sixteen new regression tests pass locally. The successful
solver executable is unchanged; no rebuild is needed. A new preflight shell
script replaces the long pasted collection block.

**Matched pair ready (2026-11-06):** The next Kestrel transcript confirms all
16 metadata regression tests passed and the four actual Co files differ only
in the proved outer `location` value. Both new cases contain the complete,
raw-identical control checkpoint and mesh. Use the existing pair
`/scratch/nwimer/sc1_conservation_pair_20261106T201549Z`; its `preflight.txt`
has been reviewed. The eight inspected dependency sources match the repository,
and the existing executable and experimental libraries have recorded hashes
and resolved dependencies. Preserve the intentional local biological/pH solver
differences. The next action is the **52-rank zero-control launch**, with a
runtime-library check for `libalphaLimitVelocity.so`; no rebuild or new pair
is needed. No conservation simulation has launched yet.

The actual case retains the saved-OUR explicit oxygen sink (`useMonadSink
false`), 200-CFD-second biology triggers, four-hour biological update size and
0.1-CFD-second pH triggers. Thus this 31.2–31.7 s comparison is a frozen-OUR
diagnostic, not a nonreacting checkpoint. No new biological update is due
during it. PB diameter drift is disabled; diffusion blending uses the default
0.3/0.7 thresholds; foam area multiplier is zero. Interpret balance residuals
with the configured six alpha subcycles, species residual-alpha stabilization
and final-iteration-only PB hooks in mind. See the
[verified configuration and launch instructions](sc1_timestep_conservation_audit.md#matched-pair-preflight-reviewed-2026-11-06).

**Zero-control audit completed (2026-11-06):** The new control in
`/scratch/nwimer/sc1_conservation_pair_20261106T201549Z/control` completed
with status 0, 10,279 steps and final time `31.6999975003429`. The final
bubble inventory and four transfer/closure diagnostics reproduce the preserved
control at displayed precision. The uploaded excerpt includes the cumulative
and final-window ledger, but not the startup/all-write records.

Over 0.500015 s the liquid source requested 248.931 mg O2 removal, transfer
added 49.047 mg, and final-iteration species correction restored 200.123 mg
(80.4% of the request). This is consistent with negative-O2 clipping in the
reviewed OF9 `normalise()` path, not physical oxygen supply. The net accounted
removal after correction is not independently validated biological uptake.
Liquid O2 residual after including that correction is only 0.00215 mg, but
bubble/headspace mass residuals are +386.461/−394.800 mg and liquid mass
residual is −1168.686 mg. Do not call the case conservative or physically
validated. The cumulative PB mass-moment correction −7.0545 g is repeated
normalization of size fractions, not measured loss of actual phase mass.

Next complete the already prepared **normal case unchanged**, using the same
audit executable, while collecting both full audit-record excerpts and old
control regression records. No rebuild, coefficient change or new biology
schedule is justified. Then target unresolved alpha/source stages,
residual-alpha species stabilization and pressure-density changes before
interpreting entrainment benefits. Read the
[control findings and next launch](sc1_timestep_conservation_audit.md#first-control-conservation-result-2026-11-06).

## 1. What we are trying to accomplish

Establish a **physically credible, mass/species-conserving, nonreacting hydrodynamic and oxygen-transfer state at 9.9 seconds using the simplest adequate phase representation**, suitable as a common baseline for subsequent **MM, EMM, and MEML** biological-model comparisons. Three fields are optional and require evidence of need. Preserve the physical liquid inventory, represent the upper-impeller/free-surface interaction and any experimentally consequential froth/entrainment, and predict credible gas holdup and O₂ transfer **without forcing a target DO by tuning arbitrary coefficients**.

The user's biological solver intentionally advances biology at an equivalent **four-hour interval per biology update** and holds OUR fixed over many CFD seconds. **Do not treat that scheduling as a bug or silently change the biological update scheme.** The immediate task is to repair/validate hydrodynamics and transfer before using any new state as a biological checkpoint.

Why this is urgent: historical reacting runs near 30 s show approximately:

| Model | Liquid inventory | Dissolved O₂ | Interpretation |
|---|---:|---:|---|
| Original three-phase MEML | ~999.5 L | ~0.0193 mg/L | Liquid inventory credible; oxygen severely depleted |
| Two-phase EMM | ~1137 L | ~4.61 mg/L | More dissolved oxygen, but unphysical excess liquid inventory |
| MAX-OUR comparison | ~760 L | ~5.92 mg/L | Not a physically comparable liquid inventory |

Those are historical diagnostics, **not matched-condition validation data**. Original MEML biological OUR is ~**0.498 g O₂/s**, versus ~0.504 g/s EMM. The oxygen limitation is a transport/physics question, not justification to weaken biology.

## 2. Confirmed original-model physics and findings

The original MEML uses `interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem` with three Eulerian phases: `liquid`, `bubbleGas` and `headspaceGas`. Its original `constant/phaseProperties` defines Henry-law O₂ composition and Higbie/spherical diffusion **only between liquid and bubbleGas**. It has **no headspaceGas–liquid Henry/diffusion pair**, and no explicit reverse `headspaceGas -> bubbleGas` entrainment mechanism. Bubble-to-headspace `phaseTransfer` deposition is present, but gas-phase deposition is **not oxygen dissolution**.

Original diagnostics: bubble-gas inventory about 7.6 L at ~30 s, with ~99% below 0.8 m and ~1.08 s inventory/flux residence proxy; sparger gas injection about **8.743 g/s**, with blended bubble-to-headspace deposition typically about **7–12 g/s**. Original bubble-to-liquid O₂ transfer ~**0.084–0.094 g/s**, far below biological OUR. Prior EMM net bubble dissolution was ~0.312 g/s. The original MEML surface forms a substantial central depression: roughly y~0.7 m inboard versus ~1.19 m near the rim. This was already evident before reactions (9.9 s).

At 9.9 s the upper impeller's area-weighted phase fractions were approximately **74.4% headspaceGas, 25.4% liquid, 0.2% bubbleGas**; the lower impeller was substantially liquid-wetted. **Do not describe the upper impeller as physically “dry.”** The user reports a foamy, bubbly upper-impeller region in experimental video; the actual video has **not been quantitatively analyzed**. A gas-rich Eulerian cell can still represent unresolved wetted froth.

The deposition term is **blended** by BiRD. The independent, unblended correlation is much larger than the actual assembled pair source. The diagnostic `orderedDepLoss_kg_s=0` is an ordered-key audit artifact and does not prove deposition is absent. Always use **signed, blended `fluid.dmdtf` and assembled phase sources** for conservation claims.

The open top `outlet` has strongly opposing gross headspace flows: near 31.3 s, gross inward ~**1.14072 kg/s**, gross outward ~**1.14847 kg/s**, signed net outward only **0.00774726 kg/s**. All three O₂ fields have `inletOutlet` composition: headspace and bubble reverse-flow O₂ mass fraction `0.23`; liquid `8e-6`. That is plausible for returning air and is **not an obvious O₂ inlet specification error**. The large bidirectional boundary flow still requires independent flow/pressure validation. A reported `headspaceBoundaryOut_kg_s` is **signed net**, not gross outward flux.

## 3. Implemented code and trial findings

### A. Separate headspace-to-liquid oxygen-transfer closure

The GitHub subtree `OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceFilmModels/` implements an optional runtime `resolvedHeadspaceFilm` model and the helper `scripts/diagnostics/prepare_sc1_headspace_o2_trial.py`. It adds conservative coupled Henry O₂ exchange across `headspaceGas`/`liquid`, retaining original `bubbleGas`/`liquid` transfer. The model uses resolved `|grad(alphaH)|` and alpha gates; optional localized foam subgrid area is exploratory. The 10 μm film thickness, foam diameter 3 mm and foam-area multiplier are **not experimentally calibrated**.

Paired independent 0.5 s runs from the same source checkpoint:

| Trial | Effective headspace–liquid area, five-write mean | Headspace-to-liquid O₂ | Bubble-to-liquid O₂ | Total modeled O₂ |
|---|---:|---:|---:|---:|
| `run_geometric_05`, foam factor 0 | ~1.72 m² | 0.003644 g/s | 0.094754 g/s | **0.098398 g/s** |
| `run_foam_05`, foam factor 1 | ~5.58 m² | 0.011627 g/s | 0.094689 g/s | **0.106316 g/s** |

These five-write arithmetic means are **not equilibrated, continuous time integrals**. Foam increased headspace O₂ transfer ~3.19× but raised total O₂ delivery by only ~0.0079 g/s. The model responds to its interface area, but that does **not** validate the area. Total supply is only ~20–21% of the ~0.498 g/s historical OUR. Do not fit area to an extrapolated “required” ~190 m² as if it were measured.

### B. Separate conservative headspace-to-bubble entrainment and PB nucleation

The GitHub subtree `OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceEntrainmentModels/` contains:

- `headspaceEntrainingDeposition.{H,C}`: preserves the original deposition correlation, plus an optional bounded countervailing headspace→bubble phase-transfer source, implemented in the *same original ordered pair* rather than adding an unexamined opposite ordered model.
- `headspaceEntrainmentNucleation.{H,C}`: assigns **local positive net** bubble-phase mass gain after deposition/blending to the existing population balance, through `popBal.eta` and a trial **2 mm** nucleation diameter.
- `scripts/diagnostics/prepare_sc1_headspace_entrainment_trial.py` and `test_sc1_headspace_entrainment_trial.py`: non-destructive scratch preparation and regression checks.

The new experimental library `libsc1HeadspaceEntrainmentModels.so` **successfully compiled and linked** on Kestrel. The diagnostic solver `putidaFoamMask01BubbleAudit` was rebuilt in scratch. All completed 52-rank runs reported status 0. The **zero-entrainment** `run_control_05` reproduced the original `run_geometric_05` at displayed precision for inventory, deposition, both O₂ transfer channels, and zero PB nucleation: the regression control **passed**.

The initial `legacyDownwardLiquid` criterion `U.liquid.y < 0` almost completely suppressed entrainment: 0–7 eligible cells among ~4,000 mixed-interface cells. Gate diagnostics showed ~0.33 m² interface in the selected upper-impeller window but only ~0–6.4×10⁻⁵ m² after the downward-liquid filter. This was an *eligibility-definition problem*, not proof that true entrainment is physically small.

The explicit opt-in `normalGasIntoLiquid` criterion instead accepts positive `(U_liquid-U_headspaceGas) · grad(alpha.headspaceGas)/|grad(alpha.headspaceGas)|`. It is a **relative-motion proxy**, not a validated breakup/detachment law. The default remains the original legacy mode. Both criteria use the same axial/radial window, phase-fraction gates and donor limitations. **The normal mode is experimental and must not be promoted into production.**

### C. Latest 0.5 s normal-inward entrainment trial — completed 2026-11-06

`run_normalactive_05`: same collated starting state, `entrainmentEfficiency 0.01`, `entrainmentVelocityCriterion normalGasIntoLiquid`, raw global cap `0.002 kg/s` (**2 g/s**), `nucleationDiameter 0.002 m`, headspace film `foamAreaMultiplier 0`. The five outputs at t=31.3–31.7 s show:

- 1,791–1,900 eligible interface cells, ~0.160–0.169 m² area;
- proposed **raw unblended entrainment** ~2.55–2.69 g/s, capped to 2.0 g/s (scale ~0.743–0.784);
- **local positive net bubble mass birth entering PB** ~0.208–0.238 g/s, five-write average ~**0.223 g/s**;
- **integrated signed bubble↔headspace exchange remains net deposition**, ~7–10 g/s. Do not call the 2 g/s cap “actual bubbles created,” and do not call signed net deposition “gross entrainment”;
- PB size-group sums report `1 1 1`; diameter remains finite, mean ~1.205 mm, minimum ~1.003 mm, maximum ~5.331 mm at the end. **BiRD clamps/renormalizes** size fractions, so this is a numerical sanity check, *not* an independent PB moment-conservation proof.

| t (s) | Bubble inventory change vs zero control | Added total O₂ transfer vs control |
|---|---:|---:|
| 31.3 | +11.09 mg | +0.282 mg/s |
| 31.4 | +24.17 mg | +0.401 mg/s |
| 31.5 | +55.87 mg | +0.558 mg/s |
| 31.6 | −4.12 mg | +0.311 mg/s |
| 31.7 | **−125.23 mg** | +0.372 mg/s |

Five-write mean modeled O₂ transfer: **0.098782756 g/s active** versus **0.098398236 g/s control**, gain only **0.000384520 g/s = 0.385 mg/s** (~0.39%). The modified case still supplies only ~19.84% of historical frozen OUR. **No sustained bubble-holdup or meaningful O₂ improvement has been demonstrated.**

### D. Unresolved conservation problem

The existing `bubbleMassAudit.H` reports **endpoint** inventories, injection, signed blended deposition, net gas–liquid transfer, and physical patch fluxes only at write times. Variable CFD timesteps and highly variable deposition make sparse 0.1 s rates insufficient to certify mass conservation. The t31.4–31.7 endpoint closure residuals are nonzero for **both** normal-active and zero-control runs. A coarse trapezoidal balance over 31.3–31.7 s misses the observed bubble inventory change by about **+0.231 g active** and **+0.382 g zero control**. This discrepancy is **not yet evidence of a real numerical leak**: the budget has not been constructed from the same timestep-discrete fluxes and source evaluations as the solver.

A further open question is whether donor O₂ is transferred and conserved correctly when headspace gas is reclassified as 2 mm bubbles, given the strong gross counterflow through the open top and the competing dissolution pathways.

## 4. File paths and working environment

All Kestrel paths below are **on the user's cluster** and are not directly available to a new ChatGPT session without connected access or user execution.

~~~bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
BIRD=$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam
BASE=/projects/microenviro/final_runs/reacting_Julia_08282026/meml/putida_sc1_final
SOURCE=/scratch/nwimer/sc1_bubble_audit_20261029/run
SOLVER=/scratch/nwimer/sc1_bubble_audit_20261029/solver
EXE=$FOAM_USER_APPBIN/putidaFoamMask01BubbleAudit

FILMROOT=/scratch/nwimer/sc1_headspace_o2_trial_20261103
FILMLIB=$FILMROOT/lib/libsc1HeadspaceFilmModels.so
ENTROOT=/scratch/nwimer/sc1_headspace_entrainment_20261104
ENTLIB=$ENTROOT/lib/libsc1HeadspaceEntrainmentModels.so
PY=/scratch/nwimer/.conda-envs/of-post/bin/python

source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc"
module load libfabric
~~~

`/scratch` may resolve to `/kfs3/scratch`. The common source checkpoint in the experimental collated `processors52` directory is **`31.199982372783058`**. Runs use **52 MPI ranks** with `srun -n 52 "$EXE" -parallel -fileHandler collated`. Current allocation from 2026-11-06: job **19721806**, node `x1006c3s5b1n1`; this is transient—**never assume a prior allocation is still active**.

Preserve existing scratch trials:
- `$FILMROOT/run_geometric_05`, `$FILMROOT/run_foam_05` (logs `log.headspaceO2.smoke`);
- `$ENTROOT/run_control_05` (log `log.entrainment.control`);
- `$ENTROOT/run_entrain_05` (log `log.entrainment.smoke`, old downward gate);
- `$ENTROOT/run_gateaudit_05` (log `log.gateaudit`);
- `$ENTROOT/run_normalgate_02` (log `log.normalgate`);
- `$ENTROOT/run_normalactive_02` (log `log.normalactive`);
- **`$ENTROOT/run_normalactive_05` (log `log.normalactive05`, latest key experiment).**

All prior named runs completed; **do not overwrite them, silently continue them, or edit their outputs**. Create fresh, uniquely named scratch clones for new experiments.

## 5. Precisely what the next session should do

**Priority 1 — Build a timestep-resolved, conservative phase/species audit in the scratch diagnostic solver.** Before writing new physics, read `OFsolvers/birdmultiphaseEulerFoam/putidaFoam/putidaFoam.C`, `bubbleMassAudit.H`, and BiRD's `PhaseTransferPhaseSystem` and `PopulationBalancePhaseSystem`. The current audit includes a post-PIMPLE, `runTime.writeTime()`-gated include; this is the likely place to extend, **but the new accumulation must run every time step and print only on selected writes**. Do not run the expensive diagnostic multiple times per PIMPLE correction.

Plan to track:
1. Whole-domain **liquid, bubbleGas, headspaceGas and combined gas** mass inventories, with start/end differences. Log physical-boundary *gross inward and gross outward* separately, especially the pressure outlet.
2. Per-timestep integrated **sparger gas injection**, correctly signed and blended gas–gas transfer, bubble–liquid and headspace–liquid O₂ transfer, and any other phase sources. Use the **actual assembled phase source and flux** and avoid double counting model terms. Accumulate source/flux integrals using the actual variable `deltaT`; distinguish an approximate endpoint/trapezoid quadrature from a provably identical solver-discrete budget.
3. Phase-resolved **O₂ inventories** `∫alpha_phase rho_phase Y_O2,phase dV`, liquid biological O₂ consumption, gas boundary O₂ in/out including inflow `inletOutlet` composition, and coupled O₂-transfer terms. Clarify whether all transported phase species are available and the role of species diffusion/advection.
4. Local positive **net** PB bubble birth and whether it is consistent with changes in size-group moments *before renormalization*, without asserting that `sum(fi)=1` proves conservation.
5. At each saved comparison time, report **mass and oxygen conservation residuals** alongside total liquid inventory, bubble holdup, transfer rates, DO and the free-surface state.

The audit must be **disabled by default**, run only in an isolated scratch diagnostic executable/case, and unit/regression tested where possible. Perform a **zero-entrainment control first**, then the matching normal-inward trial from the **same checkpoint**, using fresh scratch paths. Do not relabel the old 0.1 s endpoint closure as a solver-discrete balance.

**Priority 2 — Validate the actual free surface and impeller physics.** Re-examine the central cavity, upper-impeller gas and liquid contact, gas recirculation at the outlet, and resolved versus subgrid interface area. Obtain/inspect the actual experimental video and any measured gas holdup, bubble size, impeller wetting, free-surface shape and `k_L a`/OTR if available. Distinguish coherent headspace cavitation, dispersed bubbles and liquid-wetted froth. A positive local normal relative velocity is **not sufficient** for bubble detachment.

**Priority 3 — Decide whether the three-phase closure deserves development.** The film model is functionally active but far below OUR; the normal-inward entrainment model is numerically active but raises OTR by only ~0.385 mg/s. Do **not** blindly increase the entrainment coefficient, raw cap, 2 mm nucleus count or foam area to force DO. A different physical closure or even a revised hydrodynamic representation may be required. Preserve trustworthy features of MEML (notably liquid volume) and benchmark against EMM/MAX-OUR without mistaking their nonphysical liquid inventories for correct behavior.

**Priority 4 — Only after physics, closure and conservation checks pass:** generate/re-equilibrate the **nonreacting 9.9 s MEML checkpoint** and launch controlled, common-initial-state MM/EMM/MEML biology comparisons. Preserve the intended biology-update schedule.

**Deliverables for a successful next session:** (a) a pushed GitHub commit with a verified time-resolved audit, (b) fresh scratch-only control and active balance logs, (c) an interpretation of gas/O₂ and PB conservation, (d) an evidence-based choice of the next physical model or calibration data—not another speculative coefficient sweep.

## 6. Collaboration, execution and safety rules

- **Always commit and push every code/document change** to `codex/sc1-headspace-bubble-audit-20261029`. Use the existing GitHub connected tools or an explicitly confirmed `git push` from the user's Kestrel checkout; ensure the pushed commit is visible. The user previously specifically objected to changes without pushes.
- **Do not modify production cases** in `/projects/microenviro/final_runs`, original solver files in `$BIRD`, or already completed scratch cases. Make scratch-only source/case clones; compile experimental models into separate scratch libraries.
- **Never paste `exit`, `logout`, `exec`, top-level `set -e` / `set -euo pipefail`, or shell-exit traps into the interactive `salloc` shell.** Such commands previously ended valuable allocations. Use guarded `if` blocks and commands that return to the prompt.
- Avoid claiming you can directly inspect Kestrel files from ChatGPT. Ask the user to run **safe, read-only checks** when needed; connected GitHub source and Kestrel scratch runtime are different systems.
- Do not confuse **kg/s, g/s and mg/s**, candidate gross versus actual **net** flux, or total oxygen flux versus gas-mixture exchange.
- Do not interpret an SSH disconnect as a CFD failure if the solver already printed `End` and exit status 0. Do not mistake OpenFOAM's routine `FOAM_SIGFPE` setup line for an exception.
- OpenFOAM's `foamDictionary` can fail when reading full **binary collated field** data; the documented workaround decoded the text `boundaryField` portion read-only. Plain ASCII `constant/phaseProperties` dictionaries work with `foamDictionary`.

## 7. First five minutes in the new session

1. Read **this document**, [the detailed entrainment log](sc1_headspace_reentrainment_trial.md), [headspace O₂ design/results](sc1_headspace_oxygen_transfer_design.md) and [the historical three-phase audit](sc1_three_phase_headspace_audit.md).
2. Check the **current** GitHub branch head, pending PR and existing source. Do **not** assume the checkout/cluster job is at the same state as the last conversation.
3. Confirm the latest validated result is `run_normalactive_05` status 0 with five writes ending 31.7 s; the zero control is `run_control_05`. Do not re-run those simulations.
4. Start with the **phase and O₂ timestep-integral conservation design**, its sign conventions and scratch-only instrumentation. Push the source/document commits before asking the user to test them.
5. After tests, interpret **both** control and normal-active conservation residuals before deciding whether to modify the physical entrainment model.

**One-sentence baton pass:** We implemented the missing headspace–liquid O₂ pathway **in an isolated scratch prototype** and added capped headspace-to-bubble entrainment coupled to the PB, but the latter does not materially increase O₂ transfer or sustained bubble holdup; next we must establish trustworthy timestep-resolved phase/O₂ conservation and validate foam/entrainment physics before advancing toward a common nonreacting 9.9 s checkpoint.
