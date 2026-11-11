# SC1 isolated two-field conservation test

## Status and purpose

**All four cases pass the complete 100-step numerical screen.**
Saved-log analysis at `e7825bd` confirmed that the earlier report rejection
was caused by timestamp roundoff. The existing case root is
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T233431Z` and the verified
build is `/kfs3/scratch/nwimer/sc1_two_field_build_20261110T224950Z`.
The continuation at `56a84d7` saved the no-reference result without rerunning
CFD, then completed the gas-reference partner. The source pair subsequently
completed at `ef74cf6`, also with both screens passing. All four `result.json`
files now exist with log hashes. **This suite is complete; preserve it and do
not rerun any case or rebuild.** Next implementation is a separate paired
source test near phase-fraction bounds, followed by controlled timestep and
interface checks. No new stress-test or SC1 reactor launch is included here.
This is a numerical test harness, not a corrected SC1 reactor
case or a demonstrated solution to wall climbing, excess inventory or oxygen
clipping. The phase count remains open pending physical validation.

The prescribed-fill check showed that the archived EMM prescription is near
999 L while its 9.9 s field contains 1151 L. Even the upper cell-selection
bound is about 141 L below the observed inventory. That did not recover the
historical time-zero state or identify the missing mass-budget contributions.
The next useful step is a controlled comparison that measures those terms.

This harness uses the repository's nonreacting `multiphaseEulerFoam` main,
with observers in a **private copy** of its phase-system and pressure code.
It does not copy `putidaFoam`, rename three-field inputs or change the
intentional biological schedule. It does not use the existing three-field
observer executable. All six required BiRD libraries are rebuilt privately
from the same repository source snapshot, avoiding a mixture of new headers
and unverified installed project-library implementations.

## Four cases, run in two stages

| Pair | Domain and forcing | Initial fractions | Purpose |
|---|---|---|---|
| `closed_no_reference`, `closed_gas_reference` | Periodic 0.1 m cube; both phases initially translate at 0.1 m/s; no external/interphase sources | Two equal slabs with alphaL=0.2 and 0.8; gas is complementary | Nontrivial transport with no external mass supply or loss |
| `source_no_reference`, `source_gas_reference` | Same mesh; stationary walls and one fixed-pressure outlet; uniform liquid mass source 0.001 kg/s | Uniform alphaL=alphaG=0.5 | Known source with pressure-compatible outflow, measured phase fluxes and sources |

Within each pair the **only input difference** is `referencePhase gas;`.
Both pairs use 16 x 4 x 4 cells, fixed dt=0.001 s, 100 steps to 0.1 s, six
alpha subcycles, one alpha correction and two PIMPLE outer iterations. The
fixtures use laminar, pure isothermal phases with constant densities 1 and
1000 kg/m3, zero gravity and zero surface tension. There is no MRF, PBM,
interphase transfer or biology. Drag uses the same symmetric dispersed-pair
settings in both variants. These fixture settings are not recommendations
for the reactor.

Both measured initial phase volumes must be 0.0005 m3. The corresponding
initial masses are 0.0005 kg gas and 0.5 kg liquid. The source pair's prescribed
liquid addition is 0.0001 kg over 0.1 s; the retained liquid need not increase
by that amount because the outlet carries mass away.

The source pair cannot launch until both closed tests complete and pass the
numerical screen. A successful process exit alone is not conservation evidence.

## What is measured

`sc1PhaseFractionAudit.H` has dynamically sized phase storage. It reuses the
previously tested Euler/MULES arithmetic and hook locations, and adds a
before/after hook around the **existing** no-reference rescaling operation.
Every original phase-equation line remains in order. It reports every alpha
call, rather than combining repeated outer iterations into physical time:

- Actual MULES change, old-time reset, assembled Su/Sp, limited-flux divergence
  and the discrete equation residual, including its cellwise L1 norm.
- Original divergence-based alpha source, stored pressure dilatation, raw
  coupled dilatation, capped-source split and implicit-source response.
- Implicit correction, rescaling, reference reconstruction and final clamping.
- Actual external MULES volume flux and its difference from the integrated
  divergence, plus stored flux/subcycle averaging consistency.
- Unassigned alpha change, with volume and mass weighted by density frozen
  at alpha-solve entry. **Frozen-density terms are not boundary mass fluxes.**

A separate solver-stage observer measures actual alpha*rho inventory changes,
with the exact algebraic decomposition

`delta(alpha*rho) = rho_before*delta(alpha) + alpha_after*delta(rho)`.

Its stages telescope across all outer iterations to the accepted timestep.
The step budget uses actual stored `alphaRhoPhi` at the end of the timestep,
and source values captured from the final pressure assembly. The original
`fvModels.source(alpha,rho) & rho` expression is evaluated once, cached for
observation, then used in its original pressure-equation location. Existing
interphase source fields are observed without calling their model again.
This is end-of-step Euler source/flux quadrature; a residual remains a finding,
not something that the observer removes.

The first fixtures have constant density, so their density-change terms should
be zero. They do not validate compressible gas, interphase exchange, PBM or
species equations. No claim of numerical solution equivalence to an unobserved
binary has yet been established by execution.

## Screening and evidence integrity

The analyzer requires the complete 100-step stream, both phases, both outer
iterations, every required source/change term and each stage. It rejects
missing/duplicate records, nonfinite values, wrong source integrals, wrong
initial inventory and unsupported timestep/iteration counts. It independently
checks the main accounting identities and agreement between alpha and solver
stage observations.

For record indexing only, timestamps within 1e-12 s of the prescribed grid are
matched to that grid. Raw values remain intact, and duplicate detection runs
on the matched keys. Off-grid times and missing/extra records still fail.
Reports expose the number of adjusted keys and maximum adjustment. This
timestamp representation tolerance does not change the physical tolerances below.

The following are declared **numerical screening tolerances**, not experimental
acceptance criteria:

- Alpha within [-1e-10, 1+1e-10]; maximum internal phase-sum error <=1e-10.
- Unassigned alpha volume <1e-13 m3 and frozen-density mass <1e-10 kg.
- Closed-pair cumulative absolute mass residual <=1e-12 kg + 1e-9 times initial
  phase mass. Actual external flux must be negligible and alpha must actually
  move, so a stationary do-nothing calculation is insufficient.
- Source-pair cumulative absolute mass residual <=1e-10 kg + 0.001 times the
  larger of cumulative phase supply and net exported phase mass magnitude.

Absolute timestep residuals are summed before comparison; positive and
negative errors cannot cancel to create a pass. Reported final-outer rescaling
and source-split sums remain iteration diagnostics, not additional physical
mass added to the accepted-step budget.

Preparation hashes the source copies and installed project libraries. Make
outputs are rewritten to fresh private destinations, and build commands run
with child-only private library/app paths. Verification checks source hashes,
unchanged installed libraries, all private binaries, `ldd` dependency paths
and a mandatory observer symbol. Later runs verify the saved binary hashes.

Case preparation proves the paired input difference and saves hashes. Both
members are checked before a pair starts. Any attempted case is refused;
logs are created exclusively. Mesh generation, mesh checks and initialization
occur only in the new cases. The initialized input files are hashed before
CFD and compared afterward. Completed/failed cases are preserved.

## First Kestrel build: assignment compatibility correction

The user ran the build at commit `a0f19295d7a399de16ae2e73579850d3bcd73c46`
in `/kfs3/scratch/nwimer/sc1_two_field_build_20261110T222049Z`.
The first two library commands completed; the third, `phaseSystems`, failed
with status 2 while compiling `phaseSystemSolve.C`. The transcript identifies
this observer initialization:

```cpp
base_[i]=dgdt_[i]=requested_[i]=scalarField(phase.size(),0);
```

OpenFOAM 9's `Field::operator=` returns `void`. The assignment chain therefore
attempted to assign a `void` result to another field. This is an assistant-
introduced compile error, not evidence about the physical model or conservation.
The three fields are now initialized in separate statements. No equations,
case settings or biological scheduling were changed by the correction.

The new targeted C++ regression extracts the actual initialization statements
and compiles/runs them against a field type with that same void-assignment
contract. It also confirms that the previous chain is rejected by the compiler.
This is an API-contract regression, **not an OpenFOAM compilation**. Both new
observer headers were checked for additional chained assignments. Failure
reporting now prints the first compiler errors as well as the final 60 log
lines, preserving the full log on disk.

No CFD launched, no case was prepared, and full build/runtime verification did
not run. Preserve the failed directory and its `log.build.phaseSystems`.
The fresh corrected build described below has now succeeded. Do not patch
or rerun the failed copy.

## Corrected Kestrel build completed

The user's complete terminal transcript records a successful SSH pull to
`a5989c12737567e3afe706e537d358a1907244e1`, six successful `wmake libso`
commands, one successful `wmake`, and:

`SC1 two-field private build and runtime verification: PASS (no CFD launched)`

The verified build is
`/kfs3/scratch/nwimer/sc1_two_field_build_20261110T224950Z`.
The verifier checks the original/private source hashes, unchanged installed
project libraries, private dependency resolution, mandatory observer symbols
and executable availability. Its saved `build_verified.json` is checked again
when preparing and launching cases. The following hashes are transcribed from
the user's output; this session has not independently downloaded the binaries.

| Artifact relative to the verified build | SHA-256 |
|---|---|
| `bin/birdTwoFieldConservation` | `1483c42b8bb3aeb37be743eafbd70f1b3e59d3ed053388ccd54f7c63a044f1f5` |
| `lib/libbirdeulerianInterfacialCompositionModels.so` | `ab19f4a19b344e30829b1662f88d2c45dac65022f83193b5f4321ab8e2841649` |
| `lib/libbirdeulerianInterfacialModels.so` | `3dfd2332b13b1689fc465bee880e9df9aee4bc0d5f9fff9ac94c446b07b71232` |
| `lib/libbirdphaseSystem.so` | `5b5157c4b970e9272d012b85805f040b22366a9f14a41b511244a2c469b10558` |
| `lib/libbirdmultiphaseMomentumTransportModels.so` | `24b10084e0fb998a7297e5a9eae62157cca6681a4a0e959a264695933e9c2451` |
| `lib/libbirdmultiphaseReactions.so` | `859f19614ad51ab765b1cf8451b10803dfcb1ce7aa3525a903e17ea7deae56a9` |
| `lib/libbirdmultiphaseSystems.so` | `148f86c6cb4bd316c0da8e241d27bb75c6e11eb060f8ad21a593b45f3d14fad6` |

This establishes compilation and static runtime dependency verification, not
case-startup success or conservation. Preserve both the successful build and
the preceding failed build. No rebuild, new solver instrumentation or physics
change is needed before the closed pair.

## First closed-pair attempt: pressure dictionary correction

The first attempt created
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T231637Z`.
`blockMesh`, `checkMesh` and `setFields` completed for `closed_no_reference`.
The solver process then exited with status 1 while constructing its pressure
reference, before the timestep loop. The other closed case and both source
cases were not launched. Preserve this complete case root and its
`closed_no_reference/log.solver`; do not edit or restart it.

The assistant-generated `system/fvSolution/PIMPLE` incorrectly contained
`pMin 1000;`. Foundation OF9's
[pressureReference constructor](https://cpp.openfoam.org/v9/pressureReference_8C_source.html)
explicitly rejects this legacy entry. The case generator now retains
`pRefCell 0; pRefValue 100000;` in PIMPLE and places the pressure floor in
`system/fvConstraints` as:

```foam
limitp
{
    type limitPressure;
    min 1000;
}
```

The [OF9 constraint reader](https://cpp.openfoam.org/v9/fvConstraints_8C_source.html)
reads that system dictionary; the existing BiRD side-sparger tutorial uses the
same format. The compiled pressure equation already calls
`fvConstraints.constrain(p)`. **No solver or library rebuild is needed.**

The generator change affects only `system/fvSolution` and
`system/fvConstraints` in each of the four small cases. Each pair still differs
only in `referencePhase gas;`. A new regression parses the PIMPLE and constraint
entries for every case, rejects all six obsolete pressure/density-limit keys,
and verifies the preserved reference and the intended minimum pressure.
All 27 two-field harness tests passed at that correction. It fixed the observed startup error;
it does not establish successful timestepping or conservation.

## Second closed-pair attempt: phase-specific alpha controls

The retry at `c78199e` created
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T232902Z`.
`closed_no_reference` passed mesh generation, checking and initialization.
The solver successfully read `system/fvConstraints/limitp` with `min 1000`,
confirmed `SC1_TWO_LIBRARY version=1`, and printed initial inventories:

| Phase | Initial volume (m3) | Initial mass (kg) |
|---|---:|---:|
| gas | 0.0005 | 0.0005 |
| liquid | 0.000500000000000001 | 0.5 |

At `Time = 0.001`, PIMPLE iteration 1, it exited with status 1 because
`alpha.gas` was undefined under `system/fvSolution/solvers`. The time label
advanced, but **no timestep completed and no conservation screen passed**.
The other closed case and both source cases were not launched. Preserve this
complete root and the preceding failed root; do not edit or restart either.

The generated entry was the literal `alpha`. This satisfied
`phaseSystemSolve.C`'s `mesh_.solverDict("alpha")`, but not
`MomentumTransferPhaseSystem::implicitPhasePressure`, which requests the phase
field name before the first limiter call. The
[OF9 MULES limiter](https://cpp.openfoam.org/v9/MULESTemplates_8C_source.html)
also requests `mesh.solverDict(psi.name())`. Both consumers require entries
matching `alpha.gas` and `alpha.liquid`, even with implicit phase pressure off.

The correction uses the repository's existing quoted `"alpha.*"` pattern.
It matches `alpha`, `alpha.gas` and `alpha.liquid` with the same controls.
Only this key in `system/fvSolution` changes in each generated case; the
subcycles, corrections, audit switch, flux settings, limiter defaults and
screening tolerances are unchanged. Each pair still differs only by
`referencePhase gas;`. **Reuse the verified build; no compiled source changed.**

All 29 harness tests passed at that correction. New fixture-contract tests resolve the system and
phase field names for all four cases and reproduce the failed phase lookups
with a literal `alpha` key. These restricted Python checks are not native
OpenFOAM dictionary or CFD tests. Source review also checked the selected
pressure/momentum controls and schemes: the pure isothermal laminar fixtures
do not solve energy, species or turbulence transport equations; implicit
phase pressure defaults to false. The subsequent attempt is recorded below.

## Third attempt: solver returned, report coverage rejected

The attempt at `d7b1b2e` created
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T233431Z`.
The transcript reports `blockMesh`, `checkMesh`, `setFields`, the private
solver command, then `Incomplete or unexpected initial/step/alpha-call records`.
In the harness this error occurs only after the solver returned status 0,
initial input hashes remained unchanged, and the report found the unique
library marker and solver `End` without a fatal-error marker. At that point
the final time, record counts and budgets remained unverified. No `result.json`
was written, and the remaining three cases were not launched. The later
saved-log result below resolves the reporting failure.

Review found that record keys used exact floating-point time comparisons.
Changing a complete synthetic stream's `0.072` to `0.0720000000000001`
reproduces the same rejection. That initially demonstrated a parser defect
without proving the Kestrel cause. Record indexing now accepts only
times within 1e-12 s of the prescribed grid (one billionth of a timestep).
It preserves raw values and rejects duplicate matched keys, genuinely off-grid
times, and missing or extra records. Coverage errors now name each affected
record family, expected/observed counts and bounded missing/unexpected examples.
No solver, case input, mass-budget tolerance or required record count changed.

All 34 harness tests passed at that correction, including roundoff-equivalent streams for all four
cases with identical budgets, raw-time preservation, near-duplicate rejection,
off-grid rejection and actionable coverage diagnostics. Preserve all three
attempted roots and the verified build.

## Completed saved-log analysis: closed no-reference passes

The user ran `analyze` at `e7825bd` on the existing
`closed_no_reference/log.solver` and supplied the full JSON report:
`complete: true`, `screen_pass: true`, `errors: []`, `steps: 100`.
This is a user-supplied analysis result; the raw log remains on Kestrel.

| Quantity | Gas | Liquid |
|---|---:|---:|
| Initial mass (kg) | 0.0005 | 0.5 |
| Final mass (kg) | 0.000500000000000061 | 0.499999999999939 |
| Sum of accepted-step mass changes (kg) | 6.110443488027293e-17 | -6.11045889484697e-14 |
| Cumulative absolute physical residual (kg) | 5.769354784863862e-16 | 5.770176179219254e-13 |
| Screening limit (kg) | 1.5e-12 | 5.01e-10 |
| Final-outer rescale sum (m3) | 7.752045533271372e-20 | 1.2858637765678101e-19 |

Reported external boundary flux, option source, transfer source and
final-outer source-split contribution are all zero for both phases. The complete
screen also checks nontrivial alpha transport, bounds, phase sum, stage and
alpha-ledger identities. There were 3,360 timestamp-key adjustments, with a
maximum of `1.1102230246251565e-16` s against the `1e-12` s matching tolerance.
The successful analysis of the same saved run confirms the earlier rejection
was timestamp representation, not missing audit records.

The no-reference formulation conserves mass in this closed, constant-density
transport fixture, with negligible observed rescaling. This does not establish
conservation with sources, explain the archived EMM inventory increase, or
select the SC1 reactor phase formulation. `analyze` initially printed the
result without saving `result.json`; the completed continuation below saved it
after revalidating the evidence.

## Closed pair completed: both reference settings pass

The user ran `continue-closed` at `56a84d7` in the existing
`sc1_two_field_cases_20261110T233431Z` root. The transcript confirms
`REUSED=.../closed_no_reference`, a passing saved first report, then the four
mesh/check/initialization/solver commands for `closed_gas_reference` only.
The partner also reports `complete: true`, `screen_pass: true`, `errors: []`,
100 steps, and the runner ends with `PAIR=closed complete`.

| Quantity | No reference | Gas reference |
|---|---:|---:|
| Gas final mass (kg), from 0.0005 | 0.000500000000000061 | 0.000500000000000093 |
| Liquid final mass (kg), from 0.5 | 0.499999999999939 | 0.499999999999907 |
| Gas cumulative absolute residual (kg); limit 1.5e-12 | 5.769354784863862e-16 | 4.894737643911179e-16 |
| Liquid cumulative absolute residual (kg); limit 5.01e-10 | 5.770176179219254e-13 | 4.895160184386637e-13 |
| Gas final-outer rescale sum (m3) | 7.752045533271372e-20 | 0 |
| Liquid final-outer rescale sum (m3) | 1.2858637765678101e-19 | 0 |
| Adjusted timestamp keys | 3360 | 2128 |

Both reports have zero boundary mass export, option/transfer sources and
final-outer source-split terms. Both maximum timestamp adjustments are
`1.1102230246251565e-16` s. Gas-reference accepted-step changes are
`9.277898424502305e-17` kg gas and `-9.277435852805425e-14` kg liquid.
The differences between variants are far below the declared screening limits;
this closed fixture gives no conservation-based reason to prefer one setting.
It does not explain EMM's archived excess inventory or validate source coupling,
interphase transfer, free surfaces, oxygen, or the reactor phase count.

Recorded solver-log SHA-256 values from the supplied reports:

- `closed_no_reference`: `9793a67274cf9260750c7b822f24bc0e9cf66a00d04f4b0cd707700042c212df`
- `closed_gas_reference`: `1e563cb5179927771ffea8ccf691a89488f71cb871c500dd785117fd7e3e17ac`

These are user-supplied runtime results; the raw logs and hash-bound
`result.json` files remain on Kestrel. Preserve both cases and the earlier
failed roots.

## Source pair completed: all four baseline screens pass

The user ran the source pair at `ef74cf6` in the existing case root and
supplied both complete JSON reports and `PAIR=source complete`. Each source
case reports `complete: true`, `screen_pass: true`, `errors: []`, 100 steps.
There was no rebuild, regeneration or repeat of the closed cases.

| Liquid budget (kg) | No reference | Gas reference |
|---|---:|---:|
| Initial mass | 0.5 | 0.5 |
| Final mass | 0.500090419143175 | 0.500090419112172 |
| Sum of accepted-step mass changes | 9.041914317521618e-05 | 9.041911217155564e-05 |
| Boundary mass exported | 9.580924787716157e-06 | 9.580924769262218e-06 |
| Assembled option source | 9.999999999999999e-05 | 9.999999999999999e-05 |
| Signed physical residual | 6.796293226441216e-11 | 3.694081779173838e-11 |
| Cumulative absolute residual | 6.991205047567109e-11 | 3.8064394580038755e-11 |
| Screening limit | 1.0009999999999999e-07 | 1.0009999999999999e-07 |

The liquid balance is 0.100 g supplied, approximately 0.090419 g retained and
0.009581 g exported. Independently summing the supplied report values reproduces
each signed residual. The cumulative absolute residuals are 0.0000699121% and
0.0000380644% of the supplied liquid mass. They sum absolute timestep defects,
so cancellation is not responsible for the pass.

| Gas budget (kg) | No reference | Gas reference |
|---|---:|---:|
| Initial mass | 0.0005 | 0.0005 |
| Final mass | 0.000499909580856825 | 0.000499909580887829 |
| Sum of accepted-step mass changes | -9.041914317518895e-08 | -9.041911217155635e-08 |
| Boundary mass exported | 9.04191374331792e-08 | 9.041907543275881e-08 |
| Signed physical residual | -5.742009747825391e-15 | -3.6738797539834215e-14 |
| Cumulative absolute residual | 6.122027062821979e-15 | 3.79419149386864e-14 |
| Screening limit | 1.904191374331792e-10 | 1.9041907543275884e-10 |

Gas option sources and both phases' interphase transfer sources are zero.
The no-reference final-outer rescale sums are `3.0999895359040594e-14` m3
gas and `3.100526953236917e-14` m3 liquid; the gas-reference values are zero.
Liquid final-outer source-split sums are `5.759824041329228e-26` and
`-2.117582368135748e-26` m3, respectively, effectively zero. These are
iteration diagnostics, not additional physical source terms to add to the
accepted-step budgets. The difference between the variants' cumulative liquid
mass changes is only `3.100366054e-11` kg. The reference variant has a smaller
liquid residual and a larger gas residual; both are comfortably below their
limits. This is not a basis for declaring either setting generally superior.

Timestamp normalization remains 3360 / 2128 adjusted records, with maximum
adjustment `1.1102230246251565e-16` s in both source cases. Recorded log hashes:

- `source_no_reference`: `b3124354624ce8ea8ff56e5e7663e89d3e36c385f01e0396e723314576fa6d96`
- `source_gas_reference`: `487bec0cc44e837d2de8ae23c075511a0433344afd27736b6d8e14a3ecff7b7d`

The reports are user-supplied runtime evidence. Raw logs, initialized inputs,
manifests and hash-bound reports remain on Kestrel. All four baseline screens
are now complete. Preserve the successful root and verified build, together
with earlier failed attempts.

## Next implementation: test source splitting near phase bounds

The baseline establishes closed transport and a known external liquid-feed
balance for constant-density, isothermal two-field fixtures. It supports
continuing to evaluate the simpler representation. It does not show that
three fields are required, select a reference phase, resolve the archived EMM
excess inventory, or validate wall/free-surface behavior.

The source fixtures began at alphaL=alphaG=0.5 and produced essentially zero
source-split diagnostics. They therefore do not test the behavior that the
earlier three-field probe localized near alphaL=1. The inspected
`phaseSystemSolve.C` divides positive dilatation by `max(1-alpha,1e-4)` and
negative dilatation by `max(alpha,1e-4)`. This mechanism remains present in
the two-field solver; a mid-fraction pass cannot dismiss its effect near a
bound.

The next separate fixture should compare the same reference settings with
a near-pure-liquid initial state, a known physically feasible source and a
pressure-compatible outlet. Its specification must declare measured initial
phase inventories, expected integrated source, phase-specific tolerances
appropriate to the small gas inventory, and output proving the intended
bound regime was exercised. Retain all source-splitting, rescaling and
accepted-step ledgers. Then halve the timestep at the same physical duration
and hold other inputs fixed to distinguish temporal effects; vary outer
iterations separately if needed. Do not change floors or add source back
based only on this baseline.

The current generator and analyzer intentionally encode the baseline's
initial inventory and 100-step contract. **No near-bound launch command is
ready in this update.** Implement and validate a separate case/report contract
before giving a new CFD command; do not hand-edit completed cases or merely
change alpha in the old fixture. Interphase-transfer and resolved-interface
checks still precede a matched SC1 nonreacting reactor candidate. Existing
compiled observers can be reused if their capabilities remain sufficient.

## Source-pair launch reference (completed)

The command below has completed successfully. Do not repeat it: both source
cases are now attempted and the runner will refuse them.

The source cases were prepared with the corrected dictionaries at `d7b1b2e`
and retain the same 256-cell mesh, 100 steps to
0.1 s and numerical controls. The prescribed liquid feed is 0.001 kg/s, giving
0.0001 kg (0.1 g) total supply. The pressure outlet can export both phases;
acceptance compares inventory changes with measured boundary mass and the
assembled source, rather than requiring all supplied liquid to remain inside.
The source dictionary was checked against Foundation OF9's
[massSource implementation](https://cpp.openfoam.org/v9/massSource_8C_source.html):
it selects `phase liquid`, `rho thermo:rho.liquid`, all cells, and the existing
velocity/temperature field values. No case or compiled source changes are
were needed for this completed test. Interphase transfer and biology are absent.

The runner rechecks the private build, both closed log hashes and passing
screens, and both source input manifests before launching either source case.
It runs `source_no_reference` then `source_gas_reference`, serially. A complete
numerical failure is still reported as `screen_pass: false`; preserve that
result as a finding. If a command fails, preserve the partial root for review.
Review the two budgets before choosing the next model or numerical change.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
TWOFIELD_CASES=/scratch/nwimer/sc1_two_field_cases_20261110T233431Z

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field.py" \
        run --root "$TWOFIELD_CASES" --pair source --run
else
    echo "Source-pair setup stopped. Preserve the existing cases and paste the output."
fi
```

Both source reports and `PAIR=source complete` are recorded above. This
historical command launched no closed-case reruns and no reactor simulation.

## Closed-partner continuation reference (completed)

The following continuation has completed successfully. Do not repeat it;
the partner is now attempted, and the runner will refuse it.

`continue-closed --root CASE_ROOT --run` verifies the private build and case
manifest, the first case's attempt/build record, initialized input hashes,
prepared non-alpha inputs, mesh-check evidence, complete solver audit and
passing numerical screen. It also requires the gas-reference case to be fresh
and match its prepared input hashes. A missing first `result.json` is then
created with the saved log's SHA-256; an existing identical result is reused
without rewriting it. A conflicting result or attempted partner stops the
command. It runs only the partner's mesh/initialization/solver sequence and
saves that result normally. The existing pair runner still refuses attempted
cases. Neither completed CFD nor the source pair is launched by this command.

All 39 harness tests pass. Synthetic recovery tests verify that only the
partner's four commands run, the completed case changes only by addition of
its missing report, source cases remain unchanged, and inconsistent/partial
evidence or an attempted/modified partner blocks before mutation. These tests
mock CFD; the later Kestrel transcript above confirms the continuation worked.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
TWOFIELD_CASES=/scratch/nwimer/sc1_two_field_cases_20261110T233431Z

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field.py" \
        continue-closed --root "$TWOFIELD_CASES" --run
else
    echo "Continuation setup stopped. Preserve the existing cases and paste the output."
fi
```

The successful output is recorded above. The source pair was not automatically
launched by this continuation.

## Saved-log analysis reference (completed)

No allocation, OpenFOAM environment setup, rebuild, preparation or solver
launch is required by `analyze`. Pull the reporting change over SSH and run:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_LOG=/scratch/nwimer/sc1_two_field_cases_20261110T233431Z/closed_no_reference/log.solver

if [[ -x "$PY" && -r "$SC1_LOG" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field.py" \
        analyze --log "$SC1_LOG" --case closed_no_reference
else
    echo "Analysis setup stopped. Preserve the existing case and paste the output."
fi
```

This completed analysis command prints to the terminal without changing the
case or creating `result.json`. All four cases are now complete.

## Fresh closed-pair launch reference (not the current next action)

The command below is retained as the original launch recipe. **Do not execute
it for the completed baseline suite.**

Use the existing allocation and pull the corrected generator over SSH. Reuse
the already verified build. The command creates four fresh case dictionaries but launches only `closed_no_reference` and
`closed_gas_reference`, serially, with 256 cells and 100 timesteps each.
The two source cases remain unattempted.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
TWOFIELD_BUILD=/scratch/nwimer/sc1_two_field_build_20261110T224950Z
TWOFIELD_CASES="/scratch/nwimer/sc1_two_field_cases_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field.py" \
        prepare-cases --build-root "$TWOFIELD_BUILD" --out "$TWOFIELD_CASES"; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field.py" \
        run --root "$TWOFIELD_CASES" --pair closed --run
else
    echo "Preparation stopped. Preserve any created directory and paste the output."
fi
```

Paste the complete output, including `TWOFIELD_CASES`, both case summaries and
any error text. The command first rechecks build hashes/linkage, then generates
and checks the small meshes, initializes the fields and runs each fresh case.
It saves `log.solver` and `result.json` in each completed case. If any stage
fails, preserve the partial case root and do not repeat the command to launch
another pair before reviewing that failure.

A completed pair may still report `screen_pass: false`; process completion is
not numerical acceptance. Review both phase inventories, rescaling,
source/flux accounting, phase bounds and the recorded residuals. The later
source-pair command is `sc1_two_field.py run --root CASE_ROOT --pair source
--run`, and remains gated on successful closed-pair numerical screening.
No source-pair launch is included above.

## Validation and remaining work

All 39 local harness tests pass after the guarded continuation addition.
The preceding 208-test SC1 suite passed before these fixture/reporting/runner changes.
The harness tests cover paired fixtures,
source/hook placement, retained phase equations, private preparation, runtime
linkage rejection, attempted/modified case rejection, source-run gating and
malformed/nonconservative audit streams, void-assignment compatibility and
useful compiler-error output, the OF9 pressure-dictionary contract and
system/phase alpha-control lookup coverage, timestamp matching and strict
coverage diagnostics and evidence-preserving closed-pair continuation.
The shared standalone C++ MULES
arithmetic tests run as part of that suite. Python syntax and shell syntax
checks pass. Local tests alone did not prove OpenFOAM compatibility; the
subsequent successful Kestrel build above now supplies that compilation
evidence. The first two CFD attempts exposed corrected dictionary errors.
All four cases now pass complete audit coverage and the numerical screen.
This result-recording update changes only
documentation; the 39-test harness result above is from the preceding runner
change, not a new CFD or test-suite execution.

The passes in these small cases justify further controlled tests, not a
full reactor release. Follow with source/transfer tests near phase bounds,
timestep/outer-iteration sensitivity and explicit-interface benchmarks, then
a matched SC1 nonreacting case on a common mesh. Gas mass, liquid mass,
credible wall/free-surface behavior, and eventual oxygen/species conservation
remain separate requirements. Retain a third field only for a demonstrated
physical limitation of a simpler representation.

Primary implementation references checked alongside repository source:
[Foundation OF9 phaseSystemSolve](https://cpp.openfoam.org/v9/phaseSystemSolve_8C_source.html),
[OF9 massSource](https://cpp.openfoam.org/v9/massSource_8C_source.html), and
[OF9 fvMatrix operator interface](https://cpp.openfoam.org/v9/fvMatrix_8H_source.html).
