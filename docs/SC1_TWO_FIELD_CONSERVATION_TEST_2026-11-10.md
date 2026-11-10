# SC1 isolated two-field conservation test

## Status and purpose

**The corrected private build and runtime verification passed on Kestrel.**
Use `/kfs3/scratch/nwimer/sc1_two_field_build_20261110T224950Z`; no rebuild
is needed. The first attempt stopped on a legacy pressure-limit entry. After
that correction, the second attempt initialized both phase inventories and
entered the first timestep, but stopped on missing `alpha.gas` solver controls
before completing a step. The generator now covers both the system-level and
phase-specific alpha lookups. The next action is the closed pair in fresh cases.
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

All 29 harness tests pass. New fixture-contract tests resolve the system and
phase field names for all four cases and reproduce the failed phase lookups
with a literal `alpha` key. These restricted Python checks are not native
OpenFOAM dictionary or CFD tests. Source review also checked the selected
pressure/momentum controls and schemes: the pure isothermal laminar fixtures
do not solve energy, species or turbulence transport equations; implicit
phase pressure defaults to false. Successful timestepping still requires the
fresh Kestrel run below.

## Next Kestrel step: prepare fresh cases and run only the closed pair

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

All 29 local harness tests pass after the phase-specific alpha correction.
The preceding 208-test SC1 suite passed before these fixture-only corrections.
The harness tests cover paired fixtures,
source/hook placement, retained phase equations, private preparation, runtime
linkage rejection, attempted/modified case rejection, source-run gating and
malformed/nonconservative audit streams, void-assignment compatibility and
useful compiler-error output, the OF9 pressure-dictionary contract and
system/phase alpha-control lookup coverage.
The shared standalone C++ MULES
arithmetic tests run as part of that suite. Python syntax and shell syntax
checks pass. Local tests alone did not prove OpenFOAM compatibility; the
subsequent successful Kestrel build above now supplies that compilation
evidence. The two CFD attempts exposed the corrected dictionary errors;
successful timestepping and numerical conservation remain unproved.

A pass in these small cases would justify further controlled tests, not a
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
