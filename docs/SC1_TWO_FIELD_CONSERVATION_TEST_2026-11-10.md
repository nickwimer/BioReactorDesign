# SC1 isolated two-field conservation test

## Status and purpose

Implemented locally. The first Kestrel build stopped at the observer's field
initialization; that defect is corrected, and a **fresh build and CFD execution
remain pending**. This is a numerical test harness, not a corrected SC1 reactor
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
The next action is the corrected build in a **new timestamped directory**;
do not patch the failed copy or rerun into it.

## Next Kestrel step: fresh build only

Use the existing allocation. SSH is explicit in the pull command; no change
to `origin` and no move to another node is required by this workflow.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
TWOFIELD_BUILD="/scratch/nwimer/sc1_two_field_build_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    bash "$REPO/scripts/diagnostics/build_sc1_two_field.sh" "$TWOFIELD_BUILD"
else
    echo "Setup stopped. Preserve existing directories and paste the output."
fi
```

The command builds six libraries and one solver; it launches no CFD. It prints
`TWOFIELD_BUILD` before compilation, each build command, failure-log tails if
needed, and binary hashes after successful runtime verification. Preserve a
failed build rather than rerunning into the same destination.

After that build has been reviewed, prepare a **fresh** case root with
`sc1_two_field.py prepare-cases --build-root BUILD --out NEW_CASE_ROOT`, then
run `sc1_two_field.py run --root NEW_CASE_ROOT --pair closed --run`. The source
pair uses `--pair source` only after the closed results pass. These are the
supported next stages, not claims that they have already been run.

## Validation and remaining work

208 local SC1 tests pass, including 26 harness tests covering paired fixtures,
source/hook placement, retained phase equations, private preparation, runtime
linkage rejection, attempted/modified case rejection, source-run gating and
malformed/nonconservative audit streams, void-assignment compatibility and
useful compiler-error output. The shared standalone C++ MULES
arithmetic tests run as part of that suite. Python syntax and shell syntax
checks pass. **No OpenFOAM installation is available in this workspace, so
these results do not claim that the new OpenFOAM headers compile.**

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
