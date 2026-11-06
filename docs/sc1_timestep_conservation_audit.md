# SC1 timestep phase and O2 conservation audit

**2026-11-06: implemented and locally tested; OpenFOAM compilation and MPI
validation on Kestrel are pending. No new CFD results are claimed.**

This continues the [session handoff](SC1_MEML_SESSION_HANDOFF_2026-11-06.md),
[entrainment results](sc1_headspace_reentrainment_trial.md),
[film trials](sc1_headspace_oxygen_transfer_design.md), and
[historical audit](sc1_three_phase_headspace_audit.md).
The existing zero-control and normal-active 0.5 s runs are complete and remain
the experimental record. Their sparse outputs cannot recover missing timestep
fluxes. A fresh diagnostic repetition is justified specifically to collect
those fluxes, with unchanged physics and an identical initial state.

## Source review and what the audit can establish

The reviewed branch base is `8e10fc70f0048b1e414d7d92701781c5e6122ef0`.
Relevant code, under `OFsolvers/birdmultiphaseEulerFoam/`:

| Source | Consequence for this audit |
|---|---|
| `putidaFoam/putidaFoam.C` | `fluid.solve()` (including PB) precedes `fluid.correct()`; species solve precedes pressure/thermo updates. Capture final-iteration equation data separately from end-of-step inventories. |
| `phaseSystems/PhaseSystems/PhaseTransferPhaseSystem/PhaseTransferPhaseSystem.C` | `dmdtf(bubbleGas,headspaceGas)` is signed gain to bubbleGas. `addDmdtYfs` uses explicit donor and implicit sink terms; evaluate the stored species matrices, not simply 0.23 times net gas transfer. |
| `phaseSystems/PhaseSystems/InterfaceCompositionPhaseChangePhaseSystem/InterfaceCompositionPhaseChangePhaseSystem.C` | Liquid uptake is semi-implicit; the other phase's source contains an explicit liquid composition. At finite coupling iterations, equal-and-opposite final evaluated source actions are not guaranteed. Measure their sum before diagnosing a leak. |
| `phaseSystems/phaseModel/MultiComponentPhaseModel/MultiComponentPhaseModel.C` | The species equation includes advection, `divj`, `alpha*R`, and a residual-alpha temporal stabilizer. `correctSpecies()` normalizes species. |
| `phaseSystems/PhaseSystems/PopulationBalancePhaseSystem/PopulationBalancePhaseSystem.C` | PB solve is inside `fluid.solve()`. PB pair transfer is added to the phase-system sources. |
| `phaseSystems/populationBalanceModel/populationBalanceModel/populationBalanceModel.C` | PB source refresh has its own scheduling; it need not consume the post-PIMPLE source. |
| `phaseSystems/diameterModels/velocityGroup/velocityGroup.C` | `scale()` clamps negative groups then divides by their sum. Final sum(fi)=1 cannot independently certify conservation. |

OpenFOAM 9's [fvMatrix operator& implementation](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/finiteVolume/fvMatrices/fvMatrix/fvMatrix.C)
returns `(A*Y - source)/V`, including boundary contributions. This is the
action of a matrix appearing on either side of an equation, not necessarily a
positive source. The new ledger uses positive RHS action as gain and positive
LHS advection/diffusion as loss. A complete equation's action is its signed
unrelaxed algebraic residual.

**A potential discrepancy mechanism, not a finding of an actual leak:** the
semi-implicit receiver can consume a different O2 amount from the explicit
donor loss at finite coupling convergence. The logged sum over all three
phase species-transfer matrices tests this directly. No coefficient, transfer
equation, biological update interval, or normalization method is changed.

## Measurements and signs

Enable `PIMPLE/sc1ConservationAudit true` only in the new scratch cases. The
switch defaults to false and is read at startup. With it absent, no audit
object, field reduction, source reevaluation, or audit output is created.
The species source expression is retained once as a matrix temporary so the
enabled audit can copy the actual RHS without calling `fvModels.source` twice.

The audit starts from the loaded checkpoint **before the first advance**, so
the first interval is included. It accumulates once per global timestep after
PIMPLE, using actual `deltaT`, and prints at writes and at the configured end
time (including a final partial write interval). It does not maintain a ledger
across process restarts; use one uninterrupted process/log per trial.
The final-report predicate uses `!runTime.running()`: OpenFOAM 9 normally stops
when `t >= endTime - 0.5*deltaT`, which can precede the nominal `endTime`.
It deliberately does not call `runTime.run()`, which executes function objects.

| Record | Meaning |
|---|---|
| `SC1_BALANCE_START` | Initial time, mass, O2 mass and volume for liquid, bubbles, headspace, combined gas and total. |
| `SC1_BALANCE` | Current inventories, cumulative and write-window residuals, timestep counts, endpoint/stage rates and cumulative integrals. Liquid DO is O2 mass / liquid volume in mg/L. |
| `SC1_BOUNDARY` | Each non-coupled physical patch and phase: positive gross inward/outward mass rates and integrals, separately; face-composition O2 proxies. Processor/cyclic faces excluded. |
| `SC1_PAIR` | Signed blended bubble gain from headspace; bubble-to-liquid and headspace-to-liquid mass losses; positive-net pre-PB input proxy; donor-composition gas-gas O2 proxy. These are decompositions, **not additional ledger sources**. |
| `SC1_PB` | Volume, density-weighted mass and number moments after the final PB solve, before `fluid.correct()`; their immediate changes across `fluid.correct()` and cumulative correction impulses. |

For each phase, using kg for integrated terms:

```text
mass residual = M(t)-M(start)
              - integral(external mass + assembled dmdt + gross in - gross out)

oxygen residual = O(t)-O(start)
                - integral(external O2 + interphase O2 + reaction O2
                           - advection action - diffusion action)
                - measured pre-advance O2 jumps
                - measured final species-normalization changes
```

`oxygenExternal` is the captured `fvModels` species matrix action;
`oxygenReaction` is a same-stage reconstruction of `alpha*phase.R(O2)`.
The actual cluster implementation must establish which contains frozen OUR.
**Do not label the aggregate external source “biological consumption” without
reviewing the local source definition.** Sparger species injection can also
appear there. The liquid's two rates expose the sink without changing it.
Any direct change from scheduled `bioReact.H`/`phReact.H` between steps is
recorded as `oxygenJump`, not erroneously multiplied by CFD `deltaT`.
The intentional four-hour biology schedule remains unchanged.

For the reviewed OpenFOAM 9 `basicSpecieMixture::normalise()`, solved active
species are clipped to nonnegative values and the default species is assigned
`max(1-sum(active Yi),0)`. It does **not** rescale all active species to unit
sum on each correction. If the local frozen-OUR source requests more O2 than
is available, an explicit solve can become negative and this clipping restores
it to zero. A positive O2 correction impulse must therefore be reported as
clipping, not physical oxygen supply. Verify the actual Kestrel sink before
attributing an observed impulse to this mechanism; keep scheduled biological
updates unchanged and distinguish requested OUR from realized uptake.

The actual unrelaxed species equation and actual transfer/external matrices
are copied only for O2 in the final PIMPLE iteration. Their actions are
evaluated after all species solves, before `fluid.correctSpecies()`; this
retains explicit donor lag. Advection, diffusion and reaction matrices are
reconstructed through the same APIs used by the reviewed `YiEqn` at assembly
time, then evaluated at the solved O2. The integrated advection/diffusion
actions, **not** the patch face-Y proxies, enter the oxygen ledger. Inspect a
local `YiEqn` divergence before interpreting these reconstructed components.

The algebraic residual includes stabilization and finite solve/relaxation
error. Both signed and volume-integrated absolute (L1) residuals are reported,
so spatial cancellation cannot conceal a large equation error. These are
diagnostics only and are **not** added as another physical source.
Internal species exchange should cancel over all phases at coupling
convergence; the combined-gas source additionally includes dissolution to
liquid, so it need not be zero.

## Explicit limitations / remaining work

1. The timestep integrals are **right-endpoint quadrature**, not an exact
   solver-discrete mass balance. End-of-step pressure/density updates, alpha
   subcycling, nonlinear iteration, temporal scheme, residual-alpha
   stabilization and earlier normalization/correction stages can explain
   residuals. Capturing one complete unrelaxed species equation does not make
   the whole inventory ledger exact. `exactDiscreteBudget=0` is logged.
2. Phase `dmdts()` and density `fvModels` sources are evaluated once after
   PIMPLE; they are not intercepted at every alpha/pressure solve. Mass-source
   reevaluation should be observational for the SC1 sources, but a stateful
   user `coded` source must be checked. Compare the audited control to the
   preserved zero-control trajectory before trusting the diagnostic.
3. Gross O2 boundary values use actual patch Y (including inletOutlet inflow
   composition), weighted by phase mass flux. These are explicitly proxies
   for advective face flux. They exclude diffusion and need not equal the
   species discretization. The separately reconstructed matrix action includes
   the selected advection scheme and `divj` boundary effects in the ledger.
4. PB before/after correction changes include all work in `fluid.correct()`;
   they do not isolate `scale()` alone. A pre-solve positive-net source is a
   **proxy**, not an interception of each PB source refresh. Number changes
   are expected under breakup/coalescence; number is not generally conserved.
   Full pre-renormalization PB equation/moment closure remains outstanding.
5. Gas-liquid pair rates are total transferred mass. They are O2 rates only
   under the verified O2-only SC1 interface-composition configuration. The
   three-phase species matrices remain the authoritative O2-source diagnostic.
6. Only fixed-mesh, global-time, three-phase SC1 with active O2 in all phases
   is supported. Missing O2 hooks, moving/changing mesh or local timestepping
   stop this diagnostic; no zero/missing terms are silently substituted.
7. The first Kestrel build transcript is now available (see below), but there
   is no successful audit build or CFD runtime result yet. The existing
   Kestrel case files, raw simulation logs and experimental video are not
   mounted here. The latest simulation conclusions remain those recorded in
   the handoff, not a new independent analysis of raw cluster output.

The follow-up source review corrected an OpenFOAM 9 API compatibility error:
`Info` is a `messageStream`, so precision is set through `Info().precision()`
on its underlying stream. This source correction is not a substitute for the
pending Kestrel OpenFOAM build. A compiled regression now exercises the actual
report predicate at times below, at and beyond OpenFOAM's stopping threshold.

### First Kestrel build and diffusion-access correction (2026-11-06)

The supplied build transcript at commit `19410f4` passed all seven diagnostic
and ten analytical tests, then failed compiling `putidaFoam.C`:
`phaseModel has no member named divj`. No CFD launched. Preserve its scratch
root `/scratch/nwimer/sc1_conservation_audit_20261106T182403Z` and `build.log`.
Torch old-style-cast and biological unused-variable messages shown alongside
the error were warnings, not the reported build failure.

`divj` exists on `MovingPhaseModel`, but not the abstract `phaseModel` reference
held by the audit. The corrected audit looks up the existing phase-grouped
`thermophysicalTransportModel` using the mass-flux group and calls its `divj`.
This is the same delegate called by `MovingPhaseModel::divj`; no replacement
diffusivity, correction call, source term or library rebuild is introduced.
Missing or incompatible registered transport models fail rather than produce
a zero diffusion term. The scratch preparer appends the required header path
and direct `thermophysicalTransportModels` link flag to the copied local
`Make/options`; it preserves all original flags and libraries.

Nine local diagnostic tests now pass. New checks compile the actual lookup
block against a strict mock interface (three phase groups, missing/wrong model
rejection, no field correction) and evaluate appended build flags with GNU
Make. These tests do not compile against OpenFOAM headers or validate its ABI.
The corrected OpenFOAM build remains pending; prepare a new scratch root for
the rebuild, then review linked libraries before any CFD launch.

API evidence: [OF9 transport interface](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/thermophysicalTransportModel/thermophysicalTransportModel.H),
[grouped object registration](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/thermophysicalTransportModel/thermophysicalTransportModel.C),
and [library name](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/Make/files).

## Kestrel build gate — no CFD launch

Run on a current allocation/environment; do not assume job 19721806 still
exists. The script clones the **current scratch solver**, preserving local
biology and RHEL9 changes. It refuses changed/unrecognized hook locations,
an existing destination and overlapping paths. `Make/files` points at a
**new executable under the new scratch root/bin**, never the old executable.
Existing source dependency trees are linked for headers only; no BiRD shared
library is rebuilt. A source-hash manifest accompanies the clone.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
AUDITROOT=/scratch/nwimer/sc1_conservation_audit_20261106_v1

source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc"
module load libfabric

if git -C "$REPO" pull --ff-only; then
    if "$PY" "$REPO/scripts/diagnostics/test_sc1_conservation_audit.py"; then
        if "$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_solver.py" \
            --out "$AUDITROOT"; then
            if (cd "$AUDITROOT/solver" && wmake > "$AUDITROOT/build.log" 2>&1); then
                echo "SC1 audit build succeeded; no CFD launched."
                ldd "$AUDITROOT/bin/putidaFoamSC1ConservationAudit"
            else
                echo "SC1 audit build failed; no CFD launched."
                tail -80 "$AUDITROOT/build.log"
            fi
        fi
    fi
fi
```

If the destination exists, preserve it and choose a new suffix; do not delete
or overwrite it. Return build output before proceeding. The pure C++ ledger
test uses `g++`; if absent it explicitly reports a skip, not compilation success.

Also inspect local source before interpreting the audit. These are read-only:

```bash
BIRD="$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam"
for REL in \
  phaseSystems/PhaseSystems/PhaseTransferPhaseSystem/PhaseTransferPhaseSystem.C \
  phaseSystems/PhaseSystems/PopulationBalancePhaseSystem/PopulationBalancePhaseSystem.C \
  phaseSystems/phaseModel/MultiComponentPhaseModel/MultiComponentPhaseModel.C \
  phaseSystems/diameterModels/velocityGroup/velocityGroup.C
do
    diff -u "$REPO/OFsolvers/birdmultiphaseEulerFoam/$REL" "$BIRD/$REL"
done
sed -n '1,220p' /scratch/nwimer/sc1_bubble_audit_20261029/solver/YEqns.H
foamDictionary /scratch/nwimer/sc1_headspace_entrainment_20261104/run_control_05/system/fvSchemes \
    -entry ddtSchemes -value
```

Differences are evidence to review, not permission to replace local model
sources with the repository versions. Inspect local `fvModels` O2 source and
reaction definitions to identify the frozen biological sink explicitly.

## Matched cases — after successful build/source review

```bash
PAIRROOT=/scratch/nwimer/sc1_conservation_pair_20261106_v1
"$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_pair.py" --out "$PAIRROOT"
```

The preparer reads the explicit original checkpoint from the **completed**
`run_control_05` and `run_normalactive_05`. It hashes every initial collated
field and mesh file and requires identical inputs. It retains each existing
phaseProperties and checks that only the established efficiency/velocity-gate
switches differ; it does not regenerate or tune the physics. Constant files,
including biology, and numerical settings must match. The original zero
control's absent velocity criterion is accepted as its legacy default.
Both copies use the same 0.5 s duration, write interval and audit switch.
All data are copied, not hardlinked. Model asset files are copied if present
at case root. Review any absolute external runtime-asset paths as before.

Run **control only first**, never append to or truncate an existing log:

```bash
AUDITEXE="$AUDITROOT/bin/putidaFoamSC1ConservationAudit"
if [[ -x "$AUDITEXE" && -d "$PAIRROOT/control/processors52" && ! -e "$PAIRROOT/control/log.conservation" ]]; then
    (cd "$PAIRROOT/control" && srun -n 52 "$AUDITEXE" -parallel -fileHandler collated > log.conservation 2>&1)
    SC1_STATUS=$?
    echo "SC1 control solver status: $SC1_STATUS"
    "$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" "$PAIRROOT/control/log.conservation"
fi
```

Return the `SC1_BALANCE*`, `SC1_PAIR`, `SC1_PB`, outlet `SC1_BOUNDARY`, legacy
`BUBBLE_MASS_AUDIT` lines and solver status. Compare the control's inventory,
OTR, diameters and transfer rates to the preserved `run_control_05`, allowing
only justified numerical tolerances. Check the initial inventory against an
independent domain integration; verify the assembled interphase mass sum is
zero over all phases and examine the O2 source sum, normalization and equation
residual. No automated scientific PASS threshold is imposed.

Only after this control review, run the already prepared `normal` case with
the same executable/ranks and a new `log.conservation`. Then:

```bash
"$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" \
    "$PAIRROOT/control/log.conservation" --normal "$PAIRROOT/normal/log.conservation"
```

The reader rejects missing initial states, incomplete logs, duplicate sessions,
nonfinite output, incompatible initial inventories and unmatched report times.
It independently checks reported mass-residual arithmetic. It never integrates
old 0.1 s endpoint samples as if they were timestep data.

## Decision after the paired audit

First classify discrepancies: phase-source sum, explicit/implicit O2 imbalance,
normalization impulses, unrelaxed equation error, boundary transport or
remaining source/flux evaluation timing. If unresolved, add targeted
alpha/pressure/PB equation-stage instrumentation in a new scratch build before
changing entrainment physics. Record any timestep/convergence sensitivity as a
separate diagnostic experiment, not a modification of these matched controls.

At common saved times, combine the new balance records with the existing
surface/upper-impeller postprocessing, using new scratch VTK destinations.
The audit itself does not calculate a free-surface shape. Request the actual
experimental video and operating metadata, plus available gas holdup, bubble
size, blade wetting, free-surface and kLa/OTR measurements. Test connected cavity
versus dispersed, liquid-wetted froth explicitly; high headspace alpha alone
does not mean a dry impeller. There is no evidence yet selecting a new closure
or calibrating foam area, nucleus size or entrainment efficiency.

The nonreacting 9.9 s checkpoint and subsequent MM/EMM/MEML comparisons remain
blocked on physical validation and conservation, not on the intentional
four-hour biological updates. Neither experimental closure is production-ready.
