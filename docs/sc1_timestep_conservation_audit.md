# SC1 timestep phase and O2 conservation audit

**2026-11-06: OpenFOAM compilation and 19 checks passed on Kestrel at
`6ac7e36`. The checkpoint mismatch was proved to be four Co location values;
16 metadata checks and matched-case preparation passed at `289a606`.
Both 52-rank diagnostics completed with status 0; all five new control records
reproduce the old control at printed precision. The budgets expose large
correction terms and unresolved residuals; validation has not passed.
The next isolated stage audit compiled successfully on Kestrel at `4a1d7de`,
with all 50 checks passing and dependencies resolved. Its MPI runtime
validation is pending; reuse the `sc1_stage_audit_20261106T231907Z` build.** Read the
[paired results and current build command](sc1_conservation_pair_results_2026-11-06.md)
before using the historical commands below.

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
7. Kestrel compilation now succeeds at `6ac7e36` (see build record below),
   but there is no audit CFD runtime result yet. The existing Kestrel case
   files, raw simulation logs and experimental video are not mounted here.
   The latest simulation conclusions remain those recorded in the handoff,
   not a new independent analysis of raw cluster output.

The follow-up source review corrected an OpenFOAM 9 API compatibility error:
`Info` is a `messageStream`, so precision is set through `Info().precision()`
on its underlying stream. That API correction alone was not proof of a
successful OpenFOAM build. A compiled regression now exercises the actual
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
The corrected OpenFOAM build subsequently succeeded in a new scratch root,
as recorded below; the failed tree remains preserved.

API evidence: [OF9 transport interface](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/thermophysicalTransportModel/thermophysicalTransportModel.H),
[grouped object registration](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/thermophysicalTransportModel/thermophysicalTransportModel.C),
and [library name](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/ThermophysicalTransportModels/Make/files).

### Successful Kestrel build (2026-11-06)

The supplied follow-up transcript confirms commit
`6ac7e36b4efdd5833282d67853ed910b691a2bc0`, nine diagnostic tests and ten
analytical tests passing on Kestrel, and successful `wmake` output at:

```text
/scratch/nwimer/sc1_conservation_audit_20261106T183120Z/bin/putidaFoamSC1ConservationAudit
```

The preparer resolved this root through `/kfs3/scratch/nwimer`. The executable's
`ldd` listing has no `not found` entries: BiRD libraries resolve under
`OpenFOAM9_RHEL9/user-9`, OpenFOAM libraries under `OpenFOAM9_RHEL9/OpenFOAM-9`,
Torch under `/home/nwimer/libs/libtorch`, and MPI/libfabric resolve. The new
`libthermophysicalTransportModels.so` direct dependency is present.

No CFD launched. This validates compilation and dependency resolution only;
the existing successful build can now be used without recompilation. Still
inspect libraries loaded dynamically by `controlDict`, local source differences,
actual O2 sources, biological initialization/schedule, numerical controls,
function objects and any external asset/write paths before the control run.
Hashing binaries records identity but does not prove which source built them.

The pair preparer hashes/checks collated initial fields, mesh, constants and
numerical settings before copying. Root model/scaler assets are copied but
are not included in its matched-input validation: compare those assets across
the new cases separately. Verify each new `processors52` contains only
`constant` and `31.199982372783058`. This audit pair preserves the old reacting/
frozen-OUR configuration deliberately; do not relabel it as the nonreacting
9.9 s baseline. Any source/config mismatch is evidence to review, not an
instruction to overwrite the local source or change biology settings.

### Checkpoint mismatch before case preparation (2026-11-06)

The next supplied terminal output stopped at
`Inputs are not identical: processors52/31.199982372783058`. The collated
mesh/constant comparison preceding that check passed. This is a raw byte/file
set mismatch, not yet evidence of different physical field values. Validation
occurs before the destination is created, so this attempt created neither
matched case and launched no CFD. The successful executable is unaffected.

Use the new read-only inspection before attempting preparation again:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_BRANCH=codex/sc1-headspace-bubble-audit-20261029
if [[ "$(git -C "$REPO" branch --show-current)" == "$SC1_BRANCH" ]] \
    && git -C "$REPO" pull --ff-only origin "$SC1_BRANCH"; then
    "$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_pair.py" --inspect-inputs
else
    echo "Repository update stopped; paste the output and git status."
fi
```

The report covers both collated mesh/checkpoint trees and the constant/system
directories, continuing past differences to expose other discrepancies in
one inspection. It lists changed and one-sided files, sizes, SHA-256 hashes,
resolved paths and bounded hex/escaped-ASCII excerpts at the first byte
difference. Hashing detects ordinary size/mtime/inode changes during a read;
the report is not an atomic filesystem snapshot. Expected raw differences in
phaseProperties/controlDict are reported too; they do not bypass the preparer's
separate configuration checks.

Differences remain `classification: undetermined`: the first mismatch cannot
prove that later numerical data agree. Collated binary fields are not decoded,
normalized, stripped of headers or rewritten. Inspection creates no files;
`--out` is rejected in inspection mode. Exit status zero means inspection
completed without read errors, not that inputs are equal or scientifically
accepted. Preserve all original inputs and review the report before selecting
any common seed or changing the comparison. No solver rebuild is needed.

Eight input-inspection regression tests pass locally, covering nested missing
files, hash/size reports, differences beyond 4 MiB, an early difference with
later changed data, unequal lengths, read errors, escaped bounded output and
CLI refusal of copy/write/preparation operations. These tests do not decode
OpenFOAM fields or establish equality of the actual Kestrel checkpoints.

### Courant-only mismatch and guarded common checkpoint (2026-11-06)

The next inspection report shows:

| Input section | Raw agreement |
|---|---|
| Collated mesh | 36/36 files identical |
| Checkpoint | 55/59 files identical; only four Co files differ |
| Case constants | 50/51 identical; phaseProperties differs |
| Numerical system | 12/12 files identical |

The four files are `Co.liquid`, `Co.bubbleGas`, `Co.headspaceGas` and
`Co.mixture`. Each normal file is five bytes longer, matching the length
difference between `run_control_05` and `run_normalactive_05`. The first
difference lies in an outer-header absolute location path. This is evidence
for a metadata explanation, not yet a proof of the remaining bytes.

The opt-in `--allow-courant-location-metadata` now permits preparation only if:

1. File sets and all mesh bytes agree; every other checkpoint file is raw-identical.
2. Each changed file is exactly one of those four names and has a flat, bounded
   outer `FoamFile` header with binary format, `decomposedBlockData` class,
   matching object name and one quoted location resolving to its actual file.
3. The complete outer headers agree after replacing ONLY that location value
   in memory, and every remaining byte through EOF agrees. Embedded per-rank
   headers, binary payload, whitespace and all other metadata remain compared.
4. The proof hashes agree with the preceding raw input hashes.

There is no blanket exclusion of Co fields, broad header stripping, floating
point decoding, physical-field tolerance or fallback. The ordinary mode
retains strict byte equality. The new option does not modify the originals.
After proof, the entire control mesh/checkpoint is copied unchanged to BOTH
new cases, retaining all 59 checkpoint files. Every destination hash must
match the control's raw map. Each case retains its own already validated
physics configuration. The manifest records the common source, original Co
hashes, validated locations and full header/tail proof hashes.

Sixteen synthetic-container/regression tests pass locally, including changes
after 4 MiB, malformed headers, additional metadata or payload changes,
forbidden fields, missing files, strict default behavior and an end-to-end
59-file copy with unchanged originals. These tests do not establish that the
actual Kestrel Co tails match; the next preparation performs that proof.
No solver recompilation or CFD launch is part of this step.

### Next command: prove, prepare and collect evidence

The preflight collector below replaces the earlier long pasted block. It
writes a new `preflight.txt` only, refuses an existing report, and reads the
prepared cases, source dependency links and known experimental library paths.
It does not launch CFD. Compare those library paths against the actual
`controlDict` in its output, and review all errors, included dictionaries,
runtime model assets and external write/dependency paths. Empty source diffs
mean identical files, not proof that a shared binary was built from that source.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_BRANCH=codex/sc1-headspace-bubble-audit-20261029
AUDITROOT=/scratch/nwimer/sc1_conservation_audit_20261106T183120Z
PAIRROOT="/scratch/nwimer/sc1_conservation_pair_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ "$(git -C "$REPO" branch --show-current)" == "$SC1_BRANCH" ]] \
    && git -C "$REPO" pull --ff-only origin "$SC1_BRANCH" \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_checkpoint_metadata.py" \
    && "$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_pair.py" \
        --allow-courant-location-metadata --out "$PAIRROOT"; then
    bash "$REPO/scripts/diagnostics/collect_sc1_conservation_preflight.sh" \
        "$AUDITROOT" "$PAIRROOT"
else
    echo "Preparation stopped. Paste the output; preserve any created directory."
fi
```

Return preparation output and `preflight.txt`. If a proof/configuration check
fails, preserve its output and do not bypass the failure. Successful copying
and report collection do not establish scientific validity. The zero-control
launch follows only after input/source review.

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
"$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_pair.py" \
    --allow-courant-location-metadata --out "$PAIRROOT"
```

The preparer reads the explicit original checkpoint from the **completed**
`run_control_05` and `run_normalactive_05`. It hashes every initial collated
field and mesh file and requires identical inputs except the opt-in, fully
proved Co outer-location exception described above. It copies the complete
control state to both new cases and verifies raw target equality. It retains each existing
phaseProperties and checks that only the established efficiency/velocity-gate
switches differ; it does not regenerate or tune the physics. Constant files,
including biology, and numerical settings must match. The original zero
control's absent velocity criterion is accepted as its legacy default.
Both copies use the same 0.5 s duration, write interval and audit switch.
All data are copied, not hardlinked. Model asset files are copied if present
at case root. Review any absolute external runtime-asset paths as before.

### Matched-pair preflight reviewed (2026-11-06)

The uploaded Kestrel transcript confirms the successful preparation at
`289a606`, with all 16 metadata tests passing. The four complete Co payloads
and all other checkpoint files passed equality checks. Both cases contain
all 59 checkpoint files copied from the control, plus the identical collated
mesh. The initial `processors52` directories contain only `constant` and
`31.199982372783058`. No original run was changed and no CFD was launched.

Reuse these existing artifacts:

| Artifact | Location / SHA256 |
|---|---|
| Pair and preflight report | `/scratch/nwimer/sc1_conservation_pair_20261106T201549Z`, `preflight.txt` |
| Audit build | `/scratch/nwimer/sc1_conservation_audit_20261106T183120Z` |
| Executable | `e8b71ba803b23661fd137b82f3a42fe105e926eea4da170a74c2290edbc77dde` |
| Film library | `f823c8a3a3c83067171882cae53dfce4f9988dbc34c07bae6f6d58668c77ab4d` |
| Entrainment library | `27c8ce353108690a041454d2a38957b9c88862e2215f0f75e56fb6a3b3f847fb` |

All eight compared dependency sources match the repository: phase transfer,
interface-composition transfer, population-balance phase system,
MultiComponentPhaseModel, MovingPhaseModel, phaseSystemSolve, velocityGroup
and populationBalanceModel. `YEqns.H` also matches. Preserved local solver
differences include restart clock initialization, skipping `wellMixed.H`,
NH4OH configuration and biological/pH handling. Do not replace those sources.
Executable and experimental-library `ldd` output has no missing dependencies.
The generic `libalphaLimitVelocity.so` entry still needs its short launch-time
lookup/dependency check below. The collector did not expand `globalVars`;
print that unchanged case-local file with the launch evidence.

The actual settings resolve several previously open questions:

| Setting | Observed value and interpretation |
|---|---|
| Integration | Euler, maximum timestep `5e-5`, adaptive `maxCo 0.7`, six alpha subcycles, two outer/two pressure corrections; keep these for regression. |
| Species | Upwind; tolerance `1e-12`, relative tolerance zero; equation residual alpha `1e-5` versus phase residual alpha `1e-6`. The residual-alpha stabilizer is in the captured equation but not separately itemized in the inventory budget. |
| PB | One correction, `scale true`, final-iteration solve, source refresh every solve; no configured diameter drift. The 2 mm nucleus is exactly the f3 pivot. Final-iteration hooks do not measure every outer-iteration scaling action. |
| Free surface / transfer | Liquid–headspace compression is configured. Headspace diameter is 10 mm. No `diffusiveMassTransfer` blending override, so default 0.3/0.7 applies. Film thickness 10 micrometres and foam multiplier zero remain experimental. |
| External O2 | Sparger injects O2 mass fraction 0.23. Liquid `O2sink` uses the saved `our` field with `useMonadSink false`, an explicit source without an availability limiter. Clipping may compensate requested consumption exceeding supply; measure it. |
| Schedule | `fluid_update_time 200` CFD seconds, `reaction_update_time 4` hours per biological update, `ph_update_time 0.1` CFD seconds. No biological update is due over this 31.2–31.7 s continuation; pH updates remain active. |

Root ML asset lists are empty. The inspected loaders do not throw on missing
files at startup; model inference is inside the biological update, and the
local `wellMixed.H` call is skipped. This does not block this short diagnostic,
but required model/scaler assets must be provided and verified before any
future run crosses a biological update. Do not extend this run to test that.
The frozen-OUR comparison is not the nonreacting 9.9 s validation release.

Run **control only first**, inside the active allocation. This checks the
remaining library, uses the existing build and pair, and refuses an existing
log. Noclobber and the case working directory are confined to a subshell;
there are no interactive shell-terminating commands. No pull or rebuild is
needed for this launch because the current changes are documentation only.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
AUDITROOT=/scratch/nwimer/sc1_conservation_audit_20261106T183120Z
PAIRROOT=/scratch/nwimer/sc1_conservation_pair_20261106T201549Z
AUDITEXE="$AUDITROOT/bin/putidaFoamSC1ConservationAudit"
CONTROL="$PAIRROOT/control"

if source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" && module load libfabric; then
    SC1_LIMITER=
    IFS=: read -r -a SC1_LIBDIRS <<< "${LD_LIBRARY_PATH:-}"
    for SC1_DIR in "${SC1_LIBDIRS[@]}" "${FOAM_USER_LIBBIN:-}" "${FOAM_SITE_LIBBIN:-}" "${FOAM_LIBBIN:-}"; do
        if [[ -n "$SC1_DIR" && -f "$SC1_DIR/libalphaLimitVelocity.so" ]]; then
            SC1_LIMITER="$SC1_DIR/libalphaLimitVelocity.so"
            break
        fi
    done
    if [[ -z "${SLURM_JOB_ID:-}" || ! -x "$AUDITEXE" || ! -x "$PY" \
          || ! -d "$CONTROL/processors52/31.199982372783058" \
          || ! -f "$PAIRROOT/comparison_manifest.json" || ! -f "$PAIRROOT/preflight.txt" ]]; then
        echo "Stopped: allocation or prepared artifacts missing."
    elif [[ -e "$CONTROL/log.conservation" || -L "$CONTROL/log.conservation" ]]; then
        echo "Stopped: existing log; preserve this case."
    elif [[ -z "$SC1_LIMITER" ]]; then
        echo "Stopped: libalphaLimitVelocity.so not found in OpenFOAM library paths."
    elif ! SC1_DEPS=$(ldd "$SC1_LIMITER" 2>&1) || [[ "$SC1_DEPS" == *"not found"* ]]; then
        printf 'Stopped: limiter dependencies unresolved.\n%s\n' "$SC1_DEPS"
    else
        sha256sum "$AUDITEXE" "$SC1_LIMITER"
        cat "$CONTROL/constant/globalVars"
        (cd "$CONTROL" && set -C && \
            srun -n 52 "$AUDITEXE" -parallel -fileHandler collated > log.conservation 2>&1)
        SC1_STATUS=$?
        echo "SC1 control solver status: $SC1_STATUS"
        "$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" "$CONTROL/log.conservation"
        tail -n 60 "$CONTROL/log.conservation"
    fi
else
    echo "Stopped: OpenFOAM environment setup failed."
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

## First control conservation result (2026-11-06)

The uploaded command transcript reports **solver status 0**, `End` and
`Finalising parallel run`. The log reader accepted the complete file. The
control reached `31.6999975003429` from `31.1999823727831` in **10,279 steps**,
an actual interval of `0.50001512756` s; use this duration for rates. Reported
`ExecutionTime` is 3508.66 s. Preserve this completed control and its log.

The limiter library resolved at
`/home/nwimer/OpenFOAM9_RHEL9/user-9/platforms/linux64GccDPInt32Opt/lib/libalphaLimitVelocity.so`,
SHA256 `e63f7be7cee3503d3d1c7b9bb5131c791e5645b07bff7e061a3cca9a5ce2d5ac`.
The audit executable hash is unchanged. The final legacy bubble inventory
`0.00961191 kg`, deposition `0.00985781 kg/s`, bubble–liquid transfer
`9.5896e-05 kg/s`, headspace–liquid transfer `3.64407e-06 kg/s` and endpoint
closure `0.00170228 kg/s` match the preserved zero-control record at displayed
precision. The complete old/new trajectory comparison is still to be collected;
the submitted tail contains only the final write.

### Oxygen correction dominates the liquid budget

| Cumulative liquid O2 term | mg over the actual interval |
|---|---:|
| External source (requested sink) | −248.931443 |
| Interphase transfer | +49.046806 |
| Measured final-PIMPLE species correction | +200.123373 |
| Actual inventory change | +0.240884 |
| Remaining inventory residual, after correction | +0.002147 |

The measured correction offsets **80.393% of requested removal**. The remaining
48.808070 mg is the net accounted sink after that correction, **not a measured
biological uptake rate**. Average interphase supply is 0.09809065 g/s, versus
requested removal 0.497848 g/s. Liquid O2 rises from approximately 19.3371 to
19.5780 mg and final DO is `0.01958567 mg/L`. Direct pre-advance oxygen jumps
and captured chemistry `R(O2)` are zero in every final cumulative phase record.

The correction hook measures O2 immediately before/after `fluid.correctSpecies()`.
The reviewed [OpenFOAM 9 normalise implementation](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/src/thermophysicalModels/reactionThermo/mixtures/basicSpecieMixture/basicSpecieMixture.C)
clips active species to zero and assigns the default-species remainder; O2 is
active and the liquid default is H2O. Thus the addition is consistent with
negative-O2 clipping after the unrestricted saved-OUR sink. The installed
base-mixture source and negative-cell statistics have not been included in
this transcript, so distinguish this source-supported attribution from direct
cellwise measurement. Hooks capture final PIMPLE only. Do not count correction
as physical OTR or infer biological validity from the small corrected residual.

### Phase-mass and gas-species closure remain unresolved

| Phase | Mass residual (mg) | O2 residual (mg) |
|---|---:|---:|
| Liquid | −1168.686138 | +0.002147 |
| Bubble gas | +386.460926 | −3.598547 |
| Headspace gas | −394.800287 | +4.167649 |
| Combined gas | −8.339360 | +0.569103 |
| Total | −1177.025498 | +0.571250 |

Bubble mass increased 486.600599 mg, whereas sampled injection/transfer/flux
terms predict only 100.139673 mg. Its residual is 79.4% of the observed gain
and 8.84% of injected gas. The opposite headspace residual makes phase
allocation/source-stage consistency a leading hypothesis, not a proven cause.
Liquid gained 355.705276 mg versus expected base addition 1475.344635 mg plus
gas transfer 49.046778 mg. This constant-density liquid discrepancy deserves
an actual alpha-source/flux audit even though it is small relative to 1000 kg.
The final 0.1 s window still has bubble/headspace residuals +77.920/−79.686 mg
and liquid −197.746 mg: the discrepancy is not confined to startup.

The total sampled interphase mass integral is effectively zero
(`3.7364e-20 kg`), as expected from equal/opposite source assembly. That does
not prove the alpha equations used those same amounts. Total evaluated O2
transfer imbalance is only +0.011744 mg; it does not explain the total O2
residual +0.571250 mg. The integrated equation residual L1 is approximately
`2.36e-5 mg`, much smaller again. Headspace boundary face-Y proxy transport
agrees with the assembled advection integral to printed precision. Simply
tightening linear species tolerances is not the supported next fix.

Capture the actual alpha `Su/Sp` actions and subcycle fluxes, alpha-clamp
changes, pressure-density inventory changes and species residual-alpha
stabilization in the next targeted diagnostic revision, if the paired records
confirm these discrepancies. Do not change coefficients while attributing them.

### PB correction is not phase-mass removal

The final `SC1_PB` reports cumulative correction `−0.007054524918 kg` in the
moment `integral(alphaB*rhoB*sum(fi))`. This sums repeated normalization of
size fractions, not removals from `alphaB*rhoB`. At the final step the moment
changes from 9.612637089 g by −0.000730840 g to 9.611906249 g, close to actual
final phase mass 9.611905230 g. The later pressure-density stage explains why
these need not be exactly simultaneous quantities; it has not been separately
measured. The cumulative correction exceeds integrated deposition 4.224152 g,
so deposition alone has not explained it. PB solve/source moments, nonlinear
convergence, stabilization, absent drift and scaling stages remain to be audited.
Do not add this PB correction as a separate physical sink in the phase ledger.

### Feed and outlet interpretation

Printed `globalVars` confirms constant liquid density 1000 kg/m3, temperature
303.15 K, surface tension 0.07 N/m and a 0.41-vvm feed converted with bottom
hydrostatic gas density. Actual injection is 8.742768 g/s with O2 mass fraction
**0.23**, as independently obtained from the captured mass/species sources.
The `f_O2 1` transport-weighting entry is not evidence of pure-O2 injection.
Clarify the experimental flow-meter reference pressure/temperature later:
the same volume rate at 101325 Pa and 303.15 K gives about 7.9581 g/s, whereas
the present bottom-pressure convention is about 9.86% higher. Leave this matched
pair unchanged; neither convention is established as experimental truth here.

The final gross headspace outlet flows are 1.1272244 kg/s inward and
1.1358814 kg/s outward, with net outward 0.00865698 kg/s. The old pressure
diagnostic's `headspaceOut` is signed net, so it agrees with the new gross
ledger. Large bidirectional flow still needs physical boundary/domain-height
validation; do not disable backflow merely to reduce the gross values.

### Next: complete the unchanged normal diagnostic

The runtime and endpoint-regression evidence support completing the matched
normal diagnostic with the existing executable. This is not approval of
conservation or entrainment physics. No rebuild or coefficient/schedule change
is needed. Run once inside an active allocation, then return the summary and
the extracted records from both logs. The final loop also prints the preserved
old control's legacy records for trajectory comparison. It only reads logs.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
PAIRROOT=/scratch/nwimer/sc1_conservation_pair_20261106T201549Z
AUDITEXE=/scratch/nwimer/sc1_conservation_audit_20261106T183120Z/bin/putidaFoamSC1ConservationAudit
NORMAL="$PAIRROOT/normal"

if source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" && module load libfabric; then
    if [[ -z "${SLURM_JOB_ID:-}" || ! -x "$AUDITEXE" || ! -x "$PY" \
          || ! -d "$NORMAL/processors52/31.199982372783058" ]]; then
        echo "Stopped: allocation, executable or prepared normal case missing."
    elif [[ -e "$NORMAL/log.conservation" || -L "$NORMAL/log.conservation" ]]; then
        echo "Stopped: existing normal log; preserve this case."
    elif ! "$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" "$PAIRROOT/control/log.conservation"; then
        echo "Stopped: completed control log did not pass the reader."
    else
        (cd "$NORMAL" && set -C && \
            srun -n 52 "$AUDITEXE" -parallel -fileHandler collated > log.conservation 2>&1)
        SC1_STATUS=$?
        echo "SC1 normal solver status: $SC1_STATUS"
        "$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" \
            "$PAIRROOT/control/log.conservation" --normal "$NORMAL/log.conservation"
        for SC1_CASE in control normal; do
            printf '\nAUDIT RECORDS: %s\n' "$SC1_CASE"
            awk '/^SC1_/ {if ($1 != "SC1_BOUNDARY" || /patch=outlet /) print}
                 /^(BUBBLE_MASS_AUDIT|BUBBLE_STAGE_AUDIT|HEADSPACE_O2_AREA_AUDIT|HEADSPACE_ENTRAINMENT.*_AUDIT) / {print}' \
                "$PAIRROOT/$SC1_CASE/log.conservation"
        done
        for SC1_LOG in /scratch/nwimer/sc1_headspace_entrainment_20261104/run_control_05/log*; do
            if [[ -f "$SC1_LOG" ]]; then
                printf '\nPRESERVED CONTROL: %s\n' "$SC1_LOG"
                awk '/^BUBBLE_MASS_AUDIT / {print}' "$SC1_LOG"
            fi
        done
        tail -n 30 "$NORMAL/log.conservation"
    fi
else
    echo "Stopped: OpenFOAM environment setup failed."
fi
```

If the reader rejects report-time matching, preserve both runs and return the
records: investigate actual times rather than silently interpolating or rerunning.
Do not rerun the completed control. A later biology implementation must distinguish
requested from realized oxygen-limited uptake and preserve coupled reaction
stoichiometry; the observed clipping is not permission to weaken OUR, tune OTR
to match it, or change the four-hour update schedule.

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
