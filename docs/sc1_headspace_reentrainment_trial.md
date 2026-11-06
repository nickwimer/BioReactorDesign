# SC1: conservative headspace-to-bubble entrainment — isolated development trial

**Status (2026-11-06): isolated 52-rank numerical trials, including
a full 0.5 s normal-inward entrainment sensitivity, have completed.
The model produces positive local PB bubble birth but only a tiny OTR
gain and no persistent bubble inventory improvement. Time-integrated
gas and O₂ inventory conservation, unresolved foam morphology, and
bubble-size physical validity remain unverified. This is not an
approved production closure or common nonreacting checkpoint.**

The previous headspace-liquid Henry/film model is compiled and runtime
tested. The two 0.5-second baseline/foam trials showed total OTR only
~0.098–0.106 g/s against ~0.498 g/s frozen OUR. The top
`O2.headspaceGas` reverse-flow BC is `inletOutlet` with
`inletValue uniform 0.23`, appropriate for air.
The remaining hypothesis is that large-scale free-surface
suction entrains gas into additional dispersed bubbles at
the upper impeller. This hypothesis must be checked against
experimental gas hold-up, video morphology and bubble sizes.

## Why this requires two coupled sources, not a reverse phaseTransfer entry alone

BiRD's `PhaseTransferPhaseSystem` conserves both gas-mixture
mass and donor gas species across an ordered pair, but **its
populationBalanceModel does not automatically populate the
newly created bubbleGas size groups from a generic phase
transfer source**.

Also, `BlendedInterfacialModel::dmdtf()` uses signed
blending of the two opposite *ordered* submodels; simply
appending `(headspaceGas in bubbleGas)` beside the
existing `(bubbleGas in headspaceGas)` deposition risks
an incorrect sign or suppression of the existing deposition
correlation.

The isolated solution combines:

1. **`headspaceEntrainingDeposition`**: a scratch-only
   `phaseTransferModel` installed in the **same original
   ordered** `(bubbleGas in headspaceGas)` entry, *replacing
   only the experimental copy's* `type deposition`. It
   computes the exact original deposition correlation,
   plus a new countervailing positive bubble-production
   term. With `entrainmentEfficiency 0`, its dmdtf is
   algebraically identical to the original deposition
   term *before the unchanged BiRD blending*.
2. **`headspaceEntrainmentNucleation`**:
   a scratch-only `nucleationModel` added to the
   existing `populationBalanceCoeffs/bubbles/nucleationModels`
   list. It obtains the **signed net** bubbleGas mass
   gain from `fluid.dmdtf(bubbleGas, headspaceGas)`,
   clips negative gain (deposition), converts the
   positive remainder to volume source using `rho_bubbleGas`,
   then distributes it to the existing size classes via
   BiRD's `popBal.eta(i, V_nucleus)`.
   Default `nucleationDiameter 0.002 m` is one of the
   configured original size groups, not a measured
   entrained-bubble size.

The two models are built into **a separate library** in
`OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceEntrainmentModels`,
installed as
`/scratch/nwimer/sc1_headspace_entrainment_20261104/lib/`
`libsc1HeadspaceEntrainmentModels.so`.
The previously validated
`libsc1HeadspaceFilmModels.so` is NOT modified or
rebuilt. The original solver and production cases remain
untouched. The scratch `putidaFoamMask01BubbleAudit`
can be rebuilt independently to log both transfer channels.

## Experimental model equations and locality

The **raw candidate entrainment rate** in an eligible cell is

```text
a_resolved  = |grad(alpha.headspaceGas)| [m^-1]
g(alpha)   = clamp((alpha-alphaMin)/(alphaFull-alphaMin), 0, 1)
S_H_to_B   = rho_H * C_e * g(alpha_H) * g(alpha_L)
             * a_resolved * min(|U_L-U_H|, u_slipMax) [kg/m3/s]
```

There is also a per-cell donor limit:

```text
S_H_to_B <= rho_H * max(alpha_H-alphaMin,0) / tau_donor
```

Before inserting the source, the model sums `S_H_to_B*V`
over all MPI ranks, then applies

```text
globalScale = min(1, capTotal_kg_s / integral(S_H_to_B dV))
S_H_to_B,applied = globalScale * S_H_to_B
```

This cap applies to the **raw, unblended** entrainment.
BiRD's original phaseTransfer blending is subsequently
applied, and the actual signed transfer is the algebraic
sum of deposition and entrainment.

The trial only considers cells with:

- `0.68 <= y <= 0.96 m`, `0.05 <= r <= 0.30 m`
  about the vertical y-axis;
- appreciable local headspace and liquid phase fractions
  (`alphaMin 0.02`, `alphaFull 0.10`);
- `alpha.bubbleGas < 0.10` to restrict the source to
  a mostly not-yet-dispersed interface;
- **downward-moving liquid** `U.liquid_y < 0`;
- relative headspace/liquid slip capped at 5 m/s.

These are engineering hypotheses and exploratory local
guards, **not empirically calibrated entrainment physics**.
Their impact must be inspected on saved Eulerian fields.
No re-entrainment occurs elsewhere or from the entire
headspace volume.

Defaults for the first test:

```foam
entrainmentEfficiency 0;      // control: exactly original deposition
maxGlobalEntrainedKgS 0.002;  // experimental 2 g/s raw global ceiling
nucleationDiameter 0.002;     // m, size class, not measured
foamAreaMultiplier 0;         // existing headspace film: resolved only
```

A second test uses `entrainmentEfficiency 0.01`
with exactly the same 2 g/s raw ceiling and other fields.
This is **sensitivity analysis**, not an attempt to force
the 0.498 g/s OUR. The cap is a safety bound, not a
prescribed minimum gas re-entrainment rate.

### Mass, oxygen species and size-group accounting

The phaseTransfer model returns signed `dmdtf`
for the bubbleGas↔headspaceGas pair. Because both phases
are treated by BiRD's phase-system mass equations,
positive mass assigned to bubbleGas is taken
from headspaceGas, and donor O2 mass fraction enters
the bubbleGas species equation via
`PhaseTransferPhaseSystem::addDmdtYfs()`.
**No independent oxygen mass is created.**

For a bubble created at volume `V_nucleus`,
`nucleationRate_i = eta(i,V_nucleus) *
max(0, dmdtf_bubbleGain)/(rho_B * V_i)`
assigns the **positive net** bubble production to
the size groups. The original breakup and coalescence
models then act on those classes.

**Known limitation:** clipping to positive NET
mass transfer, rather than independently tagging
the gross entrained and deposited subpopulations,
means a cell with simultaneous bubble deposition
and re-entrainment gets a size-bin birth source only
for the net positive amount. Existing PB
normalization handles net removal, but this is
an approximation to the actual simultaneous
birth/death size distribution. Before production,
validate the size-class volume fractions and
correlations against independent measurements.

## Kestrel: build only the separate scratch library

First load environment on the current compute node if needed:

```bash
source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc"
module load libfabric
```

Do not run top-level shell termination commands in an
interactive `salloc` session. The safe commands below
build only a separate scratch library, do not launch
CFD, and never modify production cases.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
BIRD="$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam"
SRC="$REPO/OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceEntrainmentModels"
LIB=/scratch/nwimer/sc1_headspace_entrainment_20261104/lib/libsc1HeadspaceEntrainmentModels.so

if git -C "$REPO" pull --ff-only; then
    if "$PY" "$REPO/scripts/diagnostics/test_sc1_headspace_entrainment_trial.py"; then
        mkdir -p "$(dirname "$LIB")"

        if (
            cd "$SRC" &&
            SC1_BIRD_ROOT="$BIRD" \
                wmake > build.sc1Entrainment.log 2>&1
        ); then
            if [[ -f "$LIB" ]]; then
                echo "ENTRAINMENT LIBRARY BUILD SUCCESS"
                ldd "$LIB" | grep -E 'bird|not found'
            else
                echo "Build reported success but library is missing."
            fi
        else
            echo "Build failed. No CFD launched."
            grep -nE 'error:|fatal error:|ld:|collect2' \
                "$SRC/build.sc1Entrainment.log" | head -35
        fi
    else
        echo "Regression checks failed; nothing compiled."
    fi
else
    echo "Git pull failed; nothing compiled."
fi
```

Do NOT proceed to trial creation until the new
`libsc1HeadspaceEntrainmentModels.so` links cleanly,
including `libbirdeulerianInterfacialModels.so` and
`libbirdphaseSystem.so`.

### First Kestrel library build confirmed (2026-11-04)

The user updated the feature branch to `ec9136e`. On Kestrel,
`test_sc1_headspace_entrainment_trial.py` passed both
configuration/PB-injection tests and bounded-source/global-cap
checks. The separate `sc1HeadspaceEntrainmentModels` C++ library
then **built and linked successfully**, and `ldd` resolved both
`libbirdeulerianInterfacialModels.so` and
`libbirdphaseSystem.so` under the expected RHEL9
`$FOAM_USER_LIBBIN` path.

**This establishes compilation/link success only.** No runtime
trial has been executed and the current prototype's signed
transfer, PB nucleation, species transport and final bubble mass
budget have not yet been verified in a solver run.

The next gate is to prepare **two independent scratch cases**
from the same 31.199982372783058 s source checkpoint,
using `--entrainment-efficiency 0` and `0.01`,
with `--foam-multiplier 0` for both. Inspect the two
`phaseTransfer` and `nucleationModels` entries and
load both experimental libraries. Rebuild only the isolated
scratch audit solver from the updated `bubbleMassAudit.H`,
then run the zero-efficiency equivalence case FIRST.
Do not change any original solver or production files.

## Paired scratch 52-rank tests — after compile and preflight

Run independent *0.5-second* continuations from the same
unchanged source scratch checkpoint
`31.199982372783058`. Use an explicitly fresh
`--out` path for every scenario.

First ensure scratch diagnostic solver includes
the new `netBubbleGasFromHeadspace_kg_s` label
(the newly pushed `bubbleMassAudit.H` has it).
Rebuild **scratch diagnostic executable only**, as
in the existing `docs/sc1_headspace_oxygen_transfer_design.md`.

Then create two trials, with *no* headspace foam
enhancement `--foam-multiplier 0`, from the SAME source:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
ROOT=/scratch/nwimer/sc1_headspace_entrainment_20261104
SOURCE=/scratch/nwimer/sc1_bubble_audit_20261029/run
EXE="$FOAM_USER_APPBIN/putidaFoamMask01BubbleAudit"

for SPEC in "control 0" "entrain 0.01"; do
    NAME=${SPEC%% *}
    EFF=${SPEC##* }
    RUN="$ROOT/run_${NAME}_05"

    if [[ -e "$RUN" ]]; then
        echo "$RUN already exists; not overwriting."
    elif "$PY" \
        "$REPO/scripts/diagnostics/prepare_sc1_headspace_entrainment_trial.py" \
        --source "$SOURCE" --out "$RUN" \
        --entrainment-efficiency "$EFF" \
        --max-entrained-kg-s 0.002 \
        --nucleation-diameter 0.002 \
        --foam-multiplier 0 --duration 0.5
    then
        echo "Inspect the new phaseTransfer and nucleationModels first."
        foamDictionary "$RUN/constant/phaseProperties" \
            -entry phaseTransfer -value
        foamDictionary "$RUN/constant/phaseProperties" \
            -entry populationBalanceCoeffs/bubbles/nucleationModels -value

        (
            cd "$RUN" &&
            srun -n 52 "$EXE" -parallel -fileHandler collated \
                > log.entrainment.smoke 2>&1
        )
        STATUS=$?
        echo "Scenario $NAME solver status $STATUS"
        grep -E 'HEADSPACE_ENTRAINMENT|BUBBLE_MASS_AUDIT|HEADSPACE_O2_AREA_AUDIT|FOAM FATAL|^End$' \
            "$RUN/log.entrainment.smoke" | tail -25

        if [[ "$STATUS" -ne 0 ]]; then
            echo "Stop: trial failure; do not run another scenario."
            break
        fi
    else
        echo "Stop: trial preparation failed."
        break
    fi
done
```

If the `foamDictionary` path to `nucleationModels`
is not available in this particular case, print only
the relevant subsection with `grep -n -A ...`
rather than editing it.

**Important:** preserve the old,
validated `run_geometric_05` and `run_foam_05`
as independent comparisons; never overwrite them.

## Zero-efficiency control passed on Kestrel (2026-11-04)

The user rebuilt the scratch audit solver and ran
`run_control_05` from the same collated
`31.199982372783058` state as the earlier
`run_geometric_05`, with
`entrainmentEfficiency 0`, `foamAreaMultiplier 0`,
`nucleationDiameter 0.002` m and raw cap
`0.002 kg/s`.

**The trial completed normally (solver status 0).**
At all five saved t=31.3–31.7 s audits, the numerical
values **exactly match the displayed precision**
of the original geometric-only trial for
`inventory_kg`, `injection_kg_s`,
`depositionLoss_kg_s`,
`assembledBubbleDmdt_kg_s`,
`netBubbleLiquidMassLoss_kg_s`,
`netHeadspaceLiquidMassLoss_kg_s`, and
`effectiveArea_m2`. The diagnostic
`HEADSPACE_ENTRAINMENT_PB_AUDIT` confirms
`positiveNetBubbleGain_kg_s=0` at all five times,
so the companion population-balance source
adds no bubbles when entrainment is disabled.

Example at 31.7 s, in both models:

```text
inventory_kg=0.00961191
depositionLoss_kg_s=0.00985781
assembledBubbleDmdt_kg_s=-0.0099537
netBubbleLiquidMassLoss_kg_s=9.5896e-05
netHeadspaceLiquidMassLoss_kg_s=3.64407e-06
effectiveArea_m2=1.72277
endpointClosure_kg_s=0.00170228
```

The newly logged
`netBubbleGasFromHeadspace_kg_s=-0.00985781`
is the signed net bubble mass gain from that
phase pair, hence negative for normal deposition.
It is **not** a gross re-entrainment measurement.
The equivalence at display precision is an
important regression check, but it does not
validate the full solver-discrete gas or
species mass balance; the endpoint residual
remains nonzero as in the original trial.

**Next gate:** a **fresh**, independent 0.5 s
`run_entrain_05` trial at
`entrainmentEfficiency=0.01`, with
`foamAreaMultiplier=0`, the same starting
checkpoint and unchanged 2 mm PB nucleus
diameter and 2 g/s raw global cap. The
`HEADSPACE_ENTRAINMENT_CANDIDATE` output
must be checked against `positiveNetBubbleGain_kg_s`.
An unblended raw entrainment candidate
may be nonzero yet the NET bubble creation
may remain zero due to local deposition and
phaseTransfer blending. Do not equate candidate
rates with actual bubble injection or new OTR.

**Do not overwrite** `run_control_05`,
`run_geometric_05` or `run_foam_05`.
If `run_entrain_05` already exists, inspect
it instead of recreating the run.

## Nonzero entrainment pilot completed: near-zero eligibility (2026-11-04)

The user ran the separate `run_entrain_05` at
`entrainmentEfficiency=0.01`, `foamAreaMultiplier=0`,
`maxGlobalEntrainedKgS=0.002`, a 2 mm PB nucleus,
from the same `31.199982372783058` scratch checkpoint
as the validated zero-efficiency control. The 52-rank
continuation **completed successfully** (status 0)
with five audit writes at t=31.3, 31.4, 31.5, 31.6,
and 31.7 s.

| Diagnostic | t31.3 | t31.4 | t31.5 | t31.6 | t31.7 |
|---|---:|---:|---:|---:|---:|
| Eligible cells | 2 | 7 | 0 | 4 | 0 |
| Raw entrainment (kg/s) | 7.82129e-7 | 1.5131e-7 | 0 | 3.69635e-6 | 0 |
| Positive NET bubble creation to PB (kg/s) | 0 | 8.86303e-9 | 0 | 6.50491e-7 | 0 |
| Candidate interface area (m²) | 1.36809e-5 | 2.61243e-6 | 0 | 6.37885e-5 | 0 |
| Headspace-to-liquid O₂ loss (kg/s) | 3.65348e-6 | 3.62201e-6 | 3.62579e-6 | 3.65850e-6 | 3.60418e-6 |

**Critical result:** At most **7 cells** met all criteria,
and two of five output times had no eligible cells.
The maximum raw candidate was **3.69635 mg/s**
(not grams/s), with maximum positive net
bubble gain **0.650491 mg/s**. The model
never reached the 2 g/s cap: `capScale=1`
whenever the candidate was nonzero.
At t31.6, its total eligible interface area
of `6.37885e-5 m²` was roughly
**0.0037%** of the concurrent
`1.72912 m²` headspace-liquid resolved
interface area reported by the O₂ film
model. The two area definitions differ
slightly by geometric normalizations and
the entrainment restrictions, but their
large scale separation is unambiguous.

Five-sample arithmetic means, in g/s:

| Rate | Zero control | Entrainment 0.01 |
|---|---:|---:|
| bubbleGas→liquid O₂ | 0.094754 | 0.094732 |
| headspaceGas→liquid O₂ | 0.003644 | 0.003633 |
| **Total OTR** | **0.098398** | **0.098365** |

Entrainment produced **no meaningful OTR improvement**
(~-0.000033 g/s difference, below the resolution
needed to infer a physical causal signal). This
does not establish that surface gas entrainment
is physically negligible: the experimental model's
**eligibility gates** are selecting almost
none of the simulated gas-liquid interface.
Overall bubble deposition is still ~7–10 g/s,
far greater than the raw entrainment candidate.
The logged `netBubbleGasFromHeadspace_kg_s`
remains negative at every output: this
is signed **net** bubbleGas mass loss, not
a gross entrainment measurement. Rare local
positive net production is confirmed by
`HEADSPACE_ENTRAINMENT_PB_AUDIT`.

### Next gate diagnosis: no coefficient changes

Before adjusting entrainment efficiency or
the global cap, determine which eligibility
condition removes most of the candidate
interface. A **diagnostic-only** update to
`headspaceEntrainingDeposition.C` adds the
following stage counters/areas to the existing
`HEADSPACE_ENTRAINMENT_CANDIDATE` line:

- `domainMixedArea_m2`: gated `|grad(alphaH)|`
  integrated throughout the domain.
- `windowCells`: cells within the prescribed
  y/r bounding window (before phase gates).
- `windowMixedCells`, `windowMixedArea_m2`:
  window cells containing both headspace
  and liquid within the phase-fraction gates.
- `windowLowBubbleCells`,
  `windowLowBubbleArea_m2`: those remaining
  after rejecting bubble-rich cells.
- `eligibleCells`, `interfaceArea_m2`:
  those remaining after requiring liquid
  `U_y<0` (the actual entrainment eligibility).

The C++ source **changes reporting only**; the
per-cell proposal, donor cap, global mass cap,
phase-transfer signs, and PB coupling remain
unchanged. Recompile only the separate scratch
entrainment library, then prepare a **fresh 0.2 s
scratch trial** from the same original source
using `--entrainment-efficiency 0.01`.
Do **not** overwrite the completed zero/nonzero
0.5 s tests or modify production cases.

If window area is near zero, inspect the
free-surface shape and upper impeller/cavity
geometry before moving bounds. If the phase
overlap gate eliminates most area, determine
whether the Eulerian phase separation
represents a coherent cavity rather than
the observed bubbly froth. If the
downward-liquid requirement eliminates most
of the area, inspect the sign of local
interfacial normal velocity and entrainment
physics rather than blindly removing the
gate. A cap/efficiency increase cannot
repair an empty eligible set.

## Gate audit resolves the bottleneck: downward-liquid velocity (2026-11-04)

User ran the **diagnostic-only** `run_gateaudit_05`
successfully (solver status 0, 52 ranks) from the
same 31.199982372783058 s source, with
`entrainmentEfficiency=0.01`, `foamAreaMultiplier=0`,
raw cap `0.002 kg/s`. It reproduces previous
nonzero-trial rates, confirming the intermediate
gate counters did not change the applied model.

| t (s) | Domain gated area (m²) | Window gated area (m²) | Low-bubble area (m²) | After downflow (m²) | Eligible cells |
|---|---:|---:|---:|---:|---:|
| 31.3 | 1.53462 | 0.331308 | 0.331308 | 1.36809e-5 | 2 |
| 31.4 | 1.52251 | 0.327969 | 0.327969 | 2.61243e-6 | 7 |
| 31.5 | 1.52675 | 0.329672 | 0.329657 | 0 | 0 |
| 31.6 | 1.53930 | 0.337778 | 0.337778 | 6.37885e-5 | 4 |
| 31.7 | 1.52077 | 0.334136 | 0.334136 | 0 | 0 |

18,432 cells lie in the axial/radial window.
About 4,000 cells per output pass the mixed-phase
gate, with virtually all remaining after the
bubble-fraction filter. Yet only 0–7 cells pass
`requireDownwardLiquid`.

Across the five output times, **99.99518% of
the window's eligible interfacial area is
discarded by the downward-liquid velocity
requirement**, retaining only **0.00482%**
(area-weighted). The window retains
~21.73% of the domain's gated interface
area before the velocity requirement.
Neither the prescribed zone nor the
bubble-rich exclusion is the main bottleneck.
The actual 2 g/s raw source cap remained
inactive throughout.

**Interpretation:** The criterion
`U.liquid.y() < 0` does not follow from
interfacial gas transport or topological
bubble formation. It rejects nearly all
computed interface contact, possibly because
the local liquid vertical velocity is upward,
even while a nearby free surface is deforming.
This does **not** prove gas is physically
entraining: a finite positive phase-normal
relative velocity is only a more relevant
candidate proxy, still requiring breakup,
interfacial deformation, and experimental
validation.

### Next diagnostic, no physics change

Experimental source `headspaceEntrainingDeposition.C`
now logs these **additional counterfactual statistics**
inside the same upper-impeller and low-bubble
window, without applying them to phase equations:

- `downwardLiquidArea_m2` and
  `nonDownwardLiquidArea_m2`: interface
  area grouped by sign of `U.liquid.y`;
- `normalGasIntoLiquidArea_m2` and
  `normalGasAwayArea_m2`: areas grouped by
  sign of `(U.liquid-U.headspaceGas) dot n`,
  `n=grad(alpha.headspaceGas)/|grad|`;
- `meanNormalGasIntoLiquidSpeed_m_s` and
  `meanLiquidVerticalSpeed_m_s`:
  area-weighted signed velocity proxies;
- `diagnosticRawNoDownGate_kg_s`:
  hypothetical *unblended* raw entrainment
  with the same coefficient and donor limit
  but without the vertical-velocity gate;
- `diagnosticRawNormalInward_kg_s`:
  hypothetical *unblended* entrainment using
  positive normal relative speed as both
  eligibility gate and velocity magnitude.

The gas-into-liquid sign convention is
positive `(U_L-U_H) dot n` where `n`
points from liquid toward headspace.
This is a relative-velocity diagnostic,
**not** a validated criterion for gas
breakup or a permission to change the model.
The original `requireDownwardLiquid true`
source remains unchanged in this diagnostic
run. A Python regression tests the normal
velocity sign convention and shows that
downward `U_L.y` is neither necessary nor
sufficient for inward relative gas motion.

Rebuild **only the experimental entrainment
library** (do not rebuild the film library or
audit solver), prepare a fresh **0.2-second**
scratch run from the original
`31.199982372783058` checkpoint, and
read the new `HEADSPACE_ENTRAINMENT_CANDIDATE`
statistics at t≈31.3. Compare the no-down
and normal-inward hypothetical rates before
enabling any new source physics.

If normal-inward flow is also scarce, do
not bypass both kinematic checks. The
experimental video might depict a foam
interface or trapped gas, not actual
dispersed bubble breakup. Instead assess
gas-phase topology and near-impeller
vorticity/pressure before defining a
different entrainment model.

## Interface-normal gate diagnostic and opt-in trial (2026-11-04)

The user ran the diagnostic-only normal-motion comparison
`run_normalgate_02` (0.2 s from original collated
`31.199982372783058` source, 52 ranks) successfully.
The reported values at t≈31.3 s were:

| Measurement | Value |
|---|---:|
| Cells in bounded impeller window | 18,432 |
| Mixed-phase cells within window | 3,993 |
| Mixed window area after bubble exclusion | 0.331308 m² |
| Downward-liquid cells | 2 |
| Downward-liquid area | 0.0000136809 m² |
| Non-downward-liquid area | 0.331294 m² |
| **Normal gas-into-liquid cells** | **1,765** |
| **Normal gas-into-liquid area** | **0.163135 m²** |
| Normal gas-away area | 0.168172 m² |
| Area-weighted signed normal relative speed | -0.194645 m/s |
| Area-weighted liquid vertical speed | +1.33425 m/s |
| Hypothetical raw entrainment without downflow gate | 0.0189264 kg/s |
| **Hypothetical raw normal-inward entrainment** | **0.00263105 kg/s** |
| Current 0.01-efficiency downflow-gated raw entrainment | 7.82129e-7 kg/s |

This means ~49.24% of the available window mixed-interface
area has positive gas-to-liquid relative normal motion,
while only ~0.00413% survives the earlier
`U.liquid.y < 0` condition. The global *mean*
relative normal speed is negative, despite locally positive
normal motion over ~half of the area. A speed-proportional
normal-inward model would propose **2.631 g/s raw**,
then the configured 2 g/s global raw cap would scale
the proposal by ~0.760. These are diagnostics only,
not validated real bubble entrainment rates.

**New explicit scratch-only model switch**:

```foam
entrainmentVelocityCriterion legacyDownwardLiquid;   // old default
requireDownwardLiquid true;
```

Or, in an **independently prepared** scratch trial:

```foam
entrainmentVelocityCriterion normalGasIntoLiquid;    // OPT IN
requireDownwardLiquid false;
```

The code hard-rejects `normalGasIntoLiquid` combined with
`requireDownwardLiquid true`. This prevents accidental
use of two mutually exclusive kinematic gates.
The trial helper adds a
`--velocity-criterion normalGasIntoLiquid` option,
defaulting to `legacyDownwardLiquid` for all
existing workflows. The new model uses
`max(0,(U_L-U_H) dot grad(alpha_H)/|grad(alpha_H)|)`
as the effective relative penetration velocity, capped at
`maxSlipSpeed`, with unchanged phase-fraction,
spatial-window, donor-cap and global-rate limits.
The old zero-efficiency and legacy nonzero modes are
unchanged by default.

**Next numerical gate, not physics validation:**

1. Compile just the scratch `libsc1HeadspaceEntrainmentModels.so`
   library on Kestrel and run the updated Python preflight.
2. From the same original scratch checkpoint, prepare a
   **fresh** `run_normalactive_02` with
   `--entrainment-efficiency 0.01`,
   `--velocity-criterion normalGasIntoLiquid`,
   `--max-entrained-kg-s 0.002`,
   `--nucleation-diameter 0.002`,
   `--foam-multiplier 0`, and `--duration 0.2`.
   Never overwrite the completed zero control or
   other geometric/foam/entrainment trials.
3. Verify the runtime log has
   `velocityCriterion=normalGasIntoLiquid`,
   positive `eligibleCells`, raw candidate and
   capped rate with `capScale<=1`. Compare new
   PB positive net bubble gain, bubble inventory,
   bubble-liquid and headspace-liquid O₂ flux,
   deposition and endpoint residual against the
   existing zero control. Positive **gross**
   entrainment may remain locally offset by
   deposition so net PB creation need not equal
   the 2 g/s cap.
4. Confirm finite `d.bubbleGas`, normalized
   size groups and conserved species inventories
   before extending to a 0.5 s numerical comparison.
   If runtime instability or physically unacceptable
   gas holdup appears, halt investigation of
   nonzero entrainment rather than force the coefficient.

**Science caveat:** Eulerian interphase relative normal speed
is not the true moving-interface velocity, and inward motion
does not by itself establish detachment into dispersed
2 mm bubbles. The experimental foam video, gas holdup,
breakup and oxygen-transfer measurements must be used to
decide whether and where this closure is physical. This
is still a scratch-only, uncalibrated prototype; no
production configuration or solver has changed.

## First opt-in normal-direction runtime: numerical operation demonstrated (2026-11-04)

User completed the new `run_normalactive_02` successfully
(`Solver status: 0`, `End`) on 52 MPI ranks, starting
from the **same** collated checkpoint `31.199982372783058`
and ending around 31.4 s. The scratch library compiled,
the Python regression checks passed, and
`entrainmentVelocityCriterion normalGasIntoLiquid` /
`requireDownwardLiquid false` were verified in
`constant/phaseProperties`.

Configuration: `entrainmentEfficiency=0.01`,
`maxGlobalEntrainedKgS=0.002` (2 g/s raw unblended cap),
`nucleationDiameter=0.002` m,
`foamAreaMultiplier=0`. The old zero-entrainment
`run_control_05` provides matched-time controls.

| Quantity | 31.3 s normal | 31.4 s normal |
|---|---:|---:|
| Eligible normal-inward cells | 1,791 | 1,875 |
| Eligible interface area (m²) | 0.162999 | 0.168764 |
| Raw normal-inward entrainment (kg/s) | 0.00257137 | 0.00266542 |
| Global cap scale | 0.777796 | 0.750352 |
| **Capped, still UNBLENDED entrainment (kg/s)** | **0.002** | **0.002** |
| Positive **NET** local bubble gain for PB (kg/s) | 0.000227563 | 0.000207738 |
| Signed NET bubble↔headspace rate (kg/s) | -0.00713551 | -0.00753911 |
| Bubble inventory (kg) | 0.00934556 | 0.00949945 |
| BubbleGas→liquid O₂ flux (kg/s) | 9.38776e-5 | 9.45063e-5 |
| HeadspaceGas→liquid O₂ flux (kg/s) | 3.63763e-6 | 3.62483e-6 |
| Bubble-balance endpoint residual (kg/s) | N/A | 0.000429673 |

**Matched control comparison:**

| t (s) | Bubble inventory control→normal (kg) | Bubble→liquid O₂ control→normal (g/s) | Headspace→liquid O₂ control→normal (g/s) | Total OTR gain (g/s) |
|---|---|---|---|---:|
| 31.3 | 0.00933447 → 0.00934556 | 0.0935802 → 0.0938776 | 0.00365347 → 0.00363763 | +0.00028156 |
| 31.4 | 0.00947528 → 0.00949945 | 0.0941079 → 0.0945063 | 0.00362237 → 0.00362483 | +0.00040086 |

The bubble mass increases by **11.09 mg** at
t31.3 and **24.17 mg** at t31.4 relative to
the matched control. The OTR gain at these
very early outputs is only **0.28–0.40 mg O₂/s**,
far smaller than the ~0.498 **g/s** previously
measured frozen biological OUR; do not
confuse milligrams with grams. The 0.2 s smoke
test is not long enough for entrained bubbles
to establish a new residence-time distribution
or quasi-steady OTR.

**Important rate distinctions:**

- The 2 g/s cap applies to the **candidate
  mixture mass transfer before BiRD blending**.
  It is NOT 2 g/s of new 2 mm bubbles delivered
  into the PB and certainly not 2 g O₂/s.
- `HEADSPACE_ENTRAINMENT_PB_AUDIT` reports the
  sum of **local positive NET** bubble creation
  after deposition and blending: ~0.21–0.23 g/s.
  This is about 10–11% of the unblended cap.
  The source may be offset by deposition
  elsewhere, and positive net PB production
  does not by itself prove complete population
  size-moment conservation.
- Signed integrated
  `netBubbleGasFromHeadspace_kg_s` remains
  negative (~-7.1 to -7.5 g/s): deposition still
  dominates global gas exchange.
  Do **not** call this a gross entrainment rate.
  The t31.4 net deposition differs by only
  +0.02954 g/s from the zero control and its
  direction reflects changed hydrodynamics,
  not a direct measurement of 2 g/s transfer.
- The t31.4 endpoint closure residual in the
  active normal case is `0.000429673 kg/s`
  vs `0.000268565 kg/s` in the control.
  Both are instantaneous endpoint
  approximations rather than time-integrated
  discrete solver residuals. We cannot
  certify conservation from these values.
- The near-interface area-weighted liquid
  vertical velocity remains upward (~+1.35 m/s),
  while the averaged signed gas-into-liquid
  normal velocity is negative
  (~-0.15 m/s), despite positive normal
  velocity in ~half the interface. These
  Eulerian velocity proxies need independent
  validation against the real impeller/free
  surface; a net gas-film displacement
  is not necessarily bubble detachment.

**Next action:** prepare a *fresh* 0.5 s
`run_normalactive_05` from the original
source, with exactly the above configuration,
and compare all five saved times against
`run_control_05`. Do not extend by overwriting
`run_normalactive_02`. Inspect
`bubbleGas sizeGroups-sum volume fraction`
min/max and `bubbleGas Sauter mean diameter`
in the solver log, along with positive PB
bubble creation, OTR, inventory and endpoint
closure. If size-bin sum drifts substantially
outside [0,1] or diameter becomes non-finite,
do not use the model for calibration.

**Science decision:** This is an important
functional validation of controlled
headspace→bubble conversion with PB
nucleation, but not evidence that
surface entrainment can close the
~0.4 g O₂/s deficit. Validate the PB
size distribution and gas/liquid contact
physics before adjusting `C_e` or the 2 g/s
candidate limit. Production cases, original
BiRD libraries and experimental foam film
library stay untouched.

## Normal-inward PB numerical sanity gate passed (2026-11-06)

The user inspected the complete short `run_normalactive_02` and
matched `run_control_05` solver logs on a **new Kestrel
allocation** (node `x1006c3s5b1n1`; prior runs were already
finished and preserved). The last eight reported
`bubbleGas sizeGroups-sum volume fraction, min, max`
lines were **`1 1 1`** for *both* runs.
The last eight `bubbleGas Sauter mean diameter`
lines were finite and stable:

| Last reported statistic | `run_control_05` | `run_normalactive_02` |
|---|---:|---:|
| PB fraction mean / min / max | 1 / 1 / 1 | 1 / 1 / 1 |
| Sauter mean diameter (m) | 0.00120099 | 0.00120567 |
| Sauter minimum diameter (m) | 0.00100252 | 0.00100249 |
| Sauter maximum diameter (m) | 0.00532335 | 0.00536886 |

Mean Sauter diameter rises by only ~0.39% in the
new early entrainment run; the maximum rises ~0.86%.
The search for `FOAM FATAL`, `floating point exception`,
and `segmentation fault` returned only
`sigFpe : Enabling floating point exception trapping
(FOAM_SIGFPE).`, which is **normal startup output**,
not an actual exception.

**Numerical gate passed:** no evidence of
size-group normalization failure, non-finite diameters,
or fatal solver errors. However, BiRD's
`velocityGroup::scale()` applies `max(0)` to size groups
and renormalizes them on correction, so `1 1 1` is
not an independent demonstration of unmodified
pre-normalization number/volume moment conservation.
The positive net nucleation source also has not been
independently integrated against PB size-bin mass
to prove a fully conservative population moment.

**Next:** prepare a **new, non-overwriting 0.5-second**
scratch run `run_normalactive_05` from the *same*
source checkpoint `31.199982372783058`,
with `entrainmentEfficiency=0.01`,
`entrainmentVelocityCriterion normalGasIntoLiquid`,
`maxGlobalEntrainedKgS=0.002` kg/s, 2 mm
nucleation diameter, and
`foamAreaMultiplier=0`. Compare the five
t31.3–31.7 outputs with `run_control_05`.
Inspect O₂ supply, bubble inventory, net/headspace
mass exchange, PB positive production, diameter
stability, and endpoint residuals. The run is a
short-term numerical sensitivity, **not** a
validated physical closure or nonreacting 9.9 s
checkpoint.

## Full 0.5 s opt-in normal-inward trial completed; OTR remains low (2026-11-06)

**Source:** User's Kestrel 52-rank `run_normalactive_05/log.normalactive05`
output, from the same `31.199982372783058` collated source checkpoint
as `run_control_05`, using
`entrainmentVelocityCriterion normalGasIntoLiquid`,
`entrainmentEfficiency 0.01`, `nucleationDiameter 0.002` m,
`maxGlobalEntrainedKgS=0.002` kg/s, and
`foamAreaMultiplier=0`. The Python preflight passed and the trial
**completed with solver status 0** and 5 complete writes at t31.3–31.7 s.
No production files were changed.

At each of the five writes the model reported
1,791–1,900 normal-inward eligible cells, ~0.160–0.169 m²
eligible area, raw candidate 0.00255–0.00269 kg/s,
and active global cap at exactly 0.002 kg/s
(unblended candidate, NOT final dispersed bubble mass).
`capScale` varied 0.743–0.784. Positive local
**NET** bubbleGas creation allocated to the PB
was 0.000208–0.000238 kg/s (arithmetic mean
**0.00022310 kg/s = 0.2231 g/s**, ~11.2% of
the gross unblended capped candidate).
Global signed `bubbleGas↔headspaceGas`
mass exchange remained net deposition
(negative net bubble gain) of ~7–10 g/s.

**Matched-time data against the original zero-efficiency
`run_control_05`:**

| Time (s) | Bubble inventory control (kg) | Active (kg) | Active-minus-control inventory (mg) | Active total OTR (g/s) | Active-minus-control OTR (mg/s) | Positive NET PB bubble birth (g/s) |
|---|---:|---:|---:|---:|---:|---:|
| 31.3 | 0.00933447 | 0.00934556 | +11.09 | 0.09751523 | +0.28156 | 0.227563 |
| 31.4 | 0.00947528 | 0.00949945 | +24.17 | 0.09813113 | +0.40086 | 0.207738 |
| 31.5 | 0.00950284 | 0.00955871 | +55.87 | 0.09899175 | +0.55778 | 0.238131 |
| 31.6 | 0.00956277 | 0.00955865 | -4.12 | 0.09936397 | +0.31076 | 0.222635 |
| 31.7 | 0.00961191 | 0.00948668 | -125.23 | 0.09991170 | +0.37163 | 0.219445 |

Five-sample arithmetic mean total OTR:
**0.098782756 g/s** active, **0.098398236 g/s**
control, a difference **+0.000384520 g/s = +0.38452 mg/s**
(~0.391% relative OTR increase, only +0.077
percentage points of the earlier 0.498 g/s frozen OUR).
The active OTR remains ~19.84% of frozen OUR.
Headspace-liquid O₂ flux remains ~0.0036 g/s
because the foam-area multiplier was held at zero.

**Bubble holdup does not improve monotonically.**
At 31.7 s the active trial has **125.23 mg less**
bubble gas than the zero-entrainment control,
despite all five PB positive NET birth diagnostics
remaining >0. Changes in hydrodynamics,
spatial deposition and transient transport can
outweigh localized bubble production. Do NOT
claim a persistent bubble inventory gain or
that 2 g/s of gas was actually added to bubbles.

**Population balance sanity:** last five
`bubbleGas sizeGroups-sum volume fraction,
min, max` values are `1 1 1`, and the
last mean Sauter diameter is ~0.00120455 m,
minimum ~0.00100252 m and maximum
~0.00533129 m; stable and finite. These
are *post-correction normalized* groups;
not independent PB number/mass moment closure.

**Conservation remains unvalidated.** Active
sampled endpoint bubble budget residuals
(t31.4–31.7) are 0.000429673,
0.00153456, 0.00139368 and
0.000140725 kg/s, versus zero-control
0.000268565, 0.00185957, 0.00127999
and 0.00170228 kg/s. The simple
**trapezoidal approximation** to
`injection - netHeadspaceDepositionLoss -
netBubbleLiquidMassLoss - bubbleBoundaryOut`
over 31.3–31.7 predicts bubble inventory
change of approximately **-0.09007 g**
in the active case, whereas the observed
inventory change is **+0.14112 g**,
a difference of **+0.23119 g**.
The corresponding zero-control coarse
trapezoid mismatch is **+0.38246 g**.
These are **under-sampled, endpoint-rate
approximations**, not a solver-discrete
conservation defect proven by this test.
They demonstrate that write-time fluxes
alone are NOT accurate enough to certify
conservation, even for the zero-entrainment
control. A time-integrated per-timestep
source/flux budget and gas O₂ species
inventory are required.

**Scientific decision:** the new model
runs and creates localized bubble PB
births, but the 0.5 s data show no
meaningful oxygen-delivery improvement
or lasting bubble holdup gain. Do not
increase `entrainmentEfficiency` or
`maxGlobalEntrainedKgS` to force
agreement with OUR. The next priority
is independent, time-resolved phase/gas-O₂
conservation and hydrodynamic/foam
validation around the upper impeller,
before considering any production
deployment or the common nonreacting
9.9 s checkpoint.

## Validation gates

1. **Zero-efficiency equivalence:** at identical starting
   state, the control using
   `headspaceEntrainingDeposition` at efficiency
   zero must reproduce the old `deposition` run's
   bubble mass, deposition, `bubbleGasAndLiquid`
   OTR and headspace-liquid OTR to numerical tolerance.
   If it does not, treat the integration as incorrect
   and do not enable the nonzero model.
2. **Model positivity/cap:** in the nonzero test check
   `HEADSPACE_ENTRAINMENT_CANDIDATE` for eligible cells,
   raw/limited rates, and `capScale`. Actual
   phase transfer is net and blended, so do not
   interpret raw candidate as deposited or dissolved
   oxygen. The cap is applied to the **raw gas mass
   flux**, not to O₂ or the final blended source.
3. **PB group birth:** new
   `HEADSPACE_ENTRAINMENT_PB_AUDIT`
   should show positive net bubble mass production
   where appropriate. Check that `f3.bubbleGas`
   (2 mm group) is allocated nontrivially, all
   size groups remain nonnegative, their sum
   remains ~1 where bubbleGas exists, and
   `d.bubbleGas` stays finite and physically plausible.
4. **Mass and O₂:** the net headspace gas mass sink
   must equal bubble gas mass gain from entrainment;
   gas O₂ is carried at `Y_O2.headspaceGas`,
   not created from the liquid or an arbitrary
   `0.23` field source. Track total O2 inventory,
   dissolved O2, gas fluxes and both gas-liquid
   transfer channels.
5. **Geometry, liquid volume, stability:** inspect
   phase area/inventory, free-surface cavity,
   near-impeller gas holdup, MRF behavior and
   outlet pressure boundary. Reject any strong
   unphysical gas-volume loss or positive-alpha
   violations.
6. **Longer / nonreacting case:** short reacting-state
   continuation is a smoke test, NOT a
   validated free-surface equilibrium. If promising,
   re-equilibrate proper nonreacting hydrodynamics
   and create a common 9.9 s checkpoint for
   MM, EMM and MEML. Preserve intentional
   four-hour biological update and frozen OUR
   schedule; it is not a solver bug.

References for research motivation (not calibration of
this empirical model):

- Wu et al., *AIChE Journal* 2023, "Gas phase hydrodynamics
  in a surface-aerated tank with a long-short blades agitator,"
  DOI 10.1002/aic.17680: source-aerated vortex pressure
  depression and downward transport, coupled Eulerian–PBM.
- *Chemical Engineering Research and Design* 2014,
  "CFD modeling of gas entrainment in stirred tank systems,"
  DOI 10.1016/j.cherd.2013.10.025: entrainment onset
  relates to free-surface interfacial hydrodynamics.
- OpenFOAM Foundation, "Population Balance Modelling,"
  https://openfoam.org/guides/population-balance-modelling/:
  group-volume conservation and PB size-class transport.
