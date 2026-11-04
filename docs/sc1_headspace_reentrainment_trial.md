# SC1: conservative headspace-to-bubble entrainment — isolated development trial

**Status (2026-11-04): source and regression scripts are committed, but the
new C++ library and its impact on the 52-rank case have NOT been compiled
or run on Kestrel. This is NOT a production change or a validated closure.**

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
