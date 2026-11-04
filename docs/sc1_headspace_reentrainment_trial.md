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
