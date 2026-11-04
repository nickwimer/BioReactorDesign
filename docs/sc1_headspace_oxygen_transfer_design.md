# SC1 MEML: restore headspace-to-liquid oxygen transfer without false gas-volume OTR

**Status (2026-11-03): design and trial implementation committed; NOT COMPILED or VALIDATED on Kestrel yet. This is NOT a production physics change.**

## 1. Confirmed defect in the old three-phase transfer matrix

The production MEML phase system is
`interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem` with
`bubbleGas`, `headspaceGas` (reference gas phase) and `liquid`.
The actual original MEML `constant/phaseProperties` has:

- `interfaceComposition.liquid`: Henry, `species (O2)`, on **(liquid and bubbleGas)** only.
- `diffusiveMassTransfer.liquid`: Higbie and spherical models on **bubbleGas/liquid** only.
- `interfaceComposition.headspaceGas`, `interfaceComposition.bubbleGas`,
  `diffusiveMassTransfer.headspaceGas` and
  `diffusiveMassTransfer.bubbleGas`: empty.
- `phaseTransfer`: **(bubbleGas in headspaceGas)** deposition, efficiency 1.

The only saved `interfaceCompositionPhaseChange:dmdtf.*` field was
`interfaceCompositionPhaseChange:dmdtf.bubbleGasAndLiquid`.
Thus there is no explicit headspace–liquid O2 dissolution path.
Deposition transfers gas-mixture mass and species *between the two gas phases*;
it does NOT replace an interphase oxygen-transfer model.

At t=9.9 s the original MEML upper `propeller` wall has roughly
74.4% headspace gas, 25.4% liquid and 0.2% bubbles in the upper
MRF height window (area-weighted), and 66.2% of a circumferential-normal
blade-face proxy is headspace-dominated. Experimental video was described
as showing strong foaming about this impeller, but has not been
quantitatively analyzed here.

The current MEML bubble/liquid O2 net transfer near t=30–31 s is
only ~0.084–0.092 g/s versus ~0.498 g/s biological OUR.
The missing pathway is a plausible but **not yet quantified** contributor
to this supply deficit.

## 2. Implement through the existing coupled OpenFOAM species transfer

Use the existing `InterfaceCompositionPhaseChangePhaseSystem` machinery,
which already couples a Henry-law composition model, mass-transfer coefficient
`K`, gas/liquid phase mass transfer, and equal-and-opposite O2
species source terms. Source-code path:

- `phaseSystems/PhaseSystems/InterfaceCompositionPhaseChangePhaseSystem/`
  `InterfaceCompositionPhaseChangePhaseSystem.C`
- `interfacialCompositionModels/interfaceCompositionModels/Henry/Henry.C`
- `interfacialCompositionModels/diffusiveMassTransferModels/`
  `diffusiveMassTransferModel/diffusiveMassTransferModel.C`.

For dissolved species on the liquid side the code evaluates

```text
R_O2 = rho_liquid * D_O2 * K_liquid * (Y*_O2 - Y_O2,liquid)
K_liquid = a_eff / delta
k_L = D_O2 / delta
```

where `a_eff` has units m^2/m^3, `K` units m^-2,
`delta` is a film thickness in metres,
and `Y*_O2` is the existing Henry-model equilibrium liquid
mass fraction evaluated from **headspaceGas** O2 concentration.
The phase system transfers the same O2 mass from headspace gas to liquid;
no independent source is added to the biological equation.

**Critical:** do **not** use the dispersed-bubble Higbie diameter
on a continuous headspace interface. Its area expression
`6 alpha_g/d_b` is not an appropriate general area for a headspace
cavity. The older EMM/MAX-OUR simulations experienced spurious
mass transfer in gas-dominated cells.

## 3. New closure: resolved interface first, optional localized froth

Implemented in
`OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceFilmModels/`.
The new run-time model is `resolvedHeadspaceFilm`, compiled
into a **separate scratch library** `libsc1HeadspaceFilmModels.so`.

Its default geometric interface-area density is

```text
a_geom = |grad(alpha.headspaceGas)|
```

with a smooth gate requiring both headspace and liquid
fractions to be present and little dispersed bubble gas. The
geometric area is divided by `(1-minInterfaceAlpha-fullInterfaceAlpha)`
(0.88 with the trial defaults), which normalizes the integrated
area of a planar 0-to-1 liquid/headspace transition to its
resolved surface area; this is a simple calibration and not
a substitute for grid-convergence tests.

**The resolved area alone is unlikely to supply the full oxygen
deficit**: a ~1 m2 macroscopic interface cannot be assumed to
provide the area of a dense foam. The optional froth treatment is:

```text
a_foam = C_foam * 6 * alpha_headspace * alpha_liquid / d_foam
a_eff  = min( max(a_geom, a_foam), a_max )
```

The extra area is used **only** inside a specified upper-impeller
height/radius window, and is zero unless
`foamAreaMultiplier > 0`. `max`, not `sum`, is used to
avoid simply counting the same diffuse interface twice.
The `minInterfaceAlpha/fullInterfaceAlpha` gate applies to
both estimates, and bubble-rich cells are excluded.
The parameters `foamAreaMultiplier`, `foamDiameter`,
`filmThickness` and `maxAreaDensity` require calibration
and mesh-sensitivity checks. None are experimental measurements.

A reference *pilot* configuration uses:

```foam
filmThickness 1e-5;          // m, exploratory, not calibrated
foamAreaMultiplier 0;       // first test: resolved interface ONLY
foamDiameter 0.003;         // m, used only if foam enabled
foamYMin 0.65; foamYMax 1.05;
foamRMax 0.30;
minInterfaceAlpha 0.02;
fullInterfaceAlpha 0.10;
maxBubbleAlpha 0.10;
maxAreaDensity 2000;        // 1/m
```

Since `D_O2` is supplied by the actual interface-composition model,
the trial's effective liquid-side `k_L=D_O2/filmThickness`.
Do not interpret the numerical `filmThickness` as a measured
physical film without validating mass-transfer data.

The library also optionally logs:

```text
HEADSPACE_O2_AREA_AUDIT t_s=... geometricArea_m2=...
  foamCandidateArea_m2=... effectiveArea_m2=... filmThickness_m=...
```

The two `resolvedHeadspaceFilm` orientations are added to
`diffusiveMassTransfer.liquid`; this preserves the normal BiRD
blending selection across phase topology while still using
one distinct `(liquid and headspaceGas)` Henry species model.

## 4. Paired case entries and consistency requirements

The experimental case adds:

1. `interfaceComposition.liquid`: `(liquid and headspaceGas)`
   with an O2-only Henry model **cloned from the existing bubble-liquid
   Henry entry** (same `$He_O2`, `$LeLiqMix`).
2. `diffusiveMassTransfer.liquid`: ordered
   `(headspaceGas in liquid)` and `(liquid in headspaceGas)`
   `resolvedHeadspaceFilm` models, both with the same closure
   but with only one area logger active.
3. `heatTransfer.liquid` and `heatTransfer.headspaceGas`:
   matching bidirectional RanzMarshall/spherical models. The existing
   two-resistance heat-transfer phase system *requires* a model
   on both sides of every new interface-composition pair.
4. `system/controlDict`: loads only the new isolated scratch
   `libsc1HeadspaceFilmModels.so` by absolute path.

The original bubble-liquid O2 and bubble-to-headspace deposition
entries are retained. No kinetic parameters or OUR schedule are
modified, and no independent headspace-to-bubble re-entrainment
is asserted to exist.

The scratch-only installer is
`scripts/diagnostics/prepare_sc1_headspace_o2_trial.py`.
It copies the latest collated time from the
existing **scratch** bubble audit case into a **different scratch**
run folder, patches only the new scratch `phaseProperties`
and `controlDict`, and refuses an existing destination.
It does **not** silently edit any production checkpoint.

The installer now accepts `--foam-multiplier X` (0 to 10) to create
different explicit foam-area scenarios from the **same original scratch
checkpoint**, with `--out` pointing at a new, previously nonexistent
scratch case directory. Start with `--foam-multiplier 0`, and only
consider a bounded nonzero sweep **after** the baseline source field,
gas/liquid inventory balance and effective area have been validated.
Do not use such a coefficient to fit DO without independent physical
evidence for the foam area.

### First Kestrel compilation check (2026-11-03)

The new scratch library's pure-Python preflight **passed**:
dictionary insertion, interface gates, area normalization, and bounded
foam-area checks. The first C++ `wmake` then failed with five
`invalid use of incomplete type 'const class Foam::phaseSystem'`
errors at `resolvedHeadspaceFilm.C` lines 100, 102, 144, 146 and
148. Root cause: `phaseModel.H` only forward-declares
`phaseSystem`, while this model accesses
`phaseModel::fluid().phases()`. The experimental source now includes
`phaseSystem.H` explicitly after `phaseModel.H` (commit
`bb390fb`). **Recompile the same isolated scratch library; there
is no physics/configuration change and no need to rerun the
Python checks unless the installer changes.**

### Second Kestrel linking check (2026-11-03)

After fixing the missing C++ `phaseSystem.H` include, the experimental
library's source **compiled successfully**, but the final `ld.bfd`
link failed with:

```text
cannot find -lbirdeulerianInterfacialCompositionModels
cannot find -lbirdphaseSystem
```

The libraries **exist** under `$FOAM_USER_LIBBIN`. The failed link
command had only `-L$FOAM_LIBBIN`, omitting the user library
directory. The branch now adds `-L$(FOAM_USER_LIBBIN)` to
`sc1HeadspaceFilmModels/Make/options` before the `-l`
dependencies. This is a linker search-path fix only, not a
model-physics change. Re-run `wmake` and check `ldd` for
the expected full paths; do not rebuild the existing BiRD libraries.

### Third Kestrel build check: standalone library successfully linked

On 2026-11-03 the user rebuilt the dedicated
`libsc1HeadspaceFilmModels.so` successfully after adding
`-L$(FOAM_USER_LIBBIN)`. `ldd` confirmed the new library
resolves both
`libbirdeulerianInterfacialCompositionModels.so` and
`libbirdphaseSystem.so` to
`$HOME/OpenFOAM9_RHEL9/user-9/platforms/linux64GccDPInt32Opt/lib/`.

The C++ **build and link are now verified on Kestrel**. There has
still been **no test of the case dictionary additions or the new
headspace-to-liquid source at runtime**. The next action is
preparing the isolated short `31.2-to-31.4 s` trial; validate
`controlDict`, both heat-transfer sides, and new phase pairs before
launching 52-rank MPI.

The trial preparer now returns nonzero exit status on a failed
preflight (as opposed to printing an error and returning success).
It refuses an existing destination; never delete scratch or production
checkpoint directories to work around that guard. Choose a fresh
`--out` in scratch if an earlier candidate already exists.

### Fourth Kestrel gate: trial preparer now accepts omitted writeControl

The user successfully compiled the updated scratch-only
`putidaFoamMask01BubbleAudit` executable with the
`netHeadspaceLiquidMassLoss_kg_s` diagnostic. The installed
headspace film library also compiled and linked successfully.

The first attempted trial preparation stopped **before creating the
destination run directory**, reporting
`Expected exactly one writeControl in controlDict`.
The existing scratch source case simply omits the optional
`writeControl` keyword. The installer now inserts
`writeControl runTime;` into the **new scratch case only** if
absent, replaces it when explicitly present, rejects duplicate
entries, and retains strict validation of required entries.
Regression tests cover all three cases. The production and
original scratch case `controlDict` remain unchanged.

Re-run the pure-Python preflight and trial preparer; do **not**
recompile either already successful shared library or diagnostic
solver. Do not launch CFD until verifying the new
`phaseProperties` and `controlDict` entries.

### Fifth Kestrel trial-preparation failure: nested function-object controls

The user ran the revised preparation. The Python model preflight
passed, but the installer stopped before creating the scratch trial,
reporting `Expected exactly one writeControl in controlDict, found 15`.
These matches are from the source controlDict's nested function-object
settings; the previous regular expression searched all levels.

**Fix:** the installer now lexes the OpenFOAM controlDict and filters
entries to the dictionary root, respecting nested curly/round brackets,
comments, and quoted strings. A missing root `writeControl`
is inserted in the **scratch copy** even if 15 function objects each
declare `writeControl`. Root `writeInterval`, `startFrom`
and `endTime` are updated without editing nested function-object
settings. The same scope-aware handling applies to `libs`, which
may also be declared separately inside multiple function objects.
Both root `libs`-present and -missing cases have regression tests,
including a fixture with 15 nested function objects.

Re-run only the **pure Python preflight and scratch case preparation**.
The headspace film shared library and scratch diagnostic solver have
already compiled successfully. No new CFD continuation has been run.
Do not delete existing directories or manually edit production cases
to work around an installer failure.

## 5. Kestrel execution: isolated build and short smoke test

**Never paste top-level `exit`, `logout`, `exec`, `set -e`,
or shell termination traps into an interactive `salloc` shell.**
The following only builds an experimental library in scratch and
prepares a new scratch run. First restore OpenFOAM and libfabric
if needed:

```bash
source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc"
module load libfabric
```

Then:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
BIRD="$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam"
SRC="$REPO/OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceFilmModels"
LIB=/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so

if git -C "$REPO" pull --ff-only; then
    if "$PY" "$REPO/scripts/diagnostics/test_sc1_headspace_o2_trial.py"; then
        mkdir -p "$(dirname "$LIB")"
        if (cd "$SRC" && SC1_BIRD_ROOT="$BIRD" \
            wmake > build.sc1Headspace.log 2>&1); then
            if [[ -f "$LIB" ]]; then
                echo "HEADSPACE FILM LIBRARY BUILD SUCCESS"
                ldd "$LIB" | grep -E 'bird|not found'
            else
                echo "Build reported success, but scratch library is missing."
            fi
        else
            echo "Build failed. Do not launch CFD."
            grep -nE 'error:|fatal error:|undefined reference' \
                "$SRC/build.sc1Headspace.log" | head -30
        fi
    fi
fi
```

**Only after the library exists and has no missing dependencies**,
prepare the isolated short experiment:

```bash
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
RUN=/scratch/nwimer/sc1_headspace_o2_trial_20261103/run
LIB=/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so

if [[ -f "$LIB" && ! -e "$RUN" ]]; then
    "$PY" "$REPO/scripts/diagnostics/prepare_sc1_headspace_o2_trial.py"

    if [[ -f "$RUN/constant/phaseProperties" ]]; then
        foamDictionary "$RUN/constant/phaseProperties" \
            -entry interfaceComposition.liquid -value
        foamDictionary "$RUN/constant/phaseProperties" \
            -entry diffusiveMassTransfer.liquid -value
        foamDictionary "$RUN/constant/phaseProperties" \
            -entry heatTransfer.headspaceGas -value
        foamDictionary "$RUN/system/controlDict" -entry libs -value
    fi
else
    echo "Missing library or trial directory already exists; no copy attempted."
fi
```

To capture the new channel alongside the old one, copy the
latest `bubbleMassAudit.H` from the GitHub branch into the
**scratch-only** diagnostic solver, then touch its
`putidaFoam.C` and rebuild only
`putidaFoamMask01BubbleAudit`. The old production solver is
not affected. The updated audit logs
`netBubbleLiquidMassLoss_kg_s` and
`netHeadspaceLiquidMassLoss_kg_s` separately.
For this present MEML configuration, both are O2-only
interfacial transfers; after changing any species models,
recheck that assumption.

Only after validating all of the above, run the **short**
trial continuation with the same 52-rank collated parallel IO
as the existing scratch diagnostic case:

```bash
RUN=/scratch/nwimer/sc1_headspace_o2_trial_20261103/run
EXE="$FOAM_USER_APPBIN/putidaFoamMask01BubbleAudit"
LIB=/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/libsc1HeadspaceFilmModels.so

if [[ -x "$EXE" && -f "$LIB" && -d "$RUN/processors52" ]]; then
    (cd "$RUN" && srun -n 52 "$EXE" -parallel -fileHandler collated \
        > log.headspaceO2.smoke 2>&1)
    echo "Solver command status: $?"
    grep -E 'HEADSPACE_O2_AREA_AUDIT|BUBBLE_MASS_AUDIT|FOAM FATAL|End$' \
        "$RUN/log.headspaceO2.smoke" | tail -30
fi
```

**Do not use this short existing-reacting-state trial as a new
baseline for comparative MM/EMM/MEML biology.** It checks that
the code compiles, that the new pair is assembled, and that
mass/species sources remain stable. Hydrodynamics changed: any
final baseline must be regenerated from nonreacting spin-up
or re-equilibrated sufficiently with the new physics and shown
to reach stable interfacial area, phase inventory and DO.

### Scratch trial successfully prepared; ready for runtime smoke test (2026-11-03)

The user pulled commit `f6655a7` and reran
`test_sc1_headspace_o2_trial.py`: dictionary, bounding and planar
area regression tests all **passed**. The new scratch case
`/scratch/nwimer/sc1_headspace_o2_trial_20261103/run` was
prepared successfully from collated
`processors52/31.199982372783058`, with configured end time
approximately 31.4 s and `writeControl runTime`,
`writeInterval 0.1`. `controlDict` loads the existing
`libalphaLimitVelocity.so` and the new scratch
`libsc1HeadspaceFilmModels.so` using its canonical `/kfs3/scratch`
path.

The `constant/phaseProperties` trial file has:
- added `(liquid and headspaceGas)` Henry O2 entry;
- ordered `(headspaceGas in liquid)` and
  `(liquid in headspaceGas)` mass-transfer models with
  `type resolvedHeadspaceFilm`;
- `foamAreaMultiplier 0` on both added transfer orientations;
- matching heat-transfer entries; and the original bubble-liquid
  Henry and bubble-to-headspace deposition entries retained.

The upgraded scratch `putidaFoamMask01BubbleAudit` has also compiled
successfully. **Runtime behavior has NOT yet been tested**. Next
is a short 52-rank collated `srun` with the experimental scratch case,
followed by checking new `headspaceGasAndLiquid` transfer output
against the existing `bubbleGasAndLiquid` transfer and solver
audit. No production CFD case or shared library was changed.

## 5a. First 52-rank runtime result: resolved interface alone insufficient (2026-11-03)

The first isolated trial **completed successfully (solver exit status 0)**,
from collated `31.199982372783058` with
`foamAreaMultiplier=0`. The checkpoint directory
`processors52/31.299982372786378` contains
`interfaceCompositionPhaseChange:dmdtf.headspaceGasAndLiquid`.

The `HEADSPACE_O2_AREA_AUDIT` at reported `t_s=31.3` was:

```text
geometricArea_m2=1.72563
foamCandidateArea_m2=0
effectiveArea_m2=1.72563
filmThickness_m=1e-05
```

The `BUBBLE_MASS_AUDIT` at the same reported time was:

| Endpoint quantity | kg/s | g/s |
|---|---:|---:|
| Bubble-liquid O2 dissolution (net loss from bubbleGas) | `9.35802e-05` | 0.0935802 |
| **Headspace-liquid O2 dissolution (net loss from headspaceGas)** | **`3.65347e-06`** | **0.00365347** |
| **Combined modeled O2 transfer** | **`9.723367e-05`** | **0.09723367** |
| Frozen biological OUR (earlier audit, approximate) | 0.000498 | 0.498 |

At this instant, modeled supply is **19.5%** of the biological
sink, leaving approximately **0.40077 g/s** shortfall.
The restored headspace route contributes about **3.8% of
combined OTR**, and only about **0.7% of OUR**. These
are instantaneous endpoint estimates, not a time-averaged
OTR or steady biological mass balance.

The new headspace rate divided by the resolved interface area
is `0.00365347 / 1.72563 = 0.002117 g/(m2 s)`.
If film thickness, local oxygen driving force and density were
unchanged everywhere, the **total** ~0.4044 g/s headspace
delivery needed to match OUR beyond bubble transfer (an
**additional** ~0.4008 g/s above the resolved-interface result)
would imply approximately **191 m2 total effective headspace
contact area** versus the current **1.73 m2**. This ratio is only a
**nonlinear-sensitivity scale**, NOT an endorsed foam area,
film thickness, fitted parameter or evidence that such
contact physically exists. Actual foam liquid films, dissolved
O2 and gas O2 partial pressure can differ spatially.

Other signals to investigate before calibration:

- Bubble inventory `0.00933447 kg`, sparger injection
  `0.00874277 kg/s`, deposition loss `0.00728772 kg/s`;
  assembled bubble phase-transfer source
  `-0.0073813 kg/s` equals deposition plus bubble-liquid
  mass loss to within printed precision.
- `orderedDepLoss_kg_s=0` is the historical diagnostic
  **ordered-key** lookup, not evidence deposition vanished;
  the unordered `depositionLoss` and assembled mass source
  are nonzero and agree.
- `headspaceBoundaryOut_kg_s=0.00774726` is a
  **signed NET** sum, despite the name; simultaneously
  `headspaceBoundaryIn_kg_s=-1.14072` is the
  *gross inward* mass rate. Gross outward must then be
  approximately `1.14847 kg/s`. Such large counterflow
  is not by itself proof of error, but requires patch-by-patch
  inspection via the existing `HEADSPACE_PRESSURE_AUDIT`
  log lines. This diagnostic includes physical
  non-coupled patches and must not be confused with net
  sparger throughput.
- `inventoryRate_kg_s=NA` and
  `endpointClosure_kg_s=NA` are expected at the
  **first audit write**; there is no second in-process
  write interval yet. It is premature to assert
  inventory conservation under the new closure.
- Only one saved transfer-field write was reported
  before `End`, despite nominal `endTime=31.4`.
  Check `Time =` / final `End` lines and use a longer
  isolated trial from the same source checkpoint to obtain
  at least two audit writes for closure testing.

### Next controlled tests

1. **Read-only patch diagnosis:** print existing
   `HEADSPACE_PRESSURE_AUDIT` lines from
   `log.headspaceO2.smoke`, especially patch identities
   and gross `headspaceIn_kg_s`/net
   `headspaceOut_kg_s`. No recompilation needed.
2. **Independent geometric/foam comparison:** keep the
   baseline completed trial intact. Use the already
   implemented `prepare_sc1_headspace_o2_trial.py`
   with `--source /scratch/nwimer/sc1_bubble_audit_20261029/run`,
   `--out` a **fresh** scratch path,
   `--foam-multiplier 1` and `--duration 0.5` for a
   strictly exploratory upper-impeller foam-area case.
   Existing `filmThickness`, `foamDiameter` and
   geometry bounds remain unchanged. Compare
   `effectiveArea_m2`, headspace O2 rate, bubble O2
   rate, gas/liquid inventory, and pressure-patch flows
   at common output times; check solver stability.
3. **Longer zero-foam control** from the same
   `31.199982...` source checkpoint and with
   `--duration 0.5`, so both trials have multiple
   writes for bubble inventory derivatives. A
   single 0.1 s endpoint is not a steady-state result.
4. If the full candidate foam area still cannot
   supply oxygen on physically plausible `k_L`
   and geometry, investigate the distinct missing
   `headspaceGas -> bubbleGas` re-entrainment
   process before forcing `foamAreaMultiplier`
   beyond defensible values. Validate with
   experimental hydrodynamics and OTR data.

**No production CFD/kinetic code or 9.9 s reference checkpoint
has been altered.** Runtime success validates that the
new mass-transfer pair can be assembled, not its physical
parameterization or eventual oxygen supply.

## 5b. Paired 0.5 s foam-area trials (2026-11-04)

**Data:** User's Kestrel `run_geometric_05` and `run_foam_05`
scratch logs, both from collated checkpoint `31.199982372783058`
with 52 ranks, both completed to ~31.7 s (exit status 0).
The subsequent SSH "No route to host" occurred after both runs
finished; user obtained a new `salloc` allocation.

**Arithmetic mean of the five diagnostics at 31.3–31.7 s,
not a continuous time integral or equilibrium measurement:**

| Quantity | Geometric (`C_foam=0`) | Foam (`C_foam=1`) |
|---|---:|---:|
| Resolved area | 1.72215 m² | 1.72276 m² |
| Candidate foam area | 0 | 4.25672 m² |
| **Effective headspace-liquid area** | **1.72215 m²** | **5.58391 m²** |
| BubbleGas→liquid O₂ | 0.094754 g/s | 0.094689 g/s |
| **HeadspaceGas→liquid O₂** | **0.003644 g/s** | **0.011627 g/s** |
| **Total modeled OTR** | **0.098398 g/s** | **0.106316 g/s** |
| Share of prior ~0.498 g/s frozen OUR | 19.76% | 21.35% |
| Deficit against 0.498 g/s | 0.399602 g/s | 0.391684 g/s |

The foam area makes headspace O₂ transfer **3.1905 times**
larger while the effective area grows **3.2424 times**. It
increases total OTR by only **0.007918 g/s**, ~1.59 percentage
points of OUR. The existing bubble transfer remains essentially
unchanged. The headspace transfer/area ratio is about
0.0021 g/(m² s) in both cases. This is evidence the experimental
closure responds to its area input, *not* validation of that area.

The candidate foam area `4.25672 m²` at
`foamDiameter=0.003 m` implies
`integral(gate*alphaH*alphaL dV)=4.25672*(0.003)/6
=0.00212836 m³` (with bubble-rich cells excluded).
The configured cylindrical foam window
`r<=0.30 m, 0.65<=y<=1.05 m` has volume
`pi*0.30²*0.40=0.113097 m³`.
Thus the window-averaged gate-weighted `alphaH*alphaL`
is only **0.01882**, indicating that the simulated phases
remain largely segregated within this candidate foam zone.
It does not represent extensive unresolved wet froth simply
because the impeller has a headspace-rich phase fraction.

A linear area extrapolation at unchanged mass-transfer
coefficient/driving force suggests roughly **190 m²** of
headspace contact would be required to close the shortfall
by headspace transfer alone. This is *not* a physical
target or justification to increase `foamAreaMultiplier`
toward 40–50: dissolved O₂, gas composition, geometry and
interfacial mass-transfer efficiency would change.

**Mass-budget and boundary qualifications:**

- At t=31.7, geometric deposition 0.00985781 kg/s plus
  bubble-to-liquid mass loss 9.5896e-05 kg/s agrees with
  assembled bubble interphase source -0.0099537 kg/s.
  Foam deposition 0.00931156 kg/s plus bubble-to-liquid
  9.5718e-05 kg/s agrees with -0.00940727 kg/s.
  This supports transfer-source assembly consistency,
  **not** a time-integrated conservation proof.
- The sampled endpoint bubble closure residual at t=31.7
  remains 0.00170228 kg/s geometric and
  0.000408192 kg/s foam. Mean absolute residuals over
  four later writes are ~0.001278 and ~0.000984 kg/s,
  respectively. Endpoint budgets are not solver-discrete
  conservation residuals; proper time integration is needed.
- The earlier 0.2 s trial's t=31.3
  `HEADSPACE_PRESSURE_AUDIT` measured open top `outlet`
  area 0.889278 m², net outward +0.00774726 kg/s,
  **gross inward -1.14072 kg/s**, hence gross outward
  ~1.14847 kg/s. Wall and propeller mass flows were
  nearly zero. This is large counterflow through the pressure
  boundary even though its NET mass flux is small.
  Examine inlet/outlet O₂ composition BCs and gas recirculation;
  this is not automatically proof of a numerical defect.

**Scientific decision:** A headspace-liquid film model is
necessary but insufficient by itself. Next assess **surface
gas re-entrainment into the dispersed bubbleGas phase** and
bubble breakup/coalescence near the upper impeller, in light
of the user's observed experimental foam. The present
`phaseTransfer` has bubbleGas→headspaceGas deposition
but no reverse headspaceGas→bubbleGas conversion. Any
entrainment transfer must be bounded, conservative,
experimentally constrained and benchmarked against gas
holdup, oxygen transfer and liquid inventory. Do not edit
production settings or force-fit the foam coefficient.

## 6. Acceptance criteria before any 9.9 s common checkpoint

1. The new `interfaceCompositionPhaseChange:dmdtf.headspaceGasAndLiquid`
   field is written, nonzero where headspace and liquid overlap,
   and near-zero in pure headspace and pure liquid.
2. The field's volume integral agrees with the diagnostic
   `netHeadspaceLiquidMassLoss_kg_s` (allowing for evaluation-time
   differences). Check signs, not only absolute magnitude.
3. The original `bubbleGasAndLiquid` field and
   `bubbleGas in headspaceGas` deposition persist, with expected
   signs. O2 is conserved across the new gas/liquid pair and
   headspace O2 does not become negative.
4. The total modeled O2 transfer equals
   bubble-to-liquid plus headspace-to-liquid contributions;
   compare to OUR (~0.498 g/s). Track **time averages**,
   not one endpoint.
5. With `foamAreaMultiplier=0`, inspect the
   resolved-headspace contact area vs expected physical
   free-surface area (~tank cross section + cavity/waves).
   Conduct grid-sensitivity checks on integrated area, not
   just concentration fields.
6. For a foam-area sweep, vary `foamAreaMultiplier` and
   `foamDiameter` only in scratch; constrain the assumed
   foam region using the video and quantitative OTR/kLa
   or DO-response measurements. If it takes unrealistic area
   to match OUR, revisit headspace entrainment, plume
   distribution, and cavity shape rather than forcing kL.
7. Maintain stable liquid volume (~1000 L), plausible
   bubble holdup and oxygen inventory, acceptable time-step
   changes and MRF interface structure.
8. A common validated **nonreacting 9.9 s** state must feed
   all intended kinetic variants. Preserve the deliberate
   four-hour biological step and frozen-OUR cadence.

The existing MEML model has no explicit
`headspaceGas -> bubbleGas` entrainment/fragmentation, which
remains a distinct second-stage physics requirement if the
experiment shows large gas recirculation or fine bubbles
generated from the surface. Do not treat interfacial O2
diffusion as a replacement for entrainment.

## References

- OpenFOAM 9
  [InterfaceCompositionPhaseChangePhaseSystem](https://cpp.openfoam.org/v9/InterfaceCompositionPhaseChangePhaseSystem_8C_source.html):
  shows how `K*D` enters species and phase mass transfer.
- [OpenFOAM Foundation: Disperse Multiphase Flows](https://openfoam.org/guides/disperse-multiphase-flows/):
  distinguishes resolved interfaces from the distributed
  surface area of dispersed phases.
- The exact repository `Higbie.C`, `Henry.C`,
  `BlendedInterfacialModel.C` and `PhaseTransferPhaseSystem.C`
  implementations in this branch are the source of the
  SC1-specific interfaces, sign conventions and blending behavior.
