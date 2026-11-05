# SC1 three-phase headspace pressure and bubble-mass audit

> **New headspace-to-liquid oxygen fix (2026-11-03):**
> The production MEML case has **no** headspaceGas-liquid Henry/film
> oxygen-transfer pair. A separate, optional `resolvedHeadspaceFilm`
> transfer model, non-destructive scratch-case installer, and an isolated
> validation procedure have been committed to this branch.
> See [SC1 headspace oxygen-transfer design](sc1_headspace_oxygen_transfer_design.md).
> These additions have **not yet been compiled or run on Kestrel** and
> do not change any production solver or case.


## Interactive Kestrel allocation safety

**Never paste a top-level `exit`, `logout`, `exec`, `set -e` or
`set -euo pipefail` into the interactive `salloc` shell.** A top-level
`exit` in the user command examples previously relinquished an allocation.
Use non-fatal `if` / `else` checks and print warnings instead. An executable
script invoked via `bash script.sh` is a separate process, but copy-paste
interactive commands must not contain exit-on-failure logic.

Kestrel canonicalizes `/scratch/nwimer` to `/kfs3/scratch/nwimer`.
To verify cloned files, canonicalize **both** the run directory and each
checked file before comparing prefixes:

```bash
RUN=/scratch/nwimer/sc1_bubble_audit_20261029/run
RUN_REAL=$(readlink -f "$RUN")
for f in system/controlDict system/fvSolution constant/phaseProperties \
         processors52/30/alpha.bubbleGas; do
    p=$(readlink -f "$RUN/$f")
    if [[ -f "$RUN/$f" && "$p" == "$RUN_REAL/"* ]]; then
        printf 'OK: %s\n' "$f"
    else
        printf 'CHECK: %s resolves to %s\n' "$f" "$p"
    fi
done
```

No `exit` commands or fatal shell settings appear in this interactive
preflight.

**Scope:** Diagnostic-only continuation of the MEML SC1 case using the existing
`putidaFoam` source with `PIMPLE/bubbleMassAudit true`. This is **not** a
change to deposition, oxygen uptake, sparger, phase-fraction transport or pressure
boundary conditions.

## 1. Compare the real case boundary conditions before continuing

The production Kestrel cases are not part of this repository. Compare their
**saved** 30-second patch definitions; do not assume that older templates
match the simulations. Pay particular attention to the top outlet/vent patch:

- `p_rgh`: `prghPressure`, `fixedFluxPressure`, `totalPressure`,
  `fixedValue` and their values, as applicable;
- `p`: physical static-pressure boundary (generally `calculated`
  with the `p_rgh` pressure solution);
- `U.headspaceGas`, `U.bubbleGas`: outflow and backflow behavior;
- `alpha.headspaceGas`, `alpha.bubbleGas`, `alpha.liquid`: backflow
  composition and whether a phase is artificially fixed at the outlet;
- `O2.headspaceGas`, `O2.bubbleGas`, `O2.liquid` as relevant;
- pressure reference in `system/fvSolution` and any `pMin`, `pMax`
  settings.

For OpenFOAM, `p_rgh = p - rho*gh`; a static-pressure outlet must be
checked in combination with local hydrostatic head, not in isolation.
Do not confuse an impeller-depressed *physical interface* with a fixed top
boundary. Also inspect the outlet patch identity in `constant/polyMesh/boundary`.

A minimally invasive way to inspect a reconstructed snapshot:

```bash
for f in p p_rgh U.headspaceGas U.bubbleGas \
         alpha.headspaceGas alpha.bubbleGas alpha.liquid; do
    [ -f "30/$f" ] || continue
    echo "=== $f"
    foamDictionary "30/$f" -entry boundaryField -value
done
```

For legacy EMM and MAX-OUR SC1, substitute `U.gas` and `alpha.gas`.
Take care that headspace and bubbleGas are distinct phases only in the new case.

## 1a. Observed MEML outlet velocity flux at t=25–30 s

OpenFOAM `surfaceFieldValue` on the `outlet` patch reported the following
**phase velocity** fluxes, with units m^3/s. These are **not** weighted by
phase fraction or gas density, and **must not** be interpreted as true phase
mass flow:

| t (s) | headspace phi out | headspace phi in | headspace phi net |
|---|---:|---:|---:|
| 25 | 1.00346 | 0.99617 | 0.00730 |
| 27.5 | 0.81044 | 0.80429 | 0.00615 |
| 30 | 0.92693 | 0.91964 | 0.00729 |

The corresponding `phi.bubbleGas` net at 30 s was +17.05 m^3/s but
`alpha.bubbleGas` at the top outlet is nearly zero. Therefore `phi.bubbleGas`
cannot be used as a physical bubble-flow balance. The audit below instead
uses `alphaRhoPhi` for the signed mass flow and separates inflow from outflow.

These results expose substantial simultaneous headspace inflow and outflow,
not necessarily an incorrect pressure boundary. The saved MEML, EMM, and
MAX-OUR SC1 outlet `p_rgh` conditions all use `prghPressure` with
`p=101325` Pa, while MEML `alpha.liquid` uses `zeroGradient` at the
outlet compared with `inletOutlet` (zero backflow liquid) in the historical
two-phase cases.

## 1b. Building the cloned solver with local includes

The local RHEL9 `putidaFoam/Make/options` references `../include`,
`../phaseSystems/lnInclude`, `../interfacialModels/lnInclude` and
`../interfacialCompositionModels/lnInclude`. The initial preparation
script copied only `putidaFoam`, causing
`fatal error: phaseSystem.H: No such file or directory`.
The preparation script now symlinks these existing source-dependency trees
into scratch; the already prepared case can be repaired with:

```bash
SRCROOT="$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam"
DEST=/scratch/nwimer/sc1_bubble_audit_20261029
for d in include phaseSystems interfacialModels \
         interfacialCompositionModels multiphaseCompressibleMomentumTransportModels; do
    if [[ -d "$SRCROOT/$d" && ! -e "$DEST/$d" && ! -L "$DEST/$d" ]]; then
        ln -s "$SRCROOT/$d" "$DEST/$d"
    fi
done
test -f "$DEST/phaseSystems/lnInclude/phaseSystem.H"
test -f "$DEST/include/createRDeltaTf.H"
(cd "$DEST/solver" && wmake > build.bubbleAudit.log 2>&1)
```

This links already built shared-library header trees only; it does not
modify the production case. Any subsequent C++ build errors should be
examined independently rather than ignored.

### Ordered deposition key fix (2026-10-29)

The SC1 deposition model is keyed by the ordered pair
`(bubbleGas in headspaceGas)`. The phase-pair equality operator
distinguishes ordered `in` from unordered `and` keys.
An earlier draft used the unordered
`phasePairKey("bubbleGas", "headspaceGas")` in the diagnostic;
that lookup would omit the deposition model.
The corrected diagnostic now uses
`phasePairKey("bubbleGas", "headspaceGas", true)`.

Existing scratch copies and binaries made before this fix require the
**copied** `bubbleMassAudit.H` to be updated and `wmake` rerun.
Do not modify the production solver or merge the draft PR before
checking the mass-budget output.

## 2. Build the diagnostic executable on Kestrel

This branch changes `putidaFoam.C` only to include
`bubbleMassAudit.H` **when explicitly enabled**, at an ordinary
write time, after the final PIMPLE iteration. Copy the working solver source
into a separate worktree (or build a dedicated executable); do not overwrite
the production `putidaFoamMask01` binary.

Suggested configuration in **only the cloned diagnostic case**:

```
PIMPLE
{
    // Preserve existing settings here.
    bubbleMassAudit true;
}
```

Start from the 30-second checkpoint and continue over a short, explicitly
limited time interval with the same hydrodynamic boundary conditions and
`fvModels`. Set a sufficiently frequent `writeInterval` so that the
finite-difference inventory derivative is meaningful.

**Do not infer that this diagnostic compiles against the separate RHEL9
`putidaFoamMask01` source without building and testing on Kestrel.**
The GitHub code is the baseline implementation, not proof of the final
RHEL9 binary's source history.

## 2w. First normal-inward headspace entrainment runtime (2026-11-04)

52-rank SCRATCH-only trial
`run_normalactive_02` completed normally
after updating the optional entrainment model
to use positive
`(U.liquid-U.headspaceGas) dot grad(alpha.headspaceGas)/|grad(alpha.headspaceGas)|`
instead of downflowing liquid as the
kinematic criterion. At t31.3–31.4
the model reported 0.163–0.169 m²
eligible area and 1,791–1,875 eligible
cells. Raw proposed mass was
2.571–2.665 g/s and the **global
unblended 2.0 g/s cap** was correctly
enforced (scale 0.778–0.750).

This did **not** result in 2 g/s of
new dispersed bubbles: positive
**net** bubbleGas mass gain passed to
the population balance was only
0.228–0.208 g/s, while the signed
overall bubble↔headspace net
exchange remained -7.14 to
-7.54 g/s (deposition still dominates).
Compared with the zero-entrainment
matched-time case, bubble inventory
increased by 11.09–24.17 mg and
total modeled OTR improved by just
0.28–0.40 **mg O₂/s**, not g/s.

Therefore the opt-in model's
sign, cap and bubble PB source run
successfully, but oxygen delivery
is not materially repaired on this
0.2 s test. The endpoint mass-balance
residual at t31.4
(0.000429673 kg/s) is comparable
to the positive net PB bubble birth
(0.000207738 kg/s); a rigorous
time-integrated conservative
species/phase budget is still needed.
Verify PB size-fraction min/max
and diameter before extending
the test. No production edits.

Full data table and next
acceptance gates are in
[SC1 headspace re-entrainment trial]
(sc1_headspace_reentrainment_trial.md).

## 2v. Experimental gas re-entrainment and PB size injection (2026-11-04)

The coupled headspaceGas-liquid O2 film-transfer prototype compiled
and was runtime tested successfully in previous scratch experiments.
Foam-area multiplier 1 raised headspace OTR from ~0.00364 to
~0.01163 g/s, but total OTR remained only ~0.1063 g/s against
~0.498 g/s frozen biological OUR. The top oxygen
`inletOutlet` species BCs have physically plausible reverse-flow
values (gas O2 mass fraction 0.23). The experiment reportedly shows
foamy mixing at the upper impeller.

After source review, **do not add a second opposite ordered
`phaseTransfer` entry alone**, because signed BiRD blending may
combine it incorrectly with deposition. **Do not add bubble phase
mass without populating the PB size-group fractions.**

A new **separate scratch-only** library now contains
`headspaceEntrainingDeposition`, a single combined
deposition + countervailing entrainment model with algebraically
identical original deposition at `entrainmentEfficiency=0`,
and `headspaceEntrainmentNucleation`, which assigns positive
net bubble gain (after deposition and phaseTransfer blending)
to the PB's existing size bins via `popBal.eta()`.
The mixture phase transfer carries donor headspace species
(including O2) through the standard BiRD specie-transfer
equations, without creating O2 mass from nothing.
Initial test: 0 or 0.01 entrainment efficiency, raw global
entrainment cap 0.002 kg/s, 2 mm target nucleus, downward-flow
gate and bounded upper-impeller zone. The existing
headspace O2 film model is retained with `foamAreaMultiplier=0`.
The optional new library and script have **not yet been
compiled or runtime-tested on Kestrel**, and no production
cases or libraries have changed.

See [scratch entrainment/PB development and stage gates]
(sc1_headspace_reentrainment_trial.md).

## 2u. Confirmed SC1 transfer matrix: no headspace-liquid oxygen pathway

The user read the **actual production MEML**
`constant/phaseProperties` phase-qualified transfer entries.
The active phase-system type is
`interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem`.

| Model entry | Actual configured value |
|---|---|
| `interfaceComposition.bubbleGas` | empty `()` |
| `interfaceComposition.headspaceGas` | empty `()` |
| `interfaceComposition.liquid` | `(liquid and bubbleGas)`: `Henry`, `species (O2)`, `k ($He_O2)`, `Le $LeLiqMix` |
| `diffusiveMassTransfer.bubbleGas` | empty `()` |
| `diffusiveMassTransfer.headspaceGas` | empty `()` |
| `diffusiveMassTransfer.liquid` | `(bubbleGas in liquid)`: `Higbie`, `Le $LeLiqMix`; `(liquid in bubbleGas)`: `spherical`, `Le 1` |
| `phaseTransfer` | `(bubbleGas in headspaceGas)`: `deposition`, efficiency 1 |

The only existing `interfaceCompositionPhaseChange:dmdtf.*`
checkpoint field is
`interfaceCompositionPhaseChange:dmdtf.bubbleGasAndLiquid`.
**No `(headspaceGas and liquid)` interfacial composition or
diffusive transfer model is configured.** The case therefore has
explicit O2 dissolution from `bubbleGas` to `liquid`, but
**not directly from `headspaceGas` to `liquid`**.

In `PhaseTransferPhaseSystem.C`,
`addDmdtYfs()` transfers species mass along with signed
*mixture* mass transfer (using each phase's species mass
fraction), so bubble-to-headspace deposition conveys the
gas-phase mixture, including its O2 species. However the
opposite transfer `headspaceGas -> bubbleGas` has no explicit
re-entrainment model in this configuration. No O2 diffusion
across the headspace-bubble phase boundary is separately configured.

**Correction of the earlier oxygen-audit caveat:** Although the
solver-level `fluid.dmdtf(bubbleGas,liquid)` can generally sum
multiple species, **this particular case specifies only O2**
in the interfacial composition model for this pair and no
additional `phaseTransfer` model on that pair. Therefore the
measured ~0.084–0.092 g/s net bubble-liquid transfer near t=30–31 s
is appropriately interpreted as **net modeled O2 interfacial
transfer for this specific configuration**, subject to the actual
signed field, and is much less than the ~0.498 g/s biological
oxygen uptake. The updated generic diagnostic label
`netBubbleLiquidMassLoss_kg_s` remains robust against future
addition of other species.

Physical relevance: at t=9.9 s the upper-impeller
`propeller` patch is ~74.4% `headspaceGas`, 25.4%
liquid and 0.2% bubbles (area-weighted). Experimental video
reportedly shows foamy gas-liquid contact at that impeller.
Omitting headspace-liquid O2 transfer and headspace-to-bubble
re-entrainment are **plausible missing processes** in a model
intended to resolve an aerated foam, not proof of a quantified
oxygen-transfer error.

**Do not simply add a headspace-liquid Henry/Higbie pair**:
the `headspaceGas` phase includes a large connected free
surface, and the available phase diameter/interfacial area
law may generate nonphysical mass exchange in gas-dominated
cells, as seen in older EMM/MAX calculations. Need a
physically constrained interface-area and re-entrainment
representation first. Prioritize confirming the experimental
foam topology and quantifying candidate gas-liquid contact
geometry without altering production physics.

## 2t. Verify actual O2 species-transfer pairs; phase-qualified dictionary keys

A query of the original MEML `constant/phaseProperties` for
unqualified `interfaceComposition` and
`diffusiveMassTransfer` returned missing-entry errors, while
`phaseTransfer` contained the ordered
`(bubbleGas in headspaceGas)` deposition model. **This does
not demonstrate absent gas-liquid oxygen transfer.**

BiRD's `phaseSystemTemplates.C` constructs species transfer
models from **phase-qualified dictionary keys** by calling
`IOobject::groupName(modelName, phase.name())`, e.g.
`interfaceComposition.liquid` and
`diffusiveMassTransfer.bubbleGas`, rather than requiring
an unqualified parent entry. The upstream two-phase templates
use the same convention, with `interfaceComposition.liquid`
specifying a Henry model for O2/CO2 and
`diffusiveMassTransfer.liquid` a gas-liquid film model.

Check the **actual** MEML `constant/phaseProperties`:
- `type`: is the interface-composition phase system active?
- `interfaceComposition.bubbleGas`,
  `interfaceComposition.headspaceGas`,
  `interfaceComposition.liquid`: which unordered phase pair(s)
  list `O2` and `CO2` in a Henry or other composition model?
- `diffusiveMassTransfer.bubbleGas`,
  `diffusiveMassTransfer.headspaceGas`,
  `diffusiveMassTransfer.liquid`: do the paired sides have valid
  diffusive transfer models?
- List saved
  `interfaceCompositionPhaseChange:dmdtf.*` fields: these
  name actually instantiated phase pairs, although their **net
  mass-transfer rates are not O2-specific**.

**Additional audit correction**: the `bubbleMassAudit.H`
field previously logged as `oxygenDissolutionLoss_kg_s`
is computed as
`-integral(fluid.dmdtf(phasePairKey("bubbleGas", "liquid")))`
and is therefore the **net mass exchange across that pair**,
which may combine O2, CO2, water and other interfacial
species. It cannot be compared directly to biological OUR
without measuring the O2-specific contribution. The source
diagnostic was relabeled
`netBubbleLiquidMassLoss_kg_s` for future scratch
audit builds; this is **diagnostic text only, no physics change**.
Older audit outputs retain the misleading historic label.

The bubbleGas-to-headspaceGas `deposition` model
transfers **mixture mass between two gas phases**,
not gas-liquid oxygen dissolution. If the production case
has no modeled `(headspaceGas and liquid)` species-transfer
pair, that channel would be absent even where headspace-rich
froth contacts liquid. However adding a headspace-liquid
Henry/film model indiscriminately risks excessive artificial
interfacial area and gas-transfer double counting. First
verify the actual dictionary pairs, then measure species-specific
O2 dissolution, then select a physically calibrated interface
area model consistent with the observed gas-liquid morphology.

## 2s. Physical validation caveat: experimental video shows froth at upper impeller

The user reports that the **upper impeller is liquid-submerged at
startup**, and that the **experimental video shows abundant foamy
bubbles around the upper impeller during operation**. The video itself
has **not yet been inspected** here, and there is no measured vortex
height or upper-blade liquid coverage from the experiment.

**Correct interpretation:** gas-dominated MEML wall cells near the
upper impeller (`alpha.headspaceGas~0.744`, `alpha.liquid~0.254`
in the y=0.708–0.937 m patch envelope) are **not, by themselves,
evidence that the blades are dry or improperly unsubmerged**.
A churned froth or gas-ingesting upper impeller can have a high
local gas volume fraction while remaining wetted by liquid films.
Published gas-dispersion literature distinguishes impeller
*flooding*, *loading*, and *complete dispersion* and discusses
gas-filled blade cavities; these phenomena are not synonymous
with a completely dry exposed impeller.

For this case the right model-validation question is not
"should the impeller stay fully submerged?", but:
**Does the experimental video show a persistent, connected,
relatively clear headspace cavity down to the upper blade region,
or an aerated, liquid-continuous froth with many distinct bubbles?**
The optical appearance can be ambiguous and video footage alone
will not quantify local gas fractions.

**New priority physics check:** verify whether the MEML case
provides interfacial oxygen mass transfer for both
`(bubbleGas and liquid)` and `(headspaceGas and liquid)`.
The previous oxygen audit measured the former, but the
phase fraction near the upper impeller is overwhelmingly
`headspaceGas`. If the `headspaceGas-liquid` species
transfer is absent, the model could undercount oxygen transfer
in a real **gas-rich froth**, even though the bubble-deposition
law is implemented correctly. This is a hypothesis, not an
established missing source; inspect the actual active
`constant/phaseProperties` and resulting `dmdtf` output
before suggesting new transfer physics. Do **not** naively
enable a headspace model without verifying interfacial
area scaling and preventing double counting.

Maintain the existing 9.9 s checkpoint and physics until
the actual impeller geometry, real image morphology, and
oxygen-transfer model configuration have been compared.

## 2r. Three-phase identification confirms HEADSPACE GAS at MEML upper impeller (9.9 s)

The scratch-only three-field `foamToVTK` conversion returned status
zero. `sc1_meml_propeller_gas_identity.py --convert` then evaluated
`alpha.liquid`, `alpha.headspaceGas`, and `alpha.bubbleGas` on
the 9.9 s `propeller` patch. The **maximum per-face phase-sum error
was 1.03e-6**. Area-weighted values:

| Patch region | Face area m2 | alpha liquid | alpha headspace | alpha bubble | Area with headspace fraction >0.5 |
|---|---:|---:|---:|---:|---:|
| Lower window, y=0.180–0.347 m, r<0.30 | 0.23824 | 0.870 | 0.000 | 0.130 | 0.0% |
| Upper window, y=0.708–0.937 m, r<0.30 | 0.38109 | 0.254 | **0.744** | 0.002 | **73.2%** |
| Upper y=0.708–0.780 m | 0.06407 | 0.744 | 0.248 | 0.008 | 17.1% |
| Upper y=0.780–0.850 m | 0.23565 | 0.198 | **0.801** | 0.001 | **79.6%** |
| Upper y=0.850–0.920 m | 0.07963 | 0.031 | **0.969** | ~0 | **98.8%** |
| Upper r=0.10–0.175 m, y=0.780–0.850 m | 0.09796 | 0.053 | **0.946** | 0.001 | **94.9%** |
| Upper r=0.175–0.30 m, y=0.850–0.920 m | 0.04898 | 0.050 | **0.950** | ~0 | **98.1%** |
| Upper abs(n dot e_theta)>=0.8 face-normal proxy, r=0.10–0.30 m | 0.14091 | 0.311 | **0.686** | 0.003 | **66.2%** |
| Lower equivalent normal-orientation proxy | 0.12960 | 0.883 | **0.000** | 0.117 | **0.0%** |

**Key conclusion:** MEML's nonliquid upper-propeller environment is
overwhelmingly continuous headspace phase, **not dispersed bubble gas**.
The lower propeller is predominantly liquid-dominated with some
dispersed bubbles and no headspace contact. The upper impeller's
gas-contacting state is present already at the 9.9 s nonreacting
checkpoint; 20 s liquid-only diagnostics showed similar exposure.

The circumferential-normal subset is an *approximate blade-face
orientation proxy*, not an authoritative classification of blade
faces. The `propeller` patch may combine shaft, hub and blades.
Neither a wall-cell phase fraction nor `alpha.headspaceGas>0.5`
proves physical dewetting, ventilation, or an incorrect real reactor
free surface. The EMM case is also not an independent physical
reference, because its liquid volume is ~151 L too high at 9.9 s.

**Go/no-go implication for biological baselines:** Until the
intended upper blade elevation and submerged operation are compared
with the actual mesh and expected free surface, treat the existing
9.9 s MEML state as *hydrodynamically unvalidated* for a
scientifically controlled MM/EMM/MEML comparison. Do not
compensate by changing deposition efficiency, biology timing or OUR.
If the design requires both impellers submerged and the CFD
upper blades are in headspace, investigate MRF/mesh/phase-transport
and interface dynamics, then re-equilibrate a common nonreacting
checkpoint under any verified physics corrections.

## 2q. Fine-scale 9.9-second upper-propeller liquid/gas exposure

The extended `sc1_upper_impeller_coverage.py --times 9.9`
postprocessing confirmed the following *patch surface* results:

| Geometric region | MEML area (m2) | MEML mean alphaL | MEML patch area alphaL<0.1 | MEML patch area alphaL<0.5 | EMM mean alphaL |
|---|---:|---:|---:|---:|---:|
| Lower MRF y=0.180–0.347 m, r<0.30 | 0.23824 | 0.870 | 0.0% | 2.0% | 0.901 |
| Upper MRF y=0.708–0.937 m, r<0.30 | 0.38109 | 0.254 | 64.5% | 73.3% | 0.942 |
| Upper y=0.708–0.780 m, r<0.30 | 0.06407 | 0.744 | 11.4% | 17.3% | 0.948 |
| Upper y=0.780–0.850 m, r<0.30 | 0.23565 | 0.198 | 69.5% | 79.6% | 0.941 |
| Upper y=0.850–0.920 m, r<0.30 | 0.07963 | 0.031 | 91.5% | 98.8% | 0.955 |
| Upper y=0.780–0.850 m, r=0.10–0.175 | 0.09796 | 0.053 | 90.0% | 94.9% | 0.948 |
| Upper y=0.850–0.920 m, r=0.175–0.30 | 0.04898 | 0.050 | 86.2% | 98.1% | 0.946 |

All reported EMM upper-window samples have **0%** of area below
`alpha.liquid=0.5`; its lower-window case is also nearly fully
liquid-contacting. In contrast, MEML's upper window is
headspace-/bubble-dominated near the likely upper blade-bearing
locations, even though its lower rotating region is mostly liquid
contacting. Thus the central cavity intersects the **solid
propeller-wall environment**, not merely the broader MRF envelope.

**Caveat:** This does not establish which faces are real blades,
shaft or hub, and the only available exported field in this run was
`alpha.liquid`. Nonliquid fractions can include both dispersed
`bubbleGas` and connected `headspaceGas`.

The new scratch-only
`scripts/diagnostics/sc1_meml_propeller_gas_identity.py` can
convert all three liquid/headspace/bubble alpha fields at 9.9 s
on the existing scratch-only MEML mirror, and report per-face
gas-phase identity in the exact regions above. It additionally
computes a *face-normal orientation proxy* to distinguish
likely radial-blade surfaces from hub/shaft surfaces, but this
is not authoritative geometric identification. Never conclude
the real impeller is unphysically exposed without matching
the actual production geometry and operating observations.

## 2p. Upper impeller-region liquid fraction: MEML vs EMM

The read-only `sc1_upper_impeller_coverage.py` test succeeded
for both original MEML and EMM at 9.9 and 20 s, using scratch VTK.

**Upper MRF window**: y=0.70802753–0.93662247 m; r=0–0.30 m.
This is a **geometric envelope**, not an exact list of cells
assigned to the active `rotor` zone, nor a component classification
of the `propeller` patch faces.

| Case, CFD time | Volume-weighted alpha.liquid (62 L window approx.) | Propeller patch area (m2) | Area-weighted propeller alpha.liquid |
|---|---:|---:|---:|
| Original MEML, 9.9 s | 0.459 | 0.38108 | **0.254** |
| Original MEML, 20 s | 0.465 | 0.38108 | **0.256** |
| EMM, 9.9 s | 0.945 | 0.39520 | **0.942** |
| EMM, 20 s | 0.940 | 0.39520 | **0.937** |

At 9.9 s, mean MEML liquid fractions inside the window
y=0.85–0.92 m are approximately **0** at radius 0.1–0.175 m,
and **0.024** at r=0.175–0.30 m; the EMM equivalents are
**0.960** and **0.947**. The MEML upper impeller-region fluid
is **predominantly not liquid**, while EMM is liquid-dominated
at equivalent coordinates. This condition is persistent to 20 s
and therefore predates the reacting continuation.

For propeller-patch faces in this upper window at 9.9 s,
MEML area-weighted alphaL ranges from **0.003** at
r=0.05–0.10 m to **0.439** at r=0.175–0.30 m, versus
**0.931–0.952** across those bands in EMM. The sum of all
reported patch areas is 0.38108 m2 in MEML and 0.39520 m2 in EMM;
hence the area-weighted averages above.

**Major interpretative caveats:**
1. The exported `propeller` patch combines shaft/hub/blade
   surfaces; its area-weighted mean is not literal blade wetted area.
   `1-alpha.liquid` includes both `headspaceGas` and
   `bubbleGas` in MEML, not solely headspace.
2. The meshes are not identical. The MEML propeller patch
   contains **9064 faces** and spans y=**0.218–1.741 m**;
   the EMM patch has **7928 faces** spanning
   y=**0.208–1.397 m**. These different extents
   warrant verification of actual shaft and blade geometry.
3. EMM contains ~151 L excess liquid already at 9.9 s;
   its upper-impeller liquid dominance is **not a sound
   independent validation** of the true reactor condition.

The previous question of whether the upper impeller is exposed
has progressed substantially: gas-dominated **fluid** and
**propeller wall-adjacent phase fractions** are directly measured
in MEML's upper rotating-region heights. Next, establish
which patch faces are actual blade surfaces and distinguish
connected headspace from dispersed bubbles there. The diagnostic
script now emits finer axial and radial propeller-patch liquid
coverage, plus the lower-impeller span as an internal control.
It uses only already converted VTKs and can be run at
`--times 9.9` without any CFD or VTK conversion.

This **does not** by itself prove that the free-surface vortex
is too deep; compare impeller blade geometry and experimental
surface observations before changing MRF masks or chemistry.

## 2o. Exact MEML MRF axial spans and propeller-coverage follow-up

A Kestrel inspection of the production `constant/polyMesh/cellZones`,
`constant/MRFProperties` and `system/topoSetDict.mrfFromYaml` showed:

- Original and RHEL9 MEML `cellZones` files are byte-identical;
  their active `rotor` zone has **35,629 cells**.
- Historical EMM's active `rotor` has **43,680 cells**.
  The cell counts alone are not directly comparable physical volumes.
- All three cases use `MRF1`, origin (0,0,0), Y axis,
  17.38 rad/s, and empty `nonRotatingPatches`. MEML ramps
  from zero to 17.38 rad/s over the first 0.1 s;
  EMM specifies a constant angular velocity.
- MEML creates the rotating `rotor` as the **intersection of
  the original `rotorEnvelope` and two axial windows**:
  lower y=**0.179705–0.347345 m**, upper
  y=**0.70802753–0.93662247 m**.

The MEML `alpha.liquid=0.5` contour around the shaft at 9.9 s
is near **y=0.693 m**, below the upper MRF window at the axis.
At r=0.05–0.10 m the contour is around 0.745 m; at
r=0.10–0.20 m it is near 0.789 m. Thus a significant part of
the *upper MRF zone* is headspace-dominated at the same early time.

**Do not equate the MRF zone with the physical impeller blades.**
For a more decisive result, the new read-only
`scripts/diagnostics/sc1_upper_impeller_coverage.py` measures
cell-volume-weighted `alpha.liquid` for radial bands and upper
MRF height sub-bands at 9.9/20 s, using existing scratch VTKs.
It additionally inspects the actual OpenFOAM `propeller` boundary
patch VTK and estimates the face-area-weighted liquid fraction at
upper-MRF heights if `alpha.liquid` was exported on the patch.
No new VTK conversion or CFD run is needed.

Interpret any patch-area-weighted liquid fractions as a
**wetted-area proxy**, not a rigorous interfacial wetting or gas
ingestion measurement. Compare the actual production geometry
before asserting that blades are dewatered.

## 2n. MRF rotation preflight across cases

Kestrel confirmed that all three simulations use an MRF `cellZone
rotor`, an origin `(0 0 0)`, an axis `(0 1 0)`, and the
same target angular velocity **17.38 rad/s (~166 rpm)**.

- Original MEML and the RHEL9 MEML continuation use
  `omega table ((0 0) (0.1 17.38) (10 17.38) (1000 17.38))`;
  they reach target speed by 0.1 s.
- Historical EMM uses constant `omega 17.38`, with no specified ramp.
- MEML includes `system/topoSetDict.mrfFromYaml`, which
  selects an active `rotor` cell zone from a broader
  `rotorEnvelope` with exclusions. The diagnostic has **not
  yet examined that file's full selections**, nor compared the
  resulting mesh `constant/polyMesh/cellZones` to EMM.

The impeller rotation **setpoint is not different after 0.1 s**,
but MRF cell-zone geometry, nonrotating patch declarations,
phase treatment, pressure field and spin-up initialization can
differ. Since the deep MEML vortex already exists at 9.9 s,
the important next diagnostic is whether its active MRF
cell zone overlaps the depressed core/headspace cavity and
whether the EMM zone is geometrically different.

Do not assume that the `rotor` name denotes an identical
physical volume in both meshes. Read the full MRFProperties
and MEML topoSetDict, and inspect the generated cellZones.
Changing the omega value now would not address the
underlying case-to-case geometry/phase differences.

## 2m. 9.9/20-second transient surface and liquid-volume results (2026-11-03)

The corrected `sc1_surface_evolution.py --analyze-only` completed
successfully using the scratch VTK exports. These are actual measured
quantities in reconstructed 9.9/20 s root fields.

| Metric | MEML original 9.9 s | MEML original 20 s | MEML RHEL9 30 s | EMM 9.9 s | EMM 20 s | EMM 30 s |
|---|---:|---:|---:|---:|---:|---:|
| Liquid volume (L) | 999.812 | 999.824 | 999.478 | 1151.077 | 1142.971 | 1137.119 |
| Liquid volume above y=1.02 m (L) | 82.304 | 81.614 | 82.191 | 175.650 | 169.998 | 164.149 |
| alpha.liquid=0.5 at r=0–0.05 m (m) | 0.693* | 0.729* | 0.674* | 1.205 | 1.194 | 1.193 |
| alpha.liquid=0.5 at r=0.45–0.53 m (m) | 1.194 | 1.193 | 1.198 | 1.292 | 1.255 | 1.242 |

Asterisk: multiple annularly averaged crossings were detected;
the highest downward crossing was used.

MEML contains an established depressed central free surface **already
at 9.9 s**, before the scheduled biology update at 10 s. Its liquid
volume changes just +0.012 L from 9.9 to 20 s, with no meaningful
change in elevated liquid inventory (82.304 to 81.614 L).
The outer annulus stays around 1.19 m, and the central annulus
varies by a few centimeters. Therefore the deep-vortex morphology
predates the biological continuation and is intrinsic to the
pre-9.9 hydrodynamic state/configuration; evidence does not yet
show its physical validity.

EMM has **151.077 L liquid above 1000 L already at 9.9 s**;
the excess declines to 142.971 L at 20 s and 137.119 L at 30 s.
EMM's integrated liquid inventory drops 13.958 L between 9.9
and 30 s. This is a significant inventory anomaly, but to
diagnose its *cause*, compare against t=0 inventory and
integrated liquid boundary fluxes (do not silently assume a
zero-flux baseline or attribute everything to solver mass error).

The 30 s MEML point derives from a **different RHEL9 continuation**
than the original 9.9/20 s time points. Identical
`constant/phaseProperties` files do not independently verify
identical mesh, rotation, phase transport controls or runtime
boundary configurations between these cases.

**Priority before changing the solver:** identify the actual
impeller/MRF rotation rate and geometry/baffle configuration;
compare predicted vortex depth (central alphaL crossing 0.67–0.73 m,
outer around 1.19 m) to reference operating measurements.
A new 9.9 s checkpoint under unchanged physics would simply
inherit a comparable hydrodynamic condition and is not a
remedy for the deep vortex. If hydrodynamics are changed, a new
nonreacting spin-up or a validated equilibration from warm start
is required before triggering chemistry.

## 2l. VTK selector fix: boundary-patch exports were mistakenly included

The first run of `sc1_surface_evolution.py` successfully converted
the original MEML and EMM 9.9 and 20 s snapshots (both `foamToVTK`
return codes were **0**), producing the four main volume files
`VTK/MEML_BASE_9.9.vtk`, `VTK/MEML_BASE_20.vtk`,
`VTK/EMM_9.9.vtk`, and `VTK/EMM_20.vtk`.

The script then mistakenly searched `VTK/**/*.vtk` and matched
the exported boundary patches as well as the volume mesh, so no
liquid-interface results were calculated. This was a **postprocessing
selection bug**, not a failed conversion or a CFD problem.

The fixed script selects `VTK/<case_name>_<time>.vtk` directly
(or exactly one root-level fallback) and supports
`--analyze-only` to reuse the existing VTK files without a new
conversion. Invoke in Kestrel:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
if git -C "$REPO" pull --ff-only; then
    /scratch/nwimer/.conda-envs/of-post/bin/python \
        "$REPO/scripts/diagnostics/sc1_surface_evolution.py" \
        --analyze-only
fi
```

## 2k. Availability and safe conversion of 9.9/20-second snapshots

Read-only Kestrel preflight confirmed:
- Original MEML `putida_sc1_final`: reconstructed root
  `alpha.liquid`, `alpha.headspaceGas` and `alpha.bubbleGas`
  at 9.9, 10 and 20 s; no 25.5 or 30 s root time.
- RHEL9 MEML `putida_sc1_final_RHEL9_mask01_to25p5`:
  25.5 s fields in collated `processors52`; 30 s includes
  reconstructed root `alpha.liquid` and collated gas fractions.
- EMM `putida_sc1`: reconstructed `alpha.liquid`
  at 9.9, 10, 20, 25.5 and 30 s.
- Original and RHEL9 MEML `constant/phaseProperties`
  are byte-identical, but that does **not** demonstrate unchanged
  solver binaries, model setup, or hydrodynamic history.

Use `scripts/diagnostics/sc1_surface_evolution.py` for
a read-only 9.9 and 20 s MEML_BASE-versus-EMM comparison.
It creates a real scratch directory with symlinks pointing to
production case input fields, converts only `alpha.liquid`
via serial `foamToVTK`, stores VTK/logs exclusively in
scratch, and reports total liquid inventory, above-1.02 m
liquid, and annular `alpha.liquid=0.5` crossings.
No production time directories or checkpoint fields are modified.

Do **not** combine 9.9/20 s original MEML and 25.5/30 s RHEL9
trends without explicitly recognizing the continuation boundary.
Later, if needed, read collated 25.5 s fields using a separate
scratch parallel conversion; do not launch `foamToVTK` with
`-case` pointing at production unless its output is redirected.

## 2i. EMM vs MEML using the same alpha.liquid=0.5 criterion

At t=30 s, use matched radial annuli and the highest
liquid-fraction crossing from above to compare the historical 2-phase
EMM and 3-phase MEML directly. The results are **annular averages**;
the old local wall climbing or deep interface filaments may be hidden
by azimuthal averaging.

| Radius (m) | EMM alphaL 0.5 height (m) | MEML alphaL 0.5 height (m) | EMM liquid above 1.02 m (L) | MEML liquid above 1.02 m (L) |
|---|---:|---:|---:|---:|
| 0–0.05 | 1.193 | 0.674 (multiple crossings) | 1.425 | 0 |
| 0.05–0.10 | 1.198 | 0.741 | 2.609 | 0 |
| 0.10–0.20 | 1.199 | 0.788 | 15.358 | 0 |
| 0.20–0.35 | 1.201 | 0.849 | 44.688 | 0.570 |
| 0.35–0.45 | 1.214 | 1.180 | 50.198 | 36.060 |
| 0.45–0.53 | 1.242 | 1.198 | 49.871 | 45.562 |

Integrated liquid inventory: EMM **1137.119 L** versus MEML
**999.478 L**. Liquid above nominal 1.02 m: **164.149 L** vs
**82.191 L**. EMM has approximately **137.64 L** extra liquid
relative to MEML, equivalent to about **0.158 m** of depth in a
uniform unobstructed vessel of radius 0.526 m (illustrative only;
actual mesh has internals).

The annularly averaged EMM interface is relatively uniform at
1.193–1.242 m (0.049 m range), while MEML's 0.5
liquid fraction rises from 0.674 m in the core to 1.198 m
near the wall (0.524 m range). EMM's near-wall interface
is only 0.044 m above MEML's. Therefore the old case's
excess inventory is not simply an exaggeration of the near-wall
rise; it also lacks the deep central cavity resolved in MEML.
Both simulations may have different physical-model errors:
old volume conservation / gas distribution, versus new central
vortex depth / bubble residence time.

**Next validation:** resolve azimuthal variation of
`alpha.liquid=0.5` crossings in wall-side annuli, because an
annular average can conceal the historical local penetration
and wall-climbing patterns. Do not declare the EMM interface
truly flat from annular means, nor diagnose the MEML cavity as
unphysical without physical operating-condition evidence.

## 2h. Annularly averaged MEML free-surface profile at 30 s

The height where annular mean `alpha.headspaceGas` crosses its
`0.3`, `0.5` and `0.7` thresholds was measured using height
bins of 0.025 m. The uppermost low-to-high crossing was chosen;
the method is an **annular average, not a reconstructed sharp interface**.

| Radial range from axis (m) | alphaH=0.3 height (m) | alphaH=0.5 (m) | alphaH=0.7 (m) |
|---|---:|---:|---:|
| 0.00–0.05 | 0.632 | 0.687 | 0.720 |
| 0.05–0.10 | 0.704 | 0.743 | 0.761 |
| 0.10–0.20 | 0.776 | 0.788 | 0.800 |
| 0.20–0.35 | 0.832 | 0.849 | 0.864 |
| 0.35–0.53 | 1.162 | 1.190 | 1.218 |

At the 0.5 contour, the interface is approximately **0.333 m below
the nominal 1.02 m level near the shaft** and approximately
**0.170 m above nominal near the outer wall**, giving an approximate
**0.503 m central-to-wall elevation difference**. The central
0.3-to-0.7 transition is ~0.088 m thick in the vertical coordinate;
the outer annulus transition is ~0.056 m.

**Physical interpretation is provisional.** The 3-phase mixture
contains distinct `bubbleGas` and `headspaceGas`, so an
`alpha.headspaceGas=0.5` contour is not necessarily the same as
`alpha.liquid=0.5` in the bubble-rich core. The old EMM and
MAX-OUR versions had known liquid-wall climbing and interface
penetration problems; the MEML outer-annulus elevation shows
that free-surface distortion might still be present under the new
formulation. The elevated liquid volume above the nominal height
(~82 L in MEML vs ~164 L in EMM) must be interpreted together with
total conserved liquid volume (~1.00 m3 vs 1.137 m3 respectively).
A surface rise does **not** by itself prove an overall mass gain,
since central depression can offset elevation elsewhere.

**Next test:** compare EMM vs MEML
`alpha.liquid=0.5` crossings in identical radial and height bins,
and volume-integrate `alpha.liquid` above 1.02 m per radial band.
This tests both models with a common liquid-phase definition,
rather than comparing EMM's combined `gas` with MEML's
separate `headspaceGas`.

## 2g. Radial localization of blended deposition at 30 s

The 30-second MEML VTK blended-deposition field was radially binned
using radius `sqrt(x^2 + z^2)` around the vessel axis.

| Radius (m) | Bubble volume (L), entire vessel | Deposition (L/s equiv.), entire vessel | Deposition (L/s equiv.), 0.6–0.8 m | Deposition (L/s equiv.), 0.8–1.0 m |
|---|---:|---:|---:|---:|
| 0–0.05 | 1.7865 | 3.7809 | 3.7108 | 0.0131 |
| 0.05–0.10 | 2.4843 | 1.3305 | 1.2511 | 0.0782 |
| 0.10–0.20 | 2.6703 | 1.3689 | 0.6073 | 0.6353 |
| 0.20–0.35 | 0.5176 | 1.1658 | 0.1528 | 0.9598 |
| >=0.35 | 0.1359 | 0.0783 | 0.0003 | 0.0476 |
| **Total** | **7.5946** | **7.7244** | **5.7223** | **1.7339** |

At 0.6–0.8 m, **64.85% of conversion** lies within
5 cm of the axis although that radial core contains **16.6% of
bubbles in this height band**. Within 10 cm, **86.71% of conversion**
occurs where **48.3% of this band’s bubbles** reside. At 0.8–1.0 m,
**92% of conversion** occurs at radii 0.10–0.35 m instead.
The apparent outward migration of the conversion surface with increasing
height is **consistent with an impeller-depressed, funnel-shaped
headspace cavity**, but is not sufficient to validate its geometry.

The highest priority follow-up is to calculate annular radial profiles
of the local headspace fraction crossing (0.3 onset, 0.5 midpoint,
0.7 full blending) and compare against a physically credible,
experimentally validated free-surface shape. Binned averages are
diagnostic, not exact `alpha=0.5` isosurfaces.
Do not alter deposition parameters or initiate new nonreacting spin-up
until the hydrodynamic geometry has been evaluated.

## 2f. Spatial distribution of the *blended* deposition at 30 s

The existing MEML 30 s VTK fields were reanalyzed with the exact
`phaseTransfer` default linear factor,
`clamp((alpha.headspaceGas - 0.3)/0.4, 0, 1)`.

| Height (m) | Bubble inventory (L) | Raw deposition (L/s equiv.) | Blended deposition (L/s equiv.) |
|---|---:|---:|---:|
| <0.2 | 1.788 | ~0 | ~0 |
| 0.2–0.4 | 2.333 | 2.727 | 0.178 |
| 0.4–0.6 | 2.352 | 2.207 | ~0 |
| 0.6–0.8 | 1.063 | 23.880 | 5.722 |
| 0.8–1.0 | 0.057 | 3.616 | 1.734 |
| >=1.0 | 0.002 | 0.104 | 0.091 |
| **Total** | **7.595** | **32.5339** | **7.7244** |

Approximately **76.38%** of blended deposition occurs below 0.8 m,
and **74.1%** occurs within the 0.6–0.8 m band. The strong deep
deposition therefore persists after correcting for blending.
Earlier connectivity tests showed the upper headspace connects to
a central depressed cavity reaching approximately 0.6 m; this
does not by itself show the vortex depth is physically correct.

All deposition numbers in this table are **density-normalized local
volumetric equivalents**; do not directly compare their absolute
magnitudes to the nominal inlet volume flow at reference conditions.
The actual applied mass deposition is approximately 7–12 g/s
near 30–31 s, comparable to 8.743 g/s inlet gas mass.

The most informative follow-up is to bin blended deposition radially,
especially the 0.6–0.8 m band near the shaft, and then validate
the free-surface cavity geometry against physical reference data
before considering changes to the hydrodynamics or redoing the
nonreacting 9.9-second checkpoint.

## 2e. Confirmed blending/source equivalence at 31.0–31.2 s (2026-11-03)

The user ran the diagnostic-only RHEL9 continuation successfully with a
`blendedCorrelationLoss_kg_s` measurement immediately following each
`fluid.correct()`. The configured `phaseTransfer` model inherits
the `default` linear blending with headspace continuous-phase thresholds
0.3 and 0.7. Its local weight for
`(bubbleGas in headspaceGas)` is

`f = clamp((alpha.headspaceGas - 0.3)/0.4, 0, 1)`.

| CFD time (s) | raw deposition (g/s) | blended formula (g/s) | phase-system source (g/s) |
|---|---:|---:|---:|
| 31.0 | 42.8467 | 11.9344 | 11.9344 |
| 31.1 | 46.2546 | 11.2587 | 11.2587 |
| 31.2 | 35.3585 | 7.86237 | 7.86237 |

The blended formula and stored phase-system transfer **agree to the
displayed precision** at both successive PIMPLE stages. The earlier
3–4x raw-vs-applied discrepancy is therefore fully explained by
interfacial blending, not an obvious source-code mismatch, timing
difference, stale library or incorrect deposition registration.

**Correction of an earlier interpretation:** Height-resolved estimates
of ~32 L/s raw deposition and ~89% below 0.8 m omitted the blending
factor. These are *not* actual applied transfer rates or spatial
deposition fractions. Recalculate a blended height-resolved deposition
before claiming that most applied deposition occurs below the intended
surface. Static pressure outlet is 101325 Pa, and the signed net
headspace exhaust is of the same order as the ~8.743 g/s injected gas
mass, while bubble and liquid net outlet mass fluxes are negligible.
The physical adequacy of free-surface depression, bubble inventory and
interfacial area remains unvalidated.

**Next validation:** recompute spatial distribution of blended
deposition using the existing 30-second VTK
`alpha.bubbleGas`, `alpha.headspaceGas`,
`U.bubbleGas`, `U.headspaceGas`, and cell volumes,
with the same linear weighting. Since density is unavailable in that
VTK export, report this as a **density-normalized local volumetric
equivalent**, not kg/s, and do not compare it directly with the inlet
reference volume flow without pressure/density conversion.

## 2d. Critical correction: blended interfacial deposition

The locally installed RHEL9 `deposition.C` and the public GitHub
`deposition.C` were compared and found **byte-for-byte identical**.
The audit binary resolves
`libbirdeulerianInterfacialModels.so` from the expected local
RHEL9 user library, built October 21, 2026, newer than the source
timestamp March 23, 2026. Timestamps do not prove binary identity.

**The more important discovery is that the raw `deposition::dmdtf()`
correlation is NOT the full applied phase-transfer model.**

- `PhaseTransferPhaseSystem.H` defines its model table as
  `BlendedInterfacialModel<phaseTransferModel>`.
- `PhaseTransferPhaseSystem.C` invokes
  `generatePairsAndSubModels("phaseTransfer", ...)`, and its
  `correct()` calls the BLENDED `dmdtf()`.
- `BlendedInterfacialModel.C` computes
  `f1 * model1In2_->dmdtf()` for the ordered
  `(bubbleGas in headspaceGas)` model (with corresponding other-model
  terms if configured). The local blending method is read from
  `phaseProperties/blending` and chosen by `phaseTransfer` or
  `default` if no specific one exists.
- The independent audit evaluated the **raw** deposition law, missing
  the spatial blending factor. The observed applied/raw ratios of
  0.301, 0.255 and 0.258 at 30.7, 30.8 and 30.9 s respectively
  are qualitatively consistent with this factor.

The actual case `phaseProperties` blending dictionary (reported by
Kestrel, 2026-11-03) is:

```
default {
    type linear;
    minPartlyContinuousAlpha.headspaceGas 0.3;
    minFullyContinuousAlpha.headspaceGas 0.7;
    // Identical 0.3/0.7 thresholds for liquid and bubbleGas.
}
heatTransfer {
    type linear;
    minPartlyContinuousAlpha.headspaceGas 0;
    minFullyContinuousAlpha.headspaceGas 1;
    // Identical 0/1 thresholds for other phases.
}
massTransfer { $heatTransfer; }
```

Since `phaseTransfer` has no dedicated `blending` entry, BiRD
selects `default`, **not** `massTransfer`. For the ordered
`(bubbleGas in headspaceGas)` model, the local deposition factor is
`clamp((alpha.headspaceGas - 0.3)/0.4, 0, 1)`. This can reduce
the integrated applied/raw ratio to around 0.25–0.30 in the diffuse
headspace interface. The raw-to-applied difference is **expected
from interfacial blending**, rather than evidence that the deposition
library was compiled incorrectly.

The `bubbleMassTransferStageAudit.H` test has now been extended to
print `blendedCorrelationLoss_kg_s` alongside the raw formula and
the phase-system `depositionLoss_kg_s` at the same PIMPLE stage.
For the SC1 thresholds those should agree to numerical precision
if the blending interpretation is complete. Only the separate scratch
audit executable needs recompilation.

To quantitatively confirm the explanation, read the actual SC1
`constant/phaseProperties` `blending` dictionary (read-only):

```bash
RUN=/scratch/nwimer/sc1_bubble_audit_20261029/run
foamDictionary "$RUN/constant/phaseProperties" -entry blending -value
```

Do not edit `efficiency`, rebuild libraries, or claim a deposition
implementation defect based solely on the raw/blended discrepancy.
The actual applied deposition is ~7–12 g/s near 30 s and is of the
same order as the ~8.74 g/s sparger mass feed. The unblended
30–40 g/s diagnostic estimate must not be interpreted as the
simulated mass loss.

## 2c. Same-PIMPLE-stage diagnostic at 30.7–30.9 s

The scratch audit executable was rebuilt with
`bubbleMassTransferStageAudit.H` and ran successfully to 30.9 s,
producing two nearly identical PIMPLE-stage measurements per output time.

| t (s) | phase-system deposition (kg/s) | same-state formula (kg/s) | formula / phase-system |
|---|---:|---:|---:|
| 30.7 | 0.0122922 | 0.0408206 | 3.32 |
| 30.8 | 0.00938401 | 0.0367572 | 3.92 |
| 30.9 | 0.0103014 | 0.0398872 | 3.87 |

Both PIMPLE stages agree closely at a given time, so the originally
hypothesized **post-pressure-correction evaluation-time mismatch
cannot explain the 3–4x discrepancy**.

The baseline GitHub `deposition.C` computes
`1.5*sign*efficiency_*droplet.rho()*droplet*surface/surface.d()
 *mag(droplet.U()-surface.U())`.
The stage-audit independent calculation is algebraically the same for
the SC1 configuration with efficiency 1. The next diagnostic should
therefore inspect the **actual** RHEL9 local
`interfacialModels/phaseTransferModels/deposition/deposition.C`
and the loaded `libbirdeulerianInterfacialModels.so`. A locally modified
or stale library is possible; do not assert that until source and
linker output are compared.

Use read-only `diff -u`, `sed`, `stat` and `ldd` commands.
**Do not rebuild, reconfigure, or rerun CFD until inspecting this.**
This is a diagnostic-only branch; no physics changes have been made.

## 2b. Results of revised 30.0–30.6 s audit

RHEL9 run completed with status 0. The configured deposition model
`(bubbleGas in headspaceGas)` was selected, but the **unordered**
`fluid.dmdtf(phasePairKey("bubbleGas", "headspaceGas"))` returned the
active deposited source while the **ordered** query returned zero.
The aggregate `fluid.dmdts()` bubble source agreed with the sum of
deposition and O2 dissolution. The diagnostic has therefore been
corrected to label the active unordered lookup as `depositionLoss_kg_s`
and retain `orderedDepLoss_kg_s` as a control.

Representative observed values (kg/s):

| Time | Bubble injection | Phase-system deposition loss | Bubble O2 dissolution loss | Correlation evaluated at write-time |
|---|---:|---:|---:|---:|
| 30.1 | 0.00874277 | 0.00764868 | 0.00008521 | 0.0342497 |
| 30.3 | 0.00874277 | 0.00937106 | 0.0000868895 | 0.0393078 |
| 30.5 | 0.00874277 | 0.00650396 | 0.0000890494 | 0.0280942 |
| 30.6 | 0.00874277 | 0.0110743 | 0.0000901002 | 0.0388874 |

The outlet carried essentially zero bubble mass and no liquid flux;
net headspace outflow was of the same order as sparger input. However,
gross headspace outflow and inflow through the open outlet each exceeded
1 kg/s, despite net headspace outflow of ~0.008 kg/s.

The stage-sensitive issue is that the independently recomputed
deposition correlation is approximately **3.5–4.5 times larger** than
the phase-system deposited source. This could reflect the correlation
being evaluated after the final pressure/velocity correction while
the saved phase transfer field was generated earlier in PIMPLE.
It could also reflect a local RHEL9 deposition model differing from
the checked-in source. It is not yet justified to claim a physical
conversion coefficient defect or perfect mass conservation.

A non-mutating `bubbleMassTransferStageAudit.H` is now optionally
included immediately after each `fluid.correct()` call. Its
`BUBBLE_STAGE_AUDIT` output compares both rates at the same
solver state. In an existing scratch-copy solver, copy this new header,
insert the guarded include after the two `fluid.correct()` calls,
and recompile the separate audit executable. No new non-reacting
or reacting baseline simulations are needed for this audit.

## 2a. First diagnostic runtime result and follow-up

A 52-rank RHEL9 diagnostic-only continuation from 30.0 to 30.3 s
completed successfully. Its first audit showed:
- Bubble mass inventory approximately 0.0093 kg.
- Total injected bubble mass 0.00874277 kg/s.
- O2 dissolution mass loss approximately 0.000086 kg/s.
- Net bubble outflow at the vent negligible (less than 2e-9 kg/s).
- Net headspace gas outflow 0.0080–0.0088 kg/s.
- Gross inward headspace mass flow approximately 1.07–1.12 kg/s,
  with comparably large gross outward flow.
- Outlet static pressure exactly 101325 Pa; liquid and bubble
  area fractions at the vent very small, liquid vent mass flux zero.
- The original **ordered pair diagnostic** reported a zero deposition
  loss, leaving an approximately 0.008–0.010 kg/s endpoint residual.

**Do not equate zero in this diagnostic with zero applied deposition.**
The phase-pair lookup, actual phase-system transfer-source assembly and
deposition correlation must be cross-checked before interpreting
conservation. In particular the instantaneous correlation estimate
(roughly 32 L/s density-normalized) is not a mass-flow budget.

The updated header now additionally emits:
- `unorderedDepLoss_kg_s`: the same transfer query with an unordered
  phase pair.
- `assembledBubbleDmdt_kg_s`: complete phase-system interphase mass
  source for bubbleGas from `fluid.dmdts()`. Negative values indicate
  net interphase mass removal from bubbleGas.
- `assembledBubbleDmdtPresent`: whether that field exists.
- `correlationDepLoss_kg_s`: independent density-based evaluation of
  SC1's deposition correlation using saved velocity and phase fields;
  **assumes efficiency 1**.
- `assembledEndpointClosure_kg_s`: inventory-rate residual using the
  summed phase-system source instead of the pair-specific O2 and
  deposition fields.

To apply the update, fast-forward the diagnostic Git checkout, copy
the updated audit header to the **already copied** scratch solver,
recompile the separate executable, and extend the **scratch** case
end time to 30.6 s. Never alter the production checkpoints or the
original RHEL9 executable; never paste top-level `exit` or shell-fatal
settings into an interactive Kestrel allocation.

## 3. Recorded quantities

`BUBBLE_MASS_AUDIT` contains:

- `inventory_kg` — instantaneous `sum(alpha.bubbleGas*rho.bubbleGas*V)`;
- `injection_kg_s` — recomputed total instantaneous `fvModels`
  density-equation source for bubbleGas, expected to be the 15 sparger
  sources if no other density source is configured;
- `depositionLoss_kg_s` — positive loss from bubbleGas based on the
  **signed phase-system `dmdtf`** for `(bubbleGas in headspaceGas)`;
- `oxygenDissolutionLoss_kg_s` — positive bubble-phase loss computed from
  the signed `(bubbleGas and liquid)` mass transfer;
- `bubbleBoundaryOut_kg_s` — sum of outward `alphaRhoPhi.bubbleGas`
  across non-coupled physical boundary faces;
- `headspaceBoundaryOut_kg_s` — corresponding headspace phase mass flux;
- `inventoryRate_kg_s` — difference in inventory between audit write
  times; and
- `endpointClosure_kg_s` — approximate residual between the inventory
  derivative and the contemporaneous source/flux terms.

`HEADSPACE_PRESSURE_AUDIT` logs *every* non-coupled patch:
pressure and `p_rgh` area means, their runtime boundary-condition types,
the headspace-gas velocity condition, and bubble/headspace mass outflows.

**Important limitations.** The audit samples endpoint source terms, not
time-averaged or time-integrated values. The closure is not an exact discrete
finite-volume residual, especially with moving meshes, subcycling, source
updates at PIMPLE iterations, nonconservative boundary discretizations,
species normalization and phase clipping. A large residual identifies a need
for a finer, equation-level budget; it does not by itself prove an OpenFOAM
mass-conservation bug. In particular, the snapshot-derived **32 L/s**
deposition is a density-normalized rate, whereas the audit uses kg/s.
These must not be compared as though they have the same units.

## 3a. Choosing a chemistry restart after model corrections

The diagnostic-only branch **does not change physics**, so it may continue
from the saved 30 s reacting state. Do not confuse this with a validation
restart following changes to pressure/outlet BCs, deposition efficiency,
bubble transport, sparger, initial headspace or mesh.

- For chemistry-only changes, start from a **validated** non-reacting
  9.9 s state with all required reacting fields initialized consistently.
- For any changed hydrodynamics, start from 0 s for a defensible new
  non-reacting spin-up; alternatively use a 9.9 s **non-reacting warm
  start** with the new physics and demonstrate a converged interface,
  liquid volume, bubble inventory, oxygen transfer, and outlet mass flow
  before turning on biology. A 0.1 s re-equilibration is not sufficient
  evidence.
- The original reactor biology model deliberately advances **four hours
  of biology per event** and holds OUR between CFD events; preserve that
  intended cadence when diagnosing hydrodynamics.

**Important code-specific scheduling hazard:** In `putidaFoam.C`, the
initial `prvs_react_update_time` is reconstructed from
`floor(startTime)/fluid_update_time`, and `bioReact.H` is entered
only when `time - prvs_react_update_time >= fluid_update_time`. The
`fluid_update_time` value is read at startup in `createFields.H`.
Therefore starting at 9.9 s with `fluid_update_time = 200` **does
not trigger biology at 10 s**; the first trigger would be 200 s.
To trigger once at 10 s and then freeze OUR, the previously staged
setup must restart *before* 10 s with the initial interval of 10 s,
save a checkpoint immediately after the first event (e.g. 10.1 s),
then restart with the long hold interval (e.g. 200 s). Do not run the
first stage long enough to trigger another biology event at 20 s.
Confirm the actual Kestrel `microbeUpdateProperties` and log rather
than relying on historical repository templates.

## 4. Safety and validation

- Retain the original 30-second checkpoint and all production case files.
- Build/run only in a cloned scratch case and under a separate executable
  name.
- Keep the 200-second frozen OUR update interval and all source definitions
  unchanged.
- First check no diagnostic output when `bubbleMassAudit` is absent.
- Verify MPI reductions by comparing bubble inventory at a checkpoint with
  the existing OpenFOAM `volFieldValue` result.
- Check the numerical budget at multiple write times, not from one instant.

The workflow's primary purpose is to determine whether the very large
bubble-to-headspace deposition predicted by the source law is consistent with
the actual gas-phase inventory and boundary flux; a separate pressure/BC
comparison tests whether a top-outlet condition is influencing that result.
