# SC1 three-phase headspace pressure and bubble-mass audit

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
