# SC1: investigation and recommended route to a reliable simulation

**Date:** 2026-11-06. **Status:** source review, analytical checks and research;
no new Kestrel runs or experimental measurements. Source reviewed at
`d7fc94d8f2ba763ade4257093ec2b817c56ecc70`, with diagnostic fixes accompanying
this report. All experimental closures remain unapproved for production.

**Subsequent build update:** The isolated audit at `6ac7e36` successfully
compiled on Kestrel; all 19 local checks also passed there. No new CFD audit
run has occurred. See the [build record](sc1_timestep_conservation_audit.md#successful-kestrel-build-2026-11-06)
for the executable path, dependency evidence and remaining runtime checks.

## Recommendation

Retain the three-field MEML architecture provisionally: resolved coherent
headspace/liquid interface, dispersed bubbles with a PB, and conservative
exchange between them. Establish numerical conservation and validate the
hydrodynamics before developing a replacement disengagement/entrainment
closure. The most consequential new finding is that the present **net-only
gas-conversion source loses simultaneous opposing exchange of donor properties
and populations**. A second concern is the physical applicability of the
original deposition law to a macroscopic headspace interface.

There is not enough evidence to declare the vortex incorrect, choose a new
entrainment coefficient, or abandon the Eulerian framework. The best route is
a staged verification/validation program with explicit decisions below, not
another full-reactor coefficient sweep. OpenFOAM's multiphaseEulerFoam supports
both resolved interfaces and dispersed phase fractions [1]; a wholesale solver
migration is not required merely because SC1 contains both.

The authoritative prior results are in the [handoff](SC1_MEML_SESSION_HANDOFF_2026-11-06.md),
[entrainment audit](sc1_headspace_reentrainment_trial.md),
[film audit](sc1_headspace_oxygen_transfer_design.md), and
[historical hydrodynamic audit](sc1_three_phase_headspace_audit.md).
Cluster paths in those documents are not mounted in this workspace. Repository
templates are evidence of available code, not proof of current case settings.

## What the existing results establish

| Observation | Interpretation and limit |
|---|---|
| MEML liquid inventory ~1000 L; historical EMM ~1137 L near 30 s | Preserve MEML's inventory behavior. EMM's higher DO is not a matched-condition validation target. |
| Normal entrainment gives mean total OTR 0.098783 vs 0.098398 g/s | Only 0.385 mg/s improvement in five endpoint samples; no demonstrated sustained bubble-inventory improvement. |
| Historical frozen OUR ~0.498 g/s | Modeled transfer supplies ~20% of this demand. This does not imply the demand should be reduced or the biological schedule changed. |
| Gas feed 8.74277 g/s with air-like O2 fraction 0.23 | Approximate inlet O2 is 2.01084 g/s. Demand is ~24.8% of that feed; observed transfer is ~4.9%. Feed oxygen exists, but delivery efficiency is unresolved. |
| Upper impeller region is headspace-rich and liquid-containing | This does not establish dry blades. A coherent cavity and liquid-wetted froth need different interpretation. |
| 0.5 s trial at 17.38 rad/s | About 1.38 revolutions, less than the historical ~1.08 s inventory/flux residence proxy. These are short transients, not equilibrium comparisons. |

The inlet fraction calculation is a feed-capacity check, not a complete O2
balance or a prediction that the reactor can achieve the required transfer.
Offgas, liquid accumulation and all physical boundaries still matter.

## Findings ranked by consequence

### 1. Net gas conversion omits counterexchange — confirmed model limitation

`headspaceEntrainingDeposition.C` combines deposition and entrainment before
`PhaseTransferPhaseSystem::addDmdtYfs()` selects the donor by the sign of the
net source. Let E and D be **gross, post-blending** headspace-to-bubble and
bubble-to-headspace rates in the same cell, respectively. For the interpretation
of two simultaneous physical processes, the bubble O2 gain should be

```text
gross-channel O2 gain = E*Y_H - D*Y_B
current net-only gain = max(E-D,0)*Y_H + min(E-D,0)*Y_B
missing counterexchange = min(E,D)*(Y_H-Y_B)
```

For example, E=D=1 g/s, Y_H=0.23 and Y_B=0.20 gives 30 mg O2/s of
counterexchange despite zero net gas conversion. The current net model gives
zero. This is an analytical illustration, **not an estimate of SC1's error**;
actual gross blended rates and donor compositions have not been measured.
If donor compositions are equal, this O2 difference vanishes.

This mechanism can conserve total gas mass and total oxygen while representing
the wrong exchange. The same problem affects donor momentum, enthalpy and PB
turnover. Positive-net nucleation cannot replace old bubbles with new 2 mm
bubbles when simultaneous birth and removal cancel in mass.

**Recommended eventual design:** retain separate gross channels through a
single consistent transfer implementation: phase mass E-D; species EY_H-DY_B;
momentum EU_H-DU_B; corresponding donor enthalpies; PB birth E and a defined
size-resolved removal D. Apply the same blended, limited rates and evaluation
state everywhere. Do not add an opposite ordered phaseTransfer entry casually:
the existing signed blending can net the channels again. First measure gross
rates without changing the old physics, then test the redesign in closed cells.

### 2. Deposition is not a validated free-surface disengagement law

The original `interfacialModels/phaseTransferModels/deposition/deposition.H`
describes dispersed-phase deposition to a film on a third phase. Its law is
proportional to

```text
rho_B * alpha_B * alpha_surface * |U_B-U_surface| / d_surface
```

SC1 assigns headspaceGas as that surface. Thus the bubble loss depends on the
headspace diameter model, phase overlap and total slip. It is not explicitly
based on a connected free surface, normal bubble crossing or bubble rupture.
Agreement between the formula and assembled blending proves implementation
consistency; it does not validate this physical application.

**Next:** inspect actual d.headspaceGas and diameterModel, interface thickness,
local topology and slip. In a simple planar disengagement test, vary resolution
and interface thickness at fixed physical conditions. A free-surface gas-loss
prediction should not be arbitrarily controlled by a fictitious continuous-gas
diameter. If this test fails, replace the conversion closure before trying to
compensate with stronger entrainment. Do not disable physically real gas escape.

### 3. Phase and O2 residuals mix time levels — confirmed numerical structure

`phaseSystemSolve.C` constructs alpha dilatation sources from stored `divU()`
and advances alpha with MULES/subcycling. PB then solves. `fluid.correct()`
refreshes transfer. Species solve next; pressure later updates density and
`divU()`. Consequently, post-PIMPLE `dmdts()*deltaT` and stored phase fluxes
are not guaranteed to be the exact terms that generated the inventory change.

Species transfer additionally uses implicit sink and explicit donor terms.
For positive headspace-to-bubble q, the pair-source sum contains
`q*(Y_H_at_assembly-Y_H_after_solve)`. Film uptake similarly has explicit/
implicit lag. Small linear-solver residuals alone do not establish converged
phase coupling. Measure summed species-transfer actions and refine outer
coupling/timestep independently before attributing residuals to physics.

The [new audit](sc1_timestep_conservation_audit.md) is a useful first
diagnostic, explicitly labeled timestep quadrature. If residuals remain,
capture actual alpha `Su/Sp`, subcycle flux integrals, alpha correction
impulses, pressure-induced density changes, species temporal/stabilization
terms and PB source updates. Record stage and iteration identity; do not
accumulate every PIMPLE iterate as if each advanced physical time.

### 4. Depleted oxygen can expose an explicit-sink/clipping conflict

The repository `scripts/putida_stir/nonreact/constant/fvModels` template applies
`alpha_L*OUR/3600*(32/1000)*V` as an explicit oxygen sink. It has no local O2
availability factor. The actual Kestrel source may differ and must be read.
OpenFOAM 9's `basicSpecieMixture::normalise()` clips solved species below zero
and calculates the default species from the remainder [2]. It does not
proportionally rescale all active species at each correction.

If a sink drives O2 negative and clipping restores it to zero, the correction
adds oxygen numerically relative to the solved state. Distinguish requested
OUR, matrix-applied sink, actual inventory loss and clipping. This is a possible
explanation of an O2 residual, not yet an observed SC1 defect.

Keep the four-hour biological update policy intact. Once hydrodynamics are
validated, a physically oxygen-limited coupling may need availability-limited
realized uptake and an explicit unmet-demand diagnostic, consistent with
growth/products. That requires a deliberate biological-model decision; it must
not silently reduce the requested OUR to make CFD look successful.

For nonreacting validation, preventing future biology updates is insufficient
if a nonzero saved `our` still drives a coded source. Explicitly verify that
the new nonreacting test has no biological O2 sink or direct biology updates.
This test mode is separate from the preserved reacting update schedule.

### 5. PB normalization can hide removal, and its source is stage-dependent

The size-group equation transports alpha_B*f_i. Generic gas-phase deposition
is not automatically a class sink in the reviewed PB equation. Experimental
nucleation adds only positive net birth. In an isolated deposition-only,
constant-density cell with other PB sources absent, alpha_B decreases while
unnormalized alpha_B*f_i can remain constant. `velocityGroup::scale()` then
provides the apparent group-volume removal by dividing by the class sum.

Therefore a correction impulse need not prove a new leak; it can expose the
current removal approximation. But sum(f_i)=1 is not independent verification.
Inspect local driftModels, density-change and phase-change drift, external group
sources, `sourceUpdateInterval`, PB `nCorr` and `solveOnFinalIterOnly`.
The older PB reporter prints its first invocation per write, not necessarily
the final source. The new pre-solve rate is correctly only a proxy.

For the documented exact 2 mm pivot, current nucleation algebra allocates the
intended birth volume and number. Away from a pivot, its `eta_i/x_i` allocation
preserves volume but does not generally preserve target nucleus number. Do not
claim that hypothetical off-pivot limitation invalidates the present 2 mm test.
Number changes under breakup/coalescence are physical; verify their prescribed
changes and volume conservation rather than requiring constant number always.

### 6. Interface and film configuration need exact-name and area checks

Two consequential code paths must be checked against the real dictionaries:

- `phaseSystem::surfaceTension()` only applies a pair's resolved capillary
  contribution when its key exists in `interfaceCompression`/`cAlphas_`.
  A sigma entry alone does not activate that force. Key presence and compression
  strength are distinct: a zero coefficient can retain the key while removing
  sharpening. This is a configuration check, not proof SC1 omits capillarity.
- The composition system requests `diffusiveMassTransfer`. Model blending
  lookup checks that name (or its member), then `default`. A dictionary entry
  called only `massTransfer` does not select diffusive-transfer blending in
  this source. Audit actual `diffusiveMassTransfer.liquid`/generic/default
  selections, rather than assuming the historical `massTransfer` alias applies.

The film logger integrates submodel area before blending. Its reported
effectiveArea is not necessarily the assembled transfer area. Measure
`integral(K_assembled*filmThickness)` together with D, local equilibrium Y*,
driving force and OTR. The planar 0.88 gate normalization is analytically
reasonable under its stated two-phase planar assumptions, not a mesh or froth
validation. The bubble-rich cutoff alpha_B>=0.1 can omit a distinct headspace/
liquid surface; its double-counting rationale needs a topology test.
The fixed 10 micrometre film and optional foam-area law remain uncalibrated.

### 7. The open boundary requires a physical and numerical sensitivity test

Reported gross headspace inflow ~1.14072 kg/s at Y_O2~0.23 represents ~262.4 g
O2/s, while modeled liquid transfer is ~0.1 g/s. The 0.385 mg/s entrainment
signal is only ~1.5e-6 of that gross oxygen flux. Small percentages of gross
throughput or total tank inventory are inappropriate conservation tolerances.

Air-like inletOutlet composition is reasonable. Large bidirectional flow is
not automatically an error; it may represent real exchange or an artificial
nearby pressure boundary. Check outlet position, swirl, pressure reference,
backflow turbulence/phase data and absolute versus relative flux conventions.
Use physically justified headspace/plenum extent or vent representation as a
controlled sensitivity test. Do not impose one-way outflow simply to remove
an inconvenient diagnostic. Keep mesh-boundary uncertainty distinct from
species-transfer implementation error.

## Recommended sequence and decision gates

| Stage | Work | Required decision before advancing |
|---|---|---|
| 0: establish the actual system | Hash solver/source/libraries and matched starting fields; read full current phaseProperties, fvModels, fvConstraints, fvSolution, fvSchemes, MRFProperties, thermo/momentum properties and biological settings. Reuse existing logs/VTK. | Know which source/binary/settings actually ran. Resolve local/repository differences. |
| 1: small equation tests | Closed-cell one-way and opposing gas transfer, species donor contrasts, 2 mm birth, deposition-only PB, breakup/coalescence, planar film transfer; serial and small MPI comparisons. | Known conservation and source moments before correction; timestep/coupling convergence. Analytical Python tests alone do not validate OpenFOAM. |
| 2: reactor conservation | Build corrected isolated audit; fresh zero control first, then normal case from verified identical 31.199982372783058 state. Hold all physics coefficients. | Explain mass/O2 discrepancies and compare audited control with preserved control. Escalate to equation-stage capture where needed. |
| 3: nonreacting mechanics | Test liquid/headspace without sparging against measured surface/torque/flow; then prescribed sparging with a fixed measured bubble size; then PB. Include relevant surface aeration explicitly if it occurs experimentally. | Correct fill, rotation, free-surface morphology, power and gas distribution at each level of added complexity. Simplified cases are diagnostic, not the final checkpoint. |
| 4: physical gas conversion | If needed, replace deposition and net-only turnover with verified gross-channel, donor-consistent, topology-aware conversion. Use observed entrainment/bubble data. | Gas holdup, size distribution, residence and free-surface contact agree within measured/numerical uncertainty. |
| 5: oxygen-transfer validation | Nonreacting gassing-out/reoxygenation under matched medium, temperature, pressure and gas feed; measure DO response, gas compositions and sensor response. Validate K, area, Henry equilibrium and transfer jointly. | OTR/DO dynamics agree independently of biological demand. A single fitted bulk kLa cannot identify kL, area and hydrodynamics separately. |
| 6: checkpoint and biology | Demonstrate statistical stationarity and spatial/temporal convergence; checkpoint all required phase, species, PB, turbulence and flux state. Verify restart equivalence. | Release the common nonreacting state only when validated; retain intended four-hour biological updates for MM/EMM/MEML. |

The first new reactor pair repeats an already completed interval for a specific
reason: the old logs lack timestep budgets. It must use new scratch paths and
must not overwrite or silently continue completed runs. No coefficient sweep
is justified by the existing OTR comparison.

## Mechanical modelling choice

First verify blade shape/elevations, baffles, rotation sense, shaft treatment,
wall velocities and actual MRF cell-zone contents; an MRF envelope is not the
blade geometry. Check liquid versus headspace forcing and any phase-velocity
limiter. Compare torque-derived power P=torque*omega, pumping/mixing and radial
surface profiles, rather than judging a rendered surface alone.

MRF can remain the economical baseline for mean quantities. If the experiment
depends on time-dependent blade/free-surface interaction, use a matched
sliding-mesh comparison to test whether MRF misses the relevant dynamics.
Free-surface studies have compared both approaches to velocity measurements
[3], and pilot-scale VOF/LES/sliding-mesh work has validated surface motion [4].
Those results justify a targeted comparison, not a universal declaration that
sliding mesh or LES will fix SC1. Our present fixed-mesh audit does not support
a moving-mesh case: add geometric-conservation accounting before such tests.

If verified configuration and better conversion closures still cannot describe
the observed mixture of cavity and froth, evaluate a morphology-adaptive hybrid
resolved/dispersed formulation. HZDR supplies such a research solver [5], but
compatibility with this OF9/BiRD biology stack is unverified. This is a separate
port/validation project, not an immediate production upgrade. Conversely,
surface-aerated Eulerian/PB studies have obtained measured holdup, size and
transfer agreement [6]; the current architecture is not intrinsically doomed.

## Accuracy, duration and acceptance

Set tolerances against the quantity the study must resolve. For the tiny
0.385 mg/s claimed improvement, a proposed 10%-of-signal numerical error budget
is ~0.0385 mg/s, or ~0.0192 mg O2 over 0.5 s. This is an illustrative study
design criterion, not a universal OpenFOAM standard or a threshold already
achieved. If numerical/statistical uncertainty exceeds the signal, report the
entrainment effect as unresolved. Separate this stringent comparison from the
larger task of establishing ~0.1 g/s baseline transfer accurately.

Use physically relevant observables: liquid mass/volume, total and local gas
holdup, d32/distribution, actual source moments, torque/power, interface profile,
OTR and offgas. Refine timestep and outer/PB convergence separately. Use at
least three systematically related meshes where an observed-order/GCI claim is
made, and report when asymptotic convergence is not demonstrated. Conservation,
grid convergence and comparison to independent experiments answer different
questions [7]. Validate the diagnostics' own MPI reductions and restart behavior.

Do not equate a fixed timestamp with an equilibrated state. Demonstrate stable
means and fluctuations over multiple relevant rotation, gas residence and
mixing times, with sufficient sampling of slow dynamics. The desired 9.9 s
checkpoint is a release target, not proof of stationarity. If additional
nonreacting equilibration is required, retain its physical age/provenance and
explicitly validate any subsequent time remapping against biology scheduling
and old-time fields. Never silently relabel a reacting continuation as the
nonreacting baseline.

## Evidence still needed from Kestrel and the experiment

Highest-value case evidence: the current dictionaries listed in stage 0;
local source differences; startup library paths; raw completed control/normal
logs; exact initial checkpoint manifest; mesh quality and impeller/vent
geometry; source and constraints definitions, especially the actual O2 sink,
interfaceCompression, headspace diameter, diffusiveMassTransfer blending and PB
drift/source controls. None require rerunning old simulations to inspect.

Highest-value physical evidence: the actual video with scale/view/operating
conditions; actual fill, rpm, gas rate and reference conditions; impeller and
baffle dimensions; medium density, viscosity and surface tension/antifoam;
available torque/power, gas holdup, bubble size, surface profile and DO/offgas
response. Video can distinguish a connected cavity from froth but does not
uniquely determine local volume fraction or gas-liquid area.

## Changes made during this investigation

Fixed the new audit's use of the OpenFOAM Info stream API and its final partial
interval reporting predicate. Added/expanded local diagnostic tests and an
analytical transfer-accounting test. These are diagnostic corrections and
counterexamples, not changes to production physics. OpenFOAM compilation and
Kestrel runtime validation remain pending. The prior claim of local test
success was never a successful OpenFOAM build.

## Primary references

1. [OpenFOAM Foundation: Disperse Multiphase Flows](https://openfoam.org/guides/disperse-multiphase-flows/), resolved/dispersed framework capabilities.
2. [OpenFOAM 9 basicSpecieMixture.C](https://cpp.openfoam.org/v9/basicSpecieMixture_8C_source.html), especially normalise() versus startup correctMassFractions().
3. [de Lamotte et al., 2018, Chemical Engineering Science](https://discovery.ucl.ac.uk/id/eprint/10061590/), DOI 10.1016/j.ces.2018.07.024, MRF/sliding-mesh comparison in free-surface baffled stirred tanks. Evidence used here is the repository abstract.
4. [Jahoda et al., 2011, Canadian Journal of Chemical Engineering](https://onlinelibrary.wiley.com/doi/abs/10.1002/cjce.20477), pilot-scale free-surface motion. Evidence used here is the publisher abstract.
5. [HZDR Multiphase Addon, archived 2022 release](https://rodare.hzdr.de/record/1877), hybrid resolved/dispersed research formulation; not a compatibility endorsement.
6. [Wu et al., AIChE Journal 69, e17680](https://aiche.onlinelibrary.wiley.com/doi/abs/10.1002/aic.17680), experimental local holdup/size and mass-transfer comparisons for a different surface-aerated agitator. Its coefficients are not transferable calibration for SC1.
7. [NASA/NPARC verification assessment](https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html) and [grid convergence](https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html).

Implementation evidence uses the exact BiRD/experimental source in this branch:
`phaseSystem.C`, `phaseSystemSolve.C`, `phaseSystemTemplates.C`,
`PhaseTransferPhaseSystem.C`, `InterfaceCompositionPhaseChangePhaseSystem.C`,
`PopulationBalancePhaseSystem.C`, `populationBalanceModel.C`, `velocityGroup.C`,
`MultiComponentPhaseModel.C`, `putidaFoam.C`, `pU/pEqn.H`, the two experimental
libraries and the explicitly identified coded-source template. Source behavior
and unverified current-case configuration are distinguished throughout.
