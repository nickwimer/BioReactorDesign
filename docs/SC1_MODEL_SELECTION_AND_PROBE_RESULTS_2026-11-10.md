# SC1 model selection reopened; source probe completed

## User direction and corrected objective

On 2026-11-10 the user clarified that the assistant introduced the three-field
representation as a proposed remedy for conservation problems, and that it
should be retained only if needed. **Three fields are not a project requirement.**
This supersedes the earlier provisional recommendation to retain that
architecture and the handoff's former three-phase wording of the ultimate goal.

Select the simplest representation that demonstrably preserves realistic
liquid inventory, conserves phase mass and oxygen, and predicts the relevant
hydrodynamics and transfer. Establish a physically validated, nonreacting
9.9 s checkpoint for common MM/EMM/MEML biological comparisons. Preserve the
intentional four-hour biological updates and the documented CFD trigger/hold
schedule. Biological-model labels do not independently prescribe phase count.

SC1 has two physical fluids, liquid and gas. `bubbleGas` and `headspaceGas`
are two numerical representations of the gas. Splitting them can distinguish
different velocities, compositions and interfacial areas, but also introduces
conversion closures. The number of fields alone is not a conservation fix.
Returning to two fields is not yet a demonstrated fix either: it can inherit
the same source/pressure numerics. Distinguish a structural model limitation
from an initialization, boundary-condition or discretization defect before
adding a field to address it.

## What the historical comparison does and does not establish

The [historical surface audit](sc1_three_phase_headspace_audit.md) records:

| Liquid volume | MEML, three fields | EMM, two fields |
|---|---:|---:|
| 9.9 s | 999.812 L | 1151.077 L |
| 20 s | 999.824 L | 1142.971 L |
| 30 s | 999.478 L | 1137.119 L |

The MEML 30 s point comes from a different RHEL9 continuation. The comparison
has not established matched t=0 inventories, integrated liquid boundary/source
budgets, generated MRF zones or all numerical settings. The older EMM excess
already existed at 9.9 s. Therefore it does not isolate phase count as the
cause, or prove that three fields were necessary to obtain approximately 1000 L.
Similarly, the newer volume alone did not prove discrete conservation.

The diagnostic below is a 0.1 s continuation near 31.2 s. Its sub-gram budget
discrepancy must not be extrapolated to explain the historical 151 L anomaly.
The old high-DO two-field state is not an acceptable matched starting point.

## Wall climbing is a separate acceptance test

The user reiterated on 2026-11-10 that the earlier liquid interface climbed
artificially up the tank sides. Correct total liquid inventory alone does not
resolve that defect: liquid can redistribute into an incorrect wall film or
vortex while conserving its integral. Conversely, excess inventory can raise
the entire surface. Diagnose both quantities before attributing the behavior
to phase count.

At 9.9 s, the historical audit measured the annularly averaged liquid 0.5
crossing near the wall at 1.292 m for EMM and 1.194 m for MEML, versus nominal
fill height 1.02 m. MEML's central crossing was 0.693 m, with multiple crossings.
These values neither validate the newer vortex nor establish disappearance of
localized wall climbing. Annular averaging can hide films and interface
penetration. No experimental surface-profile comparison has yet accepted either
state; do not dismiss the reported artifact as normal centrifugal rise.

Prioritize the following checks before another physical-closure change:

**Follow-up:** the [completed wall preflight](sc1_wall_preflight_run_2026-11-10.md)
finds zero `walls`-to-`rotor` adjacency in both actual meshes, ruling out the
direct MRF wall-assignment mechanism below for those meshes. Sector profiles
show persistent wall-side elevation and distinct core behavior. Its evidence
bundle has now been reviewed in the [configuration findings](SC1_WALL_BUNDLE_FINDINGS_2026-11-10.md):
the cases also differ in reference-phase rescaling, resolved capillary-force
activation, outlet treatment, turbulence and initialization. Do not repeat
the preflight or treat successful collection as a physical fix.

1. **Actual MRF zone and stationary walls.** Inspect generated cell-zone
   membership at tank walls and baffles, resolved `nonRotatingPatches`, and
   phase-velocity patch values. In [Foundation OF9's MRF implementation](https://cpp.openfoam.org/v9/MRFZone_8C_source.html),
   noncoupled, nonempty boundary faces adjacent to zone cells are included
   unless excluded; `correctBoundaryVelocity` assigns included faces
   `Omega cross (Cf - origin)`. BiRD's `MovingPhaseModel::correct()` calls
   this method. An accidentally included stationary wall is therefore a
   concrete mechanism to check. The repository stirred-tank tutorial has
   `nonRotatingPatches ()`; this is not itself erroneous when the zone does not
   touch stationary boundaries, and is not evidence of the Kestrel settings.
   Inspect the actual generated zone, not just the `rotor` name or dictionary
   intent. Also inspect baffle geometry and rotating-zone/interface overlap.
2. **Wall/interface transport and pressure balance.** Read actual alpha,
   velocity, pressure and outlet/backflow conditions, gravity, surface-tension
   pairs, interface-compression settings and near-wall mesh quality. Determine
   whether a contact-angle model is active before invoking wetting as the
   explanation. The repository tutorial's liquid wall uses `zeroGradient`,
   but archived SC1 settings remain unverified. Inspect near-wall velocity and
   pressure alongside interface migration. Do not prescribe a new contact
   angle, pin wall alpha, or change compression just to suppress the symptom.
3. **Spatial evidence from existing outputs.** Reuse scratch VTK exports to
   measure azimuthal sectors in wall-side bands, retaining all liquid-alpha
   crossings and reporting missing crossings and occupied-cell coverage.
   Compare 0.1/0.5/0.9 thresholds, high-wall liquid volume, and persistence at
   available times; a single highest crossing can mistake a droplet for the
   connected surface. Relate these to the central depression and total
   inventory. Thresholds diagnose diffuse layers, not a validated foam surface.
   Compare with experimental wall height, vortex depth and upper-impeller
   coverage at the same fill and operating conditions.

If these checks leave rotation modeling unresolved, a controlled free-surface
MRF/sliding-mesh comparison is a later diagnostic, not an immediate repeat of
the completed oxygen runs. [de Lamotte et al. (2018)](https://discovery.ucl.ac.uk/id/eprint/10061590/)
compare VOF MRF and sliding-mesh predictions for baffled stirred tanks with
experiments (abstract inspected). This supports that validation approach; it
does not prove MRF caused SC1's artifact or that VOF alone resolves dispersed
bubble area and oxygen transfer.

Any two-field candidate must pass both conservation and local surface-shape
checks. Retaining three fields requires evidence that their representation
addresses a physical requirement; the historical wall-climbing symptom alone
does not supply that evidence. No boundary, mesh, solver or case change is
made by this assessment.

## Completed probe: provenance and regression

`Pasted text(20261110-204412).txt` records completion of:

`/scratch/nwimer/sc1_source_probe_control_20261110T194214Z`

The canonical path starts `/kfs3/scratch/nwimer`. The case used the verified
private build `/scratch/nwimer/sc1_alpha_source_probe_20261110T191459Z` with
the hashes documented in the [build/run record](sc1_alpha_source_probe_run_2026-11-10.md).
It ran from 31.199982372783101 to 31.299982372786399 s, 2,000 physical steps,
two outer iterations, and six alpha subcycles. Solver status is 0, followed by
`End` and parallel finalization; ExecutionTime is 749.47 s.

The launcher reports exact printed agreement of every legacy SC1 and
bubble-mass record with the baseline. It verifies all 16 probe rows. Local
analysis of the supplied transcript independently sums the four bins for each
phase/outer index and matches pressure source, source splitting, implicit
evolution, and raw source against the existing same-outer alpha ledger.
No additional CFD is needed to establish this observational regression.

The unchanged baseline retains a liquid mass-budget residual of −181.224915 mg
over the interval, actual liquid mass increase +123.513742 mg, mean OTR
96.786624 mg/s, and O2 clipping addition 40.113551 mg. These are diagnostic
continuation results, not a nonreacting release checkpoint or a physical fix.

## New source/bounds evidence

All values below are density-weighted cumulative diagnostics for **outer index
2 only**, in mg. They are not sums over outer iterations or corrected solution
inventories. The raw target retains the actual old state and limited flux.

| Liquid source-assembly alpha bin | Split contribution | Implicit evolution | Raw-target minus actual alpha, mass weighted |
|---|---:|---:|---:|
| Lower, alpha < 1e-4 | +2.291036 | +0.000182 | −2.291219 |
| Upper, alpha > 1−1e-4 | −98.729433 | −0.532106 | +99.261539 |
| Interior | ~0 | −81.563215 | +81.563215 |
| Outside [0,1] | 0 | 0 | 0 |
| Sum | −96.438397 | −82.095139 | +178.533535 |

Raw liquid source is 302.047308 mg versus stored-pressure source 304.738688 mg;
their difference is −2.691380 mg. The bin result localizes the negative split
contribution to nearly liquid-filled cells and the implicit deficit mostly
to the interior bin. It does not identify spatial coordinates.

For bubbleGas, the lower bin supplies the full +79.777346 mg split term.
Total implicit evolution is +0.150756 mg; raw-target minus actual alpha sums
to −79.928102 mg. Thus the floor changes bubble removal materially in this
diagnostic; it is not evidence that more entrainment is required.

Each phase/outer table covers 3,379,968,000 repeated cell/subcycle samples:

| Outer index / phase | Raw target above bound: samples | Raw bound excess (mg) | Stored-pressure target above bound: samples | Pressure-target bound excess (mg) |
|---|---:|---:|---:|---:|
| 1 / liquid | 101 | 0.000108112 | 2354 | 473.466958 |
| 2 / liquid | 98 | 0.000110285 | 98 | 0.000110299 |
| 1 / bubbleGas | 0 | ~2.43e-14 | 0 | ~2.43e-14 |
| 2 / bubbleGas | 0 | ~1.43e-12 | 0 | ~1.43e-12 |

No below-bound samples are classified at the 1e-12 tolerance. The excess
metric includes even sub-tolerance excursions; counts are not unique cells.
For final-outer liquid, signed identity remainder is 2.28e-9 mg and the L1
remainder is 0.000338024 mg. Bin sums match the independent alpha terms within
printed rounding. These remainders are tiny relative to the diagnosed terms.

**Interpretation:** the final-outer raw target's individual-phase bound excess
is much smaller than the source deficit; boundedness alone has not demonstrated
that the whole deficit is unavoidable. However, joint alpha sums, reference-gas
bounds, pressure/density and species/PBM consistency were not tested by this
counterfactual. The large first-outer pressure-target excess also shows why
lagged `divU` cannot simply replace the source. No addback or floor change is
accepted, and the probe provides no evidence that three fields are necessary.

## Four two-field baseline tests completed

The [isolated two-field test record](SC1_TWO_FIELD_CONSERVATION_TEST_2026-11-10.md)
now contains four completed 100-step passing screens: closed transport and a
known liquid-feed/outlet case, each without a reference phase and with gas as
reference. The closed pair completed at `56a84d7`; the source pair completed
at `ef74cf6`, in `/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T233431Z`.
The user supplied full reports and solver-log hashes. Preserve these cases;
no repeat of the baseline suite or new build is needed.

For the source pair, each 0.100 g feed is balanced by approximately 0.090419 g
retained and 0.009581 g exported liquid. Cumulative absolute liquid residuals
are `6.991205047567109e-11` and `3.8064394580038755e-11` kg, far below the
`1.001e-07` kg screen. Both phases pass, and neither reference setting is
selected by these tests. This supports continuing with the simpler candidate;
it does not prove a third field is necessary or establish full-reactor validity.

The known-source fixture starts at equal phase fractions, has constant
densities, and shows essentially zero source-split contribution. It does not
reproduce the earlier probe's near-pure-liquid regime. The [separate matched
near-bound source test](SC1_TWO_FIELD_NEAR_BOUND_TEST_2026-11-10.md) is now
implemented, with explicit inventory/tolerance contracts and matched
dt=0.001 / 0.0005 s pairs at alphaL=0.99999. It keeps the same source/outlet
balance and original solver floors, measures both budgets and accepted-step
alpha regimes, and requires an observed source-split effect. All 53 focused
tests passed locally; CFD is mocked in local runner tests. The coarse pair
then completed at `78ebbb8` in
`/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z`. Both diagnostics pass,
but both phase mass screens fail. Liquid shortfalls are about 0.900 micrograms
in either variant (8.91 times the threshold), closely matching the measured
source-split contribution. Gas-reference rescaling is zero, yet its separate
gas and liquid budgets still fail. No reference setting is selected as a
remedy. The half-timestep pair and read-only comparison completed at `75dc652`,
with both diagnostics valid and both phase screens still failing. Fine/coarse
liquid L1 ratios are 0.9998454336 / 1.0000548393; the cumulative deficit is
essentially unchanged. All four near-bound cases are complete and must be
preserved. Source algebra predicts about 9e-10 kg lost from the coupled
dilatation component at these nominal fractions, matching the measured
deficit. The [guarded private rate-split candidate](SC1_TWO_FIELD_RATE_SPLIT_TEST_2026-11-10.md)
has now built and passed all eight pressure/flux-coupled CFD regressions at
`d83bc83`, against unchanged thresholds. All original controls remain saved.
The new [MEML-mesh quiescent startup](SC1_REACTOR_STARTUP_2026-11-10.md) reuses
that verified pure-phase binary, with a measured 1000-L fill, gravity and a vent.
It completed at `e586c47`: mass/diagnostics pass, but quiescence fails, with
liquid speed 0.0963573 m/s and Co 0.623991. Source review found that its blending
nearly cancels the listed segregated drag. The comparison restoring only its
weight [diverged by 0.0009 s](SC1_REACTOR_DRAG_FAILURE_2026-11-11.md). Saved-log
collection and bundle inspection are complete; the
[face-momentum startup](SC1_REACTOR_FACE_MOMENTUM_2026-11-11.md) also failed.
Its bundle is now inspected; the [coupled-drag comparison](SC1_REACTOR_COUPLED_DRAG_2026-11-11.md)
also failed at the same first-step outer iteration. Its bundle is now verified. The
[short pressure probe](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md) has now reproduced
the pressure collapse and stopped before limiting; its bundle is now verified.
The [private pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
stopped at 0.2 ms on the liquid-speed guard; its bundle is verified.
The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md) completed; its verified bundle also reveals horizontal gas motion near the interface. The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped on the liquid-speed guard at 8.4 ms; its bundle is verified. The completed force probe identifies a gravity/pressure residual at the sampled peaks; next test the [fixed hydrostatic-reference candidate](SC1_REACTOR_VOF_HYDROSTATIC_2026-11-11.md). Physical
adequacy remains to be established.

## Next decision: continue validating the simpler candidate

The [external benchmark review](SC1_EXTERNAL_BENCHMARKS_2026-11-10.md) identifies
two-field BiRD stirred-tank implementations with oxygen transfer and PBM.
The [Foundation framework](https://openfoam.org/guides/disperse-multiphase-flows/)
supports resolved-interface and dispersed interpretations of phase fraction.
These establish a candidate route, not SC1 validation.

1. **Read-only baseline reconstruction first.** Review the archived two-field
   case's t=0 initialization, earliest available liquid inventory, physical
   boundary fluxes, sources, phase/pressure settings and generated MRF geometry.
   Include the stationary-wall and azimuthally resolved surface checks above.
   Compare with the three-field setup without conflating biological-model and
   hydrodynamic differences. Existing archive paths are recorded below.
2. **Prepare a fresh, matched nonreacting two-field candidate.** Preserve the
   experimental geometry, fill, operating conditions and free-surface domain
   wherever its interaction with the upper impeller is relevant. Define bubble
   versus resolved-interface area/drag/PBM applicability explicitly; do not
   treat headspace gas as a cloud of small bubbles. The conversion from a
   three-field checkpoint is not a trivial field rename or unweighted sum of
   compositions. Prefer a controlled common initialization unless a conservative
   mapping is separately specified and verified.
3. **Use a phase-count-independent conservation audit.** Current SC1 readers
   and instrumented stages contain hard-coded `bubbleGas`/`headspaceGas` names.
   They need an isolated generalization or a separately validated two-field
   ledger before applying the same acceptance tests. Do not launch the current
   three-field executable on renamed two-field inputs.
4. **Select on measured adequacy.** Require liquid/gas and O2 budgets, bounds,
   convergence, credible free-surface/mixing behavior, gas holdup and transfer.
   Retain a separate bubble/headspace field only if the two-field model fails
   a demonstrated physical requirement because of its representation, rather
   than one defective setup, and the extra field resolves it without
   unacceptable conservation errors. A degassing-only domain can be a reduced
   control; it is not automatically sufficient for an interface-crossing impeller.

Archived cases, read only:

- EMM: `/projects/microenviro/final_runs/reacting_Julia_08282026/emm/putida_sc1`
- MEML: `/projects/microenviro/final_runs/reacting_Julia_08282026/meml/putida_sc1_final`

These paths are not mounted in the assistant workspace. The archived t=0
inventory/budgets remain missing; later configuration evidence was reviewed
in the wall-bundle findings. The private small-case solver is built and its
four baseline screens pass, but no matched SC1 reactor case is ready to launch.
Preserve all existing runs/builds; pause further three-field closure
development while selecting the simpler benchmark. No coefficient tuning,
production change or biological schedule change is part of this decision.
