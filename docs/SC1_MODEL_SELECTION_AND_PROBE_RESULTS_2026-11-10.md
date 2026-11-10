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

## Next decision: test the simpler candidate first

The [external benchmark review](SC1_EXTERNAL_BENCHMARKS_2026-11-10.md) identifies
two-field BiRD stirred-tank implementations with oxygen transfer and PBM.
The [Foundation framework](https://openfoam.org/guides/disperse-multiphase-flows/)
supports resolved-interface and dispersed interpretations of phase fraction.
These establish a candidate route, not SC1 validation.

1. **Read-only baseline reconstruction first.** Review the archived two-field
   case's t=0 initialization, earliest available liquid inventory, physical
   boundary fluxes, sources, phase/pressure settings and generated MRF geometry.
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
inventory/budgets and complete current configurations have not been retrieved
in this review. No corrected two-field case or new solver build is ready to
launch yet. Preserve all existing runs/builds; pause further three-field closure
development while selecting the simpler benchmark. No coefficient tuning,
production change or biological schedule change is part of this decision.
