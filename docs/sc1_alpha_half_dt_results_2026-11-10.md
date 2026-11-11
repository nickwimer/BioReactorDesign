# SC1 half-timestep result and local source-feasibility implementation

**Build follow-up:** The new source-probe build compiled and passed runtime
verification at `/scratch/nwimer/sc1_alpha_source_probe_20261110T191459Z`.
Read [the verified build and fresh-control launch instructions](sc1_alpha_source_probe_run_2026-11-10.md)
for the current action. Reuse that build; the build block below is historical.

## Completed numerical sensitivity checks

`Pasted text(20261110-190007).txt` records status 0, `End` and parallel
finalization for
`/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z/half_dt`.
It completes 4,000 physical timesteps and 8,000 alpha calls, retaining two
outer iterations and six alpha subcycles. The timestep is 2.5e-5 s.
ExecutionTime is 1402.43 s. Both numerical cases in this root are now complete;
**preserve them and do not repeat their launch blocks.**

The preparation manifest, original checkpoint, baseline log and private
runtime provenance passed. All five initial phase/aggregate records match
the baseline exactly at printed precision. All three numerical comparisons
start at 31.1999823727831 and end at 31.2999823727864 s, with actual duration
0.10000000000329834 s. The reader's `normal` label is a numerical variant
with zero entrainment; its native-times warning does not imply a measured
timestamp mismatch here.

| Quantity over the matched interval | Baseline, 2 outer | 4 outer | Half timestep, 2 outer |
|---|---:|---:|---:|
| Actual liquid mass change (mg) | +123.513742 | +150.352259 | **−14.272415** |
| Legacy liquid mass residual (mg) | −181.224915 | −154.384700 | **−319.010443** |
| Legacy bubble mass residual (mg) | +77.635748 | +75.814453 | +75.702664 |
| Legacy headspace mass residual (mg) | −79.350487 | −77.571083 | −77.344707 |
| Legacy combined-gas mass residual (mg) | −1.714739 | −1.756630 | −1.642042 |
| Legacy total mass residual (mg) | −182.939654 | −156.141330 | −320.652485 |
| Mean O2 transfer to liquid (mg/s) | 96.786624 | 96.769652 | 96.780302 |
| O2 added by clipping (mg) | 40.113551 | 40.115150 | 40.112947 |
| Requested O2 removal canceled by clipping (%) | 80.573959 | 80.577167 | 80.572761 |

Halving the timestep increases the liquid shortfall by **76.0301%** relative
to baseline. Liquid inventory now decreases despite approximately +304.738 mg
of sampled net physical source and no liquid boundary flow. Mean OTR changes
only −0.0065312%. Final bubble inventory is 2.590515 mg below baseline.
These results do not establish timestep convergence; neither smaller dt nor
four outer iterations is an accepted remedy. No formal convergence order can
be inferred from this limited set of coupled transient runs.

## The changed source response accounts for the deficit

Final-outer cumulative terms in mg use
`sameCallIndexCumulativeRhoBaseMass_kg`: call 2 for baseline/half dt and call
4 for the four-outer run. They must not be added across outer iterations as
extra elapsed-time injection.

| Liquid source term (mg) | Baseline | Four outer | Half timestep |
|---|---:|---:|---:|
| Stored pressure dilatation | +304.738688 | +304.736983 | +304.738042 |
| Mixture coupling | −2.691380 | −34.863689 | −3.429184 |
| Split at source-assembly alpha | −96.438397 | −84.103473 | **−179.906786** |
| Implicit source evolution | −82.095139 | −35.417540 | **−135.674459** |
| **Effective source action** | **+123.513773** | **+150.352281** | **−14.272387** |

The half-timestep source split and implicit evolution together remove
315.581245 mg relative to the unsplit source. Final alpha clamping contributes
only −0.000027336 mg. Signed equation residual is −1.995e-8 mg; its L1 is
0.00150649 mg. The observed all-outer alpha allocation is −14.272414835005 mg,
matching the independent stage allocation within 1.34e-11 mg at printed
precision. Its difference from actual liquid mass change is about 2.80e-7 mg.
Unassigned liquid alpha allocation is −1.78e-10 mg. Gas alpha allocations
also match the independent stage allocations at printed precision.

In the last half-dt timestep, the candidate liquid alpha mass advancement is
+0.04218129 mg at outer 1 and −0.01189899 mg at outer 2. Variation between
outer iterations therefore remains visible even at the smaller timestep.
This does not establish a periodic two-cycle at every timestep. Bubble source
splitting remains +79.844682 mg, compared with +79.777346 mg in baseline.

Total legacy oxygen residual is +0.07319473 mg. Post-species inventory changes
account for +0.06935935 mg and stabilization for +0.00383501 mg; remaining
unexplained amount is −1.77e-10 mg. Liquid clipping matches the independently
integrated negative oxygen within 1.7e-13 mg. The final negative-cell count is
130,794 and minimum liquid Y is −1.24779e-11. The smaller negative Y per step
has not reduced the approximately 40.113 mg cumulative clipping correction.
The species and alpha diagnostics close, but physical conservation and
biological uptake are not thereby validated.

## Why a local feasibility check precedes replacement

The reviewed [OF9 source](https://cpp.openfoam.org/v9/MULESTemplates_8C_source.html)
and the unchanged BiRD phase equation give, for fixed mesh and unit alpha
density,

```
alpha_after = (alpha_old + dt_sub*(Su - div(alphaPhi_limited))) / (1 - dt_sub*Sp)
```

The split coefficients are assembled before subcycling. Existing ledger terms
already prove the source differs from `q + G`, where `q` is the original
mixture-divergence term and `G` is the cross-phase dilatation. Simply adding
the deficit back, changing the denominator floor, or eliminating Sp does not
prove the resulting phase fractions are admissible. With the current limited
flux held fixed, the source-consistent counterfactual is

```
rawTarget = alpha_old + dt_sub*(q + G - div(alphaPhi_limited))
pressureTarget = alpha_old + dt_sub*(stored_phase_divU - div(alphaPhi_limited))
```

An out-of-bounds `rawTarget` means one cannot simultaneously keep that old
state, fixed flux, prescribed source and phase bounds. Flux/source/pressure
coupling must change; clipping the result would reintroduce a conservation
defect. A bounded target is only a necessary local check, not proof of a valid
coupled three-phase replacement. `pressureTarget` uses stored, lagged
dilatation; it is not a freshly evaluated physical-source truth.

The previous aggregate logs cannot distinguish where these local conflicts
occur. This is the specific missing evidence addressed by the new observer;
another uninstrumented timestep or outer-count sweep would not provide it.

## Implemented isolated observer

New files:

- [`sc1_alpha_source_probe.py`](../scripts/diagnostics/sc1_alpha_source_probe.py):
  preparation/verification of a fresh private build, never a CFD launch.
- [`sc1AlphaSourceProbe.H`](../scripts/diagnostics/include/sc1AlphaSourceProbe.H):
  independent cell arithmetic for source/bounds diagnostics.
- [`build_sc1_alpha_source_probe.sh`](../scripts/diagnostics/build_sc1_alpha_source_probe.sh):
  child-process build wrapper with exclusive logs and failure handling.
- [`test_sc1_alpha_source_probe.py`](../scripts/diagnostics/test_sc1_alpha_source_probe.py):
  compiled arithmetic, observer preservation, provenance and build safeguards.

The preparer uses the established private-alpha preparation, retains local
biology, then adds observations to the copied alpha-audit header. All original
phase-equation and audit-header lines are retained. A mandatory new library
symbol prevents silently linking a build lacking the new observer. The
existing installed library, original BiRD source and completed private builds
are preserved. The executable's name remains `putidaFoamSC1AlphaAudit` within
its **new root**, never a global installed target.

The observer requires both `sc1AlphaAudit true` and the new opt-in
`sc1AlphaSourceProbe true` in alpha controls. It reports for liquid and bubble:

- source-assembly alpha bins: lower (`alpha < 1e-4`), upper
  (`alpha > 1-1e-4`), interior and outside [0,1];
- cumulative signed split/implicit/source contributions in each bin;
- raw/pressure target bound violations and density-weighted bound excess;
- raw-target correction and its independent identity remainder;
- maximum `dt_sub*abs(Sp)`, a source-stiffness indicator, not a Courant number.

Values are accumulated **separately by outer index**. Counts are repeated
phase-cell-subcycle samples, not unique cells. Persistent tables stay local
to each rank; reporting reduces a copy, sums extensive metrics and takes a
maximum for stiffness. Bounds flags use a 1e-12 diagnostic tolerance; neither
fields nor source coefficients are clipped or changed. Excess metrics retain
even sub-tolerance amounts. Headspace is reconstructed and has no independently
solved source row. Joint phase-sum, density, mass-transfer/species/PB consistency
remain required before accepting any replacement.

All **147 local SC1 tests pass**, including 16 new tests. Standalone probe
arithmetic compiled with C++11 and warnings treated as errors. Tests cover a
feasible floor correction, an impossible fixed-flux source, transport restoring
feasibility, raw-versus-pressure targets, phase maxima/tolerance, bin boundaries,
sum-versus-max accumulation and independent equation residuals. Preparation
tests preserve original sources/biology, reject reused/redirected destinations
and changed headers, verify mandatory linkage, and stop after a failed phase
build. **This is not an OpenFOAM compilation or a new CFD result.**

## Historical build instructions (completed; do not repeat)

The following was the build-only gate before the successful transcript.
It is retained as history. Use the follow-up above for the current run step.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PROBEROOT="/scratch/nwimer/sc1_alpha_source_probe_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" pull --ff-only \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    bash "$REPO/scripts/diagnostics/build_sc1_alpha_source_probe.sh" "$PROBEROOT"
else
    echo "Setup stopped. Preserve all existing directories and paste the output."
fi
```

Return `PROBEROOT` and the complete output. The wrapper uses `wmake libso` only
in the new private phase tree, then `wmake` in its private solver, and checks
both old/new mandatory symbols plus executable and loaded-library resolution.
It creates no cases and runs no CFD; build failures return to the allocation
prompt and preserve logs.

After successful build verification, prepare a fresh zero-entrainment probe
case from the original common checkpoint and baseline configuration. Its
purpose is new local source/bounds evidence, with observational regression
against the already completed baseline. That replay is justified by the new
measurements, not a repeat numerical sweep. Use those measurements to select
and test a bounded, phase-consistent source/flux correction in minimal cases
before another physics comparison. No physical closure, biological schedule,
entrainment efficiency or foam coefficient changes are made here. The intended
nonreacting 9.9 s release checkpoint and experimental hydrodynamic validation
remain subsequent gates.
