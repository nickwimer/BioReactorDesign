# SC1 matched alpha audit: completed pair and numerical convergence checks

**Numerical follow-up:** Both `four_outer` and `half_dt` are complete in
`/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z`. Neither resolves the mass
source discrepancy; half dt worsens it. Read [the latest result and new
source-feasibility build](sc1_alpha_half_dt_results_2026-11-10.md) first.
All preparation/launch blocks below are historical; do not repeat them.
The matched entrainment-pair results remain valid.

## Completed evidence

The control and normal cases in
`/scratch/nwimer/sc1_alpha_pair_20261110T155635Z` both completed successfully.
The supplied transcripts are `Pasted text(20261110-161201).txt` (control) and
`Pasted text(20261110-174718).txt` (normal). Each records 2,000 steps from
31.1999823727831 to 31.2999823727864 s, status 0, `End`, and parallel
finalization. ExecutionTime was 742.2 s and 713.84 s, respectively.
**Do not rerun either case.**

The normal preflight proved all four completed-control checkpoint mismatches
were exclusively outer Courant-header `location` changes. Every other header
byte and the complete payload matched the pristine normal copy. All 95 normal
mesh/checkpoint files matched the preparation manifest. No file was repaired,
normalized or discarded. The original refusal was a successful protection
against an unexplained mismatch; the subsequent proof resolved it.

Both runs reused the verified private build
`/scratch/nwimer/sc1_alpha_audit_20261107T001119Z`:

| Artifact | SHA256 |
|---|---|
| `bin/putidaFoamSC1AlphaAudit` | `43aeeef4984dd016ba94664bfe64fb498458237a6ab925e15a2322357a49ca96` |
| `lib/libbirdphaseSystem.so` | `bc6d9da6e63caf8392917fd318d41097f1a488826730b36be33f0e41e9b6dd15` |

Runtime checks confirmed unchanged original/private sources, installed library
and private dependency resolution. All original v1 budget fields for all five
phases/aggregates and the legacy bubble record reproduce each run's first
report in the earlier [full pair](sc1_conservation_pair_results_2026-11-06.md)
exactly at printed precision. This is observational regression, not proof of
bitwise field equality or physical validity.

## Matched results

Both endpoints and durations match. The generic reader's `--native-times`
warning describes its comparison mode; the measured endpoint offset here is
zero. OTR below is the timestep-integrated liquid O2 transfer divided by the
actual interval, not an endpoint transfer rate.

| Quantity | Zero entrainment | Normal entrainment |
|---|---:|---:|
| Mean OTR (mg/s) | 96.78662351 | 96.92974676 |
| Integrated liquid O2 transfer (mg) | 9.678662351 | 9.692974677 |
| Actual liquid mass increase (mg) | 123.5137416 | 141.1689414 |
| Legacy liquid mass residual (mg) | −181.2249153 | −163.5840279 |
| Legacy combined-gas mass residual (mg) | −1.714739 | −1.742346 |
| Legacy total mass residual (mg) | −182.9396541 | −165.3263735 |
| Final bubble mass (kg) | 0.009334465444 | 0.009345555917 |
| Requested liquid O2 removal (mg) | 49.78475862 | 49.78474527 |
| O2 added by clipping (mg) | 40.11355116 | 40.10612280 |
| Requested removal canceled by clipping (%) | 80.57396 | 80.55906 |
| Final mean liquid DO (mg/L) | 0.01935258 | 0.01935948 |

Normal entrainment adds only **0.1431233 mg/s (0.147875%)** to mean OTR,
or 0.0143123 mg O2 over the interval. Bubble inventory is 11.09047 mg higher
at this endpoint, but the earlier 0.5 s experiment ended 125.23 mg lower.
The short gain is not evidence of persistent bubble-holdup improvement.

The 17.65520 mg increase in liquid mass gain is mostly a change in the
numerical source response, not additional physical oxygen delivery.
The residuals in this table are the original sampled-source budget residuals;
the stage/alpha ledgers below explain their origin without certifying the
sampled budget as the exact discrete physical balance.

## Source accounting explains the apparent improvement

These are **final-outer** cumulative density-weighted alpha terms
(`sameCallIndexCumulativeRhoBaseMass_kg`, call 2), converted to mg.
Summing source integrals from both outer iterations would double-count
physical advancement unless the old-time resets were included.

| Liquid source term (mg) | Control | Normal | Normal minus control |
|---|---:|---:|---:|
| Stored pressure dilatation | +304.7386880 | +304.7529958 | +0.0143078 |
| Mixture-coupling difference | −2.6913799 | −2.8575212 | −0.1661412 |
| Source split at assembly alpha | −96.4383965 | −94.2928002 | +2.1455963 |
| Implicit source change as alpha evolves | −82.0951389 | −66.4337022 | +15.6614367 |
| **Effective source action** | **+123.5137726** | **+141.1689722** | **+17.6551996** |

Final liquid alpha clamping contributes only −0.00003117 mg in each run.
The signed equation residual and limited-flux divergence are below 4e-7 mg
in magnitude. Equation-residual L1 is approximately 0.0007532 mg in each run.
The complete measured alpha allocation agrees with the independent mass-stage
allocation within 5e-12 mg at printed precision. The unassigned all-outer
alpha allocation is below 5e-10 mg in magnitude.

Bubble source splitting contributes +79.7773463 mg in control and
+80.7665185 mg in normal. The legacy bubble residual is approximately
+77.63575 mg and +78.55746 mg. These are related but are not identical budgets:
gas density, pressure, flux and outer-reset contributions remain necessary.
Density-weighted divergence must not be relabeled as boundary mass flux,
and stored gas pressure dilatation must not be treated as physical transfer
alone. See [the control derivation](sc1_alpha_control_results_2026-11-10.md).

Normal's total O2 residual is +0.05112865 mg. Post-species inventory changes
contribute +0.05151231 mg, and stabilization contributes −0.0003837032 mg;
the remaining unexplained term is 3.74e-11 mg. Its clipping correction agrees
with the independent negative-O2 integral within 2.4e-13 mg. Thus the species
accounting closes, while the requested biological removal remains mostly
canceled by clipping. Accounting closure is not validation of uptake.

## Next experiment: coupling convergence before source replacement

The [OF9 PIMPLE documentation](https://doc.cfd.direct/openfoam/user-guide-v9/fvsolution)
states that `nOuterCorrectors` repeats the coupled equation system within each
timestep. A controlled outer-count change can establish sensitivity to
incomplete nonlinear coupling. It cannot alone prove convergence or repair
the source split near phase-fraction bounds.

Prepare two fresh zero-entrainment cases with
[`sc1_alpha_convergence.py`](../scripts/diagnostics/sc1_alpha_convergence.py).
Use the completed control configuration and assets, but copy the complete
original checkpoint/mesh from `comparison_manifest.json:initial_state_source`,
after verifying every original raw hash. Do not copy the completed control's
rewritten Courant headers, and do not alter either completed case.

| Case | Outer iterations | Initial/max timestep (s) | Action now |
|---|---:|---:|---|
| Completed control | 2 | 0.00005 | Reuse its result |
| `four_outer` | 4 | 0.00005 | Run first |
| `half_dt` | 2 | 0.000025 | Prepare; leave unrun pending review |

Each new case covers the same 0.1 s interval. All other configuration and
runtime assets are preserved, including six alpha subcycles, one alpha
corrector, zero entrainment, zero foam multiplier and the intentional
biological schedule. `nAlphaCorr > 1` is not supported by this audit and is
not requested. No solver/library rebuild is needed.

Review actual timestep/outer counts, phase mass residuals, call-index source
decompositions, alpha inventory identities, species coupling lag, clipping,
OTR and bubble inventory. With unchanged timestep and no early termination,
four outer iterations give 8,000 alpha calls in 2,000 steps; compare call 4
with baseline call 2, retaining all-outer reset/density terms for closure.
Do not sum iterations as extra elapsed physical time.

If increasing outer iterations reduces discrepancies, establish the remaining
coupling sensitivity before accepting settings. If substantial source defects
persist, use the half-timestep result to separate temporal sensitivity from
the phase-fraction-floor/source formulation. A replacement then needs bounded,
phase-consistent donor/receiver mass, species and PB transfer—not an arbitrary
epsilon reduction or cancellation of the implicit term.

The [current Foundation development source](https://cpp.openfoam.org/dev/phaseSystemSolve_8C_source.html),
inspected 2026-11-10, still has the corresponding positive/negative dilatation
split with a configurable `vDotResidualAlpha` floor, assembled before alpha
subcycling. It also has newer MULES options. Therefore an OpenFOAM upgrade is
not, by itself, an evidenced remedy for this case.

## Launch and verification

All **131 local SC1 tests pass**, including 12 new tests of one-variable
changes, preserved source/baseline bytes, exact checkpoint copying, input
tampering, redirected inputs, repeat-run refusal, opt-in launch and a mocked
single-case SLURM launch. These checks do not substitute for a Kestrel result.
No new numerical-variant CFD has yet run.

The launcher pins the known binary/library hashes, checks source/runtime
provenance, verifies every prepared input and baseline-log hash, requires an
active allocation and exclusively creates the selected case's log. A failed
or partial case is preserved. `run` without `--run` only checks inputs.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
NUMROOT="/scratch/nwimer/sc1_alpha_numerics_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" pull --ff-only \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && "$PY" "$REPO/scripts/diagnostics/sc1_alpha_convergence.py" prepare --out "$NUMROOT"; then
    echo "NUMROOT=$NUMROOT"
    "$PY" "$REPO/scripts/diagnostics/sc1_alpha_convergence.py" \
        run --root "$NUMROOT" --case four_outer --run
else
    echo "Preparation stopped. Preserve any created directory and paste the output."
fi
```

Return the complete printed output and `NUMROOT`. Keep `half_dt` unrun until
the first result is interpreted. These remain frozen-OUR reacting-state
diagnostics, not the nonreacting 9.9 s checkpoint. Free-surface/foam evidence,
conservative transfer/PB validation and oxygen-availability-consistent biology
remain separate gates before production or MM/EMM/MEML comparisons.
