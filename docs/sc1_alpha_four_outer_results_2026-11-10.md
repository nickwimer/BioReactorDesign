# SC1 four-outer-iteration result and half-timestep diagnostic

**Follow-up:** The prepared half-timestep case has now completed and worsened
the liquid shortfall to 319.010443 mg. Read [the result and new build-only
source-feasibility step](sc1_alpha_half_dt_results_2026-11-10.md). Both numerical
cases are complete; do not execute the historical launch block below again.

## Completed run

The user-supplied `Pasted text(20261110-182839).txt` records a successful
52-rank run in
`/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z/four_outer`.
It reports status 0, `End`, parallel finalization, 2,000 physical steps and
8,000 alpha solves. The final step contains calls 1–4, each with six MULES
subcycles for each explicitly solved phase. ExecutionTime is 1346.53 s.
**Preserve this completed case; do not prepare another copy or rerun it.**

Preparation and runtime checks passed at commit `ae04e89`. The only changed
entry relative to the completed zero-entrainment control is
`PIMPLE/nOuterCorrectors: 2 -> 4`. The original mesh/checkpoint, constants,
model assets, biology schedule and all other numerical controls were
preserved. The verified private build remains
`/scratch/nwimer/sc1_alpha_audit_20261107T001119Z`; no rebuild is needed.

All five initial phase/aggregate records exactly match the baseline. Both
runs start at 31.1999823727831 and end at 31.2999823727864 s, an actual
0.10000000000329834 s interval. The comparison reader's `normal` label means
this numerical variant, **not active entrainment**. Its native-times warning
does not imply an actual offset here: the printed endpoint offset is zero.

## More outer iterations have not resolved the discrepancy

The baseline is the completed `control` in
`/scratch/nwimer/sc1_alpha_pair_20261110T155635Z`. All mass amounts below are
integrated over the matched interval. Legacy residuals compare inventories
against the sampled physical-source/flux budget; they are not unexplained
errors in the newly instrumented alpha equation.

| Quantity | Two outer iterations | Four outer iterations |
|---|---:|---:|
| Actual liquid mass increase (mg) | 123.513742 | 150.352259 |
| Legacy liquid mass residual (mg) | −181.224915 | −154.384700 |
| Legacy bubble mass residual (mg) | +77.635748 | +75.814453 |
| Legacy headspace mass residual (mg) | −79.350487 | −77.571083 |
| Legacy combined-gas mass residual (mg) | −1.714739 | −1.756630 |
| Legacy total mass residual (mg) | −182.939654 | −156.141330 |
| Mean liquid O2 transfer (mg/s) | 96.786624 | 96.769652 |
| Integrated liquid O2 transfer (mg) | 9.678662 | 9.676965 |
| O2 added by clipping (mg) | 40.113551 | 40.115150 |
| Requested O2 removal canceled by clipping (%) | 80.573959 | 80.577167 |
| Final bubble inventory (kg) | 0.009334465444 | 0.009332789180 |
| Final liquid DO (mg/L) | 0.019352583 | 0.019352480 |

The liquid shortfall decreases by 26.840215 mg, or **14.8104%**. More than
85% remains, while the combined-gas residual magnitude increases slightly.
Mean OTR decreases by 0.0169719 mg/s (0.0175353%). Bubble inventory is
1.676264 mg lower. Thus neither oxygen delivery nor the conservation problem
is resolved by selecting four outer iterations.

## The source response is not converged

Use the final-outer `sameCallIndexCumulativeRhoBaseMass_kg` terms, converted
to mg: call 2 for the baseline and call 4 for the new run.

| Liquid source contribution (mg) | Baseline call 2 | Four-outer call 4 |
|---|---:|---:|
| Stored pressure dilatation | +304.738688 | +304.736983 |
| Mixture-coupling difference | −2.691380 | −34.863689 |
| Source split at assembly alpha | −96.438397 | −84.103473 |
| Implicit source evolution | −82.095139 | −35.417540 |
| **Effective assembled source action** | **+123.513773** | **+150.352281** |

The smaller shortfall reflects reduced source splitting/implicit evolution
partly offset by a more negative coupling term. It is not simply convergence
of every source component toward the intended pressure input.

The new run's individual outer-index records show a pronounced alternating
response. The following values exclude old-time resets and include the
limited-flux divergence, equation residual and clamp. For effectively constant
liquid density they describe the candidate physical-step alpha mass advance
at each outer index, not independent physical injections to be added together.

| Outer index in the four-outer run | Sum of candidate advances over all steps (mg) | Candidate advance in the last step (mg) |
|---|---:|---:|
| 1 | −10783.545608 | +0.145598217 |
| 2 | +133.897440 | +0.044044956 |
| 3 | +301.771022 | +0.146812260 |
| 4 | +150.352259 | +0.037562117 |

In particular, the terminal step changes substantially from outer 3 to 4;
the cumulative values also do not approach a common value. This is evidence
against accepting four iterations as a converged solution. It does not prove
a persistent two-cycle in every timestep or establish the converged limit.
Call 2 in this run is not the independent two-outer baseline: the previous
timesteps followed a different solution trajectory.

Old-time resets cancel the preceding outer candidate, as required by the
discrete identity in the [alpha audit](sc1_alpha_subcycle_audit_2026-11-07.md).
Summing all outer source terms as elapsed-time mass addition is incorrect.
Gas density varies between solves, so gas comparisons also retain the
corresponding density/reset contributions.

The measured all-outer liquid alpha allocation is 150.352259459779 mg,
agreeing with the independent mass-stage allocation within 4e-12 mg at
printed precision. It differs from actual liquid mass gain by only
3.49e-8 mg. The final-outer liquid clamp contributes −0.000021925 mg,
signed equation residual −9.89e-9 mg, and equation-residual L1
0.000753215 mg. Alpha allocation matches the independent stage allocation
at printed precision for bubble, headspace and combined gas; liquid's
unassigned alpha amount is −4.29e-10 mg. The audit itself still closes.

Bubble source splitting remains +79.269918 mg at the final outer, compared
with +79.777346 mg in baseline. Its legacy mass residual remains large.
Source splitting, density updates and the exact transport/source coupling
must still be reconciled before approving physical conservation.

## Oxygen accounting and next decision

Total legacy O2 residual rises from +0.04792068 to +0.05169483 mg. In the
four-outer run, post-species inventory changes contribute +0.04923126 mg
and stabilization contributes +0.00246352 mg. The remaining unexplained
amount is 4.61e-10 mg. The matrix-component mismatch L1 is about 1.22e-11 mg.
Clipping's measured correction agrees with the independent negative-O2
integral within 3e-13 mg. Requested liquid removal remains 49.7847615 mg;
about 80.577% is canceled by clipping. Extra outer iterations have not
validated biological uptake or eliminated the species time-level effects.

**Next run only the existing `half_dt` case** in the same numerical root.
It uses two outer iterations and initial/max timestep 2.5e-5 s, otherwise
the same baseline configuration and original checkpoint. This isolates
timestep sensitivity from the outer-count change. Compare actual interval,
step count, source split, implicit evolution, coupling, inventories and
clipping; do not assume a smaller timestep guarantees physical conservation.
Expected step count is approximately 4,000, to be verified from output.

If substantial source discrepancies persist, the next implementation pass
should address the bounded, phase-consistent source/coupling formulation with
controlled conservation tests. Neither reducing the denominator floor nor
adding more outer iterations is an established remedy. Preserve the
intentional biological schedule and keep entrainment/foam coefficients fixed.

The existing launcher already supports this run, so no code change or
rebuild is needed. Its 131-test suite passed at `ae04e89`; this follow-up
changes documentation only. The supplied baseline/new logs passed the
existing conservation reader, initial-record equality and the independent
source/allocation arithmetic checks described above. The launch block below
passes `bash -n`.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
NUMROOT=/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z

if [[ -n "${SLURM_JOB_ID:-}" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" pull --ff-only \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    "$PY" "$REPO/scripts/diagnostics/sc1_alpha_convergence.py" \
        run --root "$NUMROOT" --case half_dt --run
else
    echo "Setup stopped. Preserve all cases and paste the output."
fi
```

The launcher checks the untouched half-timestep inputs against their manifest,
baseline-log hash and private runtime provenance, then exclusively creates
that case's log. It does not compare against or rewrite the completed
four-outer checkpoint. Return the complete printed output. This remains a
frozen-OUR reacting-state diagnostic, not the nonreacting 9.9 s release
checkpoint; physical free-surface/foam and PB validation remain pending.
