# SC1 alpha-subcycle control: result and next matched experiment

**Superseded launch state:** The matched normal case has now completed
successfully. Read [the pair results and next numerical test](sc1_alpha_pair_results_2026-11-10.md).
Both alpha cases are complete; do not rerun the historical launch block below.
The control measurements and derivation remain valid.

## Completed run and observational regression

The uploaded `Pasted text(20261110-161201).txt` records a successful
52-rank control in
`/scratch/nwimer/sc1_alpha_pair_20261110T155635Z/control`.
It reports status 0, `End`, and `Finalising parallel run`, with 2,000 steps
from 31.1999823727831 to 31.2999823727864 s. ExecutionTime was 742.2 s.
The normal case in this pair remains unrun at the time of this report.

The existing build `/scratch/nwimer/sc1_alpha_audit_20261107T001119Z` was reused.
Both recorded binary hashes passed, provenance checks passed, and the
executable/film/entrainment dependency graphs resolved the private phase-system
library. Initial-state preparation proved full equality, including all four
Courant metadata exceptions before copying the complete control checkpoint to
both new cases. The checkout printed source commit `6db9259`; subsequent
commits before this run changed documentation only.

Comparison with `Pasted text(20261106-235025).txt`, using complete unique
record lines, gives exact equality at displayed precision for:

| Record | Matching records |
|---|---:|
| `SC1_BALANCE_START` | 5 |
| `SC1_BALANCE` | 5 |
| `SC1_MASS_STAGE` | 25 |
| `SC1_MASS_STAGE_TOTAL` | 5 |
| `SC1_OXYGEN_STAGE` | 5 |
| `BUBBLE_MASS_AUDIT` | 1 |

The new startup marker is version 1. Both final alpha outer-iteration blocks
are present, with solve counts 3,999 and 4,000, call indices 1 and 2, and six
MULES subcycles for each explicitly solved phase in each call. The reference
headspace phase has no independent MULES solve. This validates observational
regression of this private build over the tested interval, not all states or
physical model accuracy.

## Liquid mass discrepancy is attributed to source discretization

The final outer iteration's density-weighted alpha source budget is below.
All values are cumulative over the actual 0.1 s interval, in mg. Liquid
density is effectively constant; later-outer oldTime resets cancel the first
outer iteration's physical advancement.

| Final-outer term | Liquid mass contribution (mg) |
|---|---:|
| Stored pressure dilatation input | +304.738688007178 |
| Mixture-coupling difference | −2.691379930760 |
| Su/Sp splitting defect at source-assembly alpha | −96.438396535935 |
| Implicit source change as alpha evolves | −82.095138926709 |
| **Effective assembled source action** | **+123.513772613773** |
| Negative of limited-flux divergence | +0.000000219096 |
| MULES equation residual | +0.000000035949 |
| Final alpha bounds clamp | −0.000031166630 |
| **Predicted alpha allocation** | **+123.513741702189** |
| Directly measured alpha allocation | +123.513741702166 |

The pressure input agrees with the old sampled physical-source prediction
304.7386569 mg to about 0.0000311 mg. The measured complete liquid mass change
is 123.513741641545 mg; its approximately 6.06e-8 mg difference from alpha
allocation is consistent with the separately measured energy/density and
summation contributions. The new alpha allocation agrees with the independent
mass-stage allocation to 5e-12 mg.

Thus the approximately **181.224915 mg shortfall** is dominated by the
source split and implicit source evolution, with a smaller mixture-coupling
contribution. Final alpha clipping and liquid boundary transport are
negligible at this scale. This is distinct from the much larger *oxygen*
clipping already measured in the species equation.

The first outer iteration changes liquid alpha mass by approximately
−10,787.5868611883 mg, including its clamp. The second outer's oldTime reset is
+10,787.5868611885 mg, canceling to about 2.5e-10 mg. Adding both outer source
integrals as physical-time injection would therefore be incorrect. The
all-outer source decomposition remains available but is not substituted for
the final-outer physical advancement above.

### Why these terms arise

The source code is
[`phaseSystemSolve.C`](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseSystem/phaseSystemSolve.C).
For positive raw dilatation `G`, it uses
`Sp = -G/max(1-a,1e-4)` and adds `G/max(1-a,1e-4)` to `Su`.
For negative `G`, it adds `G/max(a,1e-4)` to `Sp`.
At assembly alpha `a`, the net dilatation source is therefore
`G*(1-a)/max(1-a,1e-4)` or `G*a/max(a,1e-4)`, respectively.
Near a phase-fraction bound, that differs from the intended raw dilatation.
The measured net splitting defect is −96.4384 mg for liquid in the final
outer iteration.

Su/Sp are assembled outside the subcycle loop. MULES applies Sp to the updated
alpha in every subcycle. The measured `Sp*(alpha_after-a)` contribution is
another −82.0951 mg for liquid. This establishes the discrete mechanism, but
does not by itself establish which replacement discretization is physically
and numerically acceptable. The discrepancy is not an introduced audit change:
the original phase equations were retained and all old records reproduce.

The two source effects sum to −178.533535463 mg, about 58.59% of stored
liquid pressure input. The remaining source difference is the measured
mixture coupling. More subcycles alone do not remove the denominator-floor
distortion; the coefficients would still be assembled outside that loop.

## Gas phase and numerical closure

The final-outer bubble source split has the opposite net effect:

| Final-outer bubble term | Contribution (mg) |
|---|---:|
| Stored pressure dilatation input | +196.892646183281 |
| Mixture-coupling difference | −0.000141207818 |
| Splitting defect at assembly alpha | +79.777346348976 |
| Implicit source change | +0.150755738604 |
| Effective assembled source action | +276.820607063043 |
| Negative of density-weighted flux divergence | −75.372040146193 |
| Net first-outer/reset carry with different frozen densities | +0.005067308952 |
| **Measured alpha allocation** | **+201.453634225800** |

The split suppresses removal for negative dilatation near vanishing phase
fraction; its net bubble contribution is positive here. The spatial shares
of positive/negative contributions have not been separately measured, so do
not assign this entire number to a particular vessel region or transfer model.

The +79.7773 mg splitting effect is comparable to the old +77.6357 mg bubble
mass residual. It is not an exact attribution of that complete residual.
Gas pressure/energy density changes, stored-pressure timing and the
density-weighted internal transport terms remain part of the full mass
balance. A weighted divergence is not the physical boundary mass flux.
In particular, the gas `divU` from
[`pEqnComps.H`](../OFsolvers/birdmultiphaseEulerFoam/putidaFoam/pEqnComps.H)
contains density/compressibility contributions in addition to external and
interphase sources. Do not equate it directly with the final sampled source.

The new alpha allocations match the independent stage allocations for bubble,
headspace and combined gas at displayed precision. Headspace changes by
−183.336301586604 mg through reference reconstruction, and combined-gas alpha
allocation is +18.1173326391958 mg. Total alpha allocation is
+141.631074341361 mg. These are alpha-stage allocations, not full-timestep
mass changes.

For the all-outer alpha equations, net residuals are approximately
6.16e-8 mg (liquid) and −2.80e-12 mg (bubble); accumulated cellwise L1 residuals
are 0.0015064 mg and 1.77e-8 mg respectively. The latter are not silently
reported as zero. Unassigned alpha inventory changes are approximately
−4.59e-10 mg for liquid and smaller for gas phases. Flux-averaging defects are
approximately 2.31e-10 mg (liquid) and −6.21e-14 mg (bubble). The optional
implicit phase-pressure correction is inactive.

The species result remains unchanged: negative-O2 correction adds
40.113551 mg and cancels 80.574% of requested removal. The successful alpha
accounting does not validate that biological uptake. These remain reacting-state
diagnostic continuations, not a nonreacting 9.9 s checkpoint.

## Next experiment and subsequent correction design

Run the **existing, unused normal case in this same pair** with the same
executable/library and checkpoint. It is the missing matched alpha audit.
Do not rerun control, create another pair, change the timestep settings,
change entrainment efficiency, or change foam area for this comparison.
The previous "leave normal unrun" gate is satisfied for this specific pair:
the control reproduces prior results and its dominant discrepancy is now
attributed. Older completed experiments and the older unrun stage-2 normal
case remain untouched.

After both logs are available, compare source splitting, implicit evolution,
gas inventory, gas/O2 residuals and OTR at their actual report times.
The numerical phase-source effects must be distinguished from the physical
entrainment increment before modifying the closure.

The next solver-development target is consistent, bounded phase-source
integration and its coupling to pressure and density. Candidate work should
separately test source-coefficient refresh within subcycles/nonlinear
iterations and a phase-coupled bounded source treatment. If transfer is limited,
the applied donor/recipient mass, species and PB birth must remain consistent.
Simply lowering the denominator floor, dropping Sp, renormalizing alpha, or
increasing an entrainment coefficient is not a validated correction.
Time-step and outer-iteration convergence will need to be measured; the current
audit explicitly supports only one alpha correction and cannot be silently
used with `nAlphaCorr > 1`.

No solver implementation or physical configuration is changed in this report.
The biological update schedule remains intentional and unchanged. Conservative
hydrodynamics, a nonreacting validation case, experimental free-surface/foam
evidence and a suitable oxygen-limited biological source treatment remain
separate acceptance requirements before the 9.9 s release checkpoint.

## Normal preflight stop and corrected launcher

The next user transcript confirms the private runtime verifier passed again.
The normal time-directory check, both binary hashes and all normal initial-state
hashes passed. The following `validate_pair(control, normal, checkpoint)` then
stopped on a raw checkpoint difference. **The normal simulation did not launch
and this failed preflight did not modify either case.**

That inline command incorrectly reused the fresh-input comparison after the
control had executed. The preparation manifest is the immutable initial-state
reference. A copied file can legitimately retain the original source's outer
location, while its completed peer can acquire a current-case location if a
diagnostic writes it again. The supplied generic error does **not** identify
which files changed or prove that only location metadata differs. A Courant
header rewrite is a hypothesis to test, not an established finding.

The replacement `scripts/diagnostics/run_sc1_alpha_normal.py` is read-only by
default. With explicit `--run`, it launches only after all checks pass:

1. Both binary hashes must equal those used by the completed control, and
   original-source integrity/private runtime dependency checks must pass.
2. Every normal mesh/checkpoint byte must still match the preparation manifest,
   including all four Courant files. No normal-file metadata exception exists.
3. The control's mesh and every non-Courant checkpoint file must still match
   that manifest. Missing/extra files also stop the launcher.
4. For changed control Courant files only, parse the outer headers and verify
   that the completed header identifies its current file and the pristine
   copy still identifies its manifest-recorded original source. Every other
   header byte and the entire binary tail must match. A changed payload,
   whitespace, type, extra header entry or unexpected location is rejected.
5. Case settings/biology must match apart from the established entrainment
   switches. The completed control log, original normal time directory, absent
   normal logs and active allocation are required.

The existing raw inspection and fresh-case preparation rules remain strict.
No source file or checkpoint is normalized or repaired. This is a Python
workflow correction; **reuse the existing compiled solver and private library**.
If the mismatch is not proved location-only, the launcher stops and names the
affected file. Do not manually bypass that failure.

All **119 local SC1 tests pass**, including 14 new prepared-pair lifecycle
and launcher checks. Tests reproduce the failed old comparison after a
synthetic location rewrite, prove the complete unchanged payload, and reject
changed physical fields, Co payload/header changes, missing/extra files,
changed biology, unexpected time directories and existing logs. A mocked
launcher verifies process-local library paths and exclusive log creation;
no OpenFOAM simulation was executed locally.

The long inline shell/Python block is superseded by this committed launcher:

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python

if [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" pull --ff-only \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    "$PY" "$REPO/scripts/diagnostics/run_sc1_alpha_normal.py" \
        --pair-root /scratch/nwimer/sc1_alpha_pair_20261110T155635Z \
        --audit-root /scratch/nwimer/sc1_alpha_audit_20261107T001119Z --run
else
    echo "Setup stopped; paste the output."
fi
```

Omit `--run` for a read-only check. With it present, a successful preflight
launches the unused normal case once, prints solver status, compares both
logs at their native times, and emits both alpha outer-iteration blocks.
Return the complete output. Any failed check or run preserves the pair;
do not overwrite or restart an attempted normal run in place.
