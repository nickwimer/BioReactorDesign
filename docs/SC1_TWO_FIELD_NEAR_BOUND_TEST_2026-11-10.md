# SC1 two-field near-bound source and timestep test

## Status

**All four near-bound cases and the read-only comparison are complete.** The
coarse pair ran at `78ebbb8`; the fine pair and comparison ran at `75dc652` in
`/kfs3/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z`. All four reports
have `complete: true`, `diagnostic_pass: true` and `screen_pass: false`; their
only errors are the gas and liquid mass thresholds. Both timestep levels
exercise the target regime and source-split effect. Halving dt leaves the
cumulative liquid deficit at about 0.900 micrograms in either reference
setting. Preserve all four cases and the verified build; do not rerun them.

All 53 focused baseline/near-bound tests passed before this runtime result,
including mocked runner execution and unchanged canonical baseline reports.
This update records the supplied results and source-algebra review. No
tolerance, recipe, solver source or binary has changed. The next implementation
is an isolated, guarded source-split candidate described below, not another
timestep/reference sweep or a reactor launch.

All four [baseline cases](SC1_TWO_FIELD_CONSERVATION_TEST_2026-11-10.md) passed
on Kestrel. Their uniform 50/50 source case had essentially zero source-split
contribution. This new fixture tests the existing `1e-4` denominator floor in
`phaseSystemSolve.C` near a phase-fraction bound without changing that floor
or adding mass back. It does not reproduce the reactor's earlier three-field
probe, establish interphase-transfer behavior, or validate the interface.

Reuse the verified build at
`/kfs3/scratch/nwimer/sc1_two_field_build_20261110T224950Z`, including executable
SHA-256 `1483c42b8bb3aeb37be743eafbd70f1b3e59d3ed053388ccd54f7c63a044f1f5`.
No compiled source, observer header, binary or original case recipe changes.
Preserve the completed baseline at
`/kfs3/scratch/nwimer/sc1_two_field_cases_20261110T233431Z` and all failed roots.

## Controlled inputs

| Setting | Coarse pair | Fine pair |
|---|---|---|
| Case names | `bound_coarse_no_reference`, `bound_coarse_gas_reference` | `bound_fine_no_reference`, `bound_fine_gas_reference` |
| Fixed timestep | 0.001 s | 0.0005 s |
| Steps / physical duration | 100 / 0.1 s | 200 / 0.1 s |
| Initial gas / liquid fractions | 0.00001 / 0.99999 | Identical |
| Gas initial volume / mass | 1e-8 m3 / 1e-8 kg | Identical |
| Liquid initial volume / mass | 0.00099999 m3 / 0.99999 kg | Identical |
| Liquid mass feed / integrated feed | 0.001 kg/s / 0.0001 kg | Identical |

Both pairs derive from the baseline source fixture: 16 x 4 x 4 cells in a
0.1 m cube, stationary walls, fixed-pressure outlet, constant densities 1 and
1000 kg/m3, laminar isothermal phases, zero gravity/surface tension, six alpha
subcycles, one alpha correction and two outer iterations. MRF, PBM, biology
and interphase transfer remain disabled. The outlet permits volume displaced
by the feed to leave; conservation is assessed from measured phase exports.

Each reference pair differs only by `referencePhase gas;`. Between timestep
levels only `deltaT`, `maxDeltaT` and the final-write interval differ. Initial
and reverse-flow outlet fractions are complementary. Gas `residualAlpha`
remains 1e-8, below the initial gas fraction. The source denominator floor
remains 1e-4. These settings are a numerical fixture, not reactor parameters.

## Reporting and acceptance

The original report entry point retains its baseline contract. A shared
validator now accepts an explicit immutable timestep/inventory/tolerance
contract; the new report adds near-bound checks. Four canonical baseline
JSON regression checks remain identical to the preceding implementation.

Each report requires complete initial, step, alpha-call, term and stage
coverage; a unique library marker; solver `End`; finite values; prescribed
initial inventories; and the measured source, flux and stage identities.
The accepted-step budget is
`deltaMass + boundaryOut - optionSource - transferSource`.
No repeated outer-iteration diagnostic is added to physical mass.

The mass screen uses the cumulative absolute accepted-step residual, with
phase-specific allowance `A + R * max(abs(totalSource), abs(totalBoundaryOut))`:

| Phase | A (kg) | R | Interpretation |
|---|---:|---:|---|
| Gas | 1e-15 | 1e-3 | Small absolute allowance relative to its 1e-8 kg initial inventory; throughput term uses measured export |
| Liquid | 1e-12 | 1e-6 | 1.01e-10 kg when the prescribed 1e-4 kg feed exceeds export magnitude |

These are predeclared numerical screening thresholds, not experimental
validation tolerances. They are deliberately separate from the baseline's
larger source-test allowance. Reports also expose residual/initial-mass and
residual/throughput ratios, maximum absolute step residual and alpha extrema.

Additional diagnostic tolerances are explicit in
`scripts/diagnostics/sc1_two_field_bounds_cases.py` (gas / liquid):

| Quantity | Tolerance |
|---|---|
| Initial volume | 1e-16 / 1e-13 m3 |
| Initial mass | 1e-16 / 1e-10 kg |
| Volume identities | 1e-16 / 1e-16 m3 |
| Mass identities | 1e-18 / 1e-12 kg |
| Frozen-density term identities | 1e-16 / 1e-11 kg |
| MULES equation L1 | 1e-16 / 1e-14 m3 |
| Budget arithmetic and per-step source | 1e-18 / 1e-12 kg |

All accepted steps must have gas alphaMax < 1e-4 and liquid alphaMin > 0.9999.
Gas alphaMin must remain above 1e-8 to keep the intended source-floor regime
separate from residual-alpha effects. The liquid final-outer source-split
absolute-volume sum must exceed 1e-18 m3. These checks provide accepted-step
regime evidence plus an observed split effect; the existing observer does
**not** provide cellwise alpha bins at each source assembly. Report this
limitation when interpreting causality. A split effect alone does not predict
the accepted mass error, which also includes the subsequent solver response.

`diagnostic_pass` checks coverage, inputs, regime and accounting. `screen_pass`
also requires both mass budgets to meet their thresholds. A complete mass-
screen failure remains evidence and can proceed to the fine pair if its
diagnostics pass. Missing evidence, wrong regimes or inconsistent accounting
block the fine launch. Do not loosen tolerances automatically after a failure.
Preserve the log and inspect the failing quantity first.

## Provenance and preservation

Preparation revalidates all four saved baseline results against their full
logs, input hashes, attempts, mesh-check evidence and verified build/runtime
libraries. It creates a fresh root with four recipes and records the baseline
log hashes, build assets and repository commit. No baseline file is written.
Every launch repeats that validation and checks the near-bound recipes.
Both target cases must be untouched before either is launched. Attempted,
modified or redirected cases are refused. There is no automatic rerun or
partial-pair recovery in this new runner.

`compare` is read-only and verifies all four completed saved reports against
their logs. It prints coarse/fine L1 residuals and ratios, final-mass and export
differences, and final-outer source-split/rescale diagnostics. Two timestep
levels establish sensitivity only, not a convergence order or a validated
numerical correction.

## Completed coarse-pair result at `78ebbb8`

The user supplied both complete reports, log hashes and
`NEAR_BOUND_LEVEL=coarse complete` for the root above. The raw logs and saved
`result.json` files remain on Kestrel; they have not been copied into this
repository. The fine runner revalidates them against the full logs, initial
inputs and build before proceeding.

Both variants covered 100 steps to 0.1 s. Initial gas mass was approximately
1e-8 kg and liquid mass 0.99999 kg. Both supplied 1e-4 kg liquid, with zero
reported gas source and zero interphase transfer. The signed mass residual
and cumulative absolute residual have equal magnitudes in these reports.

| Quantity | No reference | Gas reference |
|---|---:|---:|
| Gas signed mass residual (kg) | +8.996994139032305e-13 | +8.99699284789516e-13 |
| Gas screening tolerance (kg) | 5.739404595110341e-15 | 5.739845087698944e-15 |
| Gas residual / tolerance | 156.7583 | 156.7463 |
| Gas outlet export (kg) | 4.7394045951103415e-12 | 4.739845087698943e-12 |
| Liquid signed mass residual (kg) | -9.001408870364855e-10 | -8.996987667297552e-10 |
| Liquid screening tolerance (kg) | 1.01e-10 | 1.01e-10 |
| Liquid residual magnitude / tolerance | 8.9123 | 8.9079 |
| Liquid inventory increase (kg) | 3.839705193797733e-09 | 3.840145840428025e-09 |
| Liquid outlet export (kg) | 9.999526015391923e-05 | 9.999526015539291e-05 |
| Liquid final-outer split, frozen-density mass (kg) | -8.996993607492838e-10 | -8.996993257360911e-10 |
| Liquid final-outer rescale, frozen-density mass (kg) | -4.41516614935189e-13 | 0 |
| Liquid final-outer implicit source, frozen-density mass (kg) | -5.758221545354363e-17 | -5.758636620798509e-17 |
| Gas accepted-step alpha range | 9.99421744130655e-06 to 9.99668100770661e-06 | 9.99421038927473e-06 to 9.99668098211615e-06 |
| Liquid accepted-step alpha range | 0.999990003318992 to 0.999990005782559 | 0.999990003319018 to 0.999990005789611 |
| Near-bound accepted steps, each phase | 100 | 100 |

The liquid deficit is approximately 0.900141 / 0.899699 micrograms, or
9.00141 / 8.99699 ppm of the injected liquid. This is a shortfall in the
feed/export/inventory budget; liquid inventory still increases slightly.
The gas residual is approximately 0.900 ng, 0.0090% of its initial inventory
but about 18.98% of its small outlet export. Both are failed predeclared
screens; small absolute size does not turn them into passes.

The liquid source-split term closely matches the accepted-step deficit in
both variants. For no reference, subtracting the split and rescale terms from
the signed physical residual leaves approximately -9.67e-18 kg. For gas
reference, physical residual minus split is approximately +5.59e-16 kg.
These are comparisons of separately measured quantities, not a new physical
budget formed by adding outer-iteration terms. The result supports source
splitting as the dominant measured contribution in this fixture. It does not
quantify the historical reactor defect or justify a source addback.

Gas reference removes rescaling here but does not remove either phase's
mass-budget failure. Its gas fraction is reconstructed as the complement;
zero reported explicit gas split/rescale terms are not evidence of a closed
gas budget. In the no-reference variant the gas split is
+8.99699360748976e-13 m3, closely matching its gas residual at rho=1 kg/m3.
Phase-sum enforcement alone therefore does not establish separate phase-mass
conservation. Neither reference setting is selected as a numerical remedy.

The maximum timestamp normalization remained 1.1102230246251565e-16 s,
with 3,360 / 2,128 adjusted records. This is not a recurrence of the earlier
reporting failure: both reports are complete and their diagnostics pass.

Solver-log SHA-256:

- `bound_coarse_no_reference`: `fd9347de2238a61671094fdc1c73fdef1c3379f08bb94cdcebd8da72b07619e0`
- `bound_coarse_gas_reference`: `cbcc168865f05f6b4eb36f3681909d80be8f285a5efdb972538d4fd0093c4e82`

Both coarse diagnostics passed, allowing the planned fine comparison with
unchanged thresholds and physics. That comparison is now complete below.

## Completed fine pair and timestep comparison at `75dc652`

The uploaded `Pasted text(20261111-011306).txt` contains both full 200-step
reports, `NEAR_BOUND_LEVEL=fine complete` and the successful read-only
comparison. Its parsed JSON and provenance are preserved in
[the fine-result evidence](sc1_two_field_bounds_fine_results_2026-11-10.json).
This is a transcript-derived record, not a replacement for the raw solver
logs, initialized inputs or saved reports on Kestrel.

Both fine runs reached 0.1 s with dt=0.0005 s, valid diagnostics, zero gas
source/interphase transfer and 1e-4 kg assembled liquid feed. Their only errors
are the same phase mass-budget thresholds. Signed residual magnitudes equal
the reported L1 residuals, so cancellation is not hiding the cumulative error.

| Phase / reference | Coarse L1 residual (kg) | Fine L1 residual (kg) | Fine/coarse | Change |
|---|---:|---:|---:|---:|
| Gas / none | 8.996994139032305e-13 | 8.997591237177871e-13 | 1.0000663664 | +0.006637% |
| Gas / gas | 8.99699284789516e-13 | 8.997483498542197e-13 | 1.0000545350 | +0.005453% |
| Liquid / none | 9.001408870364855e-10 | 9.000017554931859e-10 | 0.9998454336 | -0.015457% |
| Liquid / gas | 8.996987667297552e-10 | 8.997481055502226e-10 | 1.0000548393 | +0.005484% |

Fine liquid deficits are 0.900001755 / 0.899748106 micrograms, about 8.91
times the unchanged 1.01e-10 kg threshold. Gas residuals are about 180.81 /
180.80 times their fine-case allowances (4.976220520411336e-15 /
4.976462509774087e-15 kg). The gas allowance decreased because measured export
decreased; its absolute residual stayed essentially constant. Gas residuals
are now about 22.63% of export and 0.0090% of initial gas mass.

| Fine quantity | No reference | Gas reference |
|---|---:|---:|
| Liquid final mass (kg) | 0.999990003076462 | 0.999990003076714 |
| Liquid inventory increase (kg) | 3.0764614060707575e-09 | 3.076714197192131e-09 |
| Liquid export (kg) | 9.99960235368385e-05 | 9.999602353769733e-05 |
| Liquid signed residual (kg) | -9.000017554931859e-10 | -8.997481055502226e-10 |
| Liquid split frozen-density mass (kg) | -8.997591103659746e-10 | -8.997590906558331e-10 |
| Liquid rescale frozen-density mass (kg) | -2.4191759706582156e-13 | 0 |
| Liquid implicit-source frozen-density mass (kg) | -1.5772932778359835e-17 | -1.5773566377824442e-17 |
| Gas signed residual (kg) | +8.997591237177871e-13 | +8.997483498542197e-13 |
| Gas accepted-step alpha range | 9.9956281053072e-06 to 9.99751497400528e-06 | 9.99562420733735e-06 to 9.99751495789347e-06 |
| Liquid accepted-step alpha range | 0.999990002485026 to 0.999990004371895 | 0.999990002485042 to 0.999990004375793 |
| Near-bound steps, each phase | 200 | 200 |

The liquid deficit still closely matches source splitting. Residual minus
split and rescale is approximately -7.28e-16 kg without reference, and
+1.10e-14 kg with gas reference. These differences are small against the
approximately 9e-10 kg deficit but are not asserted to be identically zero.
Final-outer diagnostics remain separate from the accepted physical budget.

The maximum per-step liquid residual approximately halves (9.4592e-12 to
4.7574e-12 kg without reference; 9.0003e-12 to 4.4996e-12 kg with gas
reference), while the number of steps doubles. The cumulative deficit does
not improve materially. Final liquid mass also changes by approximately
-7.63e-10 kg between levels, offset by increased export of similar magnitude;
the near-equal residual does not mean the full solution is timestep-converged.
Two levels do not establish a convergence order or a timestep-independent
limit, but neither half-dt nor gas reference remedies this measured failure.

Fine log SHA-256:

- `bound_fine_no_reference`: `7533ca91f3bf7bcf57d34311be3f7ed002c15794bdc9f86d2c67c56695f87046`
- `bound_fine_gas_reference`: `2d4685f0a427d380427b0a9a727d7f7ca6dfe8a121aad34eccc0820165b8fe29`

The comparison revalidated all four saved reports and reproduced their
previously recorded hashes. Fine time-key normalization has 6,720 / 4,256
adjusted records and the same maximum 1.1102230246251565e-16 s adjustment.
No new report-coverage or timestamp error occurred.

## Source algebra and next implementation

The inspected source is
`OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseSystem/phaseSystemSolve.C`.
Let `a` be the phase fraction at source assembly, `d` its raw coupled
dilatation, and `epsilon=1e-4`. Excluding the separately recorded original
divergence source, the current Su/Sp pair evaluates at that same `a` to:

```text
d > 0:  Su + Sp*a = d*(1-a)/max(1-a, epsilon)
d < 0:  Su + Sp*a = d*a/max(a, epsilon)
```

For positive liquid dilatation at a=0.99999, or negative gas dilatation at
a=0.00001, this retains approximately 10% of `d`. It suppresses approximately
90% of the **coupled dilatation component**, not 90% of the imposed feed.
This attenuation is algebraic and contains no timestep factor; refining time
alone does not remove it at fixed source-assembly fractions. Implicit evolution
can add another contribution and remains measured separately.

For this constant-density, pure-phase, fixed-mesh, no-transfer fixture,
`multiphaseEulerFoam/pEqnComps.H` supplies the phase dilatation from the
mass-source/rho term. The uniform liquid feed gives the rate estimate
`qL = (0.001 kg/s)/(1000 kg/m3 * 0.001 m3) = 0.001 s^-1`, with qG=0.
With complementary phases, `dL = aG*qL` is initially 1e-8 s^-1. Holding that
initial fraction fixed estimates the suppressed liquid mass over 0.1 s as:

```text
(1 - 1e-5/1e-4) * 1e-5 * 0.001 * 0.001 * 0.1 * 1000 = 9e-10 kg
```

This 0.900-microgram algebraic estimate agrees closely with the measured
deficit. It is a source-code calculation using nominal uniform inputs, not
a reconstruction of every stored pressure/source-assembly value.

A narrow candidate can factor the two-phase rates before dividing by a
small fraction. For two complementary phases and nonnegative frozen rates
`qi`, `qj`, the raw dilatation is `(1-a)*qi - a*qj`. The candidate is:

```text
Su_increment = qi
Sp_increment = -(qi + qj)
Su_increment + Sp_increment*a = (1-a)*qi - a*qj
```

There is no small-fraction denominator. For the source-only backward-Euler
cell update, `aNew=(aOld+dt*qi)/(1+dt*(qi+qj))`, so 0<=aNew<=1 for
0<=aOld<=1, qi,qj>=0 and dt>=0; the complementary update sums to one. This
identity and the endpoint/interior bounds were checked with exact rational
arithmetic for zero/positive rates, both pure states, near-bound and interior
fractions, and multiple timesteps. It is not a proof for the full flux-limited,
pressure-coupled solver. Su/Sp still combines with the original divergence
source and transport fluxes, and the pressure rates are lagged.

The next implementation should therefore be a **separate private experimental
variant**, initially restricted to this two-field constant-density,
nonnegative-rate, no-interphase-transfer contract. Require complementary
fractions and finite compatible rates, and reject unsupported signed-rate
states rather than silently applying this algebra. Preserve the existing
original-divergence source, pressure/source evaluation, flux limiting and
all diagnostics. Demonstrate the guard and algebra before building, then
regress closed transport, mid-fraction feed and the near-bound pairs at both
timesteps against the unchanged budgets and thresholds. A new variant needs
its own verified build and fresh cases; the current baseline build cannot be
relabeled as corrected. No new build/run command is supplied in this results
update. No production source change, floor reduction, post-solve addback,
additional timestep sweep or reactor launch follows from these results.

## Coarse-pair launch reference (completed; do not repeat)

The following command completed successfully and created the recorded root.
It is retained for provenance, not as the next action. All four cases are now
complete and must be preserved.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
BASELINE=/scratch/nwimer/sc1_two_field_cases_20261110T233431Z
NEARBOUND_CASES="/scratch/nwimer/sc1_two_field_bounds_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
        prepare --baseline-root "$BASELINE" --out "$NEARBOUND_CASES"; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
        run --root "$NEARBOUND_CASES" --level coarse --run
else
    echo "Near-bound setup stopped. Preserve any created directory and paste the output."
fi
```

The resulting reports and exact root are recorded above. No rebuild, manual
dictionary edit, baseline rerun or reactor run is needed.

## Fine-pair launch reference (completed; do not repeat)

The following command completed at `75dc652`, running both fine cases at
dt=0.0005 s to 0.1 s and printing the read-only comparison recorded above.
It is historical; the cases are now attempted and a repeated launch is refused.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
NEARBOUND_CASES=/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
        run --root "$NEARBOUND_CASES" --level fine --run; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
            compare --root "$NEARBOUND_CASES"
    else
        echo "Fine pair stopped. Preserve existing cases and paste the output."
    fi
else
    echo "Fine-pair setup stopped. Preserve existing cases and paste the output."
fi
```

Both fine reports and the comparison have been received and reviewed.

The local validation command is:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover \
    -s scripts/diagnostics -p 'test_sc1_two_field*.py' -q
```

All 53 tests passed before the coarse run, covering baseline report preservation, both timestep
contracts, missing/duplicate/off-grid records, small-inventory thresholds,
source/regime failures, build reuse, read-only evidence checks, fresh-case
guards and coarse-to-fine gating. OpenFOAM execution is mocked locally; the
existing compiled observer capability has not changed. Both timestep levels
now supply runtime evidence. This result update checks the extracted JSON,
comparison ratios, budget arithmetic, source-only algebra, shell syntax and
links without repeating unchanged CFD or harness tests.
