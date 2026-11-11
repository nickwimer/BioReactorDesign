# SC1 two-field near-bound source and timestep test

## Status

**The coarse pair completed at `78ebbb8`: both 100-step reports have
`complete: true`, `diagnostic_pass: true` and `screen_pass: false`.** Their only
errors are the gas and liquid mass-budget thresholds. The near-bound regime
and measurable source-split effect were exercised in both variants. The next
action is the already prepared half-timestep pair in
`/kfs3/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z`, followed by the
read-only comparison. Preserve both completed coarse cases; do not rerun them.

All 53 focused baseline/near-bound tests passed before this runtime result,
including mocked runner execution and unchanged canonical baseline reports.
This result update changes documentation only. No fine-pair result is available
yet, and no tolerance, recipe, solver source or binary has changed.

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

The planned fine comparison is allowed because both diagnostics pass. Keep
all thresholds and physics unchanged. Halving dt is a sensitivity test, not
an assumed correction; the earlier three-field test worsened with smaller dt.

## Coarse-pair launch reference (completed; do not repeat)

The following command completed successfully and created the recorded root.
It is retained for provenance, not as the next action. Both coarse cases are
complete; the fine cases are prepared and remain unattempted.

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

## Next Kestrel command: fine pair and read-only comparison

The coarse review is complete. Both cases pass the diagnostic prerequisite,
so run the already prepared fine pair at dt=0.0005 s, 200 steps to the same
0.1 s endpoint. This command reuses the exact root and build, verifies the
coarse evidence, runs only the two untouched fine cases, and prints a read-only
comparison if both complete. A completed failed mass screen is included in
that comparison. Invalid coarse diagnostics or attempted fine cases stop it.

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

Return both fine reports and the comparison, including any failure output.
Do not rerun an attempted case or prepare a replacement root after a failure.

The local validation command is:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover \
    -s scripts/diagnostics -p 'test_sc1_two_field*.py' -q
```

All 53 tests passed before the coarse run, covering baseline report preservation, both timestep
contracts, missing/duplicate/off-grid records, small-inventory thresholds,
source/regime failures, build reuse, read-only evidence checks, fresh-case
guards and coarse-to-fine gating. OpenFOAM execution is mocked locally; the
existing compiled observer capability has not changed. The coarse Kestrel
results above now supply runtime evidence; fine-pair behavior remains to be
measured. This documentation-only update checks its shell syntax, links and
reported arithmetic without repeating unchanged solver or harness tests.
