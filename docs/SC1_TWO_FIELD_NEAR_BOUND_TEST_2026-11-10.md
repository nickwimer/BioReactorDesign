# SC1 two-field near-bound source and timestep test

## Status

The separate fixture, reporting contract and guarded runner are implemented.
All 53 focused baseline/near-bound tests pass locally, including mocked runner
execution and unchanged canonical baseline reports. **No near-bound CFD has
run yet.** The next action is the coarse pair below, then review its reports
before the half-timestep pair. Local tests validate the harness, not CFD.

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

## Next Kestrel command: prepare and run the coarse pair

Use the active allocation. This command pulls the feature branch over SSH,
reuses the existing build, creates a fresh root and launches only the coarse
pair. The fine cases are prepared but remain unattempted.

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

Return `NEARBOUND_CASES=...`, both reports and any failure output. No rebuild,
manual dictionary edit, baseline rerun or reactor run belongs in this step.

## Fine-pair and comparison reference

After reviewing the coarse reports, use the exact root printed above. In a
new shell set `REPO`, `PY` and `NEARBOUND_CASES` explicitly to those same paths.
The command below assumes those variables are already set and the repository
still contains this implementation. It refuses incomplete or invalid coarse
evidence and attempted fine cases.

```bash
if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" && -n "${NEARBOUND_CASES:-}" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
        run --root "$NEARBOUND_CASES" --level fine --run
else
    echo "Fine-pair setup stopped. Preserve existing cases and paste the output."
fi
```

After both fine reports complete, run this read-only comparison in the same
sourced OpenFOAM/libfabric environment:

```bash
PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_bounds.py" \
    compare --root "$NEARBOUND_CASES"
```

The local validation command is:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover \
    -s scripts/diagnostics -p 'test_sc1_two_field*.py' -q
```

All 53 tests pass, covering baseline report preservation, both timestep
contracts, missing/duplicate/off-grid records, small-inventory thresholds,
source/regime failures, build reuse, read-only evidence checks, fresh-case
guards and coarse-to-fine gating. OpenFOAM execution is mocked locally; the
existing compiled observer capability has not changed, but the new fixture's
actual runtime behavior remains to be measured on Kestrel.
