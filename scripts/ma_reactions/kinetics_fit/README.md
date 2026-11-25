# Fitting the putidaFoam kinetics

One 0D model and one fit tool for both kinetics variants of `putidaFoam`:

- **EMM**: Michaelis-Menten uptake with yield coefficients and a logistic
  biomass cap (`get_rhs` in
  `OFsolvers/birdmultiphaseEulerFoam/putidaFoam/microbeModel.C`).
- **MEML**: the same uptake, with growth, H+ and muconate rates from the
  ME-model surrogate (`get_rhs_ml`).

The parameters are the ones in a case's `constant/microbeUpdateProperties`:
`q_max`, `K_o`, `K_s`, `y_xs`, `y_ms`, `y_os`, `x_max`.

putidaFoam also reads one optional entry, `cap_n`, the sharpness of the
biomass cap on growth, `1 - (X / x_max)^cap_n`:

- It is a plain number (`cap_n 5;`), not a dimensioned entry. 1 is the
  logistic cap, a larger value leaves early growth alone and stops it more
  abruptly near `x_max`, and 0 switches the cap off.
- Without the entry nothing changes: EMM keeps its logistic cap (1) and
  MEML has no cap (0).
- In the fit tool `cap_n` is fixed at 1 unless a fit frees it
  (`--free ... cap_n`). `--model meml` is MEML without a cap and
  `--model meml_cap` MEML with it.

## Files

| File | Purpose |
|---|---|
| `kinetics0d.py` | The model: rates, surrogate, 0D integration, a literal copy of the solver's `advance`, electron balance, dataset loader |
| `fit_kinetics.py` | Least-squares fit, multi-start, profiles, the old loss for comparison; command line |
| `exp17.py` | Extraction of the Exp 17 workbooks and the off-gas analysis (oxygen yield, RQ, kLa against stir rate) |
| `make_report.py` | Runs the whole fit matrix and writes `out/REPORT.md` with figures and candidate parameter sets |
| `datasets.yaml` | The runs: fixed-stir BR1-BR3 and Exp 17 reactors 49, 50, 52, 54 |
| `tests/` | Checks against solver output, the stored fits and the data |

The earlier scripts and the fit behind the production parameters are kept in
`../0dWillComp/`.

## What is not in the repository

The surrogate files, the Exp 17 workbooks and anything derived from them are
unpublished and stay local. Copy `local.example.yaml` to `local.yaml` and
point it at your copies. `local.yaml`, `local/` and `out/` are git-ignored.
Numbers that come from those inputs, such as the values some tests expect,
go in an `expected:` block of `local.yaml` and not into the code.

Extract Exp 17 once (needs pandas and openpyxl; everything else needs only
numpy, scipy, matplotlib, ruamel.yaml and, for MEML, torch):

```bash
python exp17.py extract --xlsx-dir "<Exp17 folder>" --out-dir local/exp17
```

## Use

```bash
python fit_kinetics.py --datasets fixedstir1 fixedstir2 fixedstir3 --free q_max K_o K_s y_xs y_ms x_max theta --fix y_os=2.4
```

```bash
python make_report.py
```

```bash
python -m pytest tests
```

The tests are not part of the repository's CI because torch is not a
dependency of `bird`. Tests that need the local files are skipped when those
are missing.

## How the fit works

- **Dissolved O2.** In the fixed-stir runs O2 is quasi-steady:
  `kLa (o2sat - O2) = y_os qs X`. The data then only determine
  `theta = y_os / kLa`, so `theta` is the fitted parameter and any kLa or
  `y_os` assumption follows by division. In Exp 17 the DO log is prescribed
  and the measured oxygen uptake is fitted, which determines `y_os` without
  a kLa.
- **Weights.** Each quantity is scaled by its scatter across replicates
  (`noise_levels`), so chi2 near the number of residuals means the model
  fits as well as the replicates agree with each other.
- **Search.** Least squares in log-parameters from several starting points.
  `profile` refits the other parameters along a grid of one, which shows
  what the data do and do not determine.
- **MEML.** The surrogate returns zeros outside the convex hull of its
  training inputs, as in the solver. Fits penalise leaving the hull and
  report the time spent outside.
- **Solver parity.** `advance_solver` copies `microbemodel::advance` step
  for step and reproduces a `wellmixed.dat` written by putidaFoam. One test
  also compiles the solver's `microbeModel.C` with a small driver
  (`tests/solver_driver.cpp`, no OpenFOAM needed) and compares it with the
  Python model update by update, with and without the biomass cap.
