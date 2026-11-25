# Fit behind the production putidaFoam parameters

As-run copy of the script and files that produced the kinetic parameters used
in the putidaFoam cases:

| y_xs | y_ms | y_os | x_max | q_max | K_o | K_s |
|---|---|---|---|---|---|---|
| 0.04498 | 0.3189 | 1.195 | 7.6 | 8.26 | 0.1134 | 0.4008 |

They are the rounded contents of `fixedstir1/fitted_param.dat`
(order: q_max, K_o, y_xs, y_ms, y_os, K_s).

## Settings

Hard-coded in `setmodelconstants` of `fit_bioreaction_model.py`:

- Data: `fixedstir1` only.
- kLa fixed at 336 1/h. Case files carry this value as `kLa`, some of them
  as 428; the solver reads it but does not use it.
- o2sat 0.188125 mol/m3, bio_max (x_max) 7.6 kg/m3.
- Starting guess `fixedstir1/paper_param.dat`, after Jamshidzadeh et al.
  (2025).
- Loss: sum over biomass, glucose, muconate and dissolved O2 of
  MSE / mean(obs^2). `scipy.optimize.minimize` with bounds >= 0, Radau
  integration at default tolerances.

The script was last saved on 2025-11-21 and `fitted_param.dat` was written on
2026-03-30.

## Differences from `../gDCW_conv/fit_bioreaction_model.py`

- kLa 540 -> 336 and bio_max 7.65 -> 7.6.
- `paper_param.dat` has q_max 8.26435 instead of 8.26425.
- An O2-against-glucose overlay plot was added. It needs pandas and
  `scatter_data_biomass_dilution.csv` (metabolic-model output, not included
  here). The fit and `fitted_param.dat` are written before that plot.

## Rerunning

```bash
python fit_bioreaction_model.py -d fixedstir1 -pf
```

This overwrites `fixedstir1/fitted_param.dat`. A rerun reproduces the loss
(7.1e-3) and the parameters to within a few percent.
