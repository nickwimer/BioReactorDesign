# Kestrel check of the biomass cap

The scripts that built putidaFoam with the `cap_n` change on Kestrel and ran
it, as they were run on 2026-11-25, and the script that compares the old and
new runs. They are a record of that check more than a general tool: paths
are those of the agent workspace on Kestrel
(`/scratch/nwimer/agent/accelerate`), where they live in `capn-check/`.

They depend on files of that workspace that are not in this repository:

- `runs/env.sh`: the OpenFOAM 9 environment and the run libraries
  (`build/two-phase-reset/lib`);
- `runs/rx_variant.sh`: copies a reacting case at a written time and puts
  the oxygen floor on the sink;
- `runs/stage1.sh` and each case's `Allrun`: run a case on one node;
- the reacting cases `runs/rx_emm` and `runs/rx_meml` at 60 s.

Results are not kept here. They are in `capn-check/README.txt` on Kestrel
and in the report this tool writes.

## Scripts

| Script | What it does |
|---|---|
| `build.sh` | Builds putidaFoam from a copy of the workspace checkout with the four-file change into `capn-check/bin`, leaving the checkout and its binaries alone. Also tries the same source with `wellMixed.H` switched on, compile only. |
| `setup_case.sh NAME MODE` | Copies `runs/rx_meml` at 60 s into `capn-check/cases/NAME` with one of four biologies (`ref`, `emm`: production numbers without a `cap_n` entry; `F`, `G`: the capped sets). The biology update comes on the first step. |
| `run_case.sh NAME` | Runs such a case for a few steps with the new build and prints what the solver read and the liquid averages after the update. |
| `new_case.sh NAME SRC SET` | Makes a standard reacting run with new parameters in `runs/NAME`: `rx_variant.sh` on `SRC` at 60 s with the 0.1 mg/L oxygen floor, then parameter set `F` (MEML) or `G` (EMM). |
| `run_new_batch.sh NAME END` | Batch script: runs `runs/NAME` to `END` on one node with the new build. |
| `before_after.py DOWNLOADS OUT` | Reads the `liquidO2` and `tankVolumes` monitors of the old and new runs, as downloaded from Kestrel, and writes the tank-mean dissolved oxygen against time as a plot and a table. |

## Order of use

1. `build.sh`
2. `setup_case.sh` and `run_case.sh` for `ref`, `emm`, `F` and `G`: one
   biology update each, compared with `kinetics0d.advance_solver`.
3. `new_case.sh rx_meml_setF rx_meml F` and
   `new_case.sh rx_emm_setG rx_emm G`, then `run_new_batch.sh NAME 100` for
   each.
4. `before_after.py` on the downloaded monitors.

The parameter values written by `setup_case.sh` and `new_case.sh` are sets F
and G of the report. Set F in `new_case.sh` is the one refitted with the
surrogate model file of the Kestrel cases.
