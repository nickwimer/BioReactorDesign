# Two-phase reset for the stirred pilot tank

The stirred-tank runs go back to two phases (gas and liquid). The original
two-phase setup gained about 15% liquid in ten seconds; the three-phase setup
held the liquid inventory but lost most of the oxygen transfer. The changes
here are the ones that held the inventory, without the bubble/headspace gas
split.

**Status.** Nothing in this directory has been run in OpenFOAM yet. The
scripts are tested on the two-phase cases that are in this repository
(`scripts/putida_stir/nonreact`, `tutorial_cases/stirred_tank`); the tests
check the dictionary text only. The current pilot-tank case (sparger zones,
biology inputs, `HcellProperties`) is not in git, so the first step on the
cluster is to apply these changes to a copy of it.

## Changes

`apply_fixes.py apply SRC DST` copies a case and edits the copy. The source is
never modified. Each change has its own flag; `--recommended` applies all of
them with the defaults shown.

| Flag | Edit | Reason |
|---|---|---|
| `--reference-phase gas` | `referencePhase gas;` in `phaseProperties` | With no reference phase the solver transports both fractions and rescales them to sum to one (`phaseSystemSolve.C`, "Correct the sum of the phase fractions"). In the bulk that turns any gas-volume error into liquid. With a reference phase the liquid is transported and the gas is the remainder. |
| `--mean-flux-reference` | `meanFluxReference true;` in the alpha solver controls | Limits the liquid flux against the mean flux, which keeps the fraction bounded, so the solver's clip to [0, 1] has nothing to remove. |
| `--liquid-tight-top` | top patch: `U.liquid` slip, `alpha.liquid` zeroGradient | No liquid can cross the top of the tank. Gas still leaves and enters there. |
| `--mrf-ramp 2` | `omega table ((0 0) (2 <rate>));` in `MRFProperties` | Starts the impellers over about five revolutions instead of at full speed against still liquid. |
| `--max-co 0.3 --max-delta-t 5e-4` | `controlDict` | Starting values; the old case allowed a 1 s step. |
| `--isothermal` | isothermal phase models | Temperature does not affect the result and the energy equations add cost. |
| `--monitors` | includes `tankMonitors` in `controlDict` | Liquid and gas volume, mean dissolved oxygen, interfacial transfer rate and top mass flow, every 20 steps. |

`--reference-phase liquid` is the alternative used by the lab's OpenFOAM 13
stirred-tank tutorial: the gas (which carries the bubble sizes) is transported
and the liquid is the remainder. Compare the two in the first runs.

Not automated, to be done by hand on the real case:

- **Sparger.** Keep the gas source spread over the real holes and ramped, on
  phase `gas`. A full-rate point source in a few cells is a candidate cause of
  the liquid gain.
- **Water saturation model.** The three-phase case dropped it. It moves a
  fraction of a gram of water per second into the gas; keep or drop it
  deliberately.
- **Interface compression.** Leave the table empty: compression on a single
  gas phase would also act on the bubbles.

## First runs

All nonreacting. `check_run.py` reads the liquid inventory from the solver log.

1. Old two-phase log, no run needed: `check_run.py LOG --volume 1.642`. A phase
   sum away from one confirms the rescaling.
2. Old two-phase case plus `--reference-phase gas` only, 0.5 s. Then the same
   with `--liquid-tight-top` only. This shows which change stops the gain.
3. `--recommended`, stirring without gas, a few seconds: liquid inventory
   within 0.1%, impeller power number plausible.
4. Add sparging and run to a steady state: gas holdup, bubble size, dissolved
   oxygen rise and transfer rate, from the monitors.

```bash
python apply_fixes.py report  /path/to/old_case
python apply_fixes.py apply   /path/to/old_case /scratch/.../ref_gas --reference-phase gas
python apply_fixes.py apply   /path/to/old_case /scratch/.../reset   --recommended
python check_run.py log.putidaFoam --volume 1.642 --csv inventory.csv
python -m pytest scripts/putida_stir/two_phase_reset
```

`check_run.py` exits 0 when the relative change in liquid inventory is within
`--tol` (default 1e-3). The log prints fractions to about six digits, so it
resolves changes down to about 2e-6 of the inventory; use the `phaseVolumes`
monitor for anything finer.
