# Two-phase reset for the stirred pilot tank

The stirred-tank runs go back to two phases (gas and liquid). The early
two-phase runs lost liquid; the three-phase setup held the liquid inventory
but transfers little oxygen. The changes here are the ones
that held the inventory, without the bubble/headspace gas split.

**Status.** No case with these changes has been run in OpenFOAM yet. The
scripts are tested on the two-phase cases that are in this repository
(`scripts/putida_stir/nonreact`, `tutorial_cases/stirred_tank`); the tests
check the dictionary text only. `check_run.py` and the `report` command have
been run on the production cases (below).

## What the earlier two-phase runs show (2026-11-17)

Read from the solver logs under `/projects/microenviro/final_runs` with
`check_run.py`. All are gas/liquid runs with no reference phase unless noted.

| Run | Tank top | Fill | Liquid inventory | Phase sum before rescaling |
|---|---:|---|---|---|
| `nonreact/166rpm_0.41vvm_0.23O2_outlet` (Apr) | 1.741 m | y < 1.016 m | -22% in 14 s, nearly all before 6.3 s | 0.0004 to 1.99 |
| `nonreact/..._outlet_lowlvl` (Apr) | 1.741 m | | -4.1% from 0.2 to 1.7 s, then the time step collapses | 0.0002 to 1.90 |
| `nonreacting_Julia_08282026/putida_sc1_pilot` (Aug) | 1.397 m | y < 1.19 m, 1152 L | -0.09% in 9.9 s (1152.1 to 1151.1 L) | 0.27 to 1.50 |
| `reacting_Julia_08282026/emm/putida_sc1` (Sep), continuing it | 1.397 m | | -1.7% from 9.9 to 50 s (1151.1 to 1131.4 L) | 0.0008 to 1.96 |
| `nonreact/..._legacyStock_2phase_V3` (Sep) | 1.741 m | y < 1.016 m, 1000 L | -1.6e-6 in 0.3 s | exactly 1 |
| `nonreact/..._legacyStock_2phaseRefGas_V3` (Sep), reference phase gas | 1.741 m | y < 1.016 m, 1000 L | -1.6e-6 in 0.3 s | not applicable |

- **The original failure was liquid loss.** On the full-height tank the April
  run lost 22% of its liquid in 14 s. The loss stops at 6.3 s, when the phase
  sum returns to between 0.9 and 1.
- **The production two-phase runs used a different tank and fill.** The mm,
  emm and two-phase meml cases run in a tank shortened to 55 in and filled to
  1.19 m, about 1152 L. They did not gain liquid: the "151 L excess" reported
  earlier is the difference between that fill and the 1.016 m prescription
  left in the reacting directories.
- **So the two-phase and three-phase production runs are not like for like.**
  They differ by 0.17 m in liquid level and 0.34 m in tank height, which
  affects every comparison of surface shape, impeller wetting and oxygen
  transfer between them.
- **The fixes already hold for 0.3 s.** The two September tests on the
  full-height mesh at 1000 L use the bounded flux reference, a liquid-tight
  top, an impeller ramp and a capped time step, and conserve liquid. No
  two-phase run at that fill has yet gone further without losing liquid.
- In the production two-phase run the time step falls from 6.2e-5 to 1.3e-5 s
  between 9.9 and 50 s while the maximum Courant number stays at its 0.5 limit
  and the mean stays near 0.006: a few very fast cells set the step.

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

1. Done: the earlier logs (above). The rescaling is active and large in every
   run that loses liquid.
2. Extend the September two-phase test on the full-height mesh at 1000 L from
   0.3 s to several seconds, with reference phase gas and with reference
   phase liquid. This is the first two-phase run at the real fill with the
   fixes in place.
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
