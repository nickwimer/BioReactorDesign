# Two-phase stirred tank, stage 1: stirring only

NLR pilot-plant tank (full height, 68.5 in) with 1000 L of liquid, gas above
it and both impellers turning at 166 rpm. No sparging, no bubble sizes, no
species and no mass transfer. Solver: `birdmultiphaseEulerFoam`.

The purpose is to show that the two-phase setup keeps its liquid and keeps a
believable free surface before anything else is added. Sparging comes in
stage 2, bubble sizes and oxygen transfer in stage 3.

## What the case fixes compared with the earlier two-phase runs

| Earlier two-phase setup | This case |
| --- | --- |
| Both phase fractions transported, then rescaled to sum to 1 | Only liquid is transported; gas is the remainder (`referencePhase gas`) |
| Limiter referenced to the phase flux | Limiter referenced to the mean flux (`meanFluxReference true`) |
| Liquid could leave through the top | `U.liquid` is `slip` on `outlet`: no liquid flux through the top |
| Drag for bubbles in liquid only | Drag in both directions (bubbles in liquid, drops in gas) |
| Lift, wall lubrication and turbulent dispersion active everywhere, including at the free surface | None of them at this stage; they return one at a time in stage 2 |
| Impellers at full speed from the first step (or after 0.1 s) | Speed ramped over 2 s |

No velocity limiters, no staged start-up and no interface compression are
used. If the case needs one of them, that is a finding.

Turbulence is mixture k-epsilon on both phases, as in the earlier two-phase
runs. (The three-phase template uses k-epsilon on the liquid alone, which
only works because it starts with a trace of liquid in every cell.)

## Run

```bash
# OpenFOAM 9 environment loaded, repository on PYTHONPATH
./Allrun.pre     # mesh, rotating zone, monitoring zones, 1000 L fill
./Allrun         # decompose and run with srun (NPROCS, default 52)
```

`Allrun.pre` takes `TANK_YAML`, `LIQUID_LITRES` (default 1000) and `RAMP_S`
(default 2). `Allrun` continues from the latest time when run again.

## What is recorded

Written to `postProcessing/` every 20 steps:

| Monitor | Meaning |
| --- | --- |
| `tankVolumes` | liquid and gas volume in the whole tank |
| `aboveFill`, `above5cm`, `above10cm`, `above20cm` | liquid above the initial surface and 5.6, 10.6, 20.6 cm above it |
| `wallAbove5cm`, `wallAbove10cm` | the same in the outer 6 cm of the radius: liquid climbing the wall |
| `bulkBelowSurface` | gas more than 3 cm below the initial surface: vortex and entrainment |
| `topMassFlow` | liquid and gas mass flow through the top |
| `maxSpeed` | largest liquid and gas speed |

The solver log also prints the liquid fraction and the phase sum at every
step; `scripts/putida_stir/two_phase_reset/check_run.py log.solver` reads
them.

## Pass criteria for a 10 s run

- liquid volume within 0.1 % of 1000 L throughout;
- no liquid through the top;
- time step settles instead of falling;
- liquid above the surface stays bounded (wall rise of a few cm, not a film
  that keeps growing).
