# SC1 early-velocity spatial diagnostic

## Decision and current evidence

The [pressure-coupled v5 bundle](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
is verified. Its explicit stop is `Quiescent liquid speed exceeded 1 mm/s`,
from `sc1QuiescentStep` after timestep four. The subsequent segmentation-fault
messages occur during MPI_Abort shutdown. The earlier pressure collapse does
not recur in this 0.2-ms interval: all 16 solves converge below 1e-12, pressure
stays between 99982.986886 and 111847.796701 Pa, and all 16 coupled flux checks
agree within 4.8916021e-18 m3/s. All 624 field reductions report finite values.

| Time (ms) | Maximum liquid speed where alphaLiquid >0.01 (mm/s) |
|---:|---:|
| 0.05 | 0.244979 |
| 0.10 | 0.529355 |
| 0.15 | 0.800102 |
| 0.20 | 1.059892 — guard stops run |

The maximum grows throughout the interval. Four-step mass-residual L1 is
4.8130419768e-16 kg for gas and 3.2598279932e-14 kg for liquid, but these are
partial-window observations, not a passed 200-step conservation/startup screen.
The fourth step fails quiescence. No scheduled 0.25-ms fields were written.

The current log reports velocity maxima without their locations. Its location
records belong to pressure minima, so they cannot identify whether the liquid
peak is in the transition layer, at a wall/baffle, or in bulk liquid. That
distinction is needed before changing body-force treatment, initialization or
the resolved-interface model. No momentum closure, physical coefficient or
acceptance threshold is changed by this next diagnostic.

## What the next command does

Reuse `/scratch/nwimer/sc1_pressure_coupled_build_20261111T180325Z` and its
verified executable `birdPressureCoupled`, SHA256
`e5414f2afb8192d2fe88cd4f881ac554a7605326b1938ccbfe5fc62291937dd4`.
No compilation or library change occurs. Verify the received v5 control log
`4c469e3438576ad12139cad251de18b159dd35161cdcb7aede16124b193ecaac`
and the complete existing source/build/control chain.

Create a fresh t=0 copy. The only input changes are two controlDict entries:

- `endTime 0.00015`, limiting the run to the first three timesteps.
- `writeInterval 1`, writing each of those timesteps.

All other 58 input files and every other controlDict entry remain identical.
The same timestep, eight ranks, pressure/drag equations, 1000-L fill, boundaries,
physical coefficients and synchronous pressure/flux/1-mm/s guards remain active.
The external alpha/Co/fatal monitor stays active, with a ten-minute solver cap.
The 1-mm/s guard is not bypassed; the three-step interval ends before the
already-observed fourth-step failure. This is a new spatial observation of the
early trajectory, not an attempt to pass a shortened startup.

Require normal solver completion at exactly three printed times, and compare
33 pressure, coupled-flux, phase-budget and liquid-speed records against the
same prefix of the received v5 log (absolute tolerance 1e-12 and relative 1e-9).
A failed solver is never reconstructed. After normal completion, require all
five diagnostic fields on all eight processors at all three times, reconstruct
only 0.00005–0.00015 s, then validate native field lengths, finite values,
fractions, liquid inventory and agreement with the logged liquid speed.

`result.json` contains gas/liquid peak, RMS and vertical velocities for present,
bulk, mixed-interface and boundary-adjacent cells. Peak records contain native
cell label, coordinates, radius, height relative to the initial interface,
phase fractions, both velocity vectors, slip, volume-mixture velocity,
reconstructed static pressure and its change from initialization. Boundary
adjacency uses actual face-owner labels; it is not a distance-to-wall measure.
Full reconstructed fields remain in the new case. Their hashes and the spatial
summary are included in the normal evidence bundle.

`VELOCITY_PROBE=COMPLETE` means only that this diagnostic completed. It does not
accept the failed v5 startup, produce an approved restart, or validate stirred,
aerated, oxygen-transfer or biological behavior. The peak locations will select
the next force-balance/model test; increasing drag or relaxing velocity limits
is not selected in advance.

## Next Kestrel command

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
FAILED=/scratch/nwimer/sc1_pressure_coupled_20261111T180325Z
VELOCITY_PROBE="/scratch/nwimer/sc1_velocity_probe_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_velocity_probe.py" \
        prepare --failed-root "$FAILED" --out "$VELOCITY_PROBE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_velocity_probe.py" \
            run --root "$VELOCITY_PROBE" --run
        echo "Velocity-probe status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed `BUNDLE`. Preserve all earlier roots.

## Local validation

The 122 preceding focused regressions and seven new spatial tests pass. The
combined 128-test run (before adding the boundary-addressing test) passed in
50.766 s; the final seven-test spatial suite also passes. Tests cover exact
output/duration-only changes, preserved controls and binary identities,
rejected reuse/tampering, prefix mismatch/incompleteness, no reconstruction
after solver failure, normal reconstruction of only three steps, actual
face-owner adjacency, phase masks and native pressure/velocity peak statistics.
The received v5 log yields the expected 33 prefix records.

No OpenFOAM execution is available locally. Kestrel still needs to generate
the three native snapshots; the existing compiled v5 build is preserved.
