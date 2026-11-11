# SC1 early-velocity spatial diagnostic

## Completed Kestrel diagnostic

The run at `fa0636a` completed normally, status 0, in
`/scratch/nwimer/sc1_velocity_probe_20261111T183832Z`. The
[received console result](sc1_reactor_velocity_location_result_2026-11-11.json)
reports all 33 prefix records matching the saved v5 trajectory, three
reconstructed snapshots and unchanged hashes for all 59 prepared inputs.
The solver log has an `End` marker, three paired step records and no detected
failure events. This completes the spatial diagnostic; the parent startup
still failed at 0.2 ms.

| Time (ms) | Present-liquid peak (mm/s) | Bulk-liquid peak, alpha >=0.99 (mm/s) | Gas vertical speed at liquid peak (mm/s) |
|---:|---:|---:|---:|
| 0.05 | 0.244979 downward | 0.006329 | +0.081492 |
| 0.10 | 0.529355 downward | 0.042940 | +0.280204 |
| 0.15 | 0.800102 downward | 0.070844 | +0.440398 |

All three liquid peaks occur at y=1.018662536 m, only 0.064374 mm above
the initial interface midplane, with alphaLiquid approximately 0.4948501.
The last two global peaks are not owners of boundary faces. Their wall-adjacent
counterparts have almost identical speeds in the same interface-height layer;
the first global peak is wall-adjacent. Thus a direct wall-owner-only
explanation does not cover the observations. The phase velocities are almost
entirely vertical and oppositely directed at these liquid peaks.

At 0.15 ms, bulk-liquid peak speed is about 11.3 times smaller than the
present-liquid peak. Liquid inventory remains 1000 L to printed rounding.
However, the whole-domain volume-mixture RMS reaches 1.482643 mm/s; quieter
bulk liquid does not establish a quiet gas region or an acceptable startup.
The full bundle has now been received. Its archive SHA256, all 19 bundled
file hashes and all 13 included initial dictionaries match. The 33 prefix
comparisons were independently repeated against the verified v5 log.
Native fields/mesh and compiled binaries remain on Kestrel; the fifteen
reconstructed field hashes are included in the verified result.

The complete gas summary materially extends the liquid-peak interpretation:

| Time (ms) | Gas interface peak (mm/s) | Gas bulk peak (mm/s) |
|---:|---:|---:|
| 0.05 | 17.3839 | 11.5786 |
| 0.10 | 30.0239 | 23.2463 |
| 0.15 | 40.0167 | 34.8828 |

At 0.15 ms the gas interface peak is at radius 0.270417 m and
y=1.018661409 m, alphaLiquid=0.49494025, with velocity
(0.0069859, 0.0004855, -0.0393992) m/s. The bulk-gas peak is in pure gas
at y=1.031215754 m, about 12.6 mm above the interface midplane.
Neither cell owns a boundary face. These peaks are predominantly horizontal;
the normal counterflow at the liquid maximum does not explain the full gas
disturbance. The near-vent gas peak remains only 9.26e-7 m/s.

**The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md)
stopped on the liquid-speed guard at 8.4 ms; its bundle is verified.**
The scalar-gradient comparison also failed; next use the [nonorthogonal pressure comparison](SC1_REACTOR_VOF_NONORTHOGONAL_2026-11-11.md).
It uses the same mesh, fill, initial p_rgh, gravity and physical properties, but
a shared velocity with the Foundation OF9 interFoam equations, read-only phase
budgets and the existing startup guards in a separate private application. This establishes a simpler resolved-interface comparison before
attempting a new two-fluid force closure. It is not a dispersed-bubble, PBM,
oxygen-transfer or biological replacement.

Verified evidence source:

`/scratch/nwimer/sc1_velocity_probe_20261111T183832Z_evidence/reactor_failure_evidence.tar.gz`

- Bundle SHA256: `5ba88937370743795665315d0f266522071bd592d813072a773875ada326191d`.
- Solver-log SHA256: `c0cd5cfde42651b23262ad2cc545cf8f5b66336cc64d4e95e3bc1e11a0f967e8`.

## Force-balance interpretation to test

The v5 transform retains the original phase-specific buoyancy term in
[`pUf/pEqn.H`](../OFsolvers/birdmultiphaseEulerFoam/multiphaseEulerFoam/pUf/pEqn.H):
`interpolate(phase.rho() - rho)*(g & mesh.Sf())`. It changes the coupled
momentum inverse, not the density contrast acting on each phase. The current
[Foundation face corrector](https://cpp.openfoam.org/v13/facePressureCorrector_8C_source.html)
also retains a phase/mixture density difference inside its force construction;
this is not evidence that the remaining term is an implementation mistake.

In an ideal planar continuum comparison, with mixture-hydrostatic pressure,
zero slip and other forces neglected, phase acceleration is
`a_i = (1 - rhoMixture/rho_i)*g`. The same hydrostatic pressure gradient cannot
balance the separate gravitational accelerations of both unequal-density
phases in a mixed cell without an additional interfacial force. A finite drag
coefficient produces zero drag force at exactly zero slip. This supplies a
candidate explanation for interface-normal counterflow in a resolved
transition layer, rather than demonstrating a coefficient or source error.

At the first liquid peak, the continuum liquid acceleration scale is
-4.950565 m/s2; the measured first-step vertical velocity divided by dt is
-4.899575 m/s2. This is a scale comparison only: it omits the actual face
pressure gradient, surface tension, boundary treatment, nonlinear drag and
the two outer iterations. It is not a closed discrete force budget, and it
does not predict the observed gas velocity, especially the now-measured
horizontal peaks. Removing buoyancy contrast,
adding an unverified compensating force or increasing drag is not selected
from this numerical resemblance. The next model test is the separate
single-velocity reference described above; the original source equations
and all BiRD controls remain preserved.

## Original decision based on v5 evidence

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

## What the completed command did

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

## Historical Kestrel command — completed, do not repeat

This command produced the diagnostic above. Its bundle is verified; use the
VOF-reference result and the saved-field inspection command; do not repeat its launch.

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

The console and full bundle are recorded above. Preserve the completed
diagnostic and all earlier roots.

## Local validation

The 122 preceding focused regressions and seven new spatial tests pass. The
combined 128-test run (before adding the boundary-addressing test) passed in
50.766 s; the final seven-test spatial suite also passes. Tests cover exact
output/duration-only changes, preserved controls and binary identities,
rejected reuse/tampering, prefix mismatch/incompleteness, no reconstruction
after solver failure, normal reconstruction of only three steps, actual
face-owner adjacency, phase masks and native pressure/velocity peak statistics.
The received v5 log yields the expected 33 prefix records.

No OpenFOAM execution is available locally. Kestrel has now generated the three
native snapshots, reconstructed them and completed the spatial analysis. The
update recording that console is documentation only; it does not repeat the
129 local tests or alter the diagnostic or existing compiled v5 build.
