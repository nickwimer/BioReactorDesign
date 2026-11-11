# SC1 first-step pressure probe

## Latest result: pressure collapse reproduced before limiting

The probe ran at `1dbe432` in
`/scratch/nwimer/sc1_pressure_probe_20261111T172942Z`, using the new executable
from `/scratch/nwimer/sc1_pressure_probe_build_20261111T172942Z`.
The [verified result](sc1_reactor_pressure_probe_result_2026-11-11.json) reports
all 59 initial inputs unchanged. It reaches sequence 13, the third pressure
solve (first timestep, outer 2), then emits
`SC1_PRESSURE_STOP reason=unclippedPressureBelow1000Pa` with pressure minimum
`-2499398.1871999367` Pa: exactly the printed v4 minimum.

This is the intended pre-limit guard, not a passed startup. The external monitor
reports `solver failure` at the first rank's `FOAM FATAL ERROR`, return code 139,
and wrapper status 1. The full log now confirms the synchronous pressure guard
precedes all eight rank fatal errors. MPI_Abort then reports sigSegv stack frames
and srun reports segmentation faults during shutdown; no separate pre-guard
segmentation fault is shown.
There are no completed-step records, no `End`, and no two-step completion marker.

| Quantity | Outer 1 | Outer 2 |
|---|---:|---:|
| Maximum face diagonal addition, both phases | 0.11282816668027722 | 1451184694.6572018 |
| Minimum mixture rAUf | 4.9995536786813554e-8 | 6.797243007986973e-10 |
| Unclipped pressure minimum after first correction (Pa) | 99982.986887332227 | −2499398.1871999367 |

The maximum face diagonal addition increases by about 1.286e10, while the
minimum mixture inverse coefficient falls by a factor of 73.55. These are
global extrema and need not be colocated. The first outer iteration's second
pressure solve also remains near 99983 Pa. The logged mixture Co maximum of
zero comes from pre-first-step output; it does not establish zero phase motion.

For this case with no virtual mass,
[AFfs](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.C)
adds the face drag coefficient to each phase's diagonal. The
[segregated model](../OFsolvers/birdmultiphaseEulerFoam/interfacialModels/dragModels/segregated/segregated.C)
can be rearranged as

`K = m*rhoPair*|Ur|*G/(a1*a2) + n*muAlphaI*G^2`,

where `G` is the regularized interface-gradient magnitude and `a1`, `a2` are
the residual-limited fractions. This identifies relative velocity and the
local fraction/gradient factors as quantities to inspect. It does not establish
which cells produce the measured jump or justify changing drag coefficients.

**The bundle is now received and verified.** Its archive SHA256 is
`b08f1f185061ab774fc0a16354947737b1e9578346ca21170a01587ca46f1ef0`;
all 20 bundled file hashes and all 13 included initial dictionaries match.
The log SHA256 is
`d2c5d3b55ad096ee1fc8151a68353828fbfb05db9bc1e5c65af72fa0cc5ef656`.
All 581 log lines have been inspected; the result JSON now includes all 144
pressure-probe records. The native mesh, large fields and binaries remain on
Kestrel and were not independently materialized locally.

At the end of outer 1, maximum gas speed is 0.54505945 m/s in absent-phase
cells (alpha <=1e-6), while maximum present-liquid speed is 0.000489242 m/s.
The first bad pressure is at rank 2, local cell 24771, centre
(-0.0117119085, 0.9935527213, 0.1912571129) m: alphaGas=0, alphaLiquid=1,
pre-correction gas speed 0.387578851 m/s and liquid speed 1.30908e-7 m/s.
The outer-2 liquid explicit drag-flux range dominates its predictor range.

**Next action: build and test the [private pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md).**
It eliminates the symmetric two-phase face drag before assembling both the
pressure predictor and pressure coefficient. All v4 initial inputs and original
acceptance limits remain fixed; a fresh monitored run is capped at 10 ms.
This is a source/field-supported candidate, not an established cure. Do not
rerun the completed observer or extend any failed case.

## Earlier evidence and rationale for the probe

The v4 bundle `reactor_failure_evidence.tar(2).gz` is received and verified.
Its archive hash matches `0ac651d664f0834d9ecde87cc05bc7c39989ccf1a27bec4837eb49bbb54669b3`;
all 19 bundled hashes and all 13 included initial dictionaries match. The
native mesh and large initial fields remain on Kestrel. See the
[verified result](sc1_reactor_coupled_failure_result_2026-11-11.json).

The full log confirms that coupled drag is active. At t=0.00005 s:

| Outer iteration | Pressure solve | Final residual | Minimum drag determinant | Result |
|---|---|---:|---:|---|
| 1 | 1 | 9.1420e-13 | 0.9999999999998727 | No limiter |
| 1 | 2 | 9.8065e-13 | 0.9999999999998727 | No limiter |
| 2 | 1 | 9.6392e-13 | 0.0135944981275532 | p minimum −2499398.1871999367 Pa |
| 2 | 2 | 8.7023e-13 | 0.0135944981275532 | Printed during termination |

The determinant remains positive, but its change is consistent with a much
stronger drag coupling product in outer 2. It does not establish conditioning,
the spatial location of the change, or the pressure failure's cause. The first
bad pressure is computed before that correction's drag inversion; the inversion
does not replace the pressure equation's coefficient or predicted mixture flux.
All eight printed pressure solves reach the requested residual tolerance.
All five rate-state records have zero physical source rates. First-step alpha
changes remain tiny, with zero source-split/clamp/rescale terms.

The first failure is effectively unchanged from v3. Another configuration
change was not supported by the available field evidence. The probe above
observes the first-step momentum/pressure transition with all case inputs fixed.

## What the probe changes

`sc1_reactor_pressure_probe.py build` verifies the existing rate build, copies
only its application sources to a fresh directory, adds observation hooks,
and executes one `wmake`. It links the same six verified project libraries.
The libraries, installed executables, original sources and existing runs are
unchanged. The new executable is `birdPressureProbe`; `ldd`, `nm`, source hashes
and binary hashes enforce its dependency and source identity before use.

`prepare` verifies the exact received v4 log
`d6f4e4b1345fdbe0029598c699e9be1f322479ce77152665ff3988c05052a9ea`
and the complete preserved startup/build/regression chain. It copies only `0`,
`constant`, `system` and the three mesh-preparation logs. **All 59 initial input
hashes are identical to v4**, including dt=5e-5 s, two outer iterations, two
pressure corrections, gravity, physical coefficients and boundary conditions.
The explicit command selects the new executable; the copied application entry
in `controlDict` remains unchanged. No failed time or processor fields are copied.

The observation hooks insert reductions without removing, reordering or changing
the original pressure-equation lines. They record:

- Momentum diagonal/source magnitude, old face alpha-density and face diagonal
  addition `AFf` (the drag term in this no-virtual-mass case).
- Phase and mixture inverse momentum coefficients, predicted fluxes, combined
  force fluxes, extra force fluxes and explicit drag fluxes.
- Raw pressure matrix diagonal/source, solved p_rgh and unclipped static pressure.
- Phase velocities in cells with alpha >1e-6 and in cells at/below that threshold.
- Coordinates, local cell label, phase fractions and pre-correction speeds at
  the global internal pressure minimum, with processor identity from `Pout`.

Volume-field ranges cover internal cells; surface-field ranges cover internal
and boundary faces, including duplicate processor copies. Matrix diagonal/source
ranges are raw assembled arrays, not per-volume values; boundary contributions
may be held separately by OpenFOAM. Global extrema of different fields need not
occur at the same spatial location. The fields are observed, never assigned by
the observer. The helpers use the Foundation OF9
[fvMatrix field accessors](https://cpp.openfoam.org/v9/classFoam_1_1fvMatrix.html).

After every pressure solve, **the executable itself stops before flux correction
or pressure limiting if internal static pressure is below 1000 Pa**. Nonfinite
observed values also stop it. This prevents further iterations after the first
observed invalid pressure rather than relying only on polling the solver log.
The existing external alpha/Courant/fatal-error and 30-minute monitor remains.
If no such stop occurs, the executable exits after two completed timesteps
(0.0001 s); it never extends to the original 10-ms end time.

No restart or startup acceptance is produced. There is no reconstruction step;
the probe is intended to expose the pressure failure before a scheduled write.
On either a solver stop or the two-step completion, the runner saves a compact
summary and a bundle at `<probe-root>_evidence/reactor_failure_evidence.tar.gz`.
A reproduced pressure failure returns status 1 and is useful diagnostic evidence.
The complete log contains the location and field records omitted from the compact
console summary. Rotation, aeration, oxygen and biology remain disabled.

## Historical Kestrel command — completed, do not repeat

This command built the application, prepared a fresh copy and ran the probe
above. Its bundle is now verified; do not launch a copy.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
RATE_BUILD=/scratch/nwimer/sc1_two_field_rates_build_20261111T021516Z
FAILED_COUPLED=/scratch/nwimer/sc1_reactor_coupled_20261111T165652Z
SC1_PROBE_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
PRESSURE_BUILD="/scratch/nwimer/sc1_pressure_probe_build_${SC1_PROBE_STAMP}"
PRESSURE_CASE="/scratch/nwimer/sc1_pressure_probe_${SC1_PROBE_STAMP}"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_probe.py" \
        build --base-build "$RATE_BUILD" --out "$PRESSURE_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_probe.py" \
        prepare --failed-root "$FAILED_COUPLED" --probe-build "$PRESSURE_BUILD" --out "$PRESSURE_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_probe.py" \
            run --root "$PRESSURE_CASE" --run
        echo "Pressure-probe status: $?"
    else
        echo "Build/preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed `BUNDLE`. Do not rerun or extend any
existing startup. A pressure-stop result is expected if the v4 failure is
reproduced; a two-step completion is still only diagnostic evidence.

## Verification and limits

All 114 focused baseline, near-bound, rate-split, startup, drag, failure, face,
coupled and probe tests pass. Seven new tests cover application-only compilation
orchestration, unchanged original sources/case inputs, pinned dependencies,
rejected tampering/reuse, observation ordering, pre-limit stopping hooks,
monitored failure collection and absence of reconstruction. A scalar C++ test
executes the actual finite-value range accumulator and two-step exit logic.

The full OpenFOAM application and MPI field reductions cannot be compiled/run
locally. The subsequent Kestrel execution reached the intended pressure guard,
as recorded above; its raw-log bundle is now verified. The 114-test result
is the original prelaunch result; the subsequent pressure-coupling change has
its own 122-test record in the linked candidate document.
Neither that result nor the intended diagnostic stop establishes physical
validity or a working operating simulation.
