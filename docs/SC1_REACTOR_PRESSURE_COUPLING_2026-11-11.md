# SC1 pressure-coupled face-drag candidate

## Status

**Kestrel application build passed; the startup stopped at 0.2 ms, status 1.**
The run at `408028e` is `/scratch/nwimer/sc1_pressure_coupled_20261111T180325Z`;
its build is `/scratch/nwimer/sc1_pressure_coupled_build_20261111T180325Z`.
The [console result](sc1_reactor_pressure_coupling_result_2026-11-11.json)
records executable SHA256
`e5414f2afb8192d2fe88cd4f881ac554a7605326b1938ccbfe5fc62291937dd4`
and verified linkage to the same six rate-build libraries.

**The raw bundle is now received and verified:** archive SHA256
`69e5585474cc286dc476341487f03303b7c380f23cac0129c925bd673cf00db0`,
all 19 bundled file hashes and all 13 included initial dictionaries match.
The full log SHA256 is
`4c469e3438576ad12139cad251de18b159dd35161cdcb7aede16124b193ecaac`.
All 59 inputs are reported unchanged by the HPC collector; large fields,
mesh, build binaries and external controls remain on Kestrel.

The explicit fatal message is **`Quiescent liquid speed exceeded 1 mm/s`**,
from `sc1QuiescentStep` at `sc1FaceDragCoupling.H:92`. At 0.0002 s the maximum
speed in liquid cells with alpha>0.01 is 0.0010598917525355049 m/s. MPI_Abort
and sigSegv shutdown records follow the guard. This is an intended quiescence
stop, not evidence of an independent segmentation fault preceding it.

| Time (ms) | Present-liquid maximum (mm/s) |
|---:|---:|
| 0.05 | 0.244979 |
| 0.10 | 0.529355 |
| 0.15 | 0.800102 |
| 0.20 | 1.059892 — failed |

All 16 pressure solves converge below 1e-12. Unclipped static pressures stay
between 99982.986886 and 111847.796701 Pa. All 16 mixture/phase flux checks
agree within 4.8916021e-18 m3/s; 624 field reductions contain no nonfinite values.
The previous first-step pressure collapse did not recur in this short interval.
Four-step phase-mass residual L1 is 4.8130419768e-16 kg for gas and
3.2598279932e-14 kg for liquid. These partial-window observations do not accept
the planned 200-step startup: the fourth step fails quiescence, no `End` exists,
and no 0.25-ms field write was reached.

The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md)
has now completed at `fa0636a`, using the unchanged compiled application.
All 33 compared records match this v5 prefix. The liquid peaks lie in the
approximately 50/50 interface layer with liquid moving downward and gas upward;
the bulk-liquid peak remains much smaller. The last two global peaks are not
direct wall-face owners, while wall-adjacent peaks have almost equal speeds.
Its full bundle is now verified and also shows predominantly horizontal gas
peaks of 40.0 mm/s near the interface and 34.9 mm/s in the bulk-gas layer above.
**The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md)
built, then stopped at 8.4 ms; inspect its existing evidence bundle next.**
It uses the same mesh and initial fill. The spatial diagnostic does not accept
this failed parent startup or establish a two-fluid force correction.

## Evidence and decision

The [pressure-probe bundle](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md) is received:
archive SHA256, all 20 bundled file hashes and all 13 included initial dictionary
hashes match. The full 581-line log confirms:

- The first outer iteration ends with maximum absent-gas speed 0.54505945 m/s
  (alpha <=1e-6), versus maximum present-liquid speed 0.000489242 m/s.
- In outer 2, the maximum symmetric face drag grows to 1.4511846947e9. The
  explicit liquid drag-flux range is approximately -2.389e-4 to +1.538e-4 m3/s;
  it dominates the range of the liquid predictor. These global extrema are
  not necessarily colocated.
- The pressure minimum is in a pure-liquid cell: rank 2, local cell 24771,
  centre (-0.0117119085, 0.9935527213, 0.1912571129) m. Its pre-correction gas
  speed is 0.387578851 m/s, liquid speed 1.30908e-7 m/s, alphaGas=0.
- The synchronous pressure guard fires before invalid pressure is corrected
  or limited. The subsequent MPI_Abort shutdown reports segmentation faults;
  the log does not show a separate segmentation fault preceding the guard.

The OF9 face path adds a lagged other-phase drag flux to the pressure predictor,
uses only the diagonal momentum inverse in its pressure coefficient, then calls
`partialEliminationf` after the pressure solve. The latter restores the solved
mixture flux, so it cannot repair the pressure just computed. The evidence
supports testing elimination before pressure assembly. This is a candidate
mechanism, not a proof that all quiescent-flow defects are explained.

## Numerical change

For each face, let a and b be the original phase momentum diagonals before
drag, k the common nonnegative face drag coefficient, and q the pressure-flux
gradient. The coupled equations are

```
(a+k)*phiGas - k*phiLiquid = rhsGas + cGas*q
-k*phiGas + (b+k)*phiLiquid = rhsLiquid + cLiquid*q
```

Here c is the existing residual-limited face fraction. Define the two-by-two
matrix M from the left side. Apply `inverse(M)` to both rhs and c, then assemble
the pressure equation from the alpha-weighted sum of those coupled responses.
Use those same responses to correct each phase flux after solving pressure.
There is no lagged drag RHS and no second post-pressure drag inversion or
mixture-flux restoration.

Implementation keeps the original diagonal-normalized force construction and
applies `T = inverse(M)*diag(a+k,b+k)` to it and to the diagonal pressure
response. The determinant is evaluated as `a*b+k*(a+b)`, after scaling,
avoiding subtraction of two large nearly equal products. Nonpositive a/b,
negative/asymmetric k or nonfinite coefficients stop the candidate. The
velocity boundary constraints are applied to the coupled predictor; the
existing fixedFluxPressure consistency correction follows using its coupled
mixture coefficient.

This ordering has a primary-source precedent in the Foundation
[OpenFOAM 13 face pressure corrector](https://cpp.openfoam.org/v13/facePressureCorrector_8C_source.html),
which uses a coupled inverse for force, predicted-flux and pressure responses.
This candidate is a restricted two-phase derivation for the existing OF9
equations, **not a backport of the full OF13 formulation**: it retains the OF9
drag model, residual fractions, body-force expressions and source-rate library.

The implementation only supports the exact saved quiescent inputs: two moving
pure incompressible phases, symmetric segregated drag, no virtual mass or
other interfacial momentum models, no MRF, fixed mesh and global Euler time.
The runner verifies the complete saved recipe and all 59 initial input hashes.
Do not use this executable directly on other cases or production runs.

## Isolation and acceptance

`sc1_reactor_pressure_coupling.py build` copies only the verified rate
application sources, modifies the private main and face pressure path, adds
three private headers, and runs one `wmake`. `birdPressureCoupled` links the
same six verified rate-build libraries. Source, binary, dependency and recipe
hashes are checked before use. The earlier pressure probe and its build recipe
remain unchanged and are revalidated.

`prepare` requires the exact received probe log
`d2c5d3b55ad096ee1fc8151a68353828fbfb05db9bc1e5c65af72fa0cc5ef656`
and v4 control log
`d6f4e4b1345fdbe0029598c699e9be1f322479ce77152665ff3988c05052a9ea`.
It copies only fresh t=0 inputs and mesh-preparation logs. All 59 inputs,
including fill, dt=5e-5 s, 200 steps, two outer/two pressure corrections,
gravity, surface tension, drag, boundaries and eight ranks, match v4 exactly.
The runner explicitly selects the new binary; the old controlDict application
name is retained as part of that input identity.

The executable checks unclipped pressure before every flux correction, the
mixture/phase flux identity after correction (absolute 1e-12 m3/s plus relative
1e-10), and liquid speed after every completed timestep. The existing liquid
speed limit is 1 mm/s where alphaLiquid>0.01; checking it at every step allows
earlier stopping than the saved-snapshot analysis. No speed is clamped. The
external Co/alpha/pressure-limiter/fatal-error monitor and 30-minute step cap
remain. Any numerical stop preserves the log and writes an evidence bundle;
it does not create a valid restart.

After normal completion, reconstruct all 40 snapshots (every 0.25 ms) and use
the unchanged mass, source-rate, phase-fraction, Courant, 0.01-L redistribution
and 1-mm/s quiescence limits. Require 400 pre-pressure matrix inversions,
800 valid pressure/mixture-flux checks, 200 accepted-step speed checks and
normal completion. A passing 10-ms startup remains a short quiescent screen.
Longer rest, controlled stirring, aeration and interface/wall checks must follow
before an operating checkpoint can be accepted.

## Historical Kestrel command — completed, do not repeat

This command produced the build and stopped run above. Its bundle is now
verified; use the separate spatial diagnostic instead of repeating it.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
RATE_BUILD=/scratch/nwimer/sc1_two_field_rates_build_20261111T021516Z
FAILED_COUPLED=/scratch/nwimer/sc1_reactor_coupled_20261111T165652Z
PRESSURE_PROBE=/scratch/nwimer/sc1_pressure_probe_20261111T172942Z
SC1_COUPLING_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
PRESSURE_BUILD="/scratch/nwimer/sc1_pressure_coupled_build_${SC1_COUPLING_STAMP}"
PRESSURE_CASE="/scratch/nwimer/sc1_pressure_coupled_${SC1_COUPLING_STAMP}"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_coupling.py" \
        build --base-build "$RATE_BUILD" --out "$PRESSURE_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_coupling.py" \
        prepare --failed-root "$FAILED_COUPLED" --probe-root "$PRESSURE_PROBE" \
        --candidate-build "$PRESSURE_BUILD" --out "$PRESSURE_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_pressure_coupling.py" \
            run --root "$PRESSURE_CASE" --run
        echo "Pressure-coupled startup status: $?"
    else
        echo "Build/preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

The output and verified bundle are recorded above; preserve all roots.

## Local validation

All 122 focused tests pass (55.440 s), including eight new tests. The actual
scalar C++ inverse is compiled and compared with independent 100-digit Gaussian
elimination for zero, moderate and very strong drag up to 1e25, unequal phase
inertias and invalid inputs. A separate three-equation solution checks both
phase momenta and prescribed mixture flux for absent, near-bound and mixed
phases. Other tests cover pressure/predictor consistency, private application
builds, unchanged controls, rejected tampering/reuse, failure collection without
reconstruction, and complete acceptance records without overriding an existing
failed startup screen.

The local environment has no OpenFOAM installation. These tests do not compile
the OpenFOAM field adapter or execute MPI/CFD. The subsequent Kestrel compilation passed and the actual-mesh run stopped,
as recorded above. This documentation-only update does not rerun the 122
prelaunch tests or change the candidate. No working-simulation acceptance
follows from either local tests or this incomplete run.
