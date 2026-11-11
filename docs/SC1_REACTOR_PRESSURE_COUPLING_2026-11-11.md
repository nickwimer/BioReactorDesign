# SC1 pressure-coupled face-drag candidate

## Status

Implemented and locally tested; **Kestrel compilation and CFD result pending**.
The next action is the command below. It builds only a new private application,
then runs a fresh quiescent case for at most 10 ms. Preserve all previous builds,
cases, the archived reactor and the original solver sources. No operating,
free-surface, oxygen or biology acceptance is claimed.

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

## Next Kestrel command

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

Paste the output and attach the printed `BUNDLE`. No old case is restarted.

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
the OpenFOAM field adapter or execute MPI/CFD. The command above is the pending
compilation and actual-mesh test; no claim of a working simulation follows
from the scalar or orchestration tests alone.
