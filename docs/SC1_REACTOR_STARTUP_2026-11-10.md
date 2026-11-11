# SC1 fresh actual-mesh quiescent startup

## Status and decision

The [rate-split private build and all eight numerical regressions](SC1_TWO_FIELD_RATE_SPLIT_TEST_2026-11-10.md)
passed on Kestrel at `d83bc83`. This next step runs that verified executable
on the **281,664-cell MEML vessel mesh**, with fresh, nonreacting fields.
Preparation/report/workflow checks pass locally; this actual-mesh CFD has not
yet run. No new C++ build or repeat of the eight regressions is required.

Start with an unaerated, stationary vessel. Its expected behavior is a nearly
stationary, flat liquid interface. This separates hydrostatic initialization,
vent treatment, interface numerics and mesh behavior from impeller-driven flow.
The first interval is only **10 ms: 200 fixed steps at 5e-5 s**. It is an
initial startup screen, not an equilibrium test or the 9.9-s operating checkpoint.

The completed wall preflight supplies the exact boundary addressing and
geometry. It is not rerun. The preparer copies only `constant/polyMesh` from
the archive, verifies its four previously recorded addressing hashes, records
hashes of every mesh file, and checks the new copy with `checkMesh`. No archived
fields, dictionaries, sources or schedules are copied or modified.

## Explicit physical and numerical contract

| Item | Startup setting and reason |
|---|---|
| Geometry | Archived MEML vessel, shaft and impellers, including its existing rotor cell zone; no remeshing |
| Motion | Both `walls` and `propeller` stationary/no-slip; no MRF dictionary in this case |
| Phases | Pure isothermal `gas` and `liquid`, gas reference, densities 1 and 1000 kg/m3; existing runtime guards retained |
| Initial volume | Exactly 1 m3 liquid from native OpenFOAM cell volumes; report actual mass and total gas inventory |
| Interface | Planar linear alpha transition 0.0125 m wide, midpoint solved to match 1 m3; reject midpoint outside 0.99–1.05 m |
| Pressure | Analytic hydrostatic profile integrated from that same continuous alpha/density profile; reconstructed discrete force balance is tested, not assumed |
| Gravity / vent | `(0 -9.81 0)` m/s2; vent `p_rgh=100000` Pa, both velocities `pressureInletOutletVelocity`; incoming gas fraction 1 and liquid 0 |
| Interface forces | Segregated drag `m=0.5, n=8` and surface tension 0.07 N/m from the MEML headspace/liquid setup; compression coefficient 1 |
| Transport | Laminar; 303.15 K; liquid viscosity from the archived water law, gas viscosity 1.8e-5 Pa s |
| Sources / biology | Empty fvModels, no interphase transfer, PBM or biological model |
| Numerics | Euler, upwind, two outer and two pressure correctors, six alpha subcycles, one alpha correction, audit enabled |
| Execution | Default eight MPI ranks through `srun`; separate `--run` required; serial supported with `--ranks 1` |

The gas density/viscosity are a constant-property diagnostic approximation,
not the archived compressible multicomponent gas. Gas reference avoids the
rescale branch; both reference variants passed the small regressions. No
claim is made that reference reconstruction alone guarantees conservation.
The segregated model represents this resolved headspace interface; it does
not model sparged bubbles. No contact-angle parameter is invented.

The fill is generated from `postProcess -func writeCellCentres` and
`writeCellVolumes` at t=0 in the fresh case. These are Foundation OF9 native
[cell centres](https://cpp.openfoam.org/v9/writeCellCentres_8C_source.html)
and [cell volumes](https://cpp.openfoam.org/v9/writeCellVolumes_8C_source.html),
not VTK centre/volume approximations. The preparer solves for interface height
by volume-weighted bisection and writes complementary fractions with 17 digits.

For `x=y-h`, transition width `w`, and `a=clip(1/2-x/w,0,1)`, define
`F(x)=integral_0^x a(s) ds`: `x+w/8` below the layer,
`x/2-x*x/(2*w)` inside it, and `w/8` above it. Initial pressure is
`p=100000-g*y+999*g*(w/8-F(x))`, and
`p_rgh=p+(1+999*a)*g*y`. Thus the upper gas region matches the vent value,
and `dp/dy=-rho*g`. The absolute pressure gauge is harmless for the guarded
incompressible phases. This continuum balance does not guarantee exact
discrete balance across a nonuniform mesh; the startup measures the difference.

## Acceptance and evidence

The full accepted-step gas/liquid audit and per-alpha-call decomposition are
required, with exactly 200 steps, both outer calls, all expected terms and the
candidate library/rate markers. Sources must be zero; fractions must remain
bounded/complementary. Initial inventory is checked against the independent
native-volume initialization. Reconstructed final masses must match the log.

The reactor mass threshold is defined before running as
`1e-10 kg + 1e-8*initialPhaseMass + 1e-3*abs(cumulativeBoundaryMass)` per phase.
This separate larger-domain contract does not modify any completed small-case
threshold. The L1 sum of accepted-step physical residuals must be below it;
source splitting and equation-ledger identities are checked independently.

Every 40 steps, save and reconstruct the fields. Report phase-weighted RMS and
maximum velocity in cells where that phase exceeds 1%, liquid volume, and
volume-integrated absolute alpha redistribution. The startup screen also
requires:

- Maximum logged mixture Courant number <=0.2.
- Liquid maximum present-phase speed <=1 mm/s at all five saved snapshots.
- Absolute alpha redistribution <=1e-5 m3 (0.01 L) at every snapshot.

These are prospective startup criteria, not experimental validation limits.
The field checks sample five snapshots; they do not establish continuous-time
velocity maxima or long-time equilibrium. A complete failed startup is saved
as evidence, returns status 2, and does not proceed automatically to rotation.
Setup/solver/report errors return 1. A startup pass returns 0. Every attempted
case is protected against rerun/overwrite. Read-only `analyze --root ...`
can use saved logs and reconstructed outputs without launching CFD.

The preparer and launcher revalidate all eight existing rate regressions and
their original controls, using their saved manifests and logs. The verified
build's original/private sources, libraries and executable must still match.
The archive, initialized input hashes and preparation recipe are checked again
before launch. All commands operate in the newly created scratch root.

## Kestrel command

Use the active allocation. This creates one fresh case, checks/initializes its
copied mesh, then launches the short startup with eight ranks. Keep the printed
`REACTOR_STARTUP` root and all logs even if a check fails. No `blockMesh`,
archive restart, new compilation, or eight-case CFD rerun occurs.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
RATE_CASES=/scratch/nwimer/sc1_two_field_rates_cases_20261111T021516Z
STARTUP="/scratch/nwimer/sc1_reactor_startup_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
        prepare --regression-root "$RATE_CASES" --out "$STARTUP" --ranks 8; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
            run --root "$STARTUP" --run
        echo "Reactor startup status: $?"
    else
        echo "Startup preparation stopped. Preserve the new directory and paste the output."
    fi
else
    echo "Startup setup stopped. Preserve existing files and paste the output."
fi
```

## Validation and following stage

All 77 focused checks pass locally: the 68 existing regressions plus nine new
startup tests. Local tests cover exact native-volume fill, pressure-gradient consistency,
phase complement/bounds, malformed/nonfinite field rejection, full log coverage,
zero-source budgets, saved-field/log agreement, quiescent-motion failures,
fresh-case protection, mesh preservation/change detection, MPI/reconstruction
command order and allocation/explicit-run gates. OpenFOAM commands are mocked
locally; actual mesh/solver execution requires Kestrel. The existing 68 tests
retain their original report digests and numerical thresholds.

After reviewing the startup, extend the quiescent interval to establish a
meaningful stationary-surface baseline. Then introduce the documented impeller
ramp on the same geometry/fill and check wall/azimuthal surface behavior and
mixing. Aeration requires a separately justified dispersed-versus-resolved
gas closure. Only after those checks should the nonreacting 9.9-s checkpoint,
oxygen transfer/PBM and the unchanged four-hour biological schedule be enabled.
Three fields remain optional; their physical need has not been established.
