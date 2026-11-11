# SC1 fresh actual-mesh quiescent startup

## Status and decision

The first **281,664-cell MEML-mesh startup completed at `e586c47`**. Its
mass and diagnostic screens pass, but its quiescent screen **fails**: liquid
speed reaches 0.0963573 m/s and mixture Co reaches 0.623991. Preserve this run
as the control; do not extend it or add rotation.

A configuration error has been identified in the prepared case: its linear
blending almost cancels the unordered segregated drag. The matched correction
at `0876843` **diverged by t=0.0009 s**, with alpha extrema around 1e142,
then an FPE/time-limit cancellation (wrapper status 1). Its root is
`/scratch/nwimer/sc1_reactor_drag_20261111T144728Z`. Preserve it and the v1
control; do not rerun the historical commands below. **The failure bundle is
received and verified. The monitored
[face-momentum comparison](SC1_REACTOR_FACE_MOMENTUM_2026-11-11.md) then failed
on its first timestep. Its bundle is now inspected. The subsequent
[coupled-drag comparison](SC1_REACTOR_COUPLED_DRAG_2026-11-11.md) also failed
at the same correction; its bundle is now verified. The
[short pressure probe](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md) has now reproduced
the pressure collapse and stopped before limiting; its bundle is now verified.
The [private pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
stopped at 0.2 ms on the liquid-speed guard; its bundle is verified.
The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md) completed; its verified bundle also reveals horizontal gas motion near the interface. The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped on the liquid-speed guard at 8.4 ms; its bundle is verified. The scalar-gradient comparison also failed; next use the [nonorthogonal pressure comparison](SC1_REACTOR_VOF_NONORTHOGONAL_2026-11-11.md). Do not rerun cases.**
The [private rate-split solver and eight
numerical regressions](SC1_TWO_FIELD_RATE_SPLIT_TEST_2026-11-10.md) remain valid
under their tested contract; no library rebuild or numerical-suite rerun is needed.

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

## Completed control result and diagnosed configuration error

Control root: `/kfs3/scratch/nwimer/sc1_reactor_startup_20261111T045638Z`.
Log SHA256: `987019c594cfdc0a856472d67aa1ce99885d6205d1e9ed56b8973623040fe33b`.
The [recorded console metrics](sc1_reactor_startup_result_2026-11-10.json)
retain the initial geometry, budgets and all five velocity/redistribution snapshots.
Raw logs, fields and manifests remain on Kestrel.

The native mesh volume is 1.6420022691408904 m3. Initialization measured
0.9999999999999999 m3 liquid and 0.6420022691408908 m3 gas, with interface
midheight 1.0185981622921294 m. Both phases stay in [0,1], and all source-split
volumes and assembled physical sources are zero.

| Metric | Observed | Original criterion |
|---|---:|---:|
| Gas physical residual L1 | 8.079028304339642e-14 kg | <=6.520022765429711e-9 kg |
| Liquid physical residual L1 | 3.344917609101449e-12 kg | <=1.00001e-5 kg |
| Final liquid inventory | 999.9999999999975 L | Independently checked against accepted log |
| Final absolute alpha redistribution | 7.916063403590663e-8 m3 | <=1e-5 m3 |
| Maximum mixture Courant number | 0.6239911029205629 | <=0.2: **FAIL** |
| Final liquid maximum / weighted RMS speed | 0.0963573 / 0.00802936 m/s | Maximum <=0.001 m/s: **FAIL** |
| Final gas maximum / weighted RMS speed | 2.320774 / 0.136666 m/s | Reported, no gas-speed acceptance threshold |

The phase speeds rise through all five snapshots, while the total liquid
inventory hardly changes. Conservation alone has not validated the interface
momentum balance.

The local source identifies a definite error in the startup dictionary:

1. `phaseSystemTemplates.C` selects a model-specific blend if present, otherwise
   `blending/default`. The control had only the default linear blend.
2. For gas/liquid fractions `a` and `1-a`, its 0.3/0.7 thresholds give
   `f1=clip((1-a-0.3)/(0.4+small),0,1)` and
   `f2=clip((a-0.3)/(0.4+small),0,1)`. Their sum is effectively one over [0,1].
3. `BlendedInterfacialModel.C` multiplies unordered drag by `1-f1-f2`.
   Thus the listed `(gas and liquid) segregated` model has effectively zero
   weight. There are no ordered dispersed drag models in this startup case.

Relevant source files are
[linear.C](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/BlendedInterfacialModel/blendingMethods/linear/linear.C),
[BlendedInterfacialModel.C](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/BlendedInterfacialModel/BlendedInterfacialModel.C),
[phaseSystemTemplates.C](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseSystem/phaseSystemTemplates.C)
and [noBlending.C](../OFsolvers/birdmultiphaseEulerFoam/phaseSystems/BlendedInterfacialModel/blendingMethods/noBlending/noBlending.C).
The actual local f1/f2 implementations are compiled in a focused test, with a
minimal scalar API double: over 10,001 complementary states the old weight is
below 1e-12 in magnitude, while `none/none` yields exactly one. This verifies
selection algebra, not the complete hydrodynamic remedy.

The correction adds only this override inside `blending`:

```text
drag
{
    type none;
    continuousPhase none;
}
```

Neither actual phase is named `none`, so both ordered weights are zero and
the unordered weight is one. The Foundation OF9
[damBreak4phase tutorial at 89839ae3](https://github.com/OpenFOAM/OpenFOAM-9/blob/89839ae3b8cd7b7c8d044f6898d4991177b1db41/tutorials/multiphase/multiphaseEulerFoam/laminar/damBreak4phase/constant/phaseProperties)
uses this blend for its segregated resolved-interface models. Here it is
applied only to `drag`, preserving the default blend and every other setting.
The control's `m=0.5`, `n=8`, pressure initialization, boundaries, timestep,
mesh, fill, ranks, compression and all acceptance limits remain unchanged.
No drag-coefficient tuning, velocity clamp or force alteration is introduced.

**Interpretation:** the canceled drag is established from the case and source.
Its contribution to the observed acceleration is a hypothesis tested by the
fresh comparison. Restoring it is not yet a demonstrated quiescent-flow fix.
The same default blend appeared in the earlier MEML headspace/liquid recipe;
that is a reason for a later model review, not proof of the historical reactor
wall/oxygen problem. Production dictionaries remain unchanged.

## Historical v2 command (attempted and failed; do not repeat)

`prepare-drag` revalidates the completed control, copies its original `0`,
`constant` and `system` plus the three mesh-preparation logs, and changes only
`constant/phaseProperties` as specified above. It copies no processor output,
later time, solver log or result. Hash checks enforce this one-file difference
and identical initial inventories. The saved v1 result remains readable and
unchanged; the new variant is `sc1_meml_quiescent_segregated_drag_v2`.

On a complete run, reporting requires the unique runtime line `Selecting drag blending method: none`.
It runs the same 200 steps and unchanged screens, then records
`drag_comparison.json` in the new root. The comparison also reads both saved
end fields to locate peak velocities, separate bulk/interface motion and
measure gas/liquid relative velocity and mixture velocities. These are new
read-only diagnostics, not changes to the old acceptance thresholds.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
CONTROL=/scratch/nwimer/sc1_reactor_startup_20261111T045638Z
DRAG_STARTUP="/scratch/nwimer/sc1_reactor_drag_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
        prepare-drag --baseline-root "$CONTROL" --out "$DRAG_STARTUP"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_startup.py" \
            run --root "$DRAG_STARTUP" --run
        echo "Drag-corrected startup status: $?"
    else
        echo "Preparation stopped. Preserve both roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

If CFD completes but reporting/comparison fails, preserve the whole root and
use `analyze --root ...` after resolving the reporting issue. It can print the
v2 report and comparison without rerunning CFD or overwriting saved results.
Do not repeat the launch block on an attempted destination.

## Original physical and numerical contract (retained in v2)

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

## Historical v1 Kestrel command (completed; do not repeat)

The command below produced the completed v1 control. It created one fresh case, checked/initialized its
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

The 68 numerical and nine initial startup tests passed before this control.
All 16 startup tests now pass, including seven new drag-correction checks. Local tests cover exact native-volume fill, pressure-gradient consistency,
phase complement/bounds, malformed/nonfinite field rejection, full log coverage,
zero-source budgets, saved-field/log agreement, quiescent-motion failures,
fresh-case protection, mesh preservation/change detection, MPI/reconstruction
command order and allocation/explicit-run gates. OpenFOAM commands are mocked
locally; actual mesh/solver execution requires Kestrel. The existing 68 tests
retain their original report digests and numerical thresholds.

After a corrected startup passes the unchanged screen, extend the quiescent interval to establish a
meaningful stationary-surface baseline. Then introduce the documented impeller
ramp on the same geometry/fill and check wall/azimuthal surface behavior and
mixing. Aeration requires a separately justified dispersed-versus-resolved
gas closure. Only after those checks should the nonreacting 9.9-s checkpoint,
oxygen transfer/PBM and the unchanged four-hour biological schedule be enabled.
Three fields remain optional; their physical need has not been established.
