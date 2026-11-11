# SC1 single-velocity free-surface reference

## Verified Kestrel result: liquid-speed guard stop at 8.4 ms

The corrected application compiled and passed dependency checks at `a8349f1`
in `/scratch/nwimer/sc1_vof_build_20261111T191346Z`; reported binary SHA256
`4614cd390d5586be793e260516e4dd1a0ae2facaeb59bae62a48a13a97ce3091`.
The case is `/scratch/nwimer/sc1_vof_reference_20261111T191346Z`.

**The evidence bundle is received and verified:** archive SHA256, all 16 bundled
file hashes and all 10 included initial dictionaries match. The HPC collector
reports all 49 initial inputs unchanged; the mesh/native fields and external
binaries are not included and have not been independently revalidated locally.
See the [verified result](sc1_reactor_vof_result_2026-11-11.json).

The fatal explanation is `VOF reference failed unchanged startup screens`,
from `sc1VofAudit.H:128`. Its only logged guard violation is liquid speed:
**1.003187 mm/s at 8.4 ms**, just above the unchanged 1 mm/s limit. Record 167
at 8.35 ms is 0.997218 mm/s. The guard precedes FOAM fatal/MPI_Abort and
segmentation-fault shutdown records. Return 139 does not establish an earlier
independent numerical crash. There is no End or completed 10-ms reference.

| Measured quantity | Observed prefix through 8.4 ms |
|---|---:|
| VOF step / pressure records | 168 / 672 |
| Liquid peak speed | 0.0010031870046817379 m/s |
| Gas peak speed | 0.9727328353859714 m/s |
| Pressure minimum / maximum | 99982.986897 / 111878.150896 Pa |
| Maximum logged step Co | 0.001788180146 |
| Liquid inventory | 999.999999999982 L |
| Liquid mass-residual L1 / limit | 9.29257e-11 / 1.00001e-5 kg |
| Gas mass-residual L1 / limit | 1.29325e-13 / 6.52002e-9 kg |
| Alpha redistribution L1 | 3.87786e-8 m3 (3.87786e-5 L) |

Both phase budgets independently reconcile from volume changes and the actual
transport flux; printed per-step and cumulative residuals agree. Fractions,
complement, pressure, redistribution and Co remain within their stated limits.
The final step is printed before its rejecting guard: 168 records are not
168 accepted steps and conservation does not establish quiescence.
The old collector's zero pair/alpha values were BiRD-only parser limitations;
the new VOF-aware report reads the actual audit and stock fraction lines.

Gas speed grows from 0.005837 m/s at 0.05 ms to 0.116583 at 1 ms, 0.465039 at
4 ms and 0.972733 at 8.4 ms. It is not approaching rest. Global liquid and gas
peaks need not be colocated, and their ratio is not a local force budget.
The [stock OF9 pressure correction](https://cpp.openfoam.org/v9/multiphase_2interFoam_2pEqn_8H_source.html)
combines surface tension and gravity with the pressure-equation flux in the
velocity update. The present log does not identify which spatial term creates
the observed motion; phase slip has been removed but the reference is not
quiescent. Do not relax the speed guard or extend the failed trajectory.

**Next action: [inspect the already-written native processor fields](SC1_REACTOR_VOF_SAVED_FIELDS_2026-11-11.md).**
Directory listings exist through 8.25 ms. The new read-only Python command checks
the full global-cell mapping, initial mapped fields, and selected snapshots
against their logged inventories and velocities before reporting peak locations,
directions, patch adjacency and immediate-neighbour data. No OpenFOAM executable,
reconstruction, compilation or CFD is required. The source remains untouched.

Verified archive SHA256:
`ab4e882c60e137ca99f4a646cf288c3324c2dfb6dd9348f3186f9fde3e785cb8`.
Verified solver-log SHA256:
`2811c733dc1fe00499022c7d1b447d26ff6ec45facf145261f5f5dfae9190a9b`.
Preserve every build, case and control. The operating simulation is not yet accepted.

## First build result and correction

The first Kestrel build at `9142d6a` stopped during compilation in
`/scratch/nwimer/sc1_vof_build_20261111T190909Z`. The inspected native-source
checks passed and `wmake` started, but the observer called
`Info.precision(17)`. OF9's `Info` is a `messageStream`; its explicit
`operator()()` returns the output stream that provides `precision()`.
The correction is `Info().precision(17)`, consistent with the previously
compiled SC1 audit headers and the
[OF9 messageStream API](https://cpp.openfoam.org/v9/classFoam_1_1messageStream.html).

This is a logging API correction only. No equation, case dictionary, guard,
precision value or budget definition changes. Preparation and CFD were not
reached. The subsequent fresh build at `a8349f1` passed, as recorded above.
Preserve this failed build and `log.build.solver`; its recorded header recipe
must not be rewritten. The 33 focused tests passed at the logging correction.

## Decision from the verified spatial evidence

The [velocity bundle](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md) is now verified:
archive SHA256, all 19 bundled hashes, 13 included initial dictionaries and
all 33 numerical prefix records match. The complete native-field summary adds
a disturbance missing from the compact console: gas moves mainly horizontally
near and immediately above the interface. At 0.15 ms its interface peak is
40.0167 mm/s and bulk-gas peak is 34.8828 mm/s. Those peak cells do not own
boundary faces. The liquid interface peak remains 0.800102 mm/s downward,
versus 0.070844 mm/s in bulk liquid. A vertical-counterflow-only explanation
is therefore incomplete. Pressure collapse remains absent during the observed
three-step interval, but neither the gas nor the interface is stationary.

The selected comparison is a separate single-velocity resolved-interface reference
on the actual mesh. It uses the Foundation OF9
[interFoam equations](https://cpp.openfoam.org/v9/interFoam_8C_source.html),
with read-only observation and failure guards. The pressure and momentum
equations, surface-tension calculation and MULES transport are unchanged.
This removes independent phase-velocity slip as a degree of freedom, providing
a comparison for the mesh, initialization and pressure/gravity/interface balance.
It does not prove the two-fluid formulation is generally invalid, nor identify
which term caused the horizontal gas velocities.

This is a path toward a stable free-surface hydrodynamic baseline. It is not
an aerated reactor replacement: dispersed-bubble slip, PBM, gas-liquid mass
transfer, dissolved oxygen and biology are absent. All BiRD source-rate tests,
private binaries, reactor controls and the four-hour biological update schedule
remain preserved. The eventual phase-count decision remains based on measured
physical adequacy; this reference does not require or rule out a third field.

## Matched quantities and explicit differences

| Quantity | Reference |
|---|---|
| Mesh | Unchanged 281664-cell native MEML mesh and patches |
| Initial liquid fraction, p_rgh, gravity | Byte-identical copies from the verified spatial control |
| Fill | Same measured 1000 L; same initial transition layer |
| Density | Liquid 1000, gas 1 kg/m3 |
| Viscosity | Same liquid law at 303.15 K; gas dynamic viscosity 1.8e-5 Pa s |
| Surface tension | Same 0.07 N/m |
| Walls/propeller | Stationary, no-slip; same empty patches |
| Vent | Same fixed p_rgh=100000 Pa and pressureInletOutletVelocity |
| Time/ranks | Same fixed 5e-5 s timestep, eight ranks, at most 200 steps/10 ms |
| Momentum | One mixture velocity and pressure, replacing two phase velocities/drag |
| Alpha method | Stock OF9 semi-implicit MULES, one correction and one subcycle; no previous compression correction |
| Pressure iteration | Two outer/two pressure corrections, no nonorthogonal correctors, 1e-12 pressure tolerance |
| Output | Every five steps, forty snapshots from 0.25 to 10 ms |

The VOF transport scheme is `Gauss interfaceCompression vanLeer 1`, as in the
[OF9 tutorial](https://github.com/OpenFOAM/OpenFOAM-9/blob/master/tutorials/multiphase/interFoam/laminar/damBreak/damBreak/system/fvSchemes).
This reference changes the formulation and its associated dictionaries; it is
not a one-coefficient ablation of BiRD. The preserved mesh, fill, physical
properties and initial pressure make its rest behavior a useful comparison.
The hydrostatic pressure is not claimed to be discretely balanced on this mesh.

## Private build and measured budgets

`sc1_reactor_vof.py build` verifies sixteen installed OF9 source files against
the inspected upstream source fixtures, copies them into a fresh private build,
and adds only observer calls and `sc1VofAudit.H`. The observer never writes to
alpha, U, p, a source or a flux. The Make target is private; installed OpenFOAM
and all BiRD libraries are untouched. One application `wmake` is required.
If the installed sources differ, preparation stops before compilation rather
than silently adapting unreviewed source. The executable and resolved library
hashes are recorded, and the reference rejects BiRD library dependencies.

The case preparation revalidates the existing velocity control and its fifteen
saved field hashes, the prior source/build/control chain and the new executable.
Only initial fields and mesh are copied; no failed or completed time is used as
a restart. Fresh output paths and unchanged inputs are checked before launching.

For global Euler, a single alpha subcycle and no phase change/source, the
observer records the actual limited `alphaPhi10` transport flux immediately
after each alpha solve. It retains only the final outer iteration for each
timestep budget. Liquid boundary flux is alphaPhi10; gas boundary flux is
the complement using the same pre-pressure mixture phi employed by that alpha
solve. Processor copies and empty patches are excluded. The post-pressure
mixture boundary flux is logged separately, so pressure correction is not
silently substituted for the flux that transported alpha.

Both phase masses are reconciled independently in Python from volume changes
and integrated outward transport fluxes. No clipping, inventory correction or
source compensation is added. The observer applies the existing startup limits:

- Pressure at every solve, before velocity correction: finite and >=1000 Pa.
- Fractions/complement within 1e-10; nonfinite fields stop immediately.
- Accepted-step mixture Courant <=0.2.
- Liquid speed where alphaLiquid>0.01 <=1 mm/s.
- Alpha redistribution L1 <=1e-5 m3 (0.01 L).
- Phase mass-residual L1 <=1e-10 kg + 1e-8 times initial phase mass.

Gas speed is reported at every step and in native snapshots for review; it is
not hidden behind liquid acceptance. The external failure monitor and a
30-minute solver cap remain active. A solver failure is collected without
reconstruction. After normal completion, require 200 step budgets, 800 pressure
records and all requested processor outputs, then reconstruct forty snapshots.
Native inventory and liquid speed must match the logged accepted state.

`VOF_REFERENCE=PASS` means this reference met the stated 10-ms screens. It
does not accept the failed BiRD startup or an operating checkpoint. Review gas
motion and spatial pressure results too. If the reference passes, extend rest
and establish timestep sensitivity before controlled stirring. If it fails,
the next investigation is the shared mesh/initialization/discrete force balance
within this simpler reference. Neither outcome alone selects an aeration/PBM
closure, changes biological scheduling or authorizes a production restart.

## Historical Kestrel command — already executed; do not repeat

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
CONTROL=/scratch/nwimer/sc1_velocity_probe_20261111T183832Z
SC1_VOF_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
VOF_BUILD="/scratch/nwimer/sc1_vof_build_${SC1_VOF_STAMP}"
VOF_CASE="/scratch/nwimer/sc1_vof_reference_${SC1_VOF_STAMP}"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof.py" \
        build --out "$VOF_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof.py" \
        prepare --control-root "$CONTROL" --build-root "$VOF_BUILD" --out "$VOF_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof.py" \
            run --root "$VOF_CASE" --run
        echo "VOF reference status: $?"
    else
        echo "Build/preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

The command above produced the stopped attempt recorded at the top. Its bundle
has been verified. Use the separate saved-field inspection command; preserve
every root and do not repeat the calculation.

## Local verification and remaining runtime work

The earlier 39-test suite passed for VOF-aware failure collection. The current
29 focused tests cover the reference, failure collector and saved-field reader. They check
unchanged upstream equation lines, exact source compatibility, private build
destinations, preserved controls, runtime/source tampering, boundary recipes,
known phase export/replacement budgets, rejected incomplete/nonfinite records,
native inventory/velocity checks and reconstruction only after normal solver
completion. Prior velocity, failure-collection and pressure-coupling tests pass.
The source fixture preserves upstream license notices. JSON, documentation
links, shell syntax and diff checks are also validated before push.

There is no local OpenFOAM installation. Kestrel compilation passed, but the
reference stopped before 10 ms. The new evidence tests cover individual existing
guards, exact limits, missing/duplicate initial inventories, malformed/nonfinite
and duplicate step records, stock fraction output and unchanged-case collection.
The compiled adapter and its runtime inputs are unchanged by this update.
No pass, gas equilibrium or working full-reactor simulation is claimed.
