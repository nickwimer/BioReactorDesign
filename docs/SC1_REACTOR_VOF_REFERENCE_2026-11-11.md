# SC1 single-velocity free-surface reference

## Reported Kestrel result: build passed, reference stopped

At `a8349f1`, the corrected application compiled and passed dependency checks
in `/scratch/nwimer/sc1_vof_build_20261111T191346Z`. Its reported SHA256 is
`4614cd390d5586be793e260516e4dd1a0ae2facaeb59bae62a48a13a97ce3091`.
The fresh case `/scratch/nwimer/sc1_vof_reference_20261111T191346Z` stopped
at printed time 0.0084 s, outer 2, with a generic `FOAM FATAL ERROR` banner.
There are 168 printed time labels, no `End`, solver return 139 and wrapper
status 1. These are not 168 accepted steps. The reported maximum printed
mixture Co is 0.001778043328599279; this does not establish pressure or
quiescent-speed acceptance. All 49 initial inputs are reported unchanged.

The [saved console result](sc1_reactor_vof_result_2026-11-11.json) records the
reported identities separately from independent verification. The full bundle
has not yet been received. The console omits the fatal explanation and VOF
step records, so the initiating failure is still unknown. Do not label this
as a liquid-speed failure, pressure failure or independent segmentation fault
from the generic banner/exit code alone. The 10-ms reference did not complete.

The older failure collector recognizes `SC1_TWO_STEP` and BiRD fraction lines,
not the VOF audit or stock `Phase-1 volume fraction` output. Its zero pair count
and zero largest-alpha value are parser limitations, not physical findings.
The read-only collector now adds VOF record coverage, last logged phase budgets,
pressure extrema and logged guard violations. Malformed/nonfinite records and
duplicates remain explicit; no incomplete log can gain startup acceptance.
A VOF step is printed before its guard. Its presence alone is not acceptance.
The new summary reports logged budgets; independent reconciliation remains
required when inspecting the full log. No solver, guard, build recipe or case
input changes in this reporting update.

**Next action: attach the existing archive; no new Kestrel run is needed.**

`/scratch/nwimer/sc1_vof_reference_20261111T191346Z_evidence/reactor_failure_evidence.tar.gz`

Reported archive SHA256:
`ab4e882c60e137ca99f4a646cf288c3324c2dfb6dd9348f3186f9fde3e785cb8`.
Reported solver-log SHA256:
`2811c733dc1fe00499022c7d1b447d26ff6ec45facf145261f5f5dfae9190a9b`.

Verify that archive and its included hashes, inspect the fatal context, and
reconcile the observed VOF budgets before changing a numerical or physical
choice. Preserve both VOF builds and every case/control. No reconstruction ran;
time-directory names alone do not establish usable native snapshots. Do not
rebuild, rerun, restart or recollect this attempt merely to use the updated parser.

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

The command above produced the stopped attempt recorded at the top. Attach its
existing bundle and preserve every root; no new allocation or calculation is needed.

## Local verification and remaining runtime work

Thirty-nine focused tests pass, including ten reference tests and six additional
incomplete-VOF evidence tests. They check
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
