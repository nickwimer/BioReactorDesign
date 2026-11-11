# SC1 VOF scalar-gradient comparison

## Evidence and candidate

The [verified spatial archive](SC1_REACTOR_VOF_SAVED_FIELDS_2026-11-11.md)
locates almost circumferential motion in nearly pure-phase cells immediately
above and below the initial transition layer. At 8.25 ms, bulk-gas peak speed
is 0.955494 m/s, versus 0.000115378 m/s at the interface. The peak cells do not
own boundary faces. Their neighbours shift laterally with height. Pressure
remains finite and phase inventories reconcile, so conservation alone does not
resolve this failed rest state.

The existing default is `Gauss linear`. In Foundation OF9,
[correctedSnGrad](https://cpp.openfoam.org/v9/correctedSnGrad_8C_source.html)
obtains its nonorthogonal correction from the selected `grad(field)` scheme;
the [scalar specialization](https://cpp.openfoam.org/v9/correctedSnGrads_8C_source.html)
uses that same path directly. The pressure equation and density-gradient force
therefore depend on these gradients. The
[interface normal](https://cpp.openfoam.org/v9/interfaceProperties_8C_source.html)
uses a separate named `nHat` gradient of alpha. Merely overriding
`grad(alpha.liquid)` would not change that named lookup.

The candidate selects the built-in
[least-squares gradient](https://cpp.openfoam.org/v9/leastSquaresGrad_8C_source.html)
for these scalar paths. This tests the hypothesis that gradient reconstruction
contributes to the unforced circumferential motion on the displaced cell stencil.
The stored cell-centre values are not a face-force decomposition and do not prove
that hypothesis. Because the pressure/density and interface paths change together,
even an improvement would not uniquely attribute the original error to gravity
or surface tension.

## Exact scope

Only `system/fvSchemes` differs from the control's initial inputs:

```text
gradSchemes {
    default Gauss linear;
    grad(p_rgh) leastSquares;
    grad(rho) leastSquares;
    grad(alpha.liquid) leastSquares;
    nHat leastSquares;
}
```

The default and velocity gradients retain their existing method. The mesh,
1000-L fill, pressure initialization, boundaries, viscosity, gravity, surface
tension, advection, compression, pressure iteration settings and corrected
surface-normal/laplacian schemes remain identical. No MRF, aeration or source
is introduced. The same verified `sc1InterFoamReference` application is reused;
there is no compilation. The candidate starts from a fresh copy of t=0 only.

Preparation revalidates the original control/build chain, exact received log and
manifest, the exact spatial report, and all 185 input-hash entries in that report.
After copying, the input-hash comparison must differ only at fvSchemes. The
candidate has its own manifest, logs and output. Attempted roots are never reused.

The existing limits remain: phase budgets, pressure, alpha/complement, liquid
speed <=1 mm/s, Co <=0.2 and redistribution <=0.01 L. Fixed dt remains 5e-5 s,
end time 10 ms, eight ranks, two outer/two pressure corrections, writes every
five steps and a 30-minute solver cap.

**Additional post-run requirement:** maximum gas-present speed over every logged
step must also be <=1 mm/s. This makes a good liquid result insufficient if the
headspace continues accelerating. It is checked after a normal solver exit;
the unchanged binary does not add an early gas-speed stop. A violation yields
`gradient_pass: false` and status 2. Original guards still stop the solver early.

After normal completion, require all 200 budget records, 800 pressure records
and every requested processor field before reconstruction. Both phase-present
native velocity maxima must match the audit. Compare the candidate with the
control at the same 0.25, 1, 4 and 8.25-ms snapshots. Automatic evidence collection
also runs after failure; failed/incomplete solver runs are not reconstructed.

## Next Kestrel command

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
VOF_CONTROL=/scratch/nwimer/sc1_vof_reference_20261111T191346Z
VOF_SPATIAL=/scratch/nwimer/sc1_vof_spatial_20261111T195107Z
GRADIENT_CASE="/scratch/nwimer/sc1_vof_gradients_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_gradients.py" \
        prepare --control-root "$VOF_CONTROL" --spatial-root "$VOF_SPATIAL" --out "$GRADIENT_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_gradients.py" \
            run --root "$GRADIENT_CASE" --run
        echo "VOF gradient comparison status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed `BUNDLE`. Preserve the original reference,
spatial inspection, every earlier BiRD control and both VOF build roots.

## Validation and interpretation

All 36 focused tests pass across the gradient candidate, saved-field inspector,
VOF reference and failure collector. They cover the four intended overrides,
byte-identical remaining inputs, fresh t=0 copying, control/report/native tampering,
attempt reuse, failure collection without reconstruction, all-output gates,
native/log gas agreement and rejecting excessive gas speed even when the original
liquid screen passes. Original build/case recipes are unchanged. No local
OpenFOAM runtime is available; actual CFD performance remains unmeasured.

`VOF_GRADIENTS=PASS` establishes only this 10-ms rest screen, including gas.
Longer-rest and timestep sensitivity are required before controlled stirring.
The candidate does not validate aeration, dispersed-bubble slip, PBM, oxygen
transfer or biology and does not alter the four-hour biological update schedule.
