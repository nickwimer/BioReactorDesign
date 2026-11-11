# SC1 VOF hydrostatic pressure perturbation candidate

## Received startup result

The solver compiled and completed 10 ms at `f0c220f` in
`/scratch/nwimer/sc1_vof_hydro_20261111T220400Z`. All physical screens pass:
peak liquid/gas speeds are 8.85673136e-8 / 1.95089614e-5 m/s and the fill remains
1000 L. Full acceptance is false only because 203/800 nonorthogonal correction
groups need further convergence. The attached archive and log budgets verify.
See the [full result](sc1_reactor_vof_hydrostatic_result_2026-11-11.json).

The subsequent [24-pass comparison](SC1_REACTOR_VOF_HYDROSTATIC_CONVERGENCE_2026-11-11.md)
now has a verified archive/log PASS using the same binary and limits. The correct
replacement archive is received; the earlier duplicate is historical. See the
[current handoff](SC1_CURRENT_HANDOFF_2026-11-11.md) for the remaining work.
The build/run command below is historical; do not repeat either completed case.

## Decision from the received force measurements

The two-step force probe completed at `f951b53` in
`/scratch/nwimer/sc1_vof_force_20261111T205814Z` using the private
`sc1InterFoamForces` binary. Its archive SHA-256 is
`ae87016bfac255e8fcd6ec77286452ae46204af0037e908c22da278e37961eae`.
All sixteen bundled file digests and ten included initial-input digests were
verified locally. The full [force result](sc1_reactor_vof_force_result_2026-11-11.json)
was reproduced by independently reconstructing the eight samples / sixteen peak
cells from the attached face records. The maximum component closure error is
7.13e-15 m/s; velocity-update closure is below 1.23e-18 m/s.

At the first correction (50 microseconds), the predictor is zero:

| Selected peak cell | Net increment, m/s | Gravity magnitude, m/s | Capillary magnitude, m/s | Fresh/matrix pressure difference, m/s |
| --- | ---: | ---: | ---: | ---: |
| Gas, y=1.031215754 m | 0.005836257941 | 9.904240348 | 5.55103e-8 | 8.49296e-15 |
| Liquid, y=1.006107065 m | 5.983119903e-6 | 0.009973016386 | 5.94980e-11 | 3.16774e-17 |

These are reconstructed velocity contributions including rAU, not independent
phase velocities or accelerations. Pressure opposes almost all of the large
gravity contribution; the remaining vector is predominantly horizontal.
At one nearly same-height gas-neighbour face, corrected density snGrad is
-1877.01348 kg/m4, gravity is -18988.30369 Pa/m and pressure +18872.70171 Pa/m.
The residual is -115.60095 Pa/m. This directly identifies a gravity/pressure
imbalance at the sampled peaks. It does not prove an exhaustive cause at every
mesh location. Capillarity there is five orders of magnitude smaller than the
net increment. A frozen-pressure ghf-to-interpolated-gh diagnostic changes the
gas increment by only 5.14340e-7 m/s, so that simple substitution is unsupported.

The Kestrel script compared all 83 audit records with the nonorthogonal control's
first two steps (relative 1e-9 / absolute 1e-12 tolerance). That parent log is not
in this attachment; local verification reproduced the force result using the
attached log but cannot independently repeat the external-prefix comparison.
The candidate preparation revalidates both actual logs and the complete original
build/input chain on Kestrel. No existing result is rerun or overwritten.

## Formulation implemented

Let rhoH(y), pH(y) be the fixed analytic density and physical hydrostatic pressure
of the existing 1000-L flat fill. They satisfy dpH/dy = -9.81 rhoH. Write

```text
p = pH + pPert
rho*DU/Dt = -grad(pPert) + (rho-rhoH)*g + viscous + capillary terms.
```

This continuum rearrangement retains the full evolving density in inertia and
the finite density departure in buoyancy. It follows the reference-state
pressure/density decomposition explained by Durran, *Numerical Methods for Fluid
Dynamics*, second edition, section 1.2.2, equation 1.57
([source, PDF page 39](https://www.coaps.fsu.edu/pub/eric/OCP5930/Papers/Numerical%20Methods%20for%20GFD%20by%20Durran.pdf)).
We apply that identity to this incompressible VOF candidate; that source does
not validate this implementation, mesh or application.

The private application `sc1InterFoamHydrostatic` uses the normal face force
`interpolate(rho-rhoH)*(g & Sf)/magSf` in both momentum predictor and pressure
flux assembly, with the same rAUf in the pressure equation and velocity update.
It solves for `pPert` with fixed zero vent value and flux-consistent wall pressure
conditions. At the reference state the density departure is exactly zero,
including its original density boundary extension. The reference is initialized
once and is never recomputed from an evolving density or pressure. No velocity,
alpha, physical force or curvature is reset to impose rest.

The original density-based inertia, alpha transport, surface tension, viscosity,
mesh, fill, timestep, ranks and ten nonorthogonal passes remain. Pressure reports
and writes are physical: `p=pH+pPert`, `p_rgh=p-rho*gh`. Physical outlet pressure
matches the prior vent condition in this unaerated gas-filled outlet. The analytic
pH boundary values are used on physical patches; processor boundary values retain
the coupled field representation. The pressure variable's own wall BC is derived
from its new flux equation.

This is a different discretization. Its exact equilibrium property does not
establish dynamic accuracy, energy conservation, mesh independence or resolved
bubble/transfer validity. Density perturbations still produce buoyancy by the
formula above; a physical perturbation/flow test is required after the resting
baseline succeeds. A later operating case must retain an explicit reference,
revisit its outlet conditions if liquid reaches the vent, and support restarts.
This fresh-start candidate deliberately refuses restart or moving-mesh use.

## Bounds, inputs and acceptance

`sc1_reactor_vof_hydro.py` builds from the original verified VOF application,
not from the two-step probe binary. No shared library or original executable is
changed. A new case copies the nonorthogonal control's original t=0 inputs.
The preparation independently checks the analytic alpha/pressure profiles,
zero U and exactly 1 m3 of liquid. C++ repeats the reference-profile checks.
Only these case entries are changed or added:

- Application name in controlDict, retaining 10 ms / 200 steps and 0.25-ms writes.
- Pressure solver entries use `pPert`; its gradient uses the same least-squares setting.
- `0/pPert`, initially zero, with fixed vent / fixedFluxPressure walls.
- `constant/sc1HydrostaticReference`, recording the existing interface height and width.

Physical initial fields and mesh retain the control's hashes. The solver has a
10-ms hard bound and the monitored eight-rank job has a 30-minute wall-time cap.
All previous liquid-speed, boundedness, redistribution, pressure and phase-mass
budget guards remain. Acceptance additionally requires gas speed <=1 mm/s,
all 800 final pressure passes converged under the existing criterion, normal End,
all 40 native snapshots, both phase speeds/inventories matching the log, and
physical pressure reconstructed consistently from pH, pPert and p_rgh.
A completed failed screen remains evidence. Automatic collection runs on success
or failure and bundles the log, small inputs and result when available.

## Historical Kestrel command

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
VOF_BUILD=/scratch/nwimer/sc1_vof_build_20261111T191346Z
FORCE_CONTROL=/scratch/nwimer/sc1_vof_force_20261111T205814Z
SC1_HYDRO_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
HYDRO_BUILD="/scratch/nwimer/sc1_vof_hydro_build_${SC1_HYDRO_STAMP}"
HYDRO_CASE="/scratch/nwimer/sc1_vof_hydro_${SC1_HYDRO_STAMP}"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_hydro.py" \
        build --base-build "$VOF_BUILD" --out "$HYDRO_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_hydro.py" \
        prepare --force-root "$FORCE_CONTROL" --build-root "$HYDRO_BUILD" --out "$HYDRO_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_hydro.py" \
            run --root "$HYDRO_CASE" --run
        echo "VOF hydrostatic startup status: $?"
    else
        echo "Build/preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

Paste the output and attach the printed BUNDLE. The intended next milestone is
an accepted 10-ms resting baseline, followed by a longer rest/perturbation test
and controlled flow activation. No stirred/aerated simulation or operating
checkpoint is accepted yet. The four-hour biological update schedule is unchanged.

## Local validation

All 59 focused tests pass. The test set covers the existing workflows plus the new profile,
pressure transformation, private build, reference checks, fresh preparation,
physical pressure output, missing fields, failed runs and evidence collection.
The actual profile helper compiles under C++14; its derivative matches -rhoH*g.
Independent finite-amplitude pressure/density identities and rotated-face
reconstruction tests retain nonzero buoyancy. These do not substitute for an
OpenFOAM compilation or CFD perturbation test. OpenFOAM is unavailable locally;
the private build and ten-pass run have now completed on Kestrel. See the received result above.
