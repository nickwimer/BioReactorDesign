# SC1 private two-field rate-split candidate

## Status and purpose

**Kestrel build and all eight CFD regressions passed at `d83bc83`.**
The uploaded transcript reports `RATE_REGRESSION=PASS` and status 0. Every
report is complete, with valid diagnostics, passing physical mass screens and
no errors. The original physical thresholds were retained. Preserve the
completed build and all eight cases; do not repeat the historical command below.

The first [MEML-mesh quiescent startup](SC1_REACTOR_STARTUP_2026-11-10.md)
completed with valid mass budgets but excessive velocities. Its dictionary
nearly canceled segregated drag through blending. Restoring only that weight
[diverged by 0.0009 s](SC1_REACTOR_DRAG_FAILURE_2026-11-11.md), reusing this
verified binary and its guards. The failure bundle is received and inspected;
the subsequent [face-momentum startup](SC1_REACTOR_FACE_MOMENTUM_2026-11-11.md)
also failed. Its bundle is now inspected; the
[coupled-drag comparison](SC1_REACTOR_COUPLED_DRAG_2026-11-11.md) also failed
at the same first-step outer iteration. Its bundle is now verified. The
[short pressure probe](SC1_REACTOR_PRESSURE_PROBE_2026-11-11.md) has now reproduced
the pressure collapse and stopped before limiting; its bundle is now verified.
The [private pressure-coupled candidate](SC1_REACTOR_PRESSURE_COUPLING_2026-11-11.md)
stopped at 0.2 ms on the liquid-speed guard; its bundle is verified.
The [three-step spatial diagnostic](SC1_REACTOR_VELOCITY_LOCATION_2026-11-11.md) completed; its verified bundle also reveals horizontal gas motion near the interface. The [single-velocity VOF reference](SC1_REACTOR_VOF_REFERENCE_2026-11-11.md) stopped on the liquid-speed guard at 8.4 ms; its bundle is verified. The hydrostatic startup completed 10 ms and passed all physical screens; next close its remaining [pressure-convergence gate](SC1_REACTOR_VOF_HYDROSTATIC_CONVERGENCE_2026-11-11.md) with the same binary and limits.
These regressions establish the correction for their tested source contract;
reactor wall/free-surface, aeration, oxygen and biology validation remain open.
All 68 pre-existing focused local tests passed before the Kestrel run.

The [completed control experiment](SC1_TWO_FIELD_NEAR_BOUND_TEST_2026-11-10.md)
found approximately 0.900 micrograms of liquid mass-budget deficit at both
timesteps and both reference settings. It matched the source-floor attenuation.
The completed baseline and near-bound cases and their verified build are
preserved. The new executable is `birdTwoFieldRates`, built with all six
project libraries in a new private root. No file under
`OFsolvers/birdmultiphaseEulerFoam` is edited by this implementation.

## Completed Kestrel result

- Build: `/kfs3/scratch/nwimer/sc1_two_field_rates_build_20261111T021516Z`.
- Cases: `/kfs3/scratch/nwimer/sc1_two_field_rates_cases_20261111T021516Z`.
- Executable SHA256: `34c7f5dae2739f2b1a40776dff142f9bfb4e7ebda495ca2134438202cf01bc54`.
- Phase-library SHA256: `6ae6f57fa85d8dd44b356a8de566ec991591b1e3b8ed16a34309d63facf95d8f`.

The [exact reports and comparison](sc1_two_field_rate_results_2026-11-10.json)
include all seven binary hashes, all eight candidate/control log hashes,
transcript SHA256 and source commit. These are user-supplied console results;
the raw solver logs and input/build manifests remain on Kestrel. The runtime
runner revalidated the saved original controls without rerunning CFD.

The four closed/mid-fraction cases retain their original physical results at
printed precision. The four near-bound cases now also pass:

| Candidate | Gas residual L1 (kg) | Liquid residual L1 (kg) | Liquid L1 reduction vs control |
|---|---:|---:|---:|
| `bound_coarse_no_reference` | 3.669221e-18 | 4.839702e-13 | 1,859.9× |
| `bound_coarse_gas_reference` | 2.040476e-17 | 2.007164e-14 | 44,824.4× |
| `bound_fine_no_reference` | 1.477416e-18 | 2.758261e-13 | 3,262.9× |
| `bound_fine_gas_reference` | 7.567124e-18 | 7.718720e-15 | 116,567.0× |

The near-bound liquid threshold stays `1.01e-10 kg`; gas thresholds are about
`4.98e-15` to `5.74e-15 kg`, computed by the unchanged throughput formula.
All four gas residuals are also below their thresholds. The sum of absolute
split volumes over every outer call is at most about `7.22e-24 m3`, below
`1e-18 m3`. Zero observed floor-induced splitting is the intended result.
The nonnegative-rate and complementary-fraction guards pass in every case.

Both reference choices pass. Gas reference is used for the next quiescent
startup because it avoids the rescale branch; this is a controlled implementation
choice, not evidence that reference reconstruction alone fixes conservation.
No extrapolation to signed sources, compressible gas, interphase transfer,
PBM or biological chemistry is established by these tests.

## Numerical change and limits

For two complementary phase fractions, the raw coupled dilatation for phase
i is `(1-alpha_i)*q_i - alpha_i*q_j`, where q is the stored pressure-source
rate `phase.divU()`. In the private phase-system copy, replace only its capped
Su/Sp split with:

```text
Su_increment = q_i
Sp_increment = -(q_i + q_j)
```

The original divergence source is retained, as are pressure assembly, flux
limiting, subcycles, implicit correction, reference reconstruction, rescaling,
clamping and all observers. The source is not added back after solving. The
cached option source is still evaluated once and used by the same pressure
equation. The existing source ledger measures any remaining split defect and
implicit response independently of the accepted physical mass budget.

For nonnegative frozen rates, the source-only update
`alphaNew=(alphaOld+dt*q_i)/(1+dt*(q_i+q_j))` retains bounds and complementary
fractions. This does not prove conservation of the full pressure/flux-coupled
solver; the eight CFD regressions test that interaction.

Runtime guards require:

- Exactly two moving phases named `gas` and `liquid`, in that order, in
  `basicMultiphaseSystem`, with a fixed mesh, global timestep, one alpha
  correction and the Euler phase audit enabled.
- Pure incompressible constant-density phases, with the fixture densities
  1 and 1000 kg/m3. This deliberately does not enable the candidate for an
  arbitrary reactor or biological case.
- Finite fractions in [0,1], complementary within 1e-12, and finite
  nonnegative pressure-source rates. Both rates may be unavailable together
  at startup; their initial contribution is zero. Mixed availability fails.
- Nonnegative finite option sources and zero interphase transfer. Unsupported
  states stop rather than falling back to a different numerical formula.
- Agreement between the original raw dilatation and the candidate Su/Sp
  contribution at source assembly, scaled by the rate sum, within 1e-12.

The executable requires a distinct `sc1TwoFieldRatesLibraryVersion` symbol.
The build manifest has variant `sc1_two_field_nonnegative_rates_v1`; baseline
and candidate verifiers reject each other's build manifests. Every executable
and project library is hashed and inspected with `ldd`; mandatory executable
and phase-library symbols are checked with `nm`. Installed project binaries
and original sources must remain unchanged.

## Regression sequence and criteria

The candidate uses exactly the original eight case recipes except for the
`application` name. No timestep, source, phase fraction, mesh, tolerance,
pressure setting or reference setting is retuned.

| Stage | Cases | Progression requirement |
|---|---|---|
| 1 | Closed transport, no reference and gas reference | Both full screens pass |
| 2 | Mid-fraction liquid feed, no reference and gas reference | Both full screens pass |
| 3 | Near-bound coarse pair, dt=0.001 s | Both diagnostics pass; mass failures remain findings |
| 4 | Near-bound fine pair, dt=0.0005 s | Both diagnostics pass; print full comparison |

The runner checks all eight destinations are fresh before launching the
first. Both cases in a pair are completed before its gate is evaluated. An
invalid gate stops subsequent stages and preserves attempted cases. A valid
near-bound mass-screen failure may proceed to half-dt to characterize the new
candidate, but final acceptance requires **all eight full screens to pass**.
`RATE_REGRESSION=PASS` returns status 0; a complete comparison with any failed
screen prints `RATE_REGRESSION=FAIL` and returns 2. Setup, solver or diagnostic
failures return 1. There is no automatic repeat or overwrite of attempted cases.

Preparation and launch revalidate the four completed original baseline
screens and the four completed near-bound controls against their logs,
initialization snapshots, input hashes and build assets. New reports retain
all original physical budgets and tolerances. They additionally require a
unique candidate runtime marker and one `SC1_RATE_STATE` row per alpha call,
checking coverage, source sign, complement and assembly identity.

For the candidate, a nonzero floor-induced split is no longer required. The
near-bound fractions and positive pressure-source rates must still be
exercised, and the sum of absolute `splitAtSourceAlpha` volumes across all
outer calls must be <=1e-18 m3 for each phase. This extra diagnostic is not
physical mass added to the budget. Original baseline report defaults retain
their old behavior and the completed control results remain valid.

The final read-only comparison reports each candidate phase residual alongside
its original control, its threshold, final inventory, export, split diagnostic
and both log hashes. A numerical pass does not select a physically validated
reactor model.

## Historical Kestrel command (completed; do not repeat)

This completed command is retained for provenance. It created a separate build and case
root, verifies the new runtime, then runs all eight regressions with the
progression gates above. It needs the completed controls and their original
build to remain available at their recorded paths. It launches no reactor CFD.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
CONTROL=/scratch/nwimer/sc1_two_field_bounds_20261111T004959Z
RATE_BUILD="/scratch/nwimer/sc1_two_field_rates_build_$(date -u +%Y%m%dT%H%M%SZ)"
RATE_CASES="/scratch/nwimer/sc1_two_field_rates_cases_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_rates.py" \
        build --out "$RATE_BUILD" \
        && PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_rates.py" \
        prepare --build-root "$RATE_BUILD" --control-root "$CONTROL" --out "$RATE_CASES"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_two_field_rates.py" \
            run-all --root "$RATE_CASES" --run
        echo "Rate-split regression status: $?"
    else
        echo "Build/preparation stopped. Preserve created directories and paste the output."
    fi
else
    echo "Rate-split setup stopped. Preserve existing files and paste the output."
fi
```

Return the output, including `RATE_BUILD`, `RATE_CASES`, all reports and
`RATE_REGRESSION`, or the first failure. If a solver completed but reporting
failed, preserve `log.solver`; use the read-only `analyze --log ... --case ...`
command to inspect it. Do not rerun the full block to replace completed work.
`compare --root ...` is also read-only, requires all eight saved results and
the sourced OpenFOAM/libfabric environment, and performs the same provenance
verification. Partial-run recovery is not automatic.

## Validation and route to the reactor

Local checks compile and execute the actual C++ coefficient kernel over
interior, near-bound and pure states, varied nonnegative rates and timesteps,
and invalid signed/nonfinite states. The actual runtime guard header is also
compiled against a minimal API double to exercise startup, application
coverage, raw-source identity and unsupported-model/state failures. This API
double does not substitute for OpenFOAM compilation on Kestrel.

Python tests cover preservation of every unchanged equation around the source
replacement, distinct private build identity, unchanged original sources and
baseline report output, eight case recipes, runtime evidence coverage, retained
mass failures, altered inputs/controls, progression gates, read-only comparison,
explicit launch/allocation requirements and exit status. CFD is mocked locally.

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover \
    -s scripts/diagnostics -p 'test_sc1_two_field*.py' -q
```

All eight numerical screens passed on Kestrel. The next stage is a matched
nonreacting two-field hydrodynamic candidate with the actual reactor geometry,
fill, vent and MRF treatment. The first quiescent actual-mesh startup stays within the present
pure-phase 1/1000 kg/m3 guards. Broader sources/thermodynamics require
separate evidence before extending the contract. Check liquid/gas conservation
and wall/free-surface behavior before oxygen transfer and biology. The existing
biological schedule remains unchanged, and a third phase remains optional
pending a demonstrated physical need.
