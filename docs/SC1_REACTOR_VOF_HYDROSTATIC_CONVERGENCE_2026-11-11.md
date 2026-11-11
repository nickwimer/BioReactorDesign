# SC1 hydrostatic startup: pressure convergence comparison

## Latest result: 24-pass startup PASS, archive and log verified

The run at `652de39` in
`/scratch/nwimer/sc1_vof_hydro_convergence_20261111T223454Z` reports
`VOF_HYDRO_CONVERGENCE=PASS`, status 0, normal End, all 200 steps and 800/800
converged pressure-correction groups. All original physical screens pass.
The largest final-pass initial residual is 9.98744667660074e-13, with zero
linear iterations on every final pass. The run contains 19200 pressure audits.
The existing hydrostatic binary was reused; only the prescribed pressure-pass
count changed from the ten-pass control.

Peak liquid/gas speeds remain 8.85673136e-8 / 1.95089614e-5 m/s, maximum Co
3.05069960e-8, and final native liquid inventory 999.9999999999998 L. Final
matched-time peak-speed differences from the control are approximately
6.6e-23 / -2.0e-20 m/s; liquid-volume difference is zero.

**The correct archive is now received and independently verified.** The first
attachment duplicated the old ten-pass control; the replacement closes that gap.
See [the full verified result](sc1_reactor_vof_hydro_convergence_result_2026-11-11.json).
All 17 listed file digests, 11 included initial-input hashes and ten current
recipe-file hashes match. Attempt identities agree. Re-parsing the complete log
reproduces the phase budgets, all 19200 pressure records and 800 converged groups.
The saved result exactly matches the console. The log records 282 solver wall
seconds. The old/new included fvSolution files differ only in pressure-pass count;
all 40 saved matched-time comparisons reproduce against the preserved control.

Verified archive:
`/scratch/nwimer/sc1_vof_hydro_convergence_20261111T223454Z_evidence/reactor_failure_evidence.tar.gz`,
SHA-256 `cfc5e595ba4d7a7e6085d1c458dfd3fd7f87236fb2870b9cd3ce6e9160ef5ac3`.
Solver-log SHA-256:
`704c0ccf7b59710b830efc12ad921c0e3442ec6245bff9a702320970809582e3`.
Manifest SHA-256:
`74f7662afbe39e1505b91b5d990ac37a8eb26e400bf633600ad55500638a7340`.
Saved-result SHA-256:
`fb397516c18791f0de82b625969139d8cd7c9c9a78b7f17e8f0c8befcf71660d`.

Native field arrays and runtime binaries are not bundled. Their checks remain
Kestrel-workflow evidence; local verification checks saved summaries and the log,
not a reconstruction from native arrays. This accepts the short resting startup
within that scope. Do not rerun it. See the
[current handoff](SC1_CURRENT_HANDOFF_2026-11-11.md) for longer-rest/perturbation
work and the remaining route to an operating simulation.

## Previous ten-pass control: verified result

The private hydrostatic solver compiled and completed the full 200 steps / 10 ms
at `f0c220f`, with normal End and all 40 requested reconstructed snapshots.
The new root is `/scratch/nwimer/sc1_vof_hydro_20261111T220400Z`;
the existing build is `/scratch/nwimer/sc1_vof_hydro_build_20261111T220400Z`.
Executable `sc1InterFoamHydrostatic` SHA-256:
`7da783091d8eb2d2df22ee6d657543850ce8b4c226f59fc39af8a35a61387c5f`.

All flow and mass-budget screens pass. The full acceptance flag is still false
because 203 of 800 pressure-correction groups did not meet the existing
final-pass criterion. This is a completed, numerically quiet run with an
outstanding convergence check, not a solver crash or an accepted baseline.

| Quantity | Received value | Existing limit |
| --- | ---: | ---: |
| Peak liquid speed | 8.85673136e-8 m/s | 1e-3 m/s |
| Peak gas speed | 1.95089614e-5 m/s | 1e-3 m/s |
| Maximum Courant number | 3.05069960e-8 | 0.2 |
| Final native liquid inventory | 999.9999999999998 L | Prescribed 1000 L; phase budget checked separately |
| Liquid residual L1 | 9.76996262e-12 kg | 1.00001e-5 kg |
| Gas residual L1 | 1.06581410e-14 kg | 6.52002269e-9 kg |
| Minimum physical pressure | 99982.9869167 Pa | 1000 Pa |
| Alpha redistribution | 4.24290535e-14 m3 | 1e-5 m3 |

All 17 files in the received archive match their recorded digests; all 11
included initial-input hashes match the manifest. The archive hash is
`33f5a524c9440516b0e2477c830ef7378ea64ef3e9d4b41f2c0de6f135296fad`.
The [full result](sc1_reactor_vof_hydrostatic_result_2026-11-11.json) records the
verification scope and residual histories. Log budgets and all 8000 pressure
records were independently re-parsed locally. Native snapshots and the external
build/control chain were verified by the pinned workflow on Kestrel; those
arrays and binaries are absent from the attachment and cannot be re-read locally.

## Why 24 passes

Every linear pressure solve finishes below 1e-12. The separate nonorthogonal
criterion requires the final pass to **enter** below 1e-12 and use zero linear
iterations. The pressure residual is normalized by the linear solver; a small
final residual alone is not the separate explicit-correction check
([Foundation OF9 PCG implementation](https://cpp.openfoam.org/v9/PCG_8C_source.html)).

The failed groups are distributed `[200, 3, 0, 0]` across the four corrections
per step. Thus every first correction and three second corrections need more
passes; both corrections in the final outer iteration pass on all 200 steps.
The largest final-pass initial residual is 8.01999344e-10. For the first
correction, the ratio of pass-10 to pass-9 initial residual ranges from
0.41879 to 0.54821 across the run: the iteration is still contracting.
Fourteen additional passes at the worst observed ratio would reduce the
largest residual to about 1.78e-13. This is an extrapolation motivating a test,
not a guarantee or replacement acceptance criterion.

The original nonorthogonal experiment was already converged and could not fix
the original force imbalance. That finding remains. Here the hydrostatic
formulation has passed the physical screens, and its measured pressure sequence
has a specific remaining convergence deficit. The new comparison changes only
`nNonOrthogonalCorrectors 9` to `23` (24 total passes).

## Fresh comparison and acceptance

`sc1_reactor_vof_hydro_convergence.py` uses the existing verified hydrostatic
binary. No compilation, equation change, threshold change or restart is needed.
It revalidates the received control log, manifest and result hashes, reconciles
its log budgets, and checks the original build/library/input chain. The new case
copies only original t=0 inputs; all input hashes except fvSolution must match.
The completed control and its snapshots remain untouched.

The candidate retains dt=5e-5 s, 200 steps, eight ranks, all physical inputs and
0.25-ms output cadence. It requires 19200 pressure audits, 800 correction groups,
all original phase budgets and physical limits, both phase speeds <=1 mm/s,
all 40 native snapshot checks including physical pressure reconstruction, and
the unchanged zero-iteration / <=1e-12 final-pass criterion for every group.
The report includes matched-time peak-speed and inventory differences against
the ten-pass control; agreement does not excuse a failed acceptance check.
The existing application hard stop remains 10 ms; monitored wall time is capped
at 30 minutes. The previous ten-pass run took 238 seconds of solver wall time.

## Historical Kestrel command — completed; do not rerun

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
HYDRO_CONTROL=/scratch/nwimer/sc1_vof_hydro_20261111T220400Z
CONVERGENCE_CASE="/scratch/nwimer/sc1_vof_hydro_convergence_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only git@github.com:nickwimer/BioReactorDesign.git codex/sc1-headspace-bubble-audit-20261029 \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric; then
    if PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_hydro_convergence.py" \
        prepare --control-root "$HYDRO_CONTROL" --out "$CONVERGENCE_CASE"; then
        PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/sc1_reactor_vof_hydro_convergence.py" \
            run --root "$CONVERGENCE_CASE" --run
        echo "Hydrostatic convergence status: $?"
    else
        echo "Preparation stopped. Preserve all roots and paste the output."
    fi
else
    echo "Setup stopped. Preserve existing files and paste the output."
fi
```

The user has run this command and supplied the passing console and correct
archive. Independent log/archive verification is complete as described above. Longer rest and physical perturbation testing
require new duration/restart-aware work; the current binary stops at 10 ms.
No operating checkpoint or biology result is accepted; the four-hour biological
update schedule remains unchanged.

## Validation

All 66 focused tests pass. Tests cover all 19200 pressure records, missing/duplicate passes, unchanged
initial/final convergence criteria, phase budgets and physical pressure/velocity
checks. They also exercise immutable control identities, the single dictionary
change, no-build reuse, overlap/attempt rejection, missing native fields,
matched-time comparisons and evidence collection on normal and failed runs.
The 24-pass CFD comparison has run on Kestrel and passed; its archive and log
are now independently verified. This
documentation update does not change solver code or acceptance thresholds.
