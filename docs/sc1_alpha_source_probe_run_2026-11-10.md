# SC1 source-probe build verified; fresh control run ready

**Status update:** The user reports that the probe control is running. Its
case path and completion log have not yet been supplied. Do not repeat the
launch block below. Await that result; the [external benchmark review](SC1_EXTERNAL_BENCHMARKS_2026-11-10.md)
provides references for interpreting it and selecting subsequent tests.

## Verified build, no new CFD result

`Pasted text(20261110-192420).txt` confirms successful private phase-library
and solver compilation and source-probe runtime verification at:

`/scratch/nwimer/sc1_alpha_source_probe_20261110T191459Z`

Kestrel resolves this path under `/kfs3/scratch/nwimer`. Reuse this build;
do not rebuild it or the earlier alpha-audit tree.

| Artifact | SHA256 |
|---|---|
| `bin/putidaFoamSC1AlphaAudit` | `0f38c43c6c986810c6bf4578391cac6c71fefc22d73b5b4aef05432920a1c11b` |
| Private `lib/libbirdphaseSystem.so` | `672a40b44bda7e65e43f38e32a27b3a1c2217471c13ee768cb9564914573a20a` |
| Unchanged installed `libbirdphaseSystem.so` | `222553b342e1414a857034e5613c9c5ce8202ca883764b3dd77ed001ef1828a6` |
| Film extension | `f823c8a3a3c83067171882cae53dfce4f9988dbc34c07bae6f6d58668c77ab4d` |
| Entrainment extension | `27c8ce353108690a041454d2a38957b9c88862e2215f0f75e56fb6a3b3f847fb` |
| `libalphaLimitVelocity.so` | `e63f7be7cee3503d3d1c7b9bb5131c791e5645b07bff7e061a3cca9a5ce2d5ac` |

The original 182 phase-source files, 185 prepared private phase-source files
and 34 prepared solver-source files all passed their provenance checks. The
solver, film and entrainment extensions resolve the new private phase library.
Both mandatory audit/probe symbols passed verification. The existing unused
variable warnings did not prevent linking. **No CFD has run with this build
in the supplied evidence.**

## One fresh control, unchanged equations and settings

The [half-timestep result](sc1_alpha_half_dt_results_2026-11-10.md) establishes
that neither halving dt nor using four outer iterations resolves the deficit.
The next run collects local source/bounds measurements absent from those logs.
It is a single zero-entrainment control over 0.1 s, with the baseline two outer
iterations, dt/maxDt 5e-5 s, six alpha subcycles and one alpha corrector.

[`run_sc1_alpha_source_probe.py`](../scripts/diagnostics/run_sc1_alpha_source_probe.py)
prepares this case from the completed pair
`/scratch/nwimer/sc1_alpha_pair_20261110T155635Z`:

- Copies constant/system settings and model assets from its completed control,
  checking the established control/normal configuration differences first.
- Copies the complete mesh and `31.199982372783058` checkpoint from the
  original source recorded in `comparison_manifest.json`, only after every
  initial-state hash matches. It does not copy the completed control's rewritten
  Co headers, normalize files or omit diagnostics.
- Adds only `sc1AlphaSourceProbe true` inside explicit alpha solver controls.
  All closure settings, timestep controls and biological updates remain intact.
- Pins the verified new binary and extension hashes, rechecks private linkage
  and original-source preservation, and writes a completion manifest last.
- Refuses reused, overlapping, redirected or previously attempted cases.
  Preflight verifies all copied inputs, model assets and baseline-log hash.
  Launch requires an active allocation and explicit `--run`, uses 52 ranks,
  scopes the private library path to the child process and creates the log
  exclusively. Failures preserve the case; no automatic retry occurs.

The original cases and both completed numerical variants are preserved. This
continuation still retains frozen OUR and is **not** the nonreacting 9.9 s
release checkpoint.

## Kestrel command

Run once inside the active allocation. Preparation and runtime verification
must pass before the one solver launch. The baseline took about 12 minutes;
the new observer adds work, so its runtime is not yet measured.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
PROBECASE="/scratch/nwimer/sc1_source_probe_control_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -n "${SLURM_JOB_ID:-}" && -x "$PY" ]] \
    && [[ "$(git -C "$REPO" branch --show-current)" = "codex/sc1-headspace-bubble-audit-20261029" ]] \
    && git -C "$REPO" pull --ff-only \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && "$PY" "$REPO/scripts/diagnostics/run_sc1_alpha_source_probe.py" prepare --out "$PROBECASE"; then
    echo "PROBECASE=$PROBECASE"
    "$PY" "$REPO/scripts/diagnostics/run_sc1_alpha_source_probe.py" run --case "$PROBECASE" --run
else
    echo "Setup stopped. Preserve all existing directories and paste the output."
fi
```

Return the full output, including `PROBECASE`, solver status, legacy regression
result and all `SC1_SOURCE_PROBE*` records. If reporting is interrupted after
the solver finishes, use the read-only `report --case "$PROBECASE"` subcommand;
do not rerun `prepare` or `run --run` to recover output. Omitting `--run` from
the `run` subcommand performs only preflight on an unattempted case.

## Interpretation and acceptance gates

The reporter checks a complete set of 16 probe rows: two solved phases, two
outer indices and four source-assembly alpha bins. It checks finite metrics,
terminal reporting times and the new library marker. It compares all old
`SC1_*` records and `BUBBLE_MASS_AUDIT` records against baseline at exact printed
precision, retaining record multiplicity. A mismatch is a review gate, not
proof of a physical change; preserve both logs to distinguish roundoff from
an instrumentation or setup problem. The summary's `normal` label denotes
the observer control here, with **zero** entrainment.

After observational regression, inspect source split, implicit evolution,
raw-source correction and bound excess by phase, alpha bin and outer index.
Do not add cumulative source integrals across outer indices as elapsed-time
mass. Counts represent repeated cell/subcycle samples, not unique cells or
spatial coordinates. Check bin sums against the existing same-outer alpha
ledger and the identity remainder before interpreting the local evidence.

The raw target holds the old state and limited flux fixed. The pressure target
uses the stored, lagged `divU`; it is not independently established physical
truth. Bound violations show that a source-only addback at that fixed flux
cannot retain bounds. Absence of violations is necessary, but does not prove
coupled phase-sum, pressure/density, species or PB consistency. These results
will guide a bounded coupled correction tested in minimal cases, not a blind
increase of entrainment efficiency, foam area or iteration count.

All **162 local SC1 tests pass**, including 15 new case/launch/report checks.
They cover original-state recovery after Co-header rewriting, unchanged
physics/biology/assets, modified inputs, build pins, redirected/reused cases,
allocation and exclusive-launch guards, incomplete/nonfinite records and
legacy-regression multiplicity. Local tests do not replace the pending Kestrel
observer run or physical validation. No production approval is implied.
