# SC1 matched conservation pair: results and stage-audit follow-up

Both 52-rank diagnostics completed with status 0 using the existing audit
executable (SHA256 `e8b71ba803b23661fd137b82f3a42fe105e926eea4da170a74c2290edbc77dde`).
Preserve `/scratch/nwimer/sc1_conservation_pair_20261106T201549Z/control` and
`normal`, including their logs. The user-supplied `Pasted text(10).txt` contains
both audit histories, the old control records and normal completion. The
earlier control transcript contains its completion. Neither run is a validated
nonreacting release checkpoint. All physics and the four-hour biological update
size remain unchanged.

This continues the [handoff](SC1_MEML_SESSION_HANDOFF_2026-11-06.md) and
[timestep audit](sc1_timestep_conservation_audit.md).

## Initial-state and regression checks

All five printed initial phase/aggregate dictionaries match exactly, including
mass, O2 mass, volume and start time `31.1999823727831` s. The mesh and complete
checkpoint were independently hash-verified during preparation.
All five new control `BUBBLE_MASS_AUDIT` lines match the preserved
`run_control_05/log.entrainment.control` lines **exactly at printed precision,
across every field**. This supports observational regression of the v1 audit;
it is not proof of conservation or bitwise identity of every field.

The comparison reader correctly refused a synchronized comparison: adaptive
timesteps and `writeControl runTime` produced these actual report times:

| Report | Control time (s) | Normal time (s) | Normal minus control (microseconds) |
|---|---:|---:|---:|
| 1 | 31.2999823727864 | 31.2999823727864 | 0 |
| 2 | 31.4000186180940 | 31.3999916830166 | −26.9350774 |
| 3 | 31.5000072151629 | 31.4999916830199 | −15.5321430 |
| 4 | 31.5999975003396 | 31.5999971258958 | −0.3744438 |
| 5 | 31.6999975003429 | 31.6999846530835 | −12.8472594 |

Only the first reports share exactly the same timestamp. Do not round the
others to tenths and call them simultaneous. No rerun, interpolation or
tolerance increase is needed to report the independent cumulative budgets.
The reader now offers explicit `--native-times`: each run keeps its actual
duration, and endpoint differences are labeled unsynchronized. The default
synchronized comparison still rejects the observed clock differences. A valid
run's individual summary now prints even if its peer/comparison is rejected.

## Results at each run's own endpoint

| Quantity | Control | Normal |
|---|---:|---:|
| Duration (s) | 0.5000151275598 | 0.5000022803004 |
| Steps | 10,279 | 10,767 |
| Bubble mass (g) | 9.611905230 | 9.486675010 |
| Bubble volume (L) | 7.976086855 | 7.860242576 |
| Integrated liquid O2 transfer (mg) | 49.046806487 | 49.218686042 |
| Own-duration average O2 transfer (mg/s) | 98.090645230 | 98.436923152 |
| Requested liquid O2 removal (mg) | 248.931443241 | 248.924818759 |
| Final-iteration species-correction addition (mg) | 200.123373402 | 199.888177879 |
| Correction / requested removal | 80.3930% | 80.3006% |
| Liquid O2 inventory increase (mg) | 0.240884081 | 0.184197054 |
| Final liquid DO (mg/L) | 0.019585672 | 0.019528962 |
| Cumulative PB mass-moment correction (g) | −7.054524918 | −7.332027264 |

The separate-duration mean OTR improvement is only **0.346278 mg/s (0.353%)**.
This is the average of every-timestep integrals, unlike the historical
five-write arithmetic-mean difference 0.384520 mg/s. Neither is evidence of
long-time equilibrium or experimentally validated transfer.

The native-endpoint bubble difference is **−125.230220 mg**. Its v1 ledger
decomposition (normal minus control, unequal endpoint times retained) is:

| Bubble term | Difference (mg) |
|---|---:|
| Injection integral | −0.112321 |
| Assembled transfer integral | −132.716634 |
| Boundary outflow integral | +0.000032 |
| Unresolved mass residual | +7.598766 |
| Actual inventory change | −125.230220 |

The sampled transfer integral therefore explains most of the relative
inventory decrease within this ledger. Its gas–gas signed component is
132.539444 mg more negative with normal entrainment. This is a net assembled
exchange, not a measurement of gross deposition or entrainment separately.
Normal's positive-net pre-PB input proxy integrates to 110.379056 mg; it is
sampled at another stage and cannot be subtracted from post-PIMPLE exchange
to infer a validated gross-deposition integral.

### Residuals are persistent in both runs

| Phase | Control mass residual (mg) | Normal mass residual (mg) | Control O2 residual (mg) | Normal O2 residual (mg) |
|---|---:|---:|---:|---:|
| Liquid | −1168.686138 | −1140.334210 | +0.002147 | +0.002152 |
| Bubble gas | +386.460926 | +394.059692 | −3.598547 | −3.740233 |
| Headspace gas | −394.800287 | −402.533237 | +4.167649 | +4.485964 |
| Combined gas | −8.339360 | −8.473544 | +0.569103 | +0.745730 |
| Total | −1177.025498 | −1148.807754 | +0.571250 | +0.747882 |

Each normal 0.1-second window has approximately +77 to +81 mg bubble mass
residual, with an opposing headspace residual. This is not only a restart
transient. Both liquid budgets retain the same substantial correction after
the explicit frozen-OUR sink. The correction remains consistent with negative
O2 clipping in the reviewed OF9 path; cellwise measurement is the next step.

Sampled interphase mass sources cancel to floating-point accuracy by their
equal/opposite assembly. Total evaluated O2 transfer imbalances are only
0.011744 mg (control) and 0.015558 mg (normal); they do not explain the larger
gas O2 residuals. Linear equation L1 residual integrals remain about
`2.4e-5 mg`. Tighter linear tolerances or increased entrainment efficiency are
not supported fixes. The PB correction is repeated reconciliation of
`alpha*rho*sum(fi)`, not physical removal from `alpha*rho`.

## Implemented next diagnostic revision

The next pass adds observations to a **new copy of the current local scratch
solver**, retaining its local biology/pH source and build options. It does not
replace those with repository solver physics, change shared libraries, edit
the completed pair, or change a physical coefficient.

- `SC1_OXYGEN_STAGE` retains the v1 residual definition and decomposes it using
  the species DDT matrix, residual-alpha stabilization, actual unrelaxed
  equation residual, and a reconstruction mismatch. The explicit stabilizer
  is captured at assembly time, before Y changes. Its coefficient comes from
  `solverDict("Yi")`, exactly as in `MultiComponentPhaseModel`, not the phase
  residual alpha. The decomposition also measures old-time inventory offsets,
  the DDT/inventory difference, and changes between species correction and the
  final timestep inventory. Negative O2 mass, minimum Y and cell counts are
  measured immediately before correction and compared with the observed addition.
- `SC1_MASS_STAGE` brackets alpha/PB solving, model correction, species, energy
  and pressure blocks in **every outer iteration**. For each cell it splits
  measured mass change exactly into `(alphaAfter-alphaBefore)*rhoBefore*V`
  plus `alphaAfter*(rhoAfter-rhoBefore)*V`. Both-factor changes and their cross
  term are included. The ledger reports unassigned changes and decomposition
  remainders rather than hiding them. These are inventory changes, never
  multiplied by deltaT and never presented as physical source integrals.
- The preparer records revision 2 and header hashes, and still refuses to
  patch an already instrumented source or an unrecognized include layout.
  The case preparer now supports an explicit **0.1 s** diagnostic duration;
  0.5 s remains the default. Shortening a new diagnostic does not modify any
  biological parameter or existing case.

This remains a stage-attribution audit, **not exact alpha-subcycle closure**.
It does not intercept MULES `Su/Sp`, alpha-clamp corrections or individual
pressure-corrector equations. If the discrepancy is localized inside the
alpha/PB block, the following pass must instrument those actual library stages
in a separate scratch library before changing entrainment. Small bookkeeping
remainders establish the diagnostic decomposition, not physical validity.

All **50 focused local checks pass** (10 audit, 3 oxygen-stage, 8 mass-stage,
11 reporting and 18 checkpoint/copy checks), including compiled C++ arithmetic
and assembly-contract tests. They cover accumulator identities, stage ordering and
cross terms, stabilization/clipping identities, read-only reporting at unequal
times, source-preserving hook insertion, and complete short-case copying.
**The new OpenFOAM wrapper has not yet been compiled or run on Kestrel.**

## Next command: report old logs and build only

Do not rerun the completed pair. Use the updated reader on those existing
logs and compile a fresh solver tree. Return build status/output before
preparing the short control case. This block preserves the interactive shell
and confines build/log writes to the newly prepared scratch tree.

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_BRANCH=codex/sc1-headspace-bubble-audit-20261029
PAIRROOT=/scratch/nwimer/sc1_conservation_pair_20261106T201549Z
AUDIT2ROOT="/scratch/nwimer/sc1_stage_audit_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ "$(git -C "$REPO" branch --show-current)" == "$SC1_BRANCH" ]] \
    && git -C "$REPO" pull --ff-only origin "$SC1_BRANCH" \
    && source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc" \
    && module load libfabric \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_conservation_audit.py" \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_oxygen_stage_audit.py" \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_mass_stage_audit.py" \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_conservation_reporting.py" \
    && "$PY" "$REPO/scripts/diagnostics/test_sc1_checkpoint_metadata.py" \
    && "$PY" "$REPO/scripts/diagnostics/summarize_sc1_conservation.py" \
        "$PAIRROOT/control/log.conservation" --normal "$PAIRROOT/normal/log.conservation" --native-times \
    && "$PY" "$REPO/scripts/diagnostics/prepare_sc1_conservation_solver.py" --out "$AUDIT2ROOT"; then
    (cd "$AUDIT2ROOT/solver" && set -C && wmake > "$AUDIT2ROOT/build.log" 2>&1)
    SC1_BUILD_STATUS=$?
    echo "SC1 stage-audit build status: $SC1_BUILD_STATUS"
    echo "AUDIT2ROOT=$AUDIT2ROOT"
    tail -n 90 "$AUDIT2ROOT/build.log"
    if [[ "$SC1_BUILD_STATUS" -eq 0 ]]; then
        sha256sum "$AUDIT2ROOT/bin/putidaFoamSC1ConservationAudit"
        ldd "$AUDIT2ROOT/bin/putidaFoamSC1ConservationAudit"
        cat "$AUDIT2ROOT/source_manifest.json"
    fi
else
    echo "Stopped before build. Preserve any created directory and paste the output."
fi
```

After a successful build review, prepare a new pair with
`prepare_sc1_conservation_pair.py --allow-courant-location-metadata --duration 0.1 --out NEW_DIRECTORY`
and run its control only first. Compare its first-write v1 quantities to the
existing 31.3 s control and inspect the new oxygen/mass-stage records. Do not
extend a reacting continuation or promote any closure while the budget remains
unattributed. Quantitative free-surface/foam comparison and the eventual
nonreacting 9.9 s checkpoint remain downstream of this work.
