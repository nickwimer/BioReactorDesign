# SC1 MEML oxygen-transfer investigation — new-session handoff

**Updated:** 2026-11-10

**Research system:** *Pseudomonas putida* SC1 stirred-tank simulation, OpenFOAM 9 / BiRD, NLR Kestrel  
**GitHub branch:** `codex/sc1-headspace-bubble-audit-20261029`  
**Pull request:** https://github.com/nickwimer/BioReactorDesign/pull/4  
**Status:** Experimental closures compile and run in isolated scratch cases. The oxygen deficit remains unresolved; no experimental closure has been approved for production.

**Latest state (2026-11-10):** The four-outer-iteration zero-entrainment case
`/scratch/nwimer/sc1_alpha_numerics_20261110T180204Z/four_outer` completed
successfully: 2,000 steps, 8,000 alpha calls, the same 0.1 s interval and
initial state as the completed two-outer control. Read
[the four-outer result and half-timestep launch](sc1_alpha_four_outer_results_2026-11-10.md)
first. The liquid mass-budget shortfall falls 14.8104%, from 181.224915 to
154.384700 mg, but source splitting and coupling discrepancies remain.
The liquid response varies substantially between outer indices 3 and 4,
including in the last timestep; four iterations are not a converged remedy.
Mean OTR changes only −0.0175353%, and clipping still cancels 80.5772% of
requested liquid O2 removal. The instrumented alpha/species accounting closes.

**Next action:** Run only the **already-prepared `half_dt` case** in that
same numerical root, using `sc1_alpha_convergence.py run --case half_dt --run`.
It retains two outer iterations and halves initial/max timestep to 2.5e-5 s.
Do not prepare a new numerical root or rerun `four_outer`. The existing
launcher and verified build
`/scratch/nwimer/sc1_alpha_audit_20261107T001119Z` need no changes or rebuild.
All 131 local SC1 tests passed at `ae04e89`; this follow-up changes only
results documentation. The half-timestep case has not yet run. Its result
will guide source/coupling implementation work; no coefficient tuning or
biological schedule change is authorized by these diagnostics.

Both cases in `/scratch/nwimer/sc1_alpha_pair_20261110T155635Z` also remain
complete. Their normal preflight proved all four Courant mismatches were
outer-location-only changes. The matched alpha pair reproduced earlier v1
results and found normal entrainment improved mean OTR by only 0.147875%.
See [the matched pair record](sc1_alpha_pair_results_2026-11-10.md).
Preserve all completed cases; older launch instructions are historical.
Allocation availability must be checked rather than inferred from transcripts.

The earlier [full conservation pair](sc1_conservation_pair_results_2026-11-06.md)
remains authoritative for the 0.5 s normal/control comparison: own-duration
mean OTR improves only 0.353%, and normal's native final bubble inventory is
125.23 mg lower. The diagnostic continuations retain frozen OUR; they are not
the intended nonreacting 9.9 s checkpoint. Neither experimental closure is
approved for production. Historical follow-ups below record earlier gates;
the latest state above supersedes their pending-action instructions.

> **Instruction to the next assistant:** Read this handoff and the linked audit documents before changing anything. Review the actual BiRD source and trial outputs before forming new conclusions. Make **every code or documentation change as a Git commit pushed to the GitHub feature branch**. Do not leave work as an uncommitted patch, and do not modify the production CFD runs. Do not rerun already completed experiments without a reason.

**Conservation instrumentation follow-up (2026-11-06):** An opt-in timestep
ledger, final-PIMPLE O2 equation/source capture, PB pre/post-correction moments,
fresh solver/case preparers and log comparison tool are now implemented. Nine
local regression tests pass, including compilation/execution of the independent
C++ accumulator. **The corrected audit at `6ac7e36` now compiles on Kestrel;
all nine diagnostic and ten analytical checks pass there, and executable
dependencies resolve. There are no new control/active balance logs or
conservation conclusions.**
Read [the timestep audit design and build gate](sc1_timestep_conservation_audit.md)
next. This work explicitly distinguishes timestep quadrature from fully
discrete closure; full alpha/pressure/PB equation-level closure remains open.
All existing completed experiments and the biological schedule are preserved.

**Reliability investigation follow-up (2026-11-06):** Read the
[source review and recommended validation sequence](SC1_RELIABILITY_INVESTIGATION_2026-11-06.md).
It identifies omitted simultaneous gas counterexchange in the net-only model,
the need to validate deposition as free-surface disengagement, time-level and
clipping contributions to conservation residuals, and exact configuration
checks for capillarity and diffusive-transfer blending. Ten additional
analytical counterexample tests pass; they do not validate the OpenFOAM
implementation. The diagnostic stream API and terminal-interval predicate
were corrected. No physical closure or production case was changed, and no
new Kestrel simulation or experimental comparison has been performed.

**First build follow-up (2026-11-06):** The uploaded transcript confirms the
seven then-existing diagnostic checks and ten analytical checks passed on
Kestrel, but compilation failed at audit diffusion access. The corrected
lookup uses the registered phase transport model, with its required include
and link flags appended to the preserved local build options. Two additional
local checks pass. Preserve the failed scratch tree
`/scratch/nwimer/sc1_conservation_audit_20261106T182403Z`; rebuild in a fresh
tree. See the [build record](sc1_timestep_conservation_audit.md#first-kestrel-build-and-diffusion-access-correction-2026-11-06).

**Successful build follow-up (2026-11-06):** The next uploaded transcript
confirms `wmake` success at commit `6ac7e36`. Use the existing executable
`/scratch/nwimer/sc1_conservation_audit_20261106T183120Z/bin/putidaFoamSC1ConservationAudit`;
no further rebuild is currently needed. `ldd` reports no missing executable
dependencies, including the added transport library. This does not inspect
experimental libraries loaded later by `controlDict`. Next prepare a fresh
matched pair and review case settings, runtime assets and loaded libraries
before running the zero control. The matched audit retains the original
31.199982372783058 state and biology settings; it is not yet the nonreacting
9.9 s release checkpoint.

**Checkpoint comparison follow-up (2026-11-06):** The first matched-case
preparation stopped on a raw input mismatch in
`processors52/31.199982372783058`, after the collated mesh comparison passed.
No destination/cases were created and no CFD launched. The mismatch has not
been classified as metadata, file-set or physical-state differences. Use
`prepare_sc1_conservation_pair.py --inspect-inputs` (without `--out`) to obtain
file-level hashes, sizes and bounded byte evidence. The equality requirement
remains unchanged; do not strip headers, discard fields or silently replace
either checkpoint. See the [inspection instructions](sc1_timestep_conservation_audit.md#checkpoint-mismatch-before-case-preparation-2026-11-06).

**Mismatch narrowed (2026-11-06):** The uploaded inspection reports 36/36
collated mesh files identical, 55/59 checkpoint files identical, all 12 system
files identical, and only the expected `phaseProperties` difference among
51 constant files. The four changed checkpoint files are `Co.liquid`,
`Co.bubbleGas`, `Co.headspaceGas` and `Co.mixture`; each normal file is five
bytes longer and first differs in its outer-header run-directory path.
Their remaining bytes have not yet been independently verified on Kestrel.
The new opt-in `--allow-courant-location-metadata` proves equality of the
complete files except that one validated outer-header value, or refuses
preparation. Only after proof, it copies the entire control checkpoint to
both NEW cases and verifies every raw file hash. No original is altered and
no field is omitted. Sixteen new regression tests pass locally. The successful
solver executable is unchanged; no rebuild is needed. A new preflight shell
script replaces the long pasted collection block.

**Matched pair ready (2026-11-06):** The next Kestrel transcript confirms all
16 metadata regression tests passed and the four actual Co files differ only
in the proved outer `location` value. Both new cases contain the complete,
raw-identical control checkpoint and mesh. Use the existing pair
`/scratch/nwimer/sc1_conservation_pair_20261106T201549Z`; its `preflight.txt`
has been reviewed. The eight inspected dependency sources match the repository,
and the existing executable and experimental libraries have recorded hashes
and resolved dependencies. Preserve the intentional local biological/pH solver
differences. The next action is the **52-rank zero-control launch**, with a
runtime-library check for `libalphaLimitVelocity.so`; no rebuild or new pair
is needed. No conservation simulation has launched yet.

The actual case retains the saved-OUR explicit oxygen sink (`useMonadSink
false`), 200-CFD-second biology triggers, four-hour biological update size and
0.1-CFD-second pH triggers. Thus this 31.2–31.7 s comparison is a frozen-OUR
diagnostic, not a nonreacting checkpoint. No new biological update is due
during it. PB diameter drift is disabled; diffusion blending uses the default
0.3/0.7 thresholds; foam area multiplier is zero. Interpret balance residuals
with the configured six alpha subcycles, species residual-alpha stabilization
and final-iteration-only PB hooks in mind. See the
[verified configuration and launch instructions](sc1_timestep_conservation_audit.md#matched-pair-preflight-reviewed-2026-11-06).

**Zero-control audit completed (2026-11-06):** The new control in
`/scratch/nwimer/sc1_conservation_pair_20261106T201549Z/control` completed
with status 0, 10,279 steps and final time `31.6999975003429`. The final
bubble inventory and four transfer/closure diagnostics reproduce the preserved
control at displayed precision. The uploaded excerpt includes the cumulative
and final-window ledger, but not the startup/all-write records.

Over 0.500015 s the liquid source requested 248.931 mg O2 removal, transfer
added 49.047 mg, and final-iteration species correction restored 200.123 mg
(80.4% of the request). This is consistent with negative-O2 clipping in the
reviewed OF9 `normalise()` path, not physical oxygen supply. The net accounted
removal after correction is not independently validated biological uptake.
Liquid O2 residual after including that correction is only 0.00215 mg, but
bubble/headspace mass residuals are +386.461/−394.800 mg and liquid mass
residual is −1168.686 mg. Do not call the case conservative or physically
validated. The cumulative PB mass-moment correction −7.0545 g is repeated
normalization of size fractions, not measured loss of actual phase mass.

Next complete the already prepared **normal case unchanged**, using the same
audit executable, while collecting both full audit-record excerpts and old
control regression records. No rebuild, coefficient change or new biology
schedule is justified. Then target unresolved alpha/source stages,
residual-alpha species stabilization and pressure-density changes before
interpreting entrainment benefits. Read the
[control findings and next launch](sc1_timestep_conservation_audit.md#first-control-conservation-result-2026-11-06).

## 1. What we are trying to accomplish

Establish a **physically credible, mass/species-conserving, nonreacting three-phase MEML hydrodynamic and oxygen-transfer state at 9.9 seconds**, suitable as a common baseline for subsequent **MM, EMM, and MEML** biological-model comparisons. Preserve the physical liquid inventory, resolve the upper-impeller free surface and froth/entrainment mechanisms, and predict credible bubble residence, gas holdup and O₂ transfer **without forcing a target DO by tuning arbitrary coefficients**.

The user's biological solver intentionally advances biology at an equivalent **four-hour interval per biology update** and holds OUR fixed over many CFD seconds. **Do not treat that scheduling as a bug or silently change the biological update scheme.** The immediate task is to repair/validate hydrodynamics and transfer before using any new state as a biological checkpoint.

Why this is urgent: historical reacting runs near 30 s show approximately:

| Model | Liquid inventory | Dissolved O₂ | Interpretation |
|---|---:|---:|---|
| Original three-phase MEML | ~999.5 L | ~0.0193 mg/L | Liquid inventory credible; oxygen severely depleted |
| Two-phase EMM | ~1137 L | ~4.61 mg/L | More dissolved oxygen, but unphysical excess liquid inventory |
| MAX-OUR comparison | ~760 L | ~5.92 mg/L | Not a physically comparable liquid inventory |

Those are historical diagnostics, **not matched-condition validation data**. Original MEML biological OUR is ~**0.498 g O₂/s**, versus ~0.504 g/s EMM. The oxygen limitation is a transport/physics question, not justification to weaken biology.

## 2. Confirmed original-model physics and findings

The original MEML uses `interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem` with three Eulerian phases: `liquid`, `bubbleGas` and `headspaceGas`. Its original `constant/phaseProperties` defines Henry-law O₂ composition and Higbie/spherical diffusion **only between liquid and bubbleGas**. It has **no headspaceGas–liquid Henry/diffusion pair**, and no explicit reverse `headspaceGas -> bubbleGas` entrainment mechanism. Bubble-to-headspace `phaseTransfer` deposition is present, but gas-phase deposition is **not oxygen dissolution**.

Original diagnostics: bubble-gas inventory about 7.6 L at ~30 s, with ~99% below 0.8 m and ~1.08 s inventory/flux residence proxy; sparger gas injection about **8.743 g/s**, with blended bubble-to-headspace deposition typically about **7–12 g/s**. Original bubble-to-liquid O₂ transfer ~**0.084–0.094 g/s**, far below biological OUR. Prior EMM net bubble dissolution was ~0.312 g/s. The original MEML surface forms a substantial central depression: roughly y~0.7 m inboard versus ~1.19 m near the rim. This was already evident before reactions (9.9 s).

At 9.9 s the upper impeller's area-weighted phase fractions were approximately **74.4% headspaceGas, 25.4% liquid, 0.2% bubbleGas**; the lower impeller was substantially liquid-wetted. **Do not describe the upper impeller as physically “dry.”** The user reports a foamy, bubbly upper-impeller region in experimental video; the actual video has **not been quantitatively analyzed**. A gas-rich Eulerian cell can still represent unresolved wetted froth.

The deposition term is **blended** by BiRD. The independent, unblended correlation is much larger than the actual assembled pair source. The diagnostic `orderedDepLoss_kg_s=0` is an ordered-key audit artifact and does not prove deposition is absent. Always use **signed, blended `fluid.dmdtf` and assembled phase sources** for conservation claims.

The open top `outlet` has strongly opposing gross headspace flows: near 31.3 s, gross inward ~**1.14072 kg/s**, gross outward ~**1.14847 kg/s**, signed net outward only **0.00774726 kg/s**. All three O₂ fields have `inletOutlet` composition: headspace and bubble reverse-flow O₂ mass fraction `0.23`; liquid `8e-6`. That is plausible for returning air and is **not an obvious O₂ inlet specification error**. The large bidirectional boundary flow still requires independent flow/pressure validation. A reported `headspaceBoundaryOut_kg_s` is **signed net**, not gross outward flux.

## 3. Implemented code and trial findings

### A. Separate headspace-to-liquid oxygen-transfer closure

The GitHub subtree `OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceFilmModels/` implements an optional runtime `resolvedHeadspaceFilm` model and the helper `scripts/diagnostics/prepare_sc1_headspace_o2_trial.py`. It adds conservative coupled Henry O₂ exchange across `headspaceGas`/`liquid`, retaining original `bubbleGas`/`liquid` transfer. The model uses resolved `|grad(alphaH)|` and alpha gates; optional localized foam subgrid area is exploratory. The 10 μm film thickness, foam diameter 3 mm and foam-area multiplier are **not experimentally calibrated**.

Paired independent 0.5 s runs from the same source checkpoint:

| Trial | Effective headspace–liquid area, five-write mean | Headspace-to-liquid O₂ | Bubble-to-liquid O₂ | Total modeled O₂ |
|---|---:|---:|---:|---:|
| `run_geometric_05`, foam factor 0 | ~1.72 m² | 0.003644 g/s | 0.094754 g/s | **0.098398 g/s** |
| `run_foam_05`, foam factor 1 | ~5.58 m² | 0.011627 g/s | 0.094689 g/s | **0.106316 g/s** |

These five-write arithmetic means are **not equilibrated, continuous time integrals**. Foam increased headspace O₂ transfer ~3.19× but raised total O₂ delivery by only ~0.0079 g/s. The model responds to its interface area, but that does **not** validate the area. Total supply is only ~20–21% of the ~0.498 g/s historical OUR. Do not fit area to an extrapolated “required” ~190 m² as if it were measured.

### B. Separate conservative headspace-to-bubble entrainment and PB nucleation

The GitHub subtree `OFsolvers/birdmultiphaseEulerFoam/sc1HeadspaceEntrainmentModels/` contains:

- `headspaceEntrainingDeposition.{H,C}`: preserves the original deposition correlation, plus an optional bounded countervailing headspace→bubble phase-transfer source, implemented in the *same original ordered pair* rather than adding an unexamined opposite ordered model.
- `headspaceEntrainmentNucleation.{H,C}`: assigns **local positive net** bubble-phase mass gain after deposition/blending to the existing population balance, through `popBal.eta` and a trial **2 mm** nucleation diameter.
- `scripts/diagnostics/prepare_sc1_headspace_entrainment_trial.py` and `test_sc1_headspace_entrainment_trial.py`: non-destructive scratch preparation and regression checks.

The new experimental library `libsc1HeadspaceEntrainmentModels.so` **successfully compiled and linked** on Kestrel. The diagnostic solver `putidaFoamMask01BubbleAudit` was rebuilt in scratch. All completed 52-rank runs reported status 0. The **zero-entrainment** `run_control_05` reproduced the original `run_geometric_05` at displayed precision for inventory, deposition, both O₂ transfer channels, and zero PB nucleation: the regression control **passed**.

The initial `legacyDownwardLiquid` criterion `U.liquid.y < 0` almost completely suppressed entrainment: 0–7 eligible cells among ~4,000 mixed-interface cells. Gate diagnostics showed ~0.33 m² interface in the selected upper-impeller window but only ~0–6.4×10⁻⁵ m² after the downward-liquid filter. This was an *eligibility-definition problem*, not proof that true entrainment is physically small.

The explicit opt-in `normalGasIntoLiquid` criterion instead accepts positive `(U_liquid-U_headspaceGas) · grad(alpha.headspaceGas)/|grad(alpha.headspaceGas)|`. It is a **relative-motion proxy**, not a validated breakup/detachment law. The default remains the original legacy mode. Both criteria use the same axial/radial window, phase-fraction gates and donor limitations. **The normal mode is experimental and must not be promoted into production.**

### C. Latest 0.5 s normal-inward entrainment trial — completed 2026-11-06

`run_normalactive_05`: same collated starting state, `entrainmentEfficiency 0.01`, `entrainmentVelocityCriterion normalGasIntoLiquid`, raw global cap `0.002 kg/s` (**2 g/s**), `nucleationDiameter 0.002 m`, headspace film `foamAreaMultiplier 0`. The five outputs at t=31.3–31.7 s show:

- 1,791–1,900 eligible interface cells, ~0.160–0.169 m² area;
- proposed **raw unblended entrainment** ~2.55–2.69 g/s, capped to 2.0 g/s (scale ~0.743–0.784);
- **local positive net bubble mass birth entering PB** ~0.208–0.238 g/s, five-write average ~**0.223 g/s**;
- **integrated signed bubble↔headspace exchange remains net deposition**, ~7–10 g/s. Do not call the 2 g/s cap “actual bubbles created,” and do not call signed net deposition “gross entrainment”;
- PB size-group sums report `1 1 1`; diameter remains finite, mean ~1.205 mm, minimum ~1.003 mm, maximum ~5.331 mm at the end. **BiRD clamps/renormalizes** size fractions, so this is a numerical sanity check, *not* an independent PB moment-conservation proof.

| t (s) | Bubble inventory change vs zero control | Added total O₂ transfer vs control |
|---|---:|---:|
| 31.3 | +11.09 mg | +0.282 mg/s |
| 31.4 | +24.17 mg | +0.401 mg/s |
| 31.5 | +55.87 mg | +0.558 mg/s |
| 31.6 | −4.12 mg | +0.311 mg/s |
| 31.7 | **−125.23 mg** | +0.372 mg/s |

Five-write mean modeled O₂ transfer: **0.098782756 g/s active** versus **0.098398236 g/s control**, gain only **0.000384520 g/s = 0.385 mg/s** (~0.39%). The modified case still supplies only ~19.84% of historical frozen OUR. **No sustained bubble-holdup or meaningful O₂ improvement has been demonstrated.**

### D. Unresolved conservation problem

The existing `bubbleMassAudit.H` reports **endpoint** inventories, injection, signed blended deposition, net gas–liquid transfer, and physical patch fluxes only at write times. Variable CFD timesteps and highly variable deposition make sparse 0.1 s rates insufficient to certify mass conservation. The t31.4–31.7 endpoint closure residuals are nonzero for **both** normal-active and zero-control runs. A coarse trapezoidal balance over 31.3–31.7 s misses the observed bubble inventory change by about **+0.231 g active** and **+0.382 g zero control**. This discrepancy is **not yet evidence of a real numerical leak**: the budget has not been constructed from the same timestep-discrete fluxes and source evaluations as the solver.

A further open question is whether donor O₂ is transferred and conserved correctly when headspace gas is reclassified as 2 mm bubbles, given the strong gross counterflow through the open top and the competing dissolution pathways.

## 4. File paths and working environment

All Kestrel paths below are **on the user's cluster** and are not directly available to a new ChatGPT session without connected access or user execution.

~~~bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
BIRD=$HOME/OpenFOAM9_RHEL9/BiRD/birdmultiphaseEulerFoam
BASE=/projects/microenviro/final_runs/reacting_Julia_08282026/meml/putida_sc1_final
SOURCE=/scratch/nwimer/sc1_bubble_audit_20261029/run
SOLVER=/scratch/nwimer/sc1_bubble_audit_20261029/solver
EXE=$FOAM_USER_APPBIN/putidaFoamMask01BubbleAudit

FILMROOT=/scratch/nwimer/sc1_headspace_o2_trial_20261103
FILMLIB=$FILMROOT/lib/libsc1HeadspaceFilmModels.so
ENTROOT=/scratch/nwimer/sc1_headspace_entrainment_20261104
ENTLIB=$ENTROOT/lib/libsc1HeadspaceEntrainmentModels.so
PY=/scratch/nwimer/.conda-envs/of-post/bin/python

source "$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc"
module load libfabric
~~~

`/scratch` may resolve to `/kfs3/scratch`. The common source checkpoint in the experimental collated `processors52` directory is **`31.199982372783058`**. Runs use **52 MPI ranks** with `srun -n 52 "$EXE" -parallel -fileHandler collated`. Current allocation from 2026-11-06: job **19721806**, node `x1006c3s5b1n1`; this is transient—**never assume a prior allocation is still active**.

Preserve existing scratch trials:
- `$FILMROOT/run_geometric_05`, `$FILMROOT/run_foam_05` (logs `log.headspaceO2.smoke`);
- `$ENTROOT/run_control_05` (log `log.entrainment.control`);
- `$ENTROOT/run_entrain_05` (log `log.entrainment.smoke`, old downward gate);
- `$ENTROOT/run_gateaudit_05` (log `log.gateaudit`);
- `$ENTROOT/run_normalgate_02` (log `log.normalgate`);
- `$ENTROOT/run_normalactive_02` (log `log.normalactive`);
- **`$ENTROOT/run_normalactive_05` (log `log.normalactive05`, latest key experiment).**

All prior named runs completed; **do not overwrite them, silently continue them, or edit their outputs**. Create fresh, uniquely named scratch clones for new experiments.

## 5. Precisely what the next session should do

**Priority 1 — Build a timestep-resolved, conservative phase/species audit in the scratch diagnostic solver.** Before writing new physics, read `OFsolvers/birdmultiphaseEulerFoam/putidaFoam/putidaFoam.C`, `bubbleMassAudit.H`, and BiRD's `PhaseTransferPhaseSystem` and `PopulationBalancePhaseSystem`. The current audit includes a post-PIMPLE, `runTime.writeTime()`-gated include; this is the likely place to extend, **but the new accumulation must run every time step and print only on selected writes**. Do not run the expensive diagnostic multiple times per PIMPLE correction.

Plan to track:
1. Whole-domain **liquid, bubbleGas, headspaceGas and combined gas** mass inventories, with start/end differences. Log physical-boundary *gross inward and gross outward* separately, especially the pressure outlet.
2. Per-timestep integrated **sparger gas injection**, correctly signed and blended gas–gas transfer, bubble–liquid and headspace–liquid O₂ transfer, and any other phase sources. Use the **actual assembled phase source and flux** and avoid double counting model terms. Accumulate source/flux integrals using the actual variable `deltaT`; distinguish an approximate endpoint/trapezoid quadrature from a provably identical solver-discrete budget.
3. Phase-resolved **O₂ inventories** `∫alpha_phase rho_phase Y_O2,phase dV`, liquid biological O₂ consumption, gas boundary O₂ in/out including inflow `inletOutlet` composition, and coupled O₂-transfer terms. Clarify whether all transported phase species are available and the role of species diffusion/advection.
4. Local positive **net** PB bubble birth and whether it is consistent with changes in size-group moments *before renormalization*, without asserting that `sum(fi)=1` proves conservation.
5. At each saved comparison time, report **mass and oxygen conservation residuals** alongside total liquid inventory, bubble holdup, transfer rates, DO and the free-surface state.

The audit must be **disabled by default**, run only in an isolated scratch diagnostic executable/case, and unit/regression tested where possible. Perform a **zero-entrainment control first**, then the matching normal-inward trial from the **same checkpoint**, using fresh scratch paths. Do not relabel the old 0.1 s endpoint closure as a solver-discrete balance.

**Priority 2 — Validate the actual free surface and impeller physics.** Re-examine the central cavity, upper-impeller gas and liquid contact, gas recirculation at the outlet, and resolved versus subgrid interface area. Obtain/inspect the actual experimental video and any measured gas holdup, bubble size, impeller wetting, free-surface shape and `k_L a`/OTR if available. Distinguish coherent headspace cavitation, dispersed bubbles and liquid-wetted froth. A positive local normal relative velocity is **not sufficient** for bubble detachment.

**Priority 3 — Decide whether the three-phase closure deserves development.** The film model is functionally active but far below OUR; the normal-inward entrainment model is numerically active but raises OTR by only ~0.385 mg/s. Do **not** blindly increase the entrainment coefficient, raw cap, 2 mm nucleus count or foam area to force DO. A different physical closure or even a revised hydrodynamic representation may be required. Preserve trustworthy features of MEML (notably liquid volume) and benchmark against EMM/MAX-OUR without mistaking their nonphysical liquid inventories for correct behavior.

**Priority 4 — Only after physics, closure and conservation checks pass:** generate/re-equilibrate the **nonreacting 9.9 s MEML checkpoint** and launch controlled, common-initial-state MM/EMM/MEML biology comparisons. Preserve the intended biology-update schedule.

**Deliverables for a successful next session:** (a) a pushed GitHub commit with a verified time-resolved audit, (b) fresh scratch-only control and active balance logs, (c) an interpretation of gas/O₂ and PB conservation, (d) an evidence-based choice of the next physical model or calibration data—not another speculative coefficient sweep.

## 6. Collaboration, execution and safety rules

- **Always commit and push every code/document change** to `codex/sc1-headspace-bubble-audit-20261029`. Use the existing GitHub connected tools or an explicitly confirmed `git push` from the user's Kestrel checkout; ensure the pushed commit is visible. The user previously specifically objected to changes without pushes.
- **Do not modify production cases** in `/projects/microenviro/final_runs`, original solver files in `$BIRD`, or already completed scratch cases. Make scratch-only source/case clones; compile experimental models into separate scratch libraries.
- **Never paste `exit`, `logout`, `exec`, top-level `set -e` / `set -euo pipefail`, or shell-exit traps into the interactive `salloc` shell.** Such commands previously ended valuable allocations. Use guarded `if` blocks and commands that return to the prompt.
- Avoid claiming you can directly inspect Kestrel files from ChatGPT. Ask the user to run **safe, read-only checks** when needed; connected GitHub source and Kestrel scratch runtime are different systems.
- Do not confuse **kg/s, g/s and mg/s**, candidate gross versus actual **net** flux, or total oxygen flux versus gas-mixture exchange.
- Do not interpret an SSH disconnect as a CFD failure if the solver already printed `End` and exit status 0. Do not mistake OpenFOAM's routine `FOAM_SIGFPE` setup line for an exception.
- OpenFOAM's `foamDictionary` can fail when reading full **binary collated field** data; the documented workaround decoded the text `boundaryField` portion read-only. Plain ASCII `constant/phaseProperties` dictionaries work with `foamDictionary`.

## 7. First five minutes in the new session

1. Read **this document**, [the detailed entrainment log](sc1_headspace_reentrainment_trial.md), [headspace O₂ design/results](sc1_headspace_oxygen_transfer_design.md) and [the historical three-phase audit](sc1_three_phase_headspace_audit.md).
2. Check the **current** GitHub branch head, pending PR and existing source. Do **not** assume the checkout/cluster job is at the same state as the last conversation.
3. Confirm the latest validated result is `run_normalactive_05` status 0 with five writes ending 31.7 s; the zero control is `run_control_05`. Do not re-run those simulations.
4. Start with the **phase and O₂ timestep-integral conservation design**, its sign conventions and scratch-only instrumentation. Push the source/document commits before asking the user to test them.
5. After tests, interpret **both** control and normal-active conservation residuals before deciding whether to modify the physical entrainment model.

**One-sentence baton pass:** We implemented the missing headspace–liquid O₂ pathway **in an isolated scratch prototype** and added capped headspace-to-bubble entrainment coupled to the PB, but the latter does not materially increase O₂ transfer or sustained bubble holdup; next we must establish trustworthy timestep-resolved phase/O₂ conservation and validate foam/entrainment physics before advancing toward a common nonreacting 9.9 s checkpoint.
