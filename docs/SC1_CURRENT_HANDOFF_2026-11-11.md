# SC1 current handoff — 11 November 2026

## Where we are

**The 24-pass hydrostatic candidate has passed the 10-ms resting startup on the actual MEML vessel mesh, with 1000 L of liquid and all 800 pressure-correction groups converged.** The correct archive and complete log are now received and independently verified. This is the first verified complete pass of this short startup screen. It is a useful numerical milestone, but the stirred, aerated reactor and the common 9.9-s operating checkpoint are still unfinished.

**Evidence scope:** [the verified result](sc1_reactor_vof_hydro_convergence_result_2026-11-11.json) records all 17 bundled file digests, 11 included initial-input digests, the independently reproduced phase budgets and all 19,200 pressure records. The first attachment duplicated the previous ten-pass run; the replacement is correct and closes that evidence gap. Native field arrays and external binaries are absent from this small archive, so their validation remains evidence from the pinned Kestrel workflow. No CFD rerun is needed.

This document is the current entry point. The [older session handoff](SC1_MEML_SESSION_HANDOFF_2026-11-06.md) retains the investigation history; its older “next action” paragraphs are not current launch instructions.

## Goal and constraints

Produce the simplest physically credible, mass/species-conserving, nonreacting **9.9-s common checkpoint** for MM, EMM and MEML biological comparisons. It must preserve the intended liquid inventory and support credible free-surface behavior, gas holdup and oxygen transfer. The historical oxygen deficit and artificial wall-side liquid rise remain unresolved acceptance questions.

- Three numerical fields are optional. The assistant proposed splitting bubble gas from headspace gas; the user requires that complexity only if justified by evidence.
- Do not force a target dissolved-oxygen value by arbitrary transfer, entrainment, bubble-size or area tuning.
- Preserve the intentional **four-hour biological update size** and documented CFD trigger/hold schedule. Their scheduling is not the present defect.
- Keep production solvers, original libraries and every completed/failed scratch root intact. New candidates use separate builds and cases.
- Commit and push every code or documentation change to `codex/sc1-headspace-bubble-audit-20261029`.
- There is no direct HPC access from this workspace. The user runs reviewed commands on Kestrel; an earlier allocation must not be assumed still active.

## Latest verified startup result

Run source commit: `652de393d9d12b53e115c827c7cb8a0c598ae441`. Case variant: `sc1_meml_vof_hydrostatic_convergence_v2`. Case root:

```text
/kfs3/scratch/nwimer/sc1_vof_hydro_convergence_20261111T223454Z
```

The console reports `VOF_HYDRO_CONVERGENCE=PASS` and exit status 0. The received log independently confirms normal `End`, all required records and passing physical/convergence screens. The collector reports 51 initial inputs unchanged; the 11 included initial files also verify locally. The run uses eight ranks, dt = 5e-5 s, 200 steps and 24 pressure passes per correction (`nNonOrthogonalCorrectors 23`). There are two outer iterations and two pressure corrections per outer iteration. The workflow reconstructed and checked all 40 requested native snapshots, according to its reported success; those native arrays have not been independently read here.

| Quantity | Verified log / saved-result value | Acceptance / interpretation |
| --- | ---: | --- |
| Simulated duration | 0.01 s | Short resting startup only |
| Pressure groups converged | 800 / 800 | All groups required |
| Pressure audit records | 19,200 | Complete prescribed sequence |
| Largest final-pass initial residual | 9.98744668e-13 | <=1e-12, with zero linear iterations |
| Maximum liquid speed | 8.85673136e-8 m/s | <=1e-3 m/s |
| Maximum gas speed | 1.95089614e-5 m/s | <=1e-3 m/s |
| Maximum Courant number | 3.05069960e-8 | <=0.2 |
| Minimum physical pressure | 99982.9869167 Pa | >=1000 Pa |
| Final native liquid volume | 999.9999999999998 L | Prescribed 1000 L |
| Alpha redistribution L1 | 4.24290535e-14 m3 | <=1e-5 m3 |
| Liquid mass residual L1 | 9.76996262e-12 kg | <=1.00001e-5 kg |
| Gas mass residual L1 | 1.06581410e-14 kg | <=6.52002269e-9 kg |

Final signed liquid/gas physical residuals are 2.22044605e-13 / -9.61231480e-24 kg. Logged final alpha lies in [0,1]. At 10 ms the liquid-volume difference from the ten-pass control is zero, and peak-speed differences are approximately 6.6e-23 m/s for liquid and -2.0e-20 m/s for gas. Increasing pressure passes closed the convergence requirement without materially changing the reported flow trajectory.

### Evidence received and verified

Verified bundle, now received:

```text
/kfs3/scratch/nwimer/sc1_vof_hydro_convergence_20261111T223454Z_evidence/reactor_failure_evidence.tar.gz
```

The collector uses the filename `reactor_failure_evidence.tar.gz` for this successful diagnostic too; its name is not an indication that the new run failed.

| Identity | SHA-256 |
| --- | --- |
| Correct new bundle, verified | `cfc5e595ba4d7a7e6085d1c458dfd3fd7f87236fb2870b9cd3ce6e9160ef5ac3` |
| New solver log, verified | `704c0ccf7b59710b830efc12ad921c0e3442ec6245bff9a702320970809582e3` |
| Received console transcript | `632e67b19f6adc51559fc15e3fd448705e9e9a1801c47987f490e70429c563da` |
| Earlier duplicate attachment, old ten-pass run | `33f5a524c9440516b0e2477c830ef7378ea64ef3e9d4b41f2c0de6f135296fad` |

The replacement archive identifies the correct `20261111T223454Z` root and source commit. Its paths contain only safe regular files. All 17 listed member digests, 11 included input digests and ten current recipe-file hashes match. Attempt/build identities agree with the manifest. The saved result exactly matches the earlier console; log budgets and pressure coverage reproduce independently. The log records **282 seconds of solver wall time**.

The supplied old and new manifests differ only in the `system/fvSolution` input hash. The included dictionaries confirm the sole edit is `nNonOrthogonalCorrectors 9` to `23`. All 40 saved matched-time comparisons against the old result reproduce. Their native mesh/field arrays and external runtime assets are absent, so this verifies saved summaries and log evidence, not a local recomputation from native snapshots. The original ten-pass result remains a failed convergence control.

## What changed, and what the current model represents

The current candidate is a private solver derived from Foundation OpenFOAM 9 `interFoam`, using a single velocity and a resolved gas/liquid interface. It is not the production BiRD reactor solver, and it does not yet model unresolved dispersed bubbles, PBM, oxygen transfer or biology.

The case contains 281,664 cells, total volume 1.6420022691408904 m3, liquid volume 1 m3 and gas volume 0.6420022691408908 m3. The prescribed flat interface has midpoint y=1.0185981622921294 m and transition width 0.0125 m. Densities are 1000 and 1 kg/m3; gravity is (0,-9.81,0). It is laminar, with stationary walls/propellers, no rotation, no sparging and no physical mass sources. Original viscosity and capillary terms remain active.

Earlier force measurements identified a discrete gravity/pressure imbalance at the sampled velocity peaks. The private candidate introduces a fixed analytic reference density rhoH(y) and hydrostatic pressure pH(y), with dpH/dy = -9.81 rhoH, and solves for perturbation pressure:

```text
p = pH + pPert
rho * DU/Dt = -grad(pPert) + (rho-rhoH)*g + viscous + capillary terms
p_rgh = p - rho*gh
```

The reference is frozen once; evolving density remains in inertia and the buoyancy departure. At the exact initial resting state the departure is zero. The physical vent pressure is retained. See [the formulation and source review](SC1_REACTOR_VOF_HYDROSTATIC_2026-11-11.md).

This has suppressed the observed short-time spurious motion. A continuum identity and a quiet initialized state do **not** establish a correct response to finite interface motion, rotation, aeration or restart. Those are the next validation tasks.

## Investigation findings to preserve

| Stage | Finding and implication |
| --- | --- |
| Historical inventory and wall audit | EMM has about 1151.077 L at 9.9 s; MEML about 999.812 L. Prescribed-fill estimates do not explain the EMM excess under the stated selection assumptions. Both have unresolved wall/free-surface shape concerns. |
| Actual wall/MRF preflight | Zero `walls` faces adjacent to `rotor` in both meshes rules out that direct rotating-wall assignment hypothesis. Partial `propeller` adjacency still needs physical interpretation. Do not blindly enlarge the zone or repeat this preflight. |
| Two-field conservation fixtures | Closed and mid-fraction source cases pass. Near-bound cases reveal a source-splitting conservation defect that persists under timestep halving. |
| Private rate-split correction | All eight small CFD regressions pass with original thresholds. This validates the fixtures, not the vessel or oxygen-transfer model. |
| BiRD resting-vessel variants | Initial startup conserves inventory but develops motion. Restored drag diverges; face/coupled drag collapses pressure; pressure-coupled algebra removes that immediate collapse but still violates quiescence. None is an accepted restart. |
| Original VOF and scalar-gradient comparison | Both stop at 8.4 ms on the 1 mm/s liquid-speed guard. Gas peaks approach 0.973 m/s. These are intentional guard stops, not evidence that a segmentation fault caused the physical failure. |
| Original VOF nonorthogonal test | All 672 observed pressure-correction groups converge, yet the same speed failure remains. More pressure passes alone did not fix that original force imbalance. |
| Two-step force probe | Verified face terms reconstruct eight samples / sixteen peak cells. Gravity/pressure imbalance dominates sampled motion; capillarity is much smaller there. Fresh versus matrix pressure reconstruction agrees closely. |
| Hydrostatic ten-pass control | Completes 10 ms and passes physical screens. Full acceptance fails because 203/800 pressure groups miss the final-pass criterion. Full received log/archive independently verified. |
| Hydrostatic 24-pass comparison | Correct archive/log verified; all physical screens and all 800 pressure groups pass. Short resting startup accepted within the stated evidence scope. |

The old 151-L anomaly is not quantitatively attributed by these small diagnostics. Likewise, correct inventory does not prove that wall climbing, the central depression, upper-impeller wetting or oxygen transfer is physically correct.

## Ordered work toward a working simulation

### 1. Completed: verify the short-run evidence

The correct bundle is received and verified as described above. Preserve both controls and the passing candidate. Do not repeat this run. The next implementation task is the longer-rest and perturbation validation below.

### 2. Implement and run a longer resting validation

The present executable **hard-aborts beyond 0.010000000001 s**. Its reporters also prescribe 200 steps and 40 snapshots. Editing `endTime` alone is insufficient. Prepare a separate duration-aware candidate and reporting contract, retaining the fixed reference and original physical screens, with complete log coverage, output checks and monitored stopping.

Use staged durations chosen against measured wall time and allocation; 0.02 s followed by 0.1 s are possible planning targets, not implemented commands. The ten-pass control took 238 solver wall seconds and the verified 24-pass run took 282 seconds; these are observed solver times, not guarantees for longer runs. Assess speed and kinetic-energy trends, inventory, pressure, interface drift and mesh/timestep sensitivity. Ten quiet milliseconds does not establish long-time equilibrium or stability.

### 3. Test a physical perturbation, with the reference held fixed

Demonstrate that a controlled, volume-preserving interface displacement or prescribed velocity forcing produces a credible restoring response while retaining bounded alpha and closed budgets. Define the forcing and acceptance measures before running; the resting 1 mm/s criterion is not a general operating-flow limit.

The current initializer requires t=0, zero velocity/perturbation pressure and the checked analytic flat fill. It constructs rhoH from the initial density after requiring agreement with that profile. **Do not initialize a perturbed interface and freeze that perturbed density as the reference:** that can subtract the buoyancy response the test is intended to assess. Construct the unperturbed reference independently or apply controlled forcing after the validated rest initialization. Unit algebra/profile tests cannot replace this CFD test.

### 4. Define restart support before using a checkpoint

Restart is currently unsupported: initialization rejects nonzero time, and reference rhoH/pH fields are held in memory with `NO_WRITE`. A valid continuation needs an explicit persistent reference and identity, consistent pPert/p/p_rgh, and a check that the reference is not silently regenerated from evolved density. The current 0.01-s output is evidence, not an accepted operating restart.

### 5. Add stirring, then aeration

After the resting and perturbation checks, use the verified SC1 operating speed and geometry for a controlled liquid/gas stirring case without sparging. Check actual MRF zones and nonrotating patches, torque/power, liquid inventory, boundary fluxes, Co, wall films, azimuthal free-surface profiles and upper-impeller wetting. Compare with relevant physical data and assess mesh/time sensitivity. A low mass residual alone is insufficient.

Then choose the simplest adequate aeration representation. A resolved-interface single-velocity VOF control is not automatically adequate for subgrid bubble slip, holdup, size distribution and oxygen-transfer area. Use the earlier conservative rate-split work where relevant, but validate gas/liquid/species budgets and bubble/headspace treatment in the actual model. Do not assign all headspace gas the interfacial area of dispersed bubbles. A third field remains an evidence-based option.

### 6. Establish the common 9.9-s nonreacting checkpoint

Only after hydrodynamics, aeration and oxygen-transfer checks pass should the project accept a matched nonreacting 9.9-s state for MM/EMM/MEML comparisons. Document gas/liquid/oxygen budgets, initialization, operating conditions and restart provenance. Then restore biology with the existing four-hour update and CFD trigger/hold schedule. The present short startup has not resolved the historical oxygen deficit.

## Exact working locations

All timestamps below are literal UTC directory names. `/scratch` resolves to `/kfs3/scratch` in the Kestrel manifests.

| Purpose | Location |
| --- | --- |
| Kestrel repository | `/scratch/nwimer/BioReactorDesign_bubble_audit` |
| Python | `/scratch/nwimer/.conda-envs/of-post/bin/python` |
| OpenFOAM environment | `$HOME/OpenFOAM9_RHEL9/OpenFOAM-9/etc/bashrc`, then `module load libfabric` |
| Hydrostatic build, reused by both runs | `/scratch/nwimer/sc1_vof_hydro_build_20261111T220400Z` |
| Verified ten-pass control | `/scratch/nwimer/sc1_vof_hydro_20261111T220400Z` |
| Verified passing 24-pass case | `/scratch/nwimer/sc1_vof_hydro_convergence_20261111T223454Z` |
| Verified force-probe case | `/scratch/nwimer/sc1_vof_force_20261111T205814Z` |
| Original VOF build | `/scratch/nwimer/sc1_vof_build_20261111T191346Z` |
| Original VOF nonorthogonal control | `/scratch/nwimer/sc1_vof_nonorth_20261111T203602Z` |
| Private rate-split build | `/scratch/nwimer/sc1_two_field_rates_build_20261111T021516Z` |
| Eight completed rate-split regressions | `/scratch/nwimer/sc1_two_field_rates_cases_20261111T021516Z` |

Current hydrostatic executable is `bin/sc1InterFoamHydrostatic`, SHA-256 `7da783091d8eb2d2df22ee6d657543850ce8b4c226f59fc39af8a35a61387c5f`. Both hydrostatic runs use this same binary. Preserve it and its source/build manifest. The new run's verified manifest SHA-256 is `74f7662afbe39e1505b91b5d990ac37a8eb26e400bf633600ad55500638a7340`; saved result SHA-256 is `fb397516c18791f0de82b625969139d8cd7c9c9a78b7f17e8f0c8befcf71660d`.

## Code, tests and next-session entry points

Read this handoff first, then [the 24-pass result/workflow](SC1_REACTOR_VOF_HYDROSTATIC_CONVERGENCE_2026-11-11.md), [the hydrostatic formulation](SC1_REACTOR_VOF_HYDROSTATIC_2026-11-11.md), and [the model-selection requirements](SC1_MODEL_SELECTION_AND_PROBE_RESULTS_2026-11-10.md). Full earlier evidence remains linked from the historical handoff.

Relevant implementation under `scripts/diagnostics/`:

- `sc1_reactor_vof_hydro.py`, `_build.py`, `_instrument.py`, `_report.py`: isolated hydrostatic build, preparation, instrumentation and validation.
- `include/sc1CreateHydrostatic.H`, `include/sc1HydrostaticProfile.H`: fresh-start restrictions, fixed reference and analytic profile.
- `sc1_reactor_vof_hydro_convergence.py`, `_report.py`: immutable ten-pass control verification, one-dictionary 24-pass comparison and full coverage checks.
- `sc1_reactor_watch.py`, `sc1_reactor_failure.py`: monitored stopping and evidence collection.

The latest implementation validation recorded 66 focused tests passing across hydrostatic, convergence, force, nonorthogonal, gradient, saved-field, VOF and failure-reporting modules. These include Python contracts, a standalone compiled C++ profile test and mocked workflows. The user confirmed the actual OpenFOAM 9 build on Kestrel. Local unit tests are not a CFD validation claim. This handoff update changes documentation only.

Build recipes, binaries, dependencies, inputs and parent results are pinned. Changing an old recipe helper can invalidate reuse of a preserved build; introduce a versioned candidate deliberately. Do not loosen screens merely to obtain a pass. Do not repeat completed experiments just because a historical document contains their launch command. The short-run evidence gap is closed. The immediate next action is to implement the longer-rest and fixed-reference perturbation work above.
