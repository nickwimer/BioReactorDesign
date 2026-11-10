# SC1 configuration evidence: phase count is not an isolated comparison

The full `sc1_wall_preflight_20261110T212232Z.tar.gz` has now been inspected.
Its SHA-256 is
`82486f6334397e4a4d83cbe0589c6e68544e4dc1a8c52d7d34ff2d475f4ce3af`.
It contains the completed preflight's two configuration reports, detailed
sector profiles and summary. The archived dictionaries and saved 9.9 s
boundary entries are observations; they do not pin the historical executable
or prove that every currently archived dictionary was unchanged throughout
the original spin-up. Backup files ending `~` are not active settings.

## Most consequential differences

| Setting | Historical EMM | Original MEML |
|---|---|---|
| Phase fractions | gas + liquid; no `referencePhase` | bubbleGas + headspaceGas + liquid; reference headspaceGas |
| Alpha limiting reference | `meanFluxReference` absent (false in inspected code) | true |
| Interface compression | empty list | headspaceGas/liquid pair, coefficient 1 |
| Liquid outlet velocity | pressureInletOutletVelocity | slip |
| Liquid outlet alpha | inletOutlet, incoming alphaL=0 | zeroGradient |
| Liquid wall alpha / velocity | zeroGradient / noSlip | zeroGradient / noSlip |
| Liquid turbulence | mixtureKEpsilon | kEpsilon |
| Gas turbulence | mixtureKEpsilon | both gas fields laminar |
| Thermal phase model | multiComponentPhaseModel | multiComponentIsothermalPhaseModel |
| Transfer | saturated H2O + Henry O2 | Henry O2 for bubbleGas only in this original archive |
| Gas–liquid drag | dispersed gas-in-liquid Grace | dispersed bubble/liquid closures plus segregated headspace/liquid drag |
| Gas source | point selection, full rate | 15 cell zones representing 23 holes, 0.5 s ramp |
| Initial bulk alphaL | 0.99 below y=1.016 m | 0.999 plus a partially filled next layer |
| Explicit liquid feed | 1.28287e-6 m3/s | 2.9506e-6 m3/s |
| Time-step settings | maxCo 0.5; maxDeltaT 1 s | maxCo 0.3; maxDeltaT 5e-5 s plus phase-Courant controls |
| MRF start | constant 17.38 rad/s | reaches 17.38 rad/s at 0.1 s |

The momentum/alpha advection entries in both fvSchemes use upwind, both use
Euler time integration, six alpha subcycles, two outer correctors and two
pressure correctors. Meshes and rotor masks also differ, as already recorded
in the [surface/geometry audit](sc1_three_phase_headspace_audit.md). This is not
a controlled two-field/three-field comparison. A near-1000 L MEML inventory
cannot be credited to the gas split independently of these other changes.

## 1. The older alpha branch can alter liquid inventory by rescaling

The inspected repository `phaseSystemSolve.C` and
[Foundation OF9 source](https://cpp.openfoam.org/v9/phaseSystemSolve_8C_source.html)
use different branches depending on whether `referencePhase` is declared.
Without one, the solver rescales every moving-phase alpha by
`alphaVoid/sumAlphaMoving` after transport. For two moving phases with no
stationary phase this is `alpha_i /= (alpha_gas + alpha_liquid)`.

That operation can change the liquid integral without a corresponding liquid
boundary flux or physical source. The sign and magnitude depend on the
pre-rescaling phase-sum error. The EMM dictionary selects this branch in the
inspected implementation. MEML instead obtains headspaceGas as the complement
of solved phases. **The historical 151 L anomaly has not been attributed
quantitatively to this mechanism**, and the historical executable was not
verified by this configuration collection.

A two-field model can also select a gas reference phase and avoid rescaling
the solved liquid field in this branch. That does not require a third field.
It is a candidate numerical choice, not a conservation guarantee: gas mass,
pressure coupling, reference bounds, alpha sources and species must still
close. The already measured source-floor/implicit deficits in the current
three-field continuation remain unresolved by changing a dictionary label.

## 2. Empty compression also disables this code's resolved capillary force

In repository `phaseSystem.C`, `surfaceTension()` adds the curvature force
only for pairs present in the `interfaceCompression` table. The
[Foundation OF9 implementation](https://cpp.openfoam.org/v9/phaseSystem_8C_source.html)
has the same gate. Consequently, **with the inspected implementation**, EMM's
empty table produces no resolved gas–liquid curvature-force contribution from
this function, despite defining `sigmaLiq=0.07`. Surface tension can still
enter dispersed drag/breakup/coalescence closures; do not describe all uses
of sigma as disabled.

MEML explicitly includes headspaceGas/liquid and therefore activates that
force path as well as interface sharpening. This is a concrete model-setting
difference relevant to the wall/free-surface problem, not proof that capillarity
alone explains a tank-scale elevation or that changing one coefficient fixes
it. Physical swirling flow, gravity, turbulence, mesh and boundary treatment
also control surface shape.

There is **no active liquid contact-angle wall condition** in either saved
9.9 s alphaL field: both specify zeroGradient. Do not diagnose an excessive
prescribed contact angle or introduce a fitted angle as the remedy. A broader
wall/interfacial treatment error remains possible.

For a later controlled numerical test, a pair present with coefficient zero
would retain the inspected force gate while zeroing its compression velocity;
this allows force activation and sharpening to be distinguished. It is not
an approved reactor setting. In a shared gas field, global compression and
PBM/drag/area closures must distinguish dispersed gas from a resolved headspace;
blindly compressing all gas is not the proposed final model.

## 3. Outlet and rotor evidence narrows the diagnosis

The liquid outlet is open in EMM but uses slip in MEML, constraining its normal
velocity. Changing that boundary is a separate liquid-retention intervention.
Whether an impermeable-liquid vent is a physically acceptable approximation
depends on its location, liquid carryover and the experiment; it must not be
chosen merely because it keeps volume near 1000 L. Audit actual alpha and
mass fluxes, not only the boundary type.

At the saved EMM outlet, the paired alpha.gas + alpha.liquid values range
from **0.569857 to 1**, with unweighted face mean **0.993574197**; 601 of 1920
faces differ from unity by more than 1e-6. Both fields use inletOutlet with
their respective phase phi. This is a concrete boundary-consistency finding
to examine with actual phase fluxes and correction ordering. It is not an
integrated mass residual or proof of the mechanism that produced 151 L.

Saved EMM liquid outlet Uy ranges from −96.9882 to 0.0327509 m/s, but most of
those faces have negligible/zero liquid alpha. MEML bubble outlet Uy reaches
−111.085 to 120.815 m/s while outlet alphaB is at most 8.48694e-7. These values
flag dilute-phase velocity/coupling checks; they are not credible estimates
of bulk liquid/gas flow. No area-weighted or time-integrated flux follows from
these snapshots. The inspected controlDict loads libalphaLimitVelocity in
MEML, but this bundle contains no fvConstraints file; loading a library alone
does not establish that a particular limiter was configured or active.

The direct MRF rotation of the `walls` patch remains ruled out by zero
zone adjacency. Both liquid `walls` patches use noSlip. Partial `propeller`
adjacency plus noSlip outside the active zone warrants checking actual shaft,
hub and blade treatment, not assuming all faces named propeller are blades.
MEML's two axial windows and the longer shaft patch explain why face-count
fractions alone cannot establish missing blade rotation.

## 4. The profiles do not identify a sharp wall film

In the outer r=0.50–0.53 m band at 9.9 s, EMM's highest downward crossings
span 0.9743–1.1606 m at alphaL=0.9, 1.2825–1.3216 m at 0.5, and
1.3171–1.3759 m at 0.1 (the last threshold has no crossing in one sector).
MEML's corresponding ranges are 1.0711–1.2046, 1.1577–1.2328 and
1.1864–1.3034 m, with all sectors represented.

These are separate ranges of sector averages and may involve different
crossings; subtracting their extreme endpoints is not a film-thickness
measurement. EMM has a substantial mixed wall-side region, while MEML still
has diffuse/multicrossing structure. In MEML's core no sector has a downward
0.9 crossing in the sampled height window. None of this proves connected
foam, literal blade wetting or a validated vortex.

## 5. Initial fill remains the next small evidence gap

Neither archived root listing contains 0, orig0, 0.org or 0.orig. The saved
setFields dictionaries are available:

- EMM: alphaL=0 by default; alphaL=0.99 in the box ending at y=1.016 m.
- MEML: alphaL=1e-6 by default; 0.999 in the same box; a later overlapping
  box y=1.0124–1.0249 m overrides it with alphaL=0.5814257.

Those values and a nominal level are insufficient to know the initialized
volume without integrating the actual mesh. The new
`estimate_sc1_prescribed_fill.py` applies those exact, hash-verified archived
rules to the existing 9.9 s VTK meshes. It does not call setFields or change
any field. It reports a VTK-center estimate plus conservative box-selection
bounds because VTK centers need not equal OpenFOAM finite-volume centroids.
The bounds assume valid centroids inside cell bounding boxes; VTK volume
integration uncertainty is separate. Both actual archived dictionaries parse
successfully, and six focused tests cover selection order, full 3-D selection,
uncertainty bounds, invalid inputs and VTK hash rejection.

This is a **prescribed-fill estimate, not recovered t=0 evidence**. It tests
whether the archived fill instructions could themselves account for the high
inventory. It cannot prove those instructions were actually run or replace
the historical flux/source budget. The different nonzero liquid feeds also
mean that a perfectly constant liquid inventory is not automatically the
correct physical expectation.

Run once in the existing allocation (no OpenFOAM environment/build needed):

```bash
REPO=/scratch/nwimer/BioReactorDesign_bubble_audit
PY=/scratch/nwimer/.conda-envs/of-post/bin/python
SC1_FILLROOT="/scratch/nwimer/sc1_prescribed_fill_$(date -u +%Y%m%dT%H%M%SZ)"

if [[ -z "${SLURM_JOB_ID:-}" || ! -x "$PY" ]]; then
    echo "Stopped: active allocation or Python missing."
elif [[ "$(git -C "$REPO" branch --show-current)" != "codex/sc1-headspace-bubble-audit-20261029" ]]; then
    echo "Stopped: repository is on a different branch."
elif git -C "$REPO" diff --quiet \
    && git -C "$REPO" diff --cached --quiet \
    && git -C "$REPO" pull --ff-only origin codex/sc1-headspace-bubble-audit-20261029; then
    PYTHONDONTWRITEBYTECODE=1 "$PY" "$REPO/scripts/diagnostics/estimate_sc1_prescribed_fill.py" \
        --report /scratch/nwimer/sc1_wall_preflight_20261110T212232Z/report.json \
        --out "$SC1_FILLROOT"
    echo "Prescribed-fill check status: $?"
else
    echo "Stopped before analysis. Preserve existing files and paste the output."
fi
```

## Implementation direction after the fill check

Use a common verified mesh, measured liquid fill and physically justified
vent/MRF/shaft treatment. Start a controlled two-field numerical baseline
with a gas reference phase as a candidate, and add a phase-count-independent
ledger for transport, source splitting, rescaling/reference reconstruction,
clamping and boundary fluxes. First test closed/source-free conservation and
known sources on small cases before a reactor-scale spin-up. This is a new
isolated implementation; the current three-field audit executable must not
be launched on renamed two-field inputs.

Then validate an explicitly resolved free surface without PBM/biology as a
diagnostic, followed by controlled sparging/diameter and PBM/oxygen transfer.
No-sparging tests do not validate the final aerated experiment. Keep a third
gas field only if a demonstrated dispersed/resolved-gas requirement needs it.
Do not simply restore the old EMM case, add arbitrary foam area, or tune OTR
to biological consumption. The intentional four-hour biological updates and
current trigger/hold policy remain unchanged; archived EMM/MEML biological
settings are historical differences, not instructions to change that policy.

No solver or production setting is changed by this review. The only new
executable is the read-only prescribed-fill estimate. A corrected two-field
solver/case has not yet been built or validated.
