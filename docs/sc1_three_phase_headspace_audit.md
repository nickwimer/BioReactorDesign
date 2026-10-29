# SC1 three-phase headspace pressure and bubble-mass audit

**Scope:** Diagnostic-only continuation of the MEML SC1 case using the existing
`putidaFoam` source with `PIMPLE/bubbleMassAudit true`. This is **not** a
change to deposition, oxygen uptake, sparger, phase-fraction transport or pressure
boundary conditions.

## 1. Compare the real case boundary conditions before continuing

The production Kestrel cases are not part of this repository. Compare their
**saved** 30-second patch definitions; do not assume that older templates
match the simulations. Pay particular attention to the top outlet/vent patch:

- `p_rgh`: `prghPressure`, `fixedFluxPressure`, `totalPressure`,
  `fixedValue` and their values, as applicable;
- `p`: physical static-pressure boundary (generally `calculated`
  with the `p_rgh` pressure solution);
- `U.headspaceGas`, `U.bubbleGas`: outflow and backflow behavior;
- `alpha.headspaceGas`, `alpha.bubbleGas`, `alpha.liquid`: backflow
  composition and whether a phase is artificially fixed at the outlet;
- `O2.headspaceGas`, `O2.bubbleGas`, `O2.liquid` as relevant;
- pressure reference in `system/fvSolution` and any `pMin`, `pMax`
  settings.

For OpenFOAM, `p_rgh = p - rho*gh`; a static-pressure outlet must be
checked in combination with local hydrostatic head, not in isolation.
Do not confuse an impeller-depressed *physical interface* with a fixed top
boundary. Also inspect the outlet patch identity in `constant/polyMesh/boundary`.

A minimally invasive way to inspect a reconstructed snapshot:

```bash
for f in p p_rgh U.headspaceGas U.bubbleGas \
         alpha.headspaceGas alpha.bubbleGas alpha.liquid; do
    [ -f "30/$f" ] || continue
    echo "=== $f"
    foamDictionary "30/$f" -entry boundaryField -value
done
```

For legacy EMM and MAX-OUR SC1, substitute `U.gas` and `alpha.gas`.
Take care that headspace and bubbleGas are distinct phases only in the new case.

## 2. Build the diagnostic executable on Kestrel

This branch changes `putidaFoam.C` only to include
`bubbleMassAudit.H` **when explicitly enabled**, at an ordinary
write time, after the final PIMPLE iteration. Copy the working solver source
into a separate worktree (or build a dedicated executable); do not overwrite
the production `putidaFoamMask01` binary.

Suggested configuration in **only the cloned diagnostic case**:

```
PIMPLE
{
    // Preserve existing settings here.
    bubbleMassAudit true;
}
```

Start from the 30-second checkpoint and continue over a short, explicitly
limited time interval with the same hydrodynamic boundary conditions and
`fvModels`. Set a sufficiently frequent `writeInterval` so that the
finite-difference inventory derivative is meaningful.

**Do not infer that this diagnostic compiles against the separate RHEL9
`putidaFoamMask01` source without building and testing on Kestrel.**
The GitHub code is the baseline implementation, not proof of the final
RHEL9 binary's source history.

## 3. Recorded quantities

`BUBBLE_MASS_AUDIT` contains:

- `inventory_kg` — instantaneous `sum(alpha.bubbleGas*rho.bubbleGas*V)`;
- `injection_kg_s` — recomputed total instantaneous `fvModels`
  density-equation source for bubbleGas, expected to be the 15 sparger
  sources if no other density source is configured;
- `depositionLoss_kg_s` — positive loss from bubbleGas based on the
  **signed phase-system `dmdtf`** for `(bubbleGas in headspaceGas)`;
- `oxygenDissolutionLoss_kg_s` — positive bubble-phase loss computed from
  the signed `(bubbleGas and liquid)` mass transfer;
- `bubbleBoundaryOut_kg_s` — sum of outward `alphaRhoPhi.bubbleGas`
  across non-coupled physical boundary faces;
- `headspaceBoundaryOut_kg_s` — corresponding headspace phase mass flux;
- `inventoryRate_kg_s` — difference in inventory between audit write
  times; and
- `endpointClosure_kg_s` — approximate residual between the inventory
  derivative and the contemporaneous source/flux terms.

`HEADSPACE_PRESSURE_AUDIT` logs *every* non-coupled patch:
pressure and `p_rgh` area means, their runtime boundary-condition types,
the headspace-gas velocity condition, and bubble/headspace mass outflows.

**Important limitations.** The audit samples endpoint source terms, not
time-averaged or time-integrated values. The closure is not an exact discrete
finite-volume residual, especially with moving meshes, subcycling, source
updates at PIMPLE iterations, nonconservative boundary discretizations,
species normalization and phase clipping. A large residual identifies a need
for a finer, equation-level budget; it does not by itself prove an OpenFOAM
mass-conservation bug. In particular, the snapshot-derived **32 L/s**
deposition is a density-normalized rate, whereas the audit uses kg/s.
These must not be compared as though they have the same units.

## 4. Safety and validation

- Retain the original 30-second checkpoint and all production case files.
- Build/run only in a cloned scratch case and under a separate executable
  name.
- Keep the 200-second frozen OUR update interval and all source definitions
  unchanged.
- First check no diagnostic output when `bubbleMassAudit` is absent.
- Verify MPI reductions by comparing bubble inventory at a checkpoint with
  the existing OpenFOAM `volFieldValue` result.
- Check the numerical budget at multiple write times, not from one instant.

The workflow's primary purpose is to determine whether the very large
bubble-to-headspace deposition predicted by the source law is consistent with
the actual gas-phase inventory and boundary flux; a separate pressure/BC
comparison tests whether a top-outlet condition is influencing that result.
