# External OpenFOAM benchmarks and literature for SC1

Reviewed 2026-11-10. This research supplements the
[reliability investigation](SC1_RELIABILITY_INVESTIGATION_2026-11-06.md)
and [current source probe](sc1_alpha_source_probe_run_2026-11-10.md).
The user reports that the probe is running; no completion output has been
received. No solver, case, library or biological schedule is changed here.

## Assessment

Useful code and experimental precedents exist. The most direct code comparison
is upstream BiRD's OpenFOAM 13 stirred-tank implementation. The closest
architecture reference for simultaneous headspace interfaces and dispersed
bubbles is HZDR's morphology-adaptive multifield framework. Published
stirred-tank PBM studies provide hydrodynamic and mass-transfer benchmarks.
These are complementary references, not an identified, validated replacement
for the full SC1 three-field/oxygen/PBM/biology system.

SC1's three numerical fields represent liquid, dispersed gas and headspace gas.
A paper describing gas–liquid–solid flow is not automatically a close match.
Similarly, a two-field reactor with a degassing outlet does not establish the
accuracy of an explicit headspace, gas conversion or foam representation.

## 1. Upstream BiRD: closest inspectable implementation

The [official documentation](https://natlabrockies.github.io/BioReactorDesign/)
documents OpenFOAM 9 and 13 support. The upstream repository was inspected at
commit `faf54ff57534a3f65e29fb51e1bff05dd287c6e4`, not an unpinned future main.

- [MRF stirred-tank tutorial](https://github.com/NatLabRockies/BioReactorDesign/tree/faf54ff57534a3f65e29fb51e1bff05dd287c6e4/tutorial_cases/OF13/STR_MRF)
- [Rotating-mesh stirred-tank tutorial](https://github.com/NatLabRockies/BioReactorDesign/tree/faf54ff57534a3f65e29fb51e1bff05dd287c6e4/tutorial_cases/OF13/STR_rotatingMesh)
- [Alpha solver source](https://github.com/NatLabRockies/BioReactorDesign/blob/faf54ff57534a3f65e29fb51e1bff05dd287c6e4/OFsolvers/OF13/birdmultiphaseEuler/phaseSystem/phaseSystem/phaseSystemSolve.C)

Both README files describe a 500 mL oxygen/water example at 1 VVM and 100 rpm,
including mass transfer and PBM. Both inspected `phaseProperties` files declare
gas and liquid, with liquid as reference phase. The MRF `fvSolution` enables
semi-implicit MULES. These are runnable tutorial precedents; the inspected
files do not establish experimental validation or SC1's three-field transfer.

**Proposed use:** compare source assembly, subcycling, flux accumulation,
reference-phase reconstruction and pressure coupling against our exact OF9
implementation. An isolated upstream tutorial can later provide a regression
baseline. Its numerical settings and physical coefficients are not SC1 inputs.
Compatibility with the custom putida solver and its biological schedule remains
unverified; no migration is authorized by a successful tutorial alone.

## 2. Newer Foundation alpha numerics: relevant, not a demonstrated fix

The [OF9 source](https://cpp.openfoam.org/v9/phaseSystemSolve_8C_source.html)
and [OF13 source](https://cpp.openfoam.org/v13/phaseSystemSolve_8C_source.html)
were inspected. OF13 offers an implicit predictor and limited flux corrections,
but retains sign-dependent Su/Sp splitting with denominators bounded by
`vDotResidualAlpha`. Our OF9/BiRD source uses the analogous fixed `1e-4` floor.
The inspected upstream BiRD OF13 file retains that source structure too.

**Inference:** an upgrade does not by itself prove removal of the observed
source deficit. Use the running probe to identify the operative bounds and
coupling issue, then compare minimal cases across formulations. A source-floor
change or copied MULES switch is not an established conservative correction.

## 3. HZDR: resolved interfaces and dispersed bubbles together

[Meller, Schlegel and Lucas, IJNMF 93, 748–773 (2021)](https://doi.org/10.1002/fld.4907)
develop a multifield OpenFOAM formulation with morphology-dependent coupling.
It includes rising-bubble verification and coexistence of large interfaces
with dispersed microbubbles. The paper's conclusion explicitly leaves transfer
between dispersed and continuous gas, and froth, as future work. It therefore
does not validate our entrainment/deposition or foam closures.

The older 2022 addon cited in our initial investigation has newer archived
[code](https://rodare.hzdr.de/record/4741) and
[case collections](https://rodare.hzdr.de/record/4740), released June 29, 2026,
version `13-s.1-hzdr.2` / `13-s.1-hzdr.2-cases.1`. Their records include
polydisperse bubbly-flow and MultiMorph interface cases. Archive contents were
not compiled or run in this review; the advertised capability is not a verified
SC1-compatible implementation. These releases target Foundation OpenFOAM.

**Proposed use:** borrow the verification structure and inspect interface-normal
velocity coupling, resolved-interface drag and pressure interpolation before
adapting a closure. Confirm separately what the newer implementation supports
for morphology transfer, species transport, PBM and foam. Do not infer those
capabilities from the older hydrodynamic paper.

## 4. Stirred-tank PBM and oxygen-transfer validation

[Seidel and Eibl, Processes 9, 1185 (2021)](https://doi.org/10.3390/pr9071185)
study OpenFOAM 7 Euler–Euler/PBM predictions in a 2 L working-volume stirred
bioreactor, validated against measured kLa. Bubble diameter and the kL model
strongly affect predictions; PBM improves agreement for the tested systems.
The study also documents mesh-quality checks and discretization assessment.
Its data-availability statement offers data on request, not a verified public
case archive.

**Proposed use:** validate hydrodynamics, bubble size/interfacial area and kL
independently of biological oxygen demand. A matching kLa alone can conceal
compensating errors in bubble inventory and transfer coefficient. This last
point is our validation inference, not a numerical fix supplied by the paper.

[Wang et al., Chemical Engineering Journal (2025), 161277](https://doi.org/10.1016/j.cej.2025.161277)
present a twoPhaseEulerFoam-based CFD–PBE stirred-tank framework, experimental
comparisons and a degassing boundary. The publisher's indexed abstract and
introduction were accessible; direct full-page retrieval returned 403. The
authors identify problems with applying dispersed-bubble kernels in a freeboard
where gas fraction approaches unity. Public runnable code was not verified.

**Proposed use:** a separate liquid-domain/degassing control could isolate
submerged hydrodynamics and PBM from headspace-interface errors. This is a
proposed reduced benchmark, not a substitute for SC1's eventual explicit
headspace validation. No such case has been prepared or launched.

## 5. Free-surface experimental evidence

[de Lamotte et al., Chemical Engineering Science 172, 125–142 (2017)](https://discovery.ucl.ac.uk/id/eprint/10064111/)
experimentally examine surface oxygen transfer in a baffled 20 L vessel with
two Rushton turbines, using PIV and PLIF. Their repository abstract reports
coupling to circulation and persistent vertical dissolved-gas gradients.
This review used that abstract, not the full accepted manuscript. It is an
experimental reference, not an OpenFOAM code release or a foam calibration.

**Proposed use:** validate surface transfer together with circulation and local
oxygen profiles. A single total OTR or a prescribed film thickness is not
sufficient evidence. None of the sources reviewed establishes an SC1 foam-area
multiplier; broth-specific foam evidence is still needed.

## How this changes the next implementation decision

1. Finish and interpret the current observer run, including legacy regression.
   Preserve it and all completed runs; do not launch a second probe.
2. Compare the exact source/flux/pressure sequence against the pinned upstream
   BiRD and Foundation implementations. Build an isolated minimal reproduction
   of the observed source/bounds problem before selecting a correction.
3. Require equal-and-opposite phase/species transfer, bounded fractions,
   phase-sum/flux consistency and PB moments in those small tests. Passing a
   published hydrodynamic benchmark without mass exchange does not test this.
4. Use a published stirred-tank/PBM configuration or the upstream tutorial as
   a separate reference once the numerical correction has evidence. If useful,
   add a reduced degassing control to separate submerged-flow and free-surface
   errors; do not replace the ongoing SC1 experiment with it.
5. Validate the explicit interface and gas conversion against appropriate
   benchmarks, then gas holdup, bubble size, mixing and oxygen transfer against
   SC1 evidence. Establish the nonreacting 9.9 s checkpoint before biological
   comparisons. Retain the four-hour biological update schedule throughout.

This review identifies concrete code and experiments that can reduce trial and
error. It does not claim that a published case, a new OpenFOAM version or an
uncalibrated closure already solves SC1. Changes in this commit are documentation
only; no new CFD, compilation or solver tests were performed.
