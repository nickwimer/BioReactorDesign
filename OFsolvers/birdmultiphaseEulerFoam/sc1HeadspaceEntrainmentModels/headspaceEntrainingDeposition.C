/*---------------------------------------------------------------------------*\
  SCRATCH-ONLY SC1 deposition plus spatially limited surface entrainment.

  This is a *mixture* phaseTransferModel: transfer of headspace gas to
  bubbleGas is balanced by removal from headspaceGas. The original
  deposition correlation remains unchanged at efficiencyEntr=0.
  BiRD transports donor mixture species using addDmdtYfs().
\*---------------------------------------------------------------------------*/

#include "headspaceEntrainingDeposition.H"

#include "phasePair.H"
#include "phaseSystem.H"
#include "fvcGrad.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace phaseTransferModels
{
    defineTypeNameAndDebug(headspaceEntrainingDeposition, 0);
    addToRunTimeSelectionTable
    (
        phaseTransferModel,
        headspaceEntrainingDeposition,
        dictionary
    );
}
}


Foam::phaseTransferModels::headspaceEntrainingDeposition::
headspaceEntrainingDeposition
(
    const dictionary& dict,
    const phasePair& pair
)
:
    phaseTransferModel(dict, pair),
    dropletName_(dict.lookup("droplet")),
    surfaceName_(dict.lookup("surface")),
    liquidName_(dict.lookupOrDefault<word>("liquid", "liquid")),
    depositionEfficiency_(dict.lookup<scalar>("efficiency")),
    entrainmentEfficiency_
    (
        dict.lookupOrDefault<scalar>("entrainmentEfficiency", 0)
    ),
    yMin_(dict.lookupOrDefault<scalar>("entrainmentYMin", 0.68)),
    yMax_(dict.lookupOrDefault<scalar>("entrainmentYMax", 0.96)),
    rMin_(dict.lookupOrDefault<scalar>("entrainmentRMin", 0.05)),
    rMax_(dict.lookupOrDefault<scalar>("entrainmentRMax", 0.30)),
    axisX_(dict.lookupOrDefault<scalar>("entrainmentAxisX", 0)),
    axisZ_(dict.lookupOrDefault<scalar>("entrainmentAxisZ", 0)),
    minAlpha_(dict.lookupOrDefault<scalar>("minInterfaceAlpha", 0.02)),
    fullAlpha_(dict.lookupOrDefault<scalar>("fullInterfaceAlpha", 0.10)),
    maxBubbleAlpha_(dict.lookupOrDefault<scalar>("maxBubbleAlpha", 0.10)),
    maxSlipSpeed_(dict.lookupOrDefault<scalar>("maxSlipSpeed", 5)),
    donorTimeScale_
    (
        dict.lookupOrDefault<scalar>("donorTimeScale", 0.05)
    ),
    maxGlobalEntrainedKgS_
    (
        dict.lookupOrDefault<scalar>("maxGlobalEntrainedKgS", 0.002)
    ),
    requireDownwardLiquid_
    (
        dict.lookupOrDefault<bool>("requireDownwardLiquid", true)
    ),
    entrainmentVelocityCriterion_
    (
        dict.lookupOrDefault<word>
        (
            "entrainmentVelocityCriterion", "legacyDownwardLiquid"
        )
    ),
    reportRate_(dict.lookupOrDefault<bool>("reportRate", true))
{
    if
    (
        dropletName_ == surfaceName_
     || (dropletName_ != pair.first() && dropletName_ != pair.second())
     || (surfaceName_ != pair.first() && surfaceName_ != pair.second())
     || depositionEfficiency_ < 0
     || entrainmentEfficiency_ < 0
     || yMax_ <= yMin_
     || rMin_ < 0
     || rMax_ <= rMin_
     || minAlpha_ < 0
     || fullAlpha_ <= minAlpha_
     || fullAlpha_ >= 0.5
     || maxBubbleAlpha_ < 0
     || maxBubbleAlpha_ > 1
     || maxSlipSpeed_ <= 0
     || donorTimeScale_ <= 0
     || maxGlobalEntrainedKgS_ <= 0
    )
    {
        FatalIOErrorInFunction(dict)
            << "Invalid headspaceEntrainingDeposition phase pair or "
            << "entrainment parameters."
            << exit(FatalIOError);
    }

    if
    (
        entrainmentVelocityCriterion_ != "legacyDownwardLiquid"
     && entrainmentVelocityCriterion_ != "normalGasIntoLiquid"
    )
    {
        FatalIOErrorInFunction(dict)
            << "Unknown entrainmentVelocityCriterion "
            << entrainmentVelocityCriterion_ << ". Valid options are "
            << "legacyDownwardLiquid and normalGasIntoLiquid."
            << exit(FatalIOError);
    }

    if
    (
        entrainmentVelocityCriterion_ == "normalGasIntoLiquid"
     && requireDownwardLiquid_
    )
    {
        FatalIOErrorInFunction(dict)
            << "Set requireDownwardLiquid false with "
            << "entrainmentVelocityCriterion normalGasIntoLiquid. "
            << "Both mutually exclusive velocity screens must not be active."
            << exit(FatalIOError);
    }

    if
    (
        !pair.phase1().fluid().phases().found(liquidName_)
    )
    {
        FatalIOErrorInFunction(dict)
            << "Liquid phase " << liquidName_ << " not found"
            << exit(FatalIOError);
    }
}


Foam::phaseTransferModels::headspaceEntrainingDeposition::
~headspaceEntrainingDeposition()
{}


bool Foam::phaseTransferModels::headspaceEntrainingDeposition::mixture() const
{
    return true;
}


Foam::scalar
Foam::phaseTransferModels::headspaceEntrainingDeposition::gate
(
    const scalar alpha
) const
{
    return min
    (
        max((alpha - minAlpha_)/(fullAlpha_ - minAlpha_), scalar(0)),
        scalar(1)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::headspaceEntrainingDeposition::dmdtf() const
{
    // Reproduce the original BiRD deposition correlation and sign exactly.
    // Negative bubbleGas dmdtf removes gas from the bubble phase.
    const phaseModel& droplet =
        pair_.phase1().fluid().phases()[dropletName_];
    const phaseModel& surface =
        pair_.phase1().fluid().phases()[surfaceName_];

    const scalar sign =
        dropletName_ == pair_.first() ? scalar(-1) : scalar(1);

    tmp<volScalarField> tNet
    (
        1.5*sign*depositionEfficiency_
       *droplet.rho()*droplet*surface/surface.d()
       *mag(droplet.U() - surface.U())
    );

    // Strict zero-entrainment control reproduces the original
    // deposition physics *before the existing phaseTransfer blending*.
    if (entrainmentEfficiency_ <= SMALL)
    {
        return tNet;
    }

    const phaseModel& liquid =
        pair_.phase1().fluid().phases()[liquidName_];

    const fvMesh& mesh = surface.mesh();
    const volScalarField& alphaH = surface;
    const volScalarField& alphaL = liquid;
    const volScalarField& alphaB = droplet;

    const tmp<volScalarField> tRhoH = surface.rho();
    const volScalarField& rhoH = tRhoH();

    const tmp<volVectorField> tGradH = fvc::grad(alphaH);
    const volVectorField& gradH = tGradH();

    const tmp<volVectorField> tUL = liquid.U();
    const tmp<volVectorField> tUH = surface.U();
    const volVectorField& UL = tUL();
    const volVectorField& UH = tUH();

    const volVectorField& C = mesh.C();
    const scalarField& V = mesh.V();
    scalarField entrain(mesh.nCells(), 0);

    scalar rawIntegral = 0;
    scalar interfaceArea = 0;
    label eligibleCells = 0;

    // Diagnostic gates only; the actual source formula is unchanged.
    // These count sequential losses of interfacial area through each
    // eligibility condition, not model parameter sensitivities.
    scalar wholeDomainMixedArea = 0;
    scalar windowMixedArea = 0;
    scalar windowNoBubbleArea = 0;
    label windowCells = 0;
    label windowMixedCells = 0;
    label windowNoBubbleCells = 0;

    // Diagnostic-only alternative kinematic screenings. They DO NOT
    // change the already validated downflow-gated source expression.
    scalar windowDownwardLiquidArea = 0;
    scalar windowNonDownwardLiquidArea = 0;
    scalar windowNormalGasIntoLiquidArea = 0;
    scalar windowNormalGasAwayFromLiquidArea = 0;
    scalar rawWithoutDownwardGate = 0;
    scalar rawWithNormalInwardGate = 0;
    scalar areaWeightedNormalSpeed = 0;
    scalar areaWeightedLiquidVerticalSpeed = 0;
    label windowDownwardLiquidCells = 0;
    label windowNormalGasIntoLiquidCells = 0;

    forAll(V, celli)
    {
        const scalar h = max(min(alphaH[celli], scalar(1)), scalar(0));
        const scalar l = max(min(alphaL[celli], scalar(1)), scalar(0));
        const scalar interfaceGate = gate(h)*gate(l);
        const scalar area = interfaceGate*mag(gradH[celli]);

        if (interfaceGate > SMALL)
        {
            wholeDomainMixedArea += area*V[celli];
        }

        const vector& c = C[celli];
        const scalar r = sqrt
        (
            sqr(c.x() - axisX_) + sqr(c.z() - axisZ_)
        );

        if
        (
            c.y() < yMin_ || c.y() > yMax_
         || r < rMin_ || r > rMax_
        )
        {
            continue;
        }

        ++windowCells;

        if (interfaceGate <= SMALL)
        {
            continue;
        }

        ++windowMixedCells;
        windowMixedArea += area*V[celli];

        if (alphaB[celli] >= maxBubbleAlpha_)
        {
            continue;
        }

        ++windowNoBubbleCells;
        windowNoBubbleArea += area*V[celli];

        // Interface normal n = grad(alpha.headspaceGas)/|grad| points
        // approximately from liquid towards headspace. Positive
        // (U_liquid-U_headspace).n suggests gas penetration relative
        // to liquid; it is a velocity PROXY, not proof of breakup or
        // a validated entrainment correlation.
        const scalar normalGasIntoLiquidSpeed =
            ((UL[celli] - UH[celli]) & gradH[celli])
           /max(mag(gradH[celli]), scalar(SMALL));

        areaWeightedNormalSpeed +=
            area*normalGasIntoLiquidSpeed*V[celli];
        areaWeightedLiquidVerticalSpeed +=
            area*UL[celli].y()*V[celli];

        if (UL[celli].y() < 0)
        {
            ++windowDownwardLiquidCells;
            windowDownwardLiquidArea += area*V[celli];
        }
        else
        {
            windowNonDownwardLiquidArea += area*V[celli];
        }

        if (normalGasIntoLiquidSpeed > 0)
        {
            ++windowNormalGasIntoLiquidCells;
            windowNormalGasIntoLiquidArea += area*V[celli];
        }
        else
        {
            windowNormalGasAwayFromLiquidArea += area*V[celli];
        }

        // Counterfactual rates, before original downward gate or
        // phaseTransfer blending. NOT added to phase equations.
        const scalar diagnosticSlip =
            min(mag(UL[celli] - UH[celli]), maxSlipSpeed_);
        const scalar diagnosticDonorLimit =
            rhoH[celli]*max(h - minAlpha_, scalar(0))/donorTimeScale_;
        const scalar diagnosticFullRate = min
        (
            rhoH[celli]*entrainmentEfficiency_*area*diagnosticSlip,
            diagnosticDonorLimit
        );
        rawWithoutDownwardGate += diagnosticFullRate*V[celli];

        if (normalGasIntoLiquidSpeed > 0)
        {
            const scalar diagnosticNormalRate = min
            (
                rhoH[celli]*entrainmentEfficiency_*area
               *min(normalGasIntoLiquidSpeed, maxSlipSpeed_),
                diagnosticDonorLimit
            );
            rawWithNormalInwardGate += diagnosticNormalRate*V[celli];
        }

        // Select one mutually exclusive kinematic proxy. The legacy
        // branch is algebraically identical to previous experiments.
        // The new branch requires gas to move toward the liquid
        // RELATIVE to the liquid, in the local phase-normal direction.
        // This relative-speed diagnostic does not prove interfacial
        // breakup or detachment of new dispersed bubbles.
        const bool normalMode =
            entrainmentVelocityCriterion_ == "normalGasIntoLiquid";

        if (normalMode)
        {
            if (normalGasIntoLiquidSpeed <= 0)
            {
                continue;
            }
        }
        else if (requireDownwardLiquid_ && UL[celli].y() >= 0)
        {
            continue;
        }

        const scalar slip = normalMode
          ? min(normalGasIntoLiquidSpeed, maxSlipSpeed_)
          : min(mag(UL[celli] - UH[celli]), maxSlipSpeed_);

        // kg/(m3 s): rho_H * C_e * interfacial area density * slip.
        // This is a bounded experimental closure, not a validated
        // entrainment correlation. Interfacial area and downflow should
        // be benchmarked against the experimental video.
        const scalar proposed =
            rhoH[celli]*entrainmentEfficiency_*area*slip;

        // Additional local donor-capacity guard. This is NOT a
        // numerical-time-step bound; it is a conservative rate cap.
        const scalar donorLimit =
            rhoH[celli]*max(h - minAlpha_, scalar(0))/donorTimeScale_;

        entrain[celli] = min(proposed, donorLimit);
        rawIntegral += entrain[celli]*V[celli];
        interfaceArea += area*V[celli];
        ++eligibleCells;
    }

    reduce(rawIntegral, sumOp<scalar>());
    reduce(interfaceArea, sumOp<scalar>());
    reduce(eligibleCells, sumOp<label>());
    reduce(wholeDomainMixedArea, sumOp<scalar>());
    reduce(windowMixedArea, sumOp<scalar>());
    reduce(windowNoBubbleArea, sumOp<scalar>());
    reduce(windowCells, sumOp<label>());
    reduce(windowMixedCells, sumOp<label>());
    reduce(windowNoBubbleCells, sumOp<label>());
    reduce(windowDownwardLiquidArea, sumOp<scalar>());
    reduce(windowNonDownwardLiquidArea, sumOp<scalar>());
    reduce(windowNormalGasIntoLiquidArea, sumOp<scalar>());
    reduce(windowNormalGasAwayFromLiquidArea, sumOp<scalar>());
    reduce(rawWithoutDownwardGate, sumOp<scalar>());
    reduce(rawWithNormalInwardGate, sumOp<scalar>());
    reduce(areaWeightedNormalSpeed, sumOp<scalar>());
    reduce(areaWeightedLiquidVerticalSpeed, sumOp<scalar>());
    reduce(windowDownwardLiquidCells, sumOp<label>());
    reduce(windowNormalGasIntoLiquidCells, sumOp<label>());

    const scalar globalScale =
        rawIntegral > SMALL
      ? min(scalar(1), maxGlobalEntrainedKgS_/rawIntegral)
      : scalar(0);

    volScalarField& net = tNet.ref();

    forAll(V, celli)
    {
        // The droplet sign is negative when bubbleGas is phase1.
        // Deducting a positive entrainment term from a negative
        // deposition magnitude therefore CREATES bubbleGas mass.
        // Both phase densities/species are handled by phaseSystem.
        net[celli] -= sign*globalScale*entrain[celli];
    }

    if (reportRate_ && mesh.time().writeTime())
    {
        static scalar lastReport = -GREAT;

        if (mag(mesh.time().value() - lastReport) > SMALL)
        {
            lastReport = mesh.time().value();

            if (Pstream::master())
            {
                Info<< "HEADSPACE_ENTRAINMENT_CANDIDATE"
                    << " t_s=" << mesh.time().value()
                    << " velocityCriterion="
                    << entrainmentVelocityCriterion_
                    << " interfaceArea_m2=" << interfaceArea
                    << " eligibleCells=" << eligibleCells
                    << " domainMixedArea_m2=" << wholeDomainMixedArea
                    << " windowCells=" << windowCells
                    << " windowMixedCells=" << windowMixedCells
                    << " windowMixedArea_m2=" << windowMixedArea
                    << " windowLowBubbleCells=" << windowNoBubbleCells
                    << " windowLowBubbleArea_m2=" << windowNoBubbleArea
                    << " downwardLiquidCells=" << windowDownwardLiquidCells
                    << " downwardLiquidArea_m2=" << windowDownwardLiquidArea
                    << " nonDownwardLiquidArea_m2="
                    << windowNonDownwardLiquidArea
                    << " normalGasIntoLiquidCells="
                    << windowNormalGasIntoLiquidCells
                    << " normalGasIntoLiquidArea_m2="
                    << windowNormalGasIntoLiquidArea
                    << " normalGasAwayArea_m2="
                    << windowNormalGasAwayFromLiquidArea
                    << " meanNormalGasIntoLiquidSpeed_m_s="
                    << areaWeightedNormalSpeed/max(windowNoBubbleArea, SMALL)
                    << " meanLiquidVerticalSpeed_m_s="
                    << areaWeightedLiquidVerticalSpeed
                      /max(windowNoBubbleArea, SMALL)
                    << " diagnosticRawNoDownGate_kg_s="
                    << rawWithoutDownwardGate
                    << " diagnosticRawNormalInward_kg_s="
                    << rawWithNormalInwardGate
                    << " rawEntrainment_kg_s=" << rawIntegral
                    << " capScale=" << globalScale
                    << " cappedEntrainment_kg_s="
                    << rawIntegral*globalScale
                    << " limit_kg_s=" << maxGlobalEntrainedKgS_
                    << " isUnblended=1"
                    << endl;
            }
        }
    }

    return tNet;
}

// ************************************************************************* //
