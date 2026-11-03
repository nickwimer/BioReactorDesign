/*---------------------------------------------------------------------------*\
  SC1 experimental resolved headspace-liquid interfacial mass transfer.
  Deliberately separate from the dispersed bubble-liquid Higbie model.
\*---------------------------------------------------------------------------*/

#include "resolvedHeadspaceFilm.H"

#include "phasePair.H"
#include "phaseModel.H"
#include "fvcGrad.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace diffusiveMassTransferModels
{
    defineTypeNameAndDebug(resolvedHeadspaceFilm, 0);
    addToRunTimeSelectionTable
    (
        diffusiveMassTransferModel,
        resolvedHeadspaceFilm,
        dictionary
    );
}
}


Foam::diffusiveMassTransferModels::resolvedHeadspaceFilm::
resolvedHeadspaceFilm
(
    const dictionary& dict,
    const phasePair& pair
)
:
    diffusiveMassTransferModel(dict, pair),
    headspacePhase_
    (
        dict.lookupOrDefault<word>("headspacePhase", "headspaceGas")
    ),
    liquidPhase_
    (
        dict.lookupOrDefault<word>("liquidPhase", "liquid")
    ),
    bubblePhase_
    (
        dict.lookupOrDefault<word>("bubblePhase", "bubbleGas")
    ),
    filmThickness_("filmThickness", dimLength, dict),
    foamDiameter_("foamDiameter", dimLength, dict),
    minInterfaceAlpha_
    (
        dict.lookupOrDefault<scalar>("minInterfaceAlpha", 0.02)
    ),
    fullInterfaceAlpha_
    (
        dict.lookupOrDefault<scalar>("fullInterfaceAlpha", 0.10)
    ),
    maxBubbleAlpha_
    (
        dict.lookupOrDefault<scalar>("maxBubbleAlpha", 0.10)
    ),
    foamAreaMultiplier_
    (
        dict.lookupOrDefault<scalar>("foamAreaMultiplier", 0.0)
    ),
    foamYMin_(dict.lookupOrDefault<scalar>("foamYMin", 0.65)),
    foamYMax_(dict.lookupOrDefault<scalar>("foamYMax", 1.05)),
    foamRMax_(dict.lookupOrDefault<scalar>("foamRMax", 0.30)),
    foamAxisX_(dict.lookupOrDefault<scalar>("foamAxisX", 0.0)),
    foamAxisZ_(dict.lookupOrDefault<scalar>("foamAxisZ", 0.0)),
    maxAreaDensity_
    (
        dict.lookupOrDefault<scalar>("maxAreaDensity", 2000.0)
    ),
    reportArea_(dict.lookupOrDefault<bool>("reportArea", false))
{
    if
    (
        filmThickness_.value() <= 0
     || foamDiameter_.value() <= 0
     || minInterfaceAlpha_ < 0
     || fullInterfaceAlpha_ <= minInterfaceAlpha_
     || fullInterfaceAlpha_ >= 0.5
     || maxBubbleAlpha_ < 0
     || maxBubbleAlpha_ > 1
     || foamAreaMultiplier_ < 0
     || foamYMax_ <= foamYMin_
     || foamRMax_ <= 0
     || maxAreaDensity_ <= 0
    )
    {
        FatalErrorInFunction
            << "Invalid resolvedHeadspaceFilm parameters. "
            << "Check filmThickness, foamDiameter, area limits and gates."
            << exit(FatalError);
    }

    const phaseModel& declaredHeadspace =
        pair.phase1().fluid().phases()[headspacePhase_];
    const phaseModel& declaredLiquid =
        pair.phase1().fluid().phases()[liquidPhase_];

    if
    (
        !pair.contains(declaredHeadspace)
     || !pair.contains(declaredLiquid)
    )
    {
        FatalErrorInFunction
            << "resolvedHeadspaceFilm requires the headspace-liquid pair; "
            << "got " << pair
            << exit(FatalError);
    }
}


Foam::diffusiveMassTransferModels::resolvedHeadspaceFilm::
~resolvedHeadspaceFilm()
{}


Foam::scalar
Foam::diffusiveMassTransferModels::resolvedHeadspaceFilm::
interfaceGate(const scalar alpha) const
{
    return min
    (
        max
        (
            (alpha - minInterfaceAlpha_)
           /(fullInterfaceAlpha_ - minInterfaceAlpha_),
            scalar(0)
        ),
        scalar(1)
    );
}


Foam::tmp<Foam::volScalarField>
Foam::diffusiveMassTransferModels::resolvedHeadspaceFilm::K() const
{
    const phaseModel& head =
        pair_.phase1().fluid().phases()[headspacePhase_];
    const phaseModel& liquid =
        pair_.phase1().fluid().phases()[liquidPhase_];
    const phaseModel& bubble =
        pair_.phase1().fluid().phases()[bubblePhase_];

    const volScalarField& alphaH = head;
    const volScalarField& alphaL = liquid;
    const volScalarField& alphaB = bubble;

    const fvMesh& mesh = head.mesh();

    // Geometric interface area density (surface area per cell volume):
    // a_resolved = |grad(alpha.headspaceGas)| [m^-1].
    const tmp<volVectorField> gradAlphaHTmp = fvc::grad(alphaH);
    const volVectorField& gradAlphaH = gradAlphaHTmp();

    tmp<volScalarField> tK
    (
        volScalarField::New
        (
            IOobject::groupName("resolvedHeadspaceFilm:K", pair_.name()),
            mesh,
            dimensionedScalar(dimK, 0)
        )
    );
    volScalarField& transferK = tK.ref();

    const volVectorField& centres = mesh.C();
    const scalarField& cellVolumes = mesh.V();

    scalar integralGeometricArea = 0;
    scalar integralFoamArea = 0;
    scalar integralEffectiveArea = 0;

    forAll(cellVolumes, celli)
    {
        const scalar h = min(max(alphaH[celli], scalar(0)), scalar(1));
        const scalar l = min(max(alphaL[celli], scalar(0)), scalar(1));
        const scalar b = min(max(alphaB[celli], scalar(0)), scalar(1));

        // Smoothly suppress transfer in bulk gas and bulk liquid.
        // Also exclude bubble-rich cells to avoid counting bubbleGas area
        // again in this headspaceGas-liquid mechanism.
        if (b >= maxBubbleAlpha_)
        {
            transferK[celli] = 0;
            continue;
        }

        const scalar gate = interfaceGate(h)*interfaceGate(l);

        if (gate <= SMALL)
        {
            transferK[celli] = 0;
            continue;
        }

        // For a planar 0-to-1 interface, the two smooth alpha gates
        // would otherwise discard an area fraction minAlpha+fullAlpha.
        // Normalize the resolved area to integrate to the actual plane
        // area (when alphaH+alphaL=1 and there is no bubble phase).
        const scalar geometricAreaDensity =
            gate*mag(gradAlphaH[celli])
           /(scalar(1) - minInterfaceAlpha_ - fullInterfaceAlpha_);

        // Optional extra sub-grid gas-liquid contact inside the bounded
        // impeller foam zone; does NOT operate in the entire headspace.
        // a_foam=6*alphaH*alphaL/d_foam is a model ASSUMPTION, not a
        // direct measurement of the experimental froth area.
        scalar foamAreaDensity = 0;

        if (foamAreaMultiplier_ > 0)
        {
            const vector& c = centres[celli];
            const scalar r = sqrt
            (
                sqr(c.x() - foamAxisX_) + sqr(c.z() - foamAxisZ_)
            );

            if
            (
                c.y() >= foamYMin_
             && c.y() <= foamYMax_
             && r <= foamRMax_
            )
            {
                foamAreaDensity =
                    gate*foamAreaMultiplier_*6*h*l
                   /foamDiameter_.value();
            }
        }

        // MAX avoids double counting the same resolved interface as
        // both a large-area interface and a sub-grid foamy interface.
        const scalar effectiveAreaDensity = min
        (
            max(geometricAreaDensity, foamAreaDensity),
            maxAreaDensity_
        );

        // The main solver multiplies K by the Henry-model species
        // diffusivity D: D*K=(D/delta)*a_eff=k_L*a_eff [s^-1].
        transferK[celli] =
            effectiveAreaDensity/filmThickness_.value();

        if (reportArea_)
        {
            integralGeometricArea += geometricAreaDensity*cellVolumes[celli];
            integralFoamArea += foamAreaDensity*cellVolumes[celli];
            integralEffectiveArea += effectiveAreaDensity*cellVolumes[celli];
        }
    }

    if (reportArea_ && mesh.time().writeTime())
    {
        reduce(integralGeometricArea, sumOp<scalar>());
        reduce(integralFoamArea, sumOp<scalar>());
        reduce(integralEffectiveArea, sumOp<scalar>());

        // The model may be evaluated more than once per PIMPLE timestep.
        // Emit at most one report per output time on each rank.
        static scalar lastReportedTime = -GREAT;
        if (mag(mesh.time().value() - lastReportedTime) > SMALL)
        {
            lastReportedTime = mesh.time().value();
            if (Pstream::master())
            {
                Info<< "HEADSPACE_O2_AREA_AUDIT"
                    << " t_s=" << mesh.time().value()
                    << " pair=" << pair_.name()
                    << " geometricArea_m2=" << integralGeometricArea
                    << " foamCandidateArea_m2=" << integralFoamArea
                    << " effectiveArea_m2=" << integralEffectiveArea
                    << " filmThickness_m=" << filmThickness_.value()
                    << endl;
            }
        }
    }

    return tK;
}
