/*---------------------------------------------------------------------------*\
  Coupled PB size injection for headspace -> bubbleGas mixture transfer.

  Use positive NET bubble mass gain from actual phaseSystem dmdtf,
  in the same orientation the assembled bubble-phase audit uses.
  This avoids assigning size groups to bubble gas deposited out.
\*---------------------------------------------------------------------------*/

#include "headspaceEntrainmentNucleation.H"

#include "phasePairKey.H"
#include "phaseSystem.H"
#include "velocityGroup.H"
#include "sizeGroup.H"
#include "mathematicalConstants.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{
    defineTypeNameAndDebug(headspaceEntrainmentNucleation, 0);
    addToRunTimeSelectionTable
    (
        nucleationModel,
        headspaceEntrainmentNucleation,
        dictionary
    );
}
}
}

using Foam::constant::mathematical::pi;


Foam::diameterModels::nucleationModels::headspaceEntrainmentNucleation::
headspaceEntrainmentNucleation
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    nucleationModel(popBal, dict),
    dNuc_("nucleationDiameter", dimLength, dict),
    velGroup_
    (
        refCast<const velocityGroup>
        (
            popBal.mesh().lookupObject<phaseModel>
            (
                IOobject::groupName
                (
                    "alpha",
                    dict.lookup("velocityGroup")
                )
            ).dPtr()()
        )
    ),
    headspaceName_(dict.lookup("headspacePhase"))
{
    if
    (
        dNuc_.value() < velGroup_.sizeGroups().first().dSph().value()
     || dNuc_.value() > velGroup_.sizeGroups().last().dSph().value()
    )
    {
        FatalIOErrorInFunction(dict)
            << "headspaceEntrainmentNucleation diameter "
            << dNuc_.value() << " m outside configured size groups"
            << exit(FatalIOError);
    }

    if (velGroup_.popBalName() != popBal.name())
    {
        FatalIOErrorInFunction(dict)
            << "velocityGroup belongs to a different population balance"
            << exit(FatalIOError);
    }

    if (headspaceName_ == velGroup_.phase().name())
    {
        FatalIOErrorInFunction(dict)
            << "headspacePhase must differ from velocityGroup"
            << exit(FatalIOError);
    }
}


Foam::diameterModels::nucleationModels::headspaceEntrainmentNucleation::
~headspaceEntrainmentNucleation()
{}


void
Foam::diameterModels::nucleationModels::headspaceEntrainmentNucleation::
addToNucleationRate
(
    volScalarField& nucleationRate,
    const label i
)
{
    const sizeGroup& fi = popBal_.sizeGroups()[i];

    // The model's velocityGroup must match the PB size group being
    // populated. This particular SC1 PB has one velocityGroup.
    if (&fi.phase() != &velGroup_.phase())
    {
        return;
    }

    // This is the same bubble-first signed pair used by the
    // BUBBLE_MASS_AUDIT. Positive values create bubbleGas; negative
    // values remove bubbleGas and do not nucleate new bubbles.
    // By accessing fluid.dmdtf(), we include BiRD's actual phase
    // blending and netting against deposition, not raw candidate
    // entrainment from the unblended submodel.
    const tmp<volScalarField> tBubbleGain
    (
        popBal_.fluid().dmdtf
        (
            phasePairKey(velGroup_.phase().name(), headspaceName_)
        )
    );
    const volScalarField& bubbleGain = tBubbleGain();

    const tmp<volScalarField> tRho = fi.phase().rho();
    const volScalarField& rho = tRho();

    const dimensionedScalar zeroRate
    (
        "zeroRate", dimDensity/dimTime, 0
    );

    const dimensionedScalar groupAllocation
    (
        popBal_.eta(i, pi/6*pow3(dNuc_))
    );

    const tmp<volScalarField> tPositiveMass
    (
        max(bubbleGain, zeroRate)
    );
    const volScalarField& positiveMass = tPositiveMass();

    nucleationRate +=
        groupAllocation*positiveMass/rho/fi.x();

    // Report net positive bubble creation once per write. All size
    // groups share this same amount, allocated by eta(i,x_nuc).
    if (i == 0 && popBal_.mesh().time().writeTime())
    {
        static scalar lastReport = -GREAT;
        if (mag(popBal_.mesh().time().value() - lastReport) > SMALL)
        {
            lastReport = popBal_.mesh().time().value();

            scalar positiveBubbleMass = 0;
            const scalarField& V = popBal_.mesh().V();
            forAll(V, celli)
            {
                positiveBubbleMass += positiveMass[celli]*V[celli];
            }
            reduce(positiveBubbleMass, sumOp<scalar>());

            if (Pstream::master())
            {
                Info<< "HEADSPACE_ENTRAINMENT_PB_AUDIT"
                    << " t_s=" << popBal_.mesh().time().value()
                    << " positiveNetBubbleGain_kg_s="
                    << positiveBubbleMass
                    << " nucleationDiameter_m=" << dNuc_.value()
                    << endl;
            }
        }
    }
}

// ************************************************************************* //
