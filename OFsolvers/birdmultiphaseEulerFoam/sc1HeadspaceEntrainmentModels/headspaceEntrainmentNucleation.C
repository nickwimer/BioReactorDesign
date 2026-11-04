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
        "zeroRate", phaseTransferModel::dimDmdt, 0
    );

    const dimensionedScalar groupAllocation
    (
        popBal_.eta(i, pi/6*pow3(dNuc_))
    );

    nucleationRate +=
        groupAllocation*max(bubbleGain, zeroRate)/rho/fi.x();
}

// ************************************************************************* //
