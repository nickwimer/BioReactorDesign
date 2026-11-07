/*---------------------------------------------------------------------------*\
  Scratch-only size-group injection for re-entrained headspace gas.
\*---------------------------------------------------------------------------*/

#ifndef headspaceEntrainmentNucleation_H
#define headspaceEntrainmentNucleation_H

#include "nucleationModel.H"

namespace Foam
{
namespace diameterModels
{
namespace nucleationModels
{

class headspaceEntrainmentNucleation
:
    public nucleationModel
{
    const dimensionedScalar dNuc_;
    const velocityGroup& velGroup_;
    const word headspaceName_;

public:

    TypeName("headspaceEntrainmentNucleation");

    headspaceEntrainmentNucleation
    (
        const populationBalanceModel& popBal,
        const dictionary& dict
    );

    virtual ~headspaceEntrainmentNucleation();

    virtual void addToNucleationRate
    (
        volScalarField& nucleationRate,
        const label i
    );
};

}
}
}

#endif
