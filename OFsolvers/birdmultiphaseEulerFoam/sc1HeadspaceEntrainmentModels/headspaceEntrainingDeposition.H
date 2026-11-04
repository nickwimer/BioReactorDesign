/*---------------------------------------------------------------------------*\
  Experimental SC1 free-surface gas re-entrainment, optional scratch model.
  Retains original deposition exactly when entrainmentEfficiency = 0.
\*---------------------------------------------------------------------------*/

#ifndef headspaceEntrainingDeposition_H
#define headspaceEntrainingDeposition_H

#include "phaseTransferModel.H"

namespace Foam
{
namespace phaseTransferModels
{

class headspaceEntrainingDeposition
:
    public phaseTransferModel
{
    const word dropletName_;
    const word surfaceName_;
    const word liquidName_;
    const scalar depositionEfficiency_;
    const scalar entrainmentEfficiency_;

    const scalar yMin_;
    const scalar yMax_;
    const scalar rMin_;
    const scalar rMax_;
    const scalar axisX_;
    const scalar axisZ_;

    const scalar minAlpha_;
    const scalar fullAlpha_;
    const scalar maxBubbleAlpha_;
    const scalar maxSlipSpeed_;
    const scalar donorTimeScale_;

    // Maximum sum of UNBLENDED positive entrainment rates, kg/s.
    const scalar maxGlobalEntrainedKgS_;
    const bool requireDownwardLiquid_;

    // legacyDownwardLiquid preserves the original experiment exactly.
    // normalGasIntoLiquid uses positive phase-normal relative velocity.
    const word entrainmentVelocityCriterion_;

    const bool reportRate_;

    scalar gate(const scalar alpha) const;

public:

    TypeName("headspaceEntrainingDeposition");

    headspaceEntrainingDeposition
    (
        const dictionary& dict,
        const phasePair& pair
    );

    virtual ~headspaceEntrainingDeposition();

    virtual bool mixture() const;

    virtual tmp<volScalarField> dmdtf() const;
};

}
}

#endif
