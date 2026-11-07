/*---------------------------------------------------------------------------*\
  SC1 experimental headspace-to-liquid oxygen-transfer closure.
  Diagnostic branch: not production validated.
\*---------------------------------------------------------------------------*/

#ifndef resolvedHeadspaceFilm_H
#define resolvedHeadspaceFilm_H

#include "diffusiveMassTransferModel.H"

namespace Foam
{
namespace diffusiveMassTransferModels
{

class resolvedHeadspaceFilm
:
    public diffusiveMassTransferModel
{
    // This class returns K [m^-2]. The interface-composition system
    // multiplies it by D_O2 [m^2/s] and rho_liquid [kg/m^3].
    // K = a_eff / filmThickness, with a_eff [m^-1].

    word headspacePhase_;
    word liquidPhase_;
    word bubblePhase_;

    dimensionedScalar filmThickness_;
    dimensionedScalar foamDiameter_;

    scalar minInterfaceAlpha_;
    scalar fullInterfaceAlpha_;
    scalar maxBubbleAlpha_;

    // Optional sub-grid foam area: all zeros unless explicitly enabled.
    scalar foamAreaMultiplier_;
    scalar foamYMin_;
    scalar foamYMax_;
    scalar foamRMax_;
    scalar foamAxisX_;
    scalar foamAxisZ_;

    // Cap in inverse metres, applied to the area per volume.
    scalar maxAreaDensity_;
    bool reportArea_;

    scalar interfaceGate(const scalar alpha) const;

public:

    TypeName("resolvedHeadspaceFilm");

    resolvedHeadspaceFilm
    (
        const dictionary& dict,
        const phasePair& pair
    );

    virtual ~resolvedHeadspaceFilm();

    virtual tmp<volScalarField> K() const;
};

} // namespace diffusiveMassTransferModels
} // namespace Foam

#endif
