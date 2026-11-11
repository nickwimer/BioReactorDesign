#ifndef sc1FaceDragCoupling_H
#define sc1FaceDragCoupling_H
#include "sc1FaceDragAlgebra.H"

namespace Foam
{
// Restricted candidate: two moving incompressible phases, symmetric face drag,
// no virtual mass, fixed mesh, and the exact saved quiescent dictionaries.
// The Python runner enforces those dictionaries before launching this binary.
class sc1FaceDragCoupling
{
    surfaceScalarField t00_,t01_,t10_,t11_;

    static label coefficients
    (
        const UList<scalar>& a, const UList<scalar>& b,
        const UList<scalar>& k0, const UList<scalar>& k1,
        UList<scalar>& t00, UList<scalar>& t01,
        UList<scalar>& t10, UList<scalar>& t11
    )
    {
        label bad=0;
        forAll(a,i)
        {
            if (k0[i]!=k1[i] || !sc1FaceDragTransform
                (a[i],b[i],k0[i],t00[i],t01[i],t10[i],t11[i])) ++bad;
        }
        return bad;
    }

public:
    sc1FaceDragCoupling
    (
        const PtrList<surfaceScalarField>& a,
        const PtrList<surfaceScalarField>& k,
        const PtrList<surfaceScalarField>& rAUf
    )
    : t00_("sc1DragT00",rAUf[0]/rAUf[0]),
      t01_("sc1DragT01",rAUf[0]/rAUf[0]),
      t10_("sc1DragT10",rAUf[0]/rAUf[0]),
      t11_("sc1DragT11",rAUf[0]/rAUf[0])
    {
        label bad=coefficients
        (
            a[0].primitiveField(),a[1].primitiveField(),
            k[0].primitiveField(),k[1].primitiveField(),
            t00_.primitiveFieldRef(),t01_.primitiveFieldRef(),
            t10_.primitiveFieldRef(),t11_.primitiveFieldRef()
        );
        forAll(t00_.boundaryField(),p)
            bad+=coefficients
            (
                a[0].boundaryField()[p],a[1].boundaryField()[p],
                k[0].boundaryField()[p],k[1].boundaryField()[p],
                t00_.boundaryFieldRef()[p],t01_.boundaryFieldRef()[p],
                t10_.boundaryFieldRef()[p],t11_.boundaryFieldRef()[p]
            );
        reduce(bad,sumOp<label>());
        if (bad)
            FatalErrorInFunction << "Invalid two-phase face-drag matrix on "
                << bad << " faces including processor copies" << exit(FatalError);
        Info<< "SC1_FACE_COUPLING version=1 order=beforePressure" << endl;
    }

    void apply(surfaceScalarField& f0, surfaceScalarField& f1) const
    {
        const surfaceScalarField old0("sc1Uncoupled0",f0);
        const surfaceScalarField old1("sc1Uncoupled1",f1);
        f0=t00_*old0+t01_*old1;
        f1=t10_*old0+t11_*old1;
    }
};

// Check accepted timestep velocities, using the existing startup acceptance
// threshold and phase-presence mask. Stop before another timestep on failure.
inline void sc1QuiescentStep(const phaseSystem& fluid)
{
    const phaseModel& liquid=fluid.phases()[1];
    scalar peak=0;
    label bad=0;
    forAll(liquid,i)
    {
        const scalar speed=mag(liquid.U()()[i]);
        if (!std::isfinite(speed)) ++bad;
        else if (liquid[i]>.01) peak=max(peak,speed);
    }
    reduce(peak,maxOp<scalar>());
    reduce(bad,sumOp<label>());
    Info<< "SC1_QUIESCENT_STEP t_s=" << fluid.mesh().time().value()
        << " maxLiquidPresent_m_s=" << peak << endl;
    if (bad || peak>1e-3)
        FatalErrorInFunction << "Quiescent liquid speed exceeded 1 mm/s"
            << exit(FatalError);
}

inline void sc1CoupledFluxCheck
(
    const surfaceScalarField& mixture, const surfaceScalarField& difference
)
{
    scalar error=0, scale=0;
    label bad=0;
    auto inspect=[&](const UList<scalar>& f,const UList<scalar>& d)
    {
        forAll(f,i)
        {
            if (!std::isfinite(f[i]) || !std::isfinite(d[i])) ++bad;
            else { error=max(error,mag(d[i])); scale=max(scale,mag(f[i])); }
        }
    };
    inspect(mixture.primitiveField(),difference.primitiveField());
    forAll(mixture.boundaryField(),p)
        inspect(mixture.boundaryField()[p],difference.boundaryField()[p]);
    reduce(error,maxOp<scalar>()); reduce(scale,maxOp<scalar>());
    reduce(bad,sumOp<label>());
    Info<< "SC1_COUPLED_FLUX maxAbsError_m3_s=" << error
        << " maxAbsMixture_m3_s=" << scale << endl;
    if (bad || error>1e-12+1e-10*scale)
        FatalErrorInFunction << "Coupled phase fluxes disagree with pressure mixture flux"
            << exit(FatalError);
}
}
#endif
