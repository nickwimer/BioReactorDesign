// Private experiment only: two moving, pure, constant-density gas/liquid phases.
#ifndef sc1TwoFieldRateSplit_H
#define sc1TwoFieldRateSplit_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include "sc1TwoFieldRateAlgebra.H"
#include <algorithm>
namespace Foam
{
class sc1TwoFieldRateSplit
{
    const phaseSystem& fluid_;
    scalarField rates_[2];
    bool valid_=false, applied_[2]={false,false};
    scalar sumError_=0, identityScaled_=0;
    label call_;
    static void require(bool ok, const char* reason)
    {
        if (!ok) FatalErrorInFunction << "SC1 rate split: " << reason << exit(FatalError);
    }
public:
    sc1TwoFieldRateSplit(const phaseSystem& fluid, label call, bool LTS, label corrections, bool audit)
    : fluid_(fluid), call_(call)
    {
        require(audit && !LTS && corrections==1 && fluid.phases().size()==2
            && !fluid.stationaryPhases().size() && !fluid.mesh().moving()
            && !fluid.mesh().changing(), "requires two moving phases, fixed mesh, Euler audit and one alpha correction");
        require(fluid.lookup<word>("type")=="basicMultiphaseSystem",
            "unsupported phase-system type");
        require(fluid.phases()[0].name()=="gas" && fluid.phases()[1].name()=="liquid",
            "requires phases (gas liquid)");
        const tmp<volScalarField> q0=fluid.phases()[0].divU();
        const tmp<volScalarField> q1=fluid.phases()[1].divU();
        require(q0.valid()==q1.valid(), "inconsistent pressure-rate availability");
        valid_=q0.valid();
        for (label i=0;i<2;++i)
        {
            const phaseModel& phase=fluid.phases()[i];
            require(phase.index()==i && phase.pure() && phase.isochoric() && phase.incompressible(),
                "requires pure incompressible constant-density phases");
            const tmp<volScalarField> density=phase.rho();
            forAll(phase,c)
                require(std::isfinite(density()[c]) && mag(density()[c]-(i==0 ? 1 : 1000))<1e-12,
                    "unexpected phase density");
            rates_[i]=scalarField(phase.size(),0);
        }
        if (valid_)
        {
            rates_[0]=q0().primitiveField();
            rates_[1]=q1().primitiveField();
        }
        forAll(fluid.mesh().V(),c)
        {
            const scalar a=fluid.phases()[0][c], b=fluid.phases()[1][c];
            sc1RateSource::Coefficients coefficients;
            require(sc1RateSource::coefficients(a,b,rates_[0][c],rates_[1][c],coefficients),
                "nonfinite/negative rate or noncomplementary/out-of-bounds fractions");
            sumError_=max(sumError_,mag(a+b-1));
        }
    }
    void apply(const phaseModel& phase, volScalarField::Internal& sp,
        volScalarField::Internal& su, const volScalarField::Internal& raw)
    {
        const label i=phase.index(), j=1-i;
        require(valid_ && !applied_[i], "missing/repeated source assembly");
        forAll(phase,c)
        {
            sc1RateSource::Coefficients coefficients;
            require(sc1RateSource::coefficients(phase[c],fluid_.phases()[j][c],
                rates_[i][c],rates_[j][c],coefficients), "unsupported source state");
            const scalar error=mag(coefficients.su+coefficients.sp*phase[c]-raw[c])
                /(1e-30+rates_[i][c]+rates_[j][c]);
            require(std::isfinite(error) && error<=1e-12, "raw coupled source identity failed");
            identityScaled_=max(identityScaled_,error);
            // Preserve the existing divergence source and every later equation.
            su[c]+=coefficients.su;
            sp[c]+=coefficients.sp;
        }
        applied_[i]=true;
    }
    void report() const
    {
        require(applied_[0]==valid_ && applied_[1]==valid_, "incomplete rate assembly");
        scalar sumError=sumError_, identity=identityScaled_;
        reduce(sumError,maxOp<scalar>()); reduce(identity,maxOp<scalar>());
        const unsigned precision=Info().precision(15);
        Info<< "SC1_RATE_STATE t_s=" << fluid_.mesh().time().value()
            << " call=" << call_ << " valid=" << label(valid_)
            << " applied=" << label(applied_[0])+label(applied_[1])
            << " qGasMin=" << gMin(rates_[0]) << " qGasMax=" << gMax(rates_[0])
            << " qLiquidMin=" << gMin(rates_[1]) << " qLiquidMax=" << gMax(rates_[1])
            << " alphaSumMaxError=" << sumError << " rawIdentityScaled=" << identity << endl;
        Info().precision(precision);
    }
    static void checkOptionSource(const volScalarField& source)
    {
        forAll(source,c) require(std::isfinite(source[c]) && source[c]>=0,
            "negative/nonfinite option source is unsupported");
    }
    static void checkTransfer(const volScalarField& source)
    {
        forAll(source,c) require(source[c]==0, "interphase transfer is unsupported");
    }
};
}
#endif
