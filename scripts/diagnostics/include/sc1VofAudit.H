#ifndef sc1VofAudit_H
#define sc1VofAudit_H
#include <cmath>

namespace Foam
{
// Read-only observer for this fixed-mesh, Euler, no-source VOF reference.
class sc1VofAudit
{
    const fvMesh& mesh_;
    scalarField initialAlpha_;
    scalar previousL_, previousG_, initialL_, initialG_;
    scalar liquidOut_=0, gasOut_=0, l1L_=0, l1G_=0;
    label calls_=0;

    scalar volume(const volScalarField& a) const
    {
        return gSum(a.primitiveField()*mesh_.V().field());
    }

    scalar boundary(const surfaceScalarField& f) const
    {
        scalar result=0;
        label bad=0;
        forAll(f.boundaryField(), patchi)
        {
            if (f.boundaryField()[patchi].coupled()
             || mesh_.boundary()[patchi].type()=="empty") continue;
            forAll(f.boundaryField()[patchi], facei)
            {
                const scalar x=f.boundaryField()[patchi][facei];
                if (!std::isfinite(x)) ++bad;
                else result+=x;
            }
        }
        reduce(bad,sumOp<label>()); reduce(result,sumOp<scalar>());
        if (bad) FatalErrorInFunction << "Nonfinite VOF boundary flux" << exit(FatalError);
        return result;
    }

public:
    sc1VofAudit(const volScalarField& liquid, const volScalarField& gas)
    : mesh_(liquid.mesh()), initialAlpha_(liquid.primitiveField()),
      previousL_(volume(liquid)), previousG_(volume(gas)),
      initialL_(previousL_), initialG_(previousG_)
    {
        Info.precision(17);
        Info<< "SC1_VOF_INITIAL liquid_m3=" << initialL_
            << " gas_m3=" << initialG_ << " version=1" << endl;
    }

    void alphaFlux(const surfaceScalarField& alphaPhi, const surfaceScalarField& phi)
    {
        // Last outer iteration only: alphaPhi is the actual accepted Euler
        // transport flux. phi here is the flux used by that alpha solve,
        // before the subsequent pressure correction. No flux is modified.
        liquidOut_=boundary(alphaPhi);
        gasOut_=boundary(phi)-liquidOut_;
        ++calls_;
    }

    void pressure(const volScalarField& prgh, const volScalarField& rho,
                  const volScalarField& gh) const
    {
        scalar low=GREAT, high=-GREAT; label bad=0;
        forAll(prgh,i)
        {
            const scalar p=prgh[i]+rho[i]*gh[i];
            if (!std::isfinite(p)) ++bad;
            else { low=min(low,p); high=max(high,p); }
        }
        reduce(low,minOp<scalar>()); reduce(high,maxOp<scalar>());
        reduce(bad,sumOp<label>());
        Info<< "SC1_VOF_PRESSURE t_s=" << mesh_.time().value()
            << " min_Pa=" << low << " max_Pa=" << high << " nonfinite=" << bad << endl;
        if (bad || low<1000)
            FatalErrorInFunction << "Invalid VOF unclipped pressure" << exit(FatalError);
    }

    void step(const volScalarField& liquid, const volScalarField& gas,
              const volVectorField& U, const surfaceScalarField& phi)
    {
        scalar low=GREAT,high=-GREAT,complement=0,speedL=0,speedG=0,change=0;
        label bad=0;
        forAll(liquid,i)
        {
            const scalar a=liquid[i], b=gas[i], speed=mag(U[i]);
            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(speed)) ++bad;
            else
            {
                low=min(low,min(a,b)); high=max(high,max(a,b));
                complement=max(complement,mag(a+b-1));
                if (a>.01) speedL=max(speedL,speed);
                if (b>.01) speedG=max(speedG,speed);
                change+=mag(a-initialAlpha_[i])*mesh_.V()[i];
            }
        }
        forAll(phi,i) if (!std::isfinite(phi[i])) ++bad;
        forAll(phi.boundaryField(),patchi)
            forAll(phi.boundaryField()[patchi],facei)
                if (!std::isfinite(phi.boundaryField()[patchi][facei])) ++bad;
        reduce(bad,sumOp<label>()); reduce(low,minOp<scalar>()); reduce(high,maxOp<scalar>());
        reduce(complement,maxOp<scalar>()); reduce(speedL,maxOp<scalar>());
        reduce(speedG,maxOp<scalar>()); reduce(change,sumOp<scalar>());
        if (bad) FatalErrorInFunction << "Nonfinite VOF accepted state" << exit(FatalError);
        const scalar dt=mesh_.time().deltaTValue();
        const scalar vl=volume(liquid), vg=volume(gas);
        const scalar rl=1000*(vl-previousL_+dt*liquidOut_);
        const scalar rg=vg-previousG_+dt*gasOut_;
        l1L_+=mag(rl); l1G_+=mag(rg);
        const scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
        const scalar co=.5*dt*gMax(sumPhi/mesh_.V().field());
        const scalar postOut=boundary(phi);
        Info<< "SC1_VOF_STEP t_s=" << mesh_.time().value()
            << " dt_s=" << dt << " liquid_m3=" << vl << " gas_m3=" << vg
            << " liquidOut_m3_s=" << liquidOut_ << " gasOut_m3_s=" << gasOut_
            << " postPressureOut_m3_s=" << postOut
            << " liquidResidual_kg=" << rl << " gasResidual_kg=" << rg
            << " liquidL1_kg=" << l1L_ << " gasL1_kg=" << l1G_
            << " alphaMin=" << low << " alphaMax=" << high << " complement=" << complement
            << " liquidMax_m_s=" << speedL << " gasMax_m_s=" << speedG
            << " alphaChangeL1_m3=" << change << " Co=" << co
            << " alphaCalls=" << calls_ << endl;
        previousL_=vl; previousG_=vg;
        if (calls_!=2 || low< -1e-10 || high>1+1e-10 || complement>1e-10
         || !std::isfinite(co) || co>.2 || speedL>1e-3 || change>1e-5
         || l1L_>1e-10+1e-8*1000*initialL_ || l1G_>1e-10+1e-8*initialG_)
            FatalErrorInFunction << "VOF reference failed unchanged startup screens" << exit(FatalError);
        calls_=0;
    }
};
}
#endif
