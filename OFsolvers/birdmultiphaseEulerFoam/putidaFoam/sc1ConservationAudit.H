// Opt-in instrumentation for an isolated SC1 executable. No model corrections
// are called here. See docs/sc1_timestep_conservation_audit.md for limitations.
#ifndef sc1ConservationAudit_H
#define sc1ConservationAudit_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include "thermophysicalTransportModel.H"
#include "velocityGroup.H"
#include "sizeGroup.H"
#include "sc1BudgetLedger.H"
#include <vector>

namespace Foam
{
class sc1ConservationAudit
{
    phaseSystem& fluid_;
    const fvMesh& mesh_;
    wordList names_;
    std::vector<sc1Budget::Ledger> ledgers_;
    std::array<sc1Budget::Rates, 5> rates_{};
    std::array<scalar, 5> previousO_{}, jump_{}, normalization_{};
    std::array<bool, 3> speciesSeen_{};
    PtrList<fvScalarMatrix> equations_, transfer_, external_, advection_;
    PtrList<fvScalarMatrix> diffusion_, reaction_;
    std::array<scalar, 3> beforeNormalization_{};
    std::array<scalar, 3> pbBefore_{}, pbCorrection_{}, pbCorrectionSum_{};
    std::array<scalar, 5> pairIntegral_{};
    scalar pbInput_ = 0;
    bool pbSeen_ = false;
    // Per phase/patch: positive gross mass in/out and O2 in/out proxy integrals.
    std::vector<std::array<scalar, 4>> patchIntegral_;

    scalar integral(const tmp<volScalarField>& f) const
    { return fvc::domainIntegrate(f()).value(); }

    label index(const word& phase) const
    {
        for (label i=0; i<3; ++i) if (names_[i] == phase) return i;
        FatalErrorInFunction << "Unexpected SC1 phase " << phase
            << exit(FatalError);
        return -1;
    }

    void inventory(std::array<scalar, 5>& mass,
                   std::array<scalar, 5>& oxygen,
                   std::array<scalar, 5>& volume) const
    {
        mass.fill(0); oxygen.fill(0); volume.fill(0);
        for (label i=0; i<3; ++i)
        {
            const phaseModel& p = fluid_.phases()[names_[i]];
            const tmp<volScalarField> rho = p.rho();
            mass[i] = integral(p*rho());
            oxygen[i] = integral(p*rho()*p.Y("O2"));
            const volScalarField& alpha = p;
            volume[i] = fvc::domainIntegrate(alpha).value();
        }
        for (label i=0; i<3; ++i)
        {
            mass[4]+=mass[i]; oxygen[4]+=oxygen[i]; volume[4]+=volume[i];
            if (i>0)
            { mass[3]+=mass[i]; oxygen[3]+=oxygen[i]; volume[3]+=volume[i]; }
        }
    }

    std::array<scalar, 3> pbMoments() const
    {
        const phaseModel& p = fluid_.phases()["bubbleGas"];
        if (!isA<diameterModels::velocityGroup>(p.dPtr()()))
            FatalErrorInFunction << "SC1 audit requires bubble velocityGroup"
                << exit(FatalError);
        const auto& vg = refCast<const diameterModels::velocityGroup>(p.dPtr()());
        const auto& groups = vg.sizeGroups();
        const tmp<volScalarField> rho = p.rho();
        std::array<scalar, 3> result{};
        forAll(groups, i)
        {
            const auto& fi = groups[i];
            forAll(mesh_.V(), c)
            {
                const scalar v = p[c]*fi[c]*mesh_.V()[c];
                result[0]+=v; result[1]+=v*rho()[c];
                result[2]+=v/fi.x().value();
            }
        }
        for (auto& x : result) reduce(x, sumOp<scalar>());
        return result;
    }

public:
    sc1ConservationAudit(phaseSystem& fluid, bool LTS)
    : fluid_(fluid), mesh_(fluid.mesh()), names_(5),
      equations_(3), transfer_(3), external_(3), advection_(3),
      diffusion_(3), reaction_(3), patchIntegral_(3*mesh_.boundary().size())
    {
        if (LTS || mesh_.moving() || fluid.phases().size()!=3)
            FatalErrorInFunction << "SC1 audit requires global dt, fixed mesh, three phases"
                << exit(FatalError);
        names_[0]="liquid"; names_[1]="bubbleGas"; names_[2]="headspaceGas";
        names_[3]="gas"; names_[4]="total";
        std::array<scalar,5> m{}, o{}, v{};
        inventory(m,o,v);
        previousO_=o;
        const unsigned oldPrecision = Info().precision(15);
        for (label i=0; i<5; ++i)
        {
            ledgers_.emplace_back(mesh_.time().value(),m[i],o[i]);
            Info<< "SC1_BALANCE_START phase=" << names_[i]
                << " t_s=" << mesh_.time().value()
                << " mass_kg=" << m[i] << " oxygen_kg=" << o[i]
                << " volume_m3=" << v[i] << endl;
        }
        Info<< "SC1_BALANCE_METHOD version=1 quadrature=rightEndpointEveryStep"
            << " exactDiscreteBudget=0 speciesStage=finalPimpleBeforeNormalization"
            << " boundaryOxygen=faceYProxy PBSource=preSolveProxy" << endl;
        Info().precision(oldPrecision);
    }

    // After scheduled biology/pH updates, immediately before runTime++.
    // Record impulses; never change the intentional biology schedule.
    void beforeAdvance()
    {
        std::array<scalar,5> m{},o{},v{};
        inventory(m,o,v);
        for (label i=0; i<5; ++i)
        {
            jump_[i]=o[i]-previousO_[i];
            normalization_[i]=0; rates_[i].fill(0);
        }
        speciesSeen_.fill(false); pbSeen_=false; pbInput_=0;
    }

    // Final PIMPLE only, immediately before/after fluid.solve(), then
    // immediately after fluid.correct(). No repeated PB solve or scaling.
    void beforePB()
    {
        const tmp<volScalarField> s = fluid_.dmdtf
            (phasePairKey("bubbleGas", "headspaceGas"));
        pbInput_=0;
        forAll(mesh_.V(), c) pbInput_+=max(s()[c],scalar(0))*mesh_.V()[c];
        reduce(pbInput_,sumOp<scalar>());
    }
    void afterPBSolve() { pbBefore_=pbMoments(); pbSeen_=true; }
    void afterPBCorrect()
    {
        const auto after=pbMoments();
        for (label i=0; i<3; ++i) pbCorrection_[i]=after[i]-pbBefore_[i];
    }

    // Copy the ACTUAL unrelaxed equation and RHS matrices before solve.
    // Auxiliary advection/diffusion/reaction matrices are reconstructed at
    // the same state, using the same API as MultiComponentPhaseModel::YiEqn.
    void equation(phaseModel& p, volScalarField& y, const fvScalarMatrix& eq,
                  const fvScalarMatrix& transfer, const fvScalarMatrix& source)
    {
        const label i=index(p.name());
        equations_.set(i,new fvScalarMatrix(eq));
        transfer_.set(i,new fvScalarMatrix(transfer));
        external_.set(i,new fvScalarMatrix(source));
        const tmp<surfaceScalarField> flux=p.alphaRhoPhi();
        advection_.set(i,fvm::div(flux(),y,"div("+flux().name()+",Yi)").ptr());
        // OF9 phaseModel has no divj interface. MovingPhaseModel delegates
        // to this same registered, phase-grouped transport model. A missing
        // object/type must fail rather than silently assuming zero diffusion.
        const thermophysicalTransportModel& transport =
            mesh_.lookupObject<thermophysicalTransportModel>
            (
                IOobject::groupName
                (
                    thermophysicalTransportModel::typeName, flux().group()
                )
            );
        diffusion_.set(i,transport.divj(y).ptr());
        const volScalarField& alpha = p;
        tmp<fvScalarMatrix> reaction(alpha*p.R(y));
        reaction_.set(i,reaction.ptr());
        speciesSeen_[i]=true;
    }

    void beforeSpeciesCorrection()
    {
        for (label i=0; i<3; ++i)
        {
            if (!speciesSeen_[i])
                FatalErrorInFunction << "Missing final O2 equation for " << names_[i]
                    << exit(FatalError);
            const phaseModel& p=fluid_.phases()[names_[i]];
            const volScalarField& y=p.Y("O2");
            auto& r=rates_[i];
            // operator& gives (A*y-source)/V: positive RHS action is gain.
            r[sc1Budget::oxygenTransfer]=integral(transfer_[i]&y);
            r[sc1Budget::oxygenExternal]=integral(external_[i]&y);
            r[sc1Budget::oxygenAdvection]=integral(advection_[i]&y);
            r[sc1Budget::oxygenDiffusion]=integral(diffusion_[i]&y);
            r[sc1Budget::oxygenReaction]=integral(reaction_[i]&y);
            const tmp<volScalarField> residual(equations_[i]&y);
            r[sc1Budget::equationResidual]=integral(residual);
            r[sc1Budget::equationResidualL1]=integral(mag(residual()));
            beforeNormalization_[i]=integral(p*p.rho()*y);
        }
    }
    void afterSpeciesCorrection()
    {
        for (label i=0; i<3; ++i)
        {
            const phaseModel& p=fluid_.phases()[names_[i]];
            normalization_[i]=integral(p*p.rho()*p.Y("O2"))-beforeNormalization_[i];
        }
    }

    void sample(bool report)
    {
        if (mesh_.moving() || mesh_.changing())
            FatalErrorInFunction << "Moving/changing mesh unsupported by SC1 audit"
                << exit(FatalError);
        const scalar t=mesh_.time().value(), dt=mesh_.time().deltaTValue();
        const unsigned oldPrecision=Info().precision(15);
        std::array<scalar,5> m{},o{},v{};
        inventory(m,o,v);
        const PtrList<volScalarField> sources(fluid_.dmdts());
        for (label i=0; i<3; ++i)
        {
            if (!speciesSeen_[i] || !pbSeen_)
                FatalErrorInFunction << "Incomplete final-PIMPLE SC1 audit hooks"
                    << exit(FatalError);
            const phaseModel& p=fluid_.phases()[names_[i]];
            const volScalarField& alpha=p;
            const tmp<volScalarField> rho=p.rho();
            auto& r=rates_[i];
            if (sources.set(p.index()))
                r[sc1Budget::massTransfer]=fvc::domainIntegrate(sources[p.index()]).value();
            if (fluid_.fvModels().addsSupToField(rho().name()))
                r[sc1Budget::massExternal]=integral(fluid_.fvModels().source(alpha,rho())&rho());
            const tmp<surfaceScalarField> flux=p.alphaRhoPhi();
            const auto& y=p.Y("O2");
            forAll(mesh_.boundary(), patchi)
            {
                const auto& f=flux().boundaryField()[patchi];
                if (f.coupled()) continue;
                std::array<scalar,4> q{};
                forAll(f, facei)
                {
                    const scalar yf=y.boundaryField()[patchi][facei];
                    q[0]-=min(f[facei],scalar(0)); q[1]+=max(f[facei],scalar(0));
                    q[2]-=min(f[facei],scalar(0))*yf;
                    q[3]+=max(f[facei],scalar(0))*yf;
                }
                auto& cumulative=patchIntegral_[i*mesh_.boundary().size()+patchi];
                for (label k=0;k<4;++k)
                { reduce(q[k],sumOp<scalar>()); cumulative[k]+=dt*q[k]; }
                r[sc1Budget::massIn]+=q[0]; r[sc1Budget::massOut]+=q[1];
                r[sc1Budget::oxygenInProxy]+=q[2]; r[sc1Budget::oxygenOutProxy]+=q[3];
                if (report)
                    Info<< "SC1_BOUNDARY phase=" << names_[i]
                        << " patch=" << mesh_.boundary()[patchi].name() << " t_s=" << t
                        << " in_kg_s=" << q[0] << " out_kg_s=" << q[1]
                        << " oxygenInProxy_kg_s=" << q[2] << " oxygenOutProxy_kg_s=" << q[3]
                        << " cumIn_kg=" << cumulative[0] << " cumOut_kg=" << cumulative[1]
                        << " cumOxygenInProxy_kg=" << cumulative[2]
                        << " cumOxygenOutProxy_kg=" << cumulative[3] << endl;
            }
            for (label k=0;k<sc1Budget::nTerms;++k)
            {
                rates_[4][k]+=r[k]; if(i>0) rates_[3][k]+=r[k];
            }
            normalization_[4]+=normalization_[i];
            if(i>0) normalization_[3]+=normalization_[i];
        }

        const tmp<volScalarField> bh=fluid_.dmdtf(phasePairKey("bubbleGas","headspaceGas"));
        const auto& b=fluid_.phases()["bubbleGas"]; const auto& h=fluid_.phases()["headspaceGas"];
        scalar donor=0;
        forAll(mesh_.V(),c)
            donor+=(max(bh()[c],scalar(0))*h.Y("O2")[c]
                  +min(bh()[c],scalar(0))*b.Y("O2")[c])*mesh_.V()[c];
        reduce(donor,sumOp<scalar>());
        const std::array<scalar,5> pairRates = {{integral(bh),
            -integral(fluid_.dmdtf(phasePairKey("bubbleGas","liquid"))),
            -integral(fluid_.dmdtf(phasePairKey("headspaceGas","liquid"))),pbInput_,donor}};
        for(label k=0;k<5;++k) pairIntegral_[k]+=dt*pairRates[k];
        for(label k=0;k<3;++k) pbCorrectionSum_[k]+=pbCorrection_[k];

        static const char* labels[sc1Budget::nTerms]={
            "massExternal","massTransfer","massIn","massOut","oxygenExternal",
            "oxygenTransfer","oxygenReaction","oxygenAdvection","oxygenDiffusion",
            "oxygenInProxy","oxygenOutProxy","equationResidual","equationResidualL1"};
        for(label i=0;i<5;++i)
        {
            auto& l=ledgers_[i];
            l.advance(t,dt,rates_[i],jump_[i],normalization_[i]);
            if(report)
            {
                Info<< "SC1_BALANCE phase=" << names_[i] << " t_s=" << t
                    << " start_s=" << l.start << " windowStart_s=" << l.windowStart
                    << " steps=" << l.steps << " windowSteps=" << l.windowSteps
                    << " mass_kg=" << m[i] << " oxygen_kg=" << o[i]
                    << " volume_m3=" << v[i]
                    << " massDelta_kg=" << m[i]-l.initialMass
                    << " oxygenDelta_kg=" << o[i]-l.initialOxygen
                    << " massResidual_kg=" << m[i]-l.initialMass-sc1Budget::Ledger::massExpected(l.integral)
                    << " oxygenResidual_kg=" << o[i]-l.initialOxygen-sc1Budget::Ledger::oxygenExpected(l.integral)-l.oxygenJump-l.normalization
                    << " windowMassResidual_kg=" << m[i]-l.windowMass-sc1Budget::Ledger::massExpected(l.window)
                    << " windowOxygenResidual_kg=" << o[i]-l.windowOxygen-sc1Budget::Ledger::oxygenExpected(l.window)-l.windowJump-l.windowNormalization
                    << " oxygenJump_kg=" << l.oxygenJump
                    << " normalization_kg=" << l.normalization;
                for(label k=0;k<sc1Budget::nTerms;++k)
                    Info<< " " << labels[k] << "_kg_s=" << rates_[i][k]
                        << " " << labels[k] << "Integral_kg=" << l.integral[k];
                if(i==0) Info<< " DO_mg_L=" << 1000*o[i]/max(v[i],SMALL);
                Info<< endl;
                l.resetWindow(t,m[i],o[i]);
            }
        }
        if(report)
        {
            Info<< "SC1_PAIR t_s=" << t
                << " bubbleFromHeadspaceIntegral_kg=" << pairIntegral_[0]
                << " bubbleToLiquidMassIntegral_kg=" << pairIntegral_[1]
                << " headspaceToLiquidMassIntegral_kg=" << pairIntegral_[2]
                << " positiveNetPBInputProxyIntegral_kg=" << pairIntegral_[3]
                << " donorOxygenToBubbleProxyIntegral_kg=" << pairIntegral_[4] << endl;
            Info<< "SC1_PB t_s=" << t << " positiveNetInputProxy_kg_s=" << pbInput_
                << " preCorrectVolume_m3=" << pbBefore_[0]
                << " preCorrectMass_kg=" << pbBefore_[1]
                << " preCorrectNumber=" << pbBefore_[2]
                << " correctionVolume_m3=" << pbCorrection_[0]
                << " correctionMass_kg=" << pbCorrection_[1]
                << " correctionNumber=" << pbCorrection_[2]
                << " correctionVolumeSum_m3=" << pbCorrectionSum_[0]
                << " correctionMassSum_kg=" << pbCorrectionSum_[1]
                << " correctionNumberSum=" << pbCorrectionSum_[2] << endl;
        }
        previousO_=o;
        Info().precision(oldPrecision);
    }
};
}
#endif
