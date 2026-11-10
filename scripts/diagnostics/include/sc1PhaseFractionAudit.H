// Observer for a private phaseSystem build. No phase names/count are hard-coded.
// Fixed mesh, global Euler, one alpha correction; observations only.
#ifndef sc1PhaseFractionAudit_H
#define sc1PhaseFractionAudit_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include "fvcSurfaceIntegrate.H"
#include "sc1AlphaSubcycleLedger.H"
#include <vector>

namespace Foam { namespace sc1PhaseDiagnostic {
namespace budget {
using namespace sc1AlphaBudget;
const label rescale = sc1AlphaBudget::nTerms;
const label nTerms = rescale + 1;
}
using Pair = std::array<scalar,2>;
using Terms = std::array<Pair,budget::nTerms>;
struct Totals { label timeIndex=-1, call=0; };
class Audit
{
    const phaseSystem& fluid_;
    const fvMesh& mesh_;
    const label n_, expectedMules_;
    const word reference_;
    scalar time_, dt_;
    label call_;
    std::vector<Terms> terms_;
    std::vector<scalarField> entry_, rho_, sourceAlpha_, base_, dgdt_, requested_, before_, old_;
    std::vector<label> event_, mules_;
    std::vector<bool> baseSeen_;
    std::vector<scalar> boundary_;
    void check(bool ok, const char* reason) const
    {
        if (!ok) FatalErrorInFunction << reason << exit(FatalError);
    }
    label index(const phaseModel& phase) const { return phase.index(); }
    void add(label i, label term, label c, scalar amount)
    {
        const scalar v=amount*mesh_.V()[c];
        terms_[i][term][0]+=v; terms_[i][term][1]+=v*rho_[i][c];
    }
public:
    Audit(const phaseSystem& fluid, Totals& totals, bool LTS, label nAlphaCorr)
    : fluid_(fluid), mesh_(fluid.mesh()), n_(fluid.phases().size()),
      expectedMules_(mesh_.solverDict("alpha").lookup<label>("nAlphaSubCycles")),
      reference_(fluid.lookupOrDefault<word>("referencePhase",word::null)),
      time_(mesh_.time().value()), dt_(mesh_.time().deltaTValue()),
      terms_(n_), entry_(n_), rho_(n_), sourceAlpha_(n_), base_(n_), dgdt_(n_),
      requested_(n_), before_(n_), old_(n_), event_(n_,-1), mules_(n_,0),
      baseSeen_(n_,false), boundary_(n_,0)
    {
        check(!LTS && nAlphaCorr==1 && expectedMules_>0 && n_>=2
            && !fluid.stationaryPhases().size() && !mesh_.moving() && !mesh_.changing(),
            "Phase audit needs fixed mesh, global Euler, one alpha correction and moving phases only");
        const IOdictionary schemes(IOobject("fvSchemes",mesh_.time().system(),mesh_,
            IOobject::MUST_READ,IOobject::NO_WRITE,false));
        const dictionary& ddt=schemes.subDict("ddtSchemes");
        check(ddt.size()==1 && ddt.lookupOrDefault<word>("default","none")=="Euler",
            "Phase audit needs default-only Euler");
        if (totals.timeIndex!=mesh_.time().timeIndex())
        { totals.timeIndex=mesh_.time().timeIndex(); totals.call=0; }
        call_=++totals.call;
        for (label i=0;i<n_;++i)
        {
            const phaseModel& phase=fluid.phases()[i];
            check(phase.index()==i,"Noncontiguous phase indices");
            entry_[i]=phase.primitiveField(); rho_[i]=phase.rho()().primitiveField();
            sourceAlpha_[i]=entry_[i];
            // OpenFOAM Field::operator= returns void: do not chain assignments.
            base_[i]=scalarField(phase.size(),0);
            dgdt_[i]=scalarField(phase.size(),0);
            requested_[i]=scalarField(phase.size(),0);
        }
    }
    void originalSource(const phaseModel& phase, const volScalarField::Internal& su)
    {
        const label i=index(phase); base_[i]=su.field(); sourceAlpha_[i]=phase.primitiveField();
        const tmp<volScalarField> requested=phase.divU();
        if (requested.valid()) requested_[i]=requested().primitiveField();
        baseSeen_[i]=true;
    }
    void rawDilatation(const phaseModel& phase, const volScalarField::Internal& d)
    { dgdt_[index(phase)]=d.field(); }
    void beforeChange(const phaseModel& phase, label event)
    {
        const label i=index(phase); check(event_[i]==-1,"Overlapping phase audit hook");
        before_[i]=phase.primitiveField(); event_[i]=event;
    }
    void afterChange(const phaseModel& phase, label event)
    {
        const label i=index(phase); check(event_[i]==event,"Mismatched phase audit hook");
        forAll(mesh_.V(),c) add(i,event,c,phase[c]-before_[i][c]);
        event_[i]=-1;
    }
    void beforeMules(const phaseModel& phase)
    {
        beforeChange(phase,sc1AlphaBudget::mulesChange);
        const label i=index(phase); check(baseSeen_[i],"Missing assembled alpha source");
        const volScalarField& alpha=phase;
        old_[i]=alpha.oldTime().primitiveField();
    }
    void afterMules(const phaseModel& phase, const surfaceScalarField& flux,
        const volScalarField::Internal& sp, const volScalarField::Internal& su)
    {
        const label i=index(phase);
        check(event_[i]==sc1AlphaBudget::mulesChange,"Mismatched MULES hook");
        scalarField div(phase.size(),0); fvc::surfaceIntegrate(div,flux);
        const scalar subDt=mesh_.time().deltaTValue();
        forAll(mesh_.V(),c)
        {
            const auto cell=sc1AlphaBudget::mulesCell(before_[i][c],old_[i][c],phase[c],subDt,
                div[c],su[c],sp[c],sourceAlpha_[i][c],base_[i][c],dgdt_[i][c],requested_[i][c]);
            for (label j=0;j<sc1AlphaBudget::nTerms;++j) add(i,j,c,cell[j]);
        }
        forAll(flux.boundaryField(),p)
            if (!flux.boundaryField()[p].coupled())
                forAll(flux.boundaryField()[p],f) boundary_[i]+=subDt*flux.boundaryField()[p][f];
        ++mules_[i]; event_[i]=-1;
    }
    void afterImplicit(const phaseModel& phase, const fvScalarMatrix& equation)
    {
        using namespace sc1AlphaBudget;
        const label i=index(phase); afterChange(phase,implicitCorrection);
        const tmp<surfaceScalarField> flux=equation.flux();
        scalarField div(phase.size(),0); fvc::surfaceIntegrate(div,flux());
        forAll(mesh_.V(),c)
        {
            const scalar d=mesh_.time().deltaTValue()*div[c];
            const scalar r=phase[c]-before_[i][c]+d;
            add(i,implicitDivergence,c,d); add(i,implicitEquationResidual,c,r);
            add(i,implicitEquationResidualL1,c,mag(r));
        }
    }
    void finish()
    {
        using namespace sc1AlphaBudget;
        static const char* labels[budget::nTerms]={"mulesChange","oldTimeReset","explicitSu",
            "implicitSp","divergence","equationResidual","equationResidualL1","originalSu",
            "rawDgdt","pressureDilatation","couplingDifference","splitAtSourceAlpha",
            "implicitSourceChange","sourceDecompositionRemainder","implicitCorrection",
            "referenceSubcycle","clamp","referenceFinal","implicitDivergence",
            "implicitEquationResidual","implicitEquationResidualL1","finalTransportDivergence",
            "fluxAveragingDefect","rescale"};
        const unsigned precision=Info().precision(15);
        for (label i=0;i<n_;++i)
        {
            const phaseModel& phase=fluid_.phases()[i];
            const bool solved=phase.name()!=reference_;
            check(event_[i]==-1 && mules_[i]==(solved ? expectedMules_ : 0),
                "Incomplete phase audit hooks or unsupported subcycling");
            if (solved)
            {
                const tmp<surfaceScalarField> flux=phase.alphaPhi();
                scalarField div(phase.size(),0); fvc::surfaceIntegrate(div,flux());
                forAll(mesh_.V(),c) add(i,finalTransportDivergence,c,dt_*div[c]);
                for (label k=0;k<2;++k) terms_[i][fluxAveragingDefect][k]=
                    terms_[i][finalTransportDivergence][k]-terms_[i][divergence][k]
                    -terms_[i][implicitDivergence][k];
            }
            scalarField packed(2*budget::nTerms+4,0);
            for (label j=0;j<budget::nTerms;++j)
                for (label k=0;k<2;++k) packed[2*j+k]=terms_[i][j][k];
            forAll(mesh_.V(),c)
            {
                const scalar d=(phase[c]-entry_[i][c])*mesh_.V()[c];
                packed[2*budget::nTerms]+=d;
                packed[2*budget::nTerms+1]+=d*rho_[i][c];
                packed[2*budget::nTerms+3]+=mag(d);
            }
            packed[2*budget::nTerms+2]=boundary_[i];
            Pstream::listCombineGather(packed,plusEqOp<scalar>());
            Pstream::listCombineScatter(packed);
            Pair unassigned{{packed[2*budget::nTerms],packed[2*budget::nTerms+1]}};
            const label changes[]={mulesChange,implicitCorrection,referenceSubcycle,
                clamp,referenceFinal,budget::rescale};
            for (const label j:changes)
                for (label k=0;k<2;++k) unassigned[k]-=packed[2*j+k];
            Info<< "SC1_PHASE_SOLVE phase=" << phase.name() << " t_s=" << time_
                << " dt_s=" << dt_ << " call=" << call_ << " solved=" << solved
                << " mulesCalls=" << mules_[i]
                << " deltaVolume_m3=" << packed[2*budget::nTerms]
                << " deltaFrozenRhoMass_kg=" << packed[2*budget::nTerms+1]
                << " changeL1Volume_m3=" << packed[2*budget::nTerms+3]
                << " unassignedVolume_m3=" << unassigned[0]
                << " unassignedFrozenRhoMass_kg=" << unassigned[1]
                << " mulesBoundaryVolume_m3=" << packed[2*budget::nTerms+2]
                << " divBoundaryDifference_m3=" << packed[2*divergence]-packed[2*budget::nTerms+2]
                << endl;
            for (label j=0;j<budget::nTerms;++j)
            {
                // Reconstructed phases do not have an explicit-source budget.
                if (!solved && j!=referenceSubcycle && j!=referenceFinal) continue;
                Info<< "SC1_PHASE_TERM phase=" << phase.name() << " t_s=" << time_
                    << " call=" << call_ << " term=" << labels[j]
                    << " volume_m3=" << packed[2*j] << " frozenRhoMass_kg=" << packed[2*j+1] << endl;
            }
        }
        Info().precision(precision);
    }
};
}}
#endif
