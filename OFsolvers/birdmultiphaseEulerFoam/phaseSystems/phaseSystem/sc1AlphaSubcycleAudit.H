// Included ONLY by the instrumented scratch phaseSystemSolve.C.
// All observations use already assembled alpha sources and actual limited flux.
// No phase source/model/correction is reevaluated. No phaseSystem ABI changes.
#ifndef sc1AlphaSubcycleAudit_H
#define sc1AlphaSubcycleAudit_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include "fvcSurfaceIntegrate.H"
#include "sc1AlphaSubcycleLedger.H"
#include <map>

namespace Foam
{
namespace sc1AlphaDiagnostic
{
using Pair = std::array<scalar,2>; // volume [m3], mass weighted by entry rho [kg]
using PhaseTerms = std::array<Pair,sc1AlphaBudget::nTerms>;
using Table = std::array<PhaseTerms,5>;
using Inventories = std::array<Pair,5>;

struct Totals
{
    // Primitive values only: no registered objects or mesh lifetime dependency.
    Table all{};
    Inventories measured{}, unassigned{};
    std::array<scalar,3> boundaryVolume{};
    std::map<label,Table> byCallInStep;
    label lastTimeIndex = -1, callInStep = 0;
    unsigned long solveCalls = 0;
};

class Audit
{
    const phaseSystem& fluid_;
    const fvMesh& mesh_;
    Totals& totals_;
    wordList names_;
    Table terms_{};
    std::array<scalarField,3> entryAlpha_, entryRho_, sourceAlpha_;
    std::array<scalarField,3> originalSu_, dgdt_, requestedDilatation_;
    std::array<scalarField,3> before_, old_;
    std::array<bool,3> baseSeen_{}, sourceSeen_{}, open_{};
    std::array<label,3> event_{};
    std::array<unsigned long,3> mulesCalls_{};
    std::array<scalar,3> boundaryVolume_{};
    scalar time_, dt_;
    label callInStep_;

    label index(const phaseModel& phase) const
    {
        for (label i=0;i<3;++i) if (phase.name()==names_[i]) return i;
        FatalErrorInFunction << "Unexpected phase in SC1 alpha audit: " << phase.name()
            << exit(FatalError);
        return -1;
    }

    void fixedMesh() const
    {
        if (mesh_.moving() || mesh_.changing())
            FatalErrorInFunction << "SC1 alpha audit requires fixed mesh"
                << exit(FatalError);
    }

    void add(label phase, label term, label cell, scalar amount)
    {
        const scalar volume = amount*mesh_.V()[cell];
        terms_[phase][term][0] += volume;
        terms_[phase][term][1] += volume*entryRho_[phase][cell];
    }

    template<class T> static void aggregate(std::array<T,5>& data)
    {
        for (unsigned j=0;j<data[0].size();++j)
        {
            data[3][j]=data[1][j]+data[2][j];
            data[4][j]=data[0][j]+data[3][j];
        }
    }

    static void aggregateTable(Table& data)
    {
        for (label j=0;j<sc1AlphaBudget::nTerms;++j)
            for (label k=0;k<2;++k)
            {
                data[3][j][k]=data[1][j][k]+data[2][j][k];
                data[4][j][k]=data[0][j][k]+data[3][j][k];
            }
    }

public:
    Audit(const phaseSystem& fluid, Totals& totals, bool LTS, label nAlphaCorr)
    : fluid_(fluid), mesh_(fluid.mesh()), totals_(totals), names_(5),
      time_(mesh_.time().value()), dt_(mesh_.time().deltaTValue())
    {
        fixedMesh();
        if (LTS || nAlphaCorr!=1 || fluid.phases().size()!=3
            || fluid.stationaryPhases().size())
            FatalErrorInFunction
                << "SC1 alpha audit requires global dt, one alpha correction,"
                << " three moving phases and no stationary phases"
                << exit(FatalError);
        if (!totals_.solveCalls)
        {
            const IOdictionary schemes
            (
                IOobject("fvSchemes",mesh_.time().system(),mesh_,
                    IOobject::MUST_READ,IOobject::NO_WRITE,false)
            );
            const dictionary& ddt=schemes.subDict("ddtSchemes");
            if (ddt.size()!=1 || ddt.lookupOrDefault<word>("default","none")!="Euler")
                FatalErrorInFunction << "SC1 alpha audit requires default-only Euler"
                    << exit(FatalError);
        }
        if (fluid.lookupOrDefault<word>("referencePhase","none")!="headspaceGas")
            FatalErrorInFunction << "SC1 alpha audit requires headspaceGas reference phase"
                << exit(FatalError);
        names_[0]="liquid"; names_[1]="bubbleGas"; names_[2]="headspaceGas";
        names_[3]="gas"; names_[4]="total";
        if (totals_.lastTimeIndex!=mesh_.time().timeIndex())
        {
            totals_.lastTimeIndex=mesh_.time().timeIndex();
            totals_.callInStep=0;
        }
        callInStep_=++totals_.callInStep;
        ++totals_.solveCalls;
        for (label i=0;i<3;++i)
        {
            const phaseModel& phase=fluid_.phases()[names_[i]];
            const tmp<volScalarField> rho=phase.rho();
            entryAlpha_[i]=phase.primitiveField();
            entryRho_[i]=rho().primitiveField();
            sourceAlpha_[i]=entryAlpha_[i];
            originalSu_[i]=scalarField(phase.size(),0);
            dgdt_[i]=scalarField(phase.size(),0);
            requestedDilatation_[i]=scalarField(phase.size(),0);
        }
    }

    void originalSource(const phaseModel& phase, const volScalarField::Internal& su)
    {
        const label i=index(phase);
        originalSu_[i]=su.field();
        sourceAlpha_[i]=phase.primitiveField();
        const tmp<volScalarField> requested=phase.divU();
        if (requested.valid()) requestedDilatation_[i]=requested().primitiveField();
        baseSeen_[i]=true;
    }

    void rawDilatation(const phaseModel& phase, const volScalarField::Internal& dgdt)
    { dgdt_[index(phase)]=dgdt.field(); }

    void beforeMules(const phaseModel& phase)
    {
        const label i=index(phase);
        if (!baseSeen_[i] || open_[i])
            FatalErrorInFunction << "SC1 alpha audit: missing source or nested MULES hook"
                << exit(FatalError);
        before_[i]=phase.primitiveField();
        const volScalarField& alpha=phase;
        // MULES::limit has already accessed this existing level; const access
        // does not create an additional old-old-time field.
        old_[i]=alpha.oldTime().primitiveField();
        open_[i]=true;
        event_[i]=sc1AlphaBudget::mulesChange;
    }

    void afterMules(const phaseModel& phase, const surfaceScalarField& alphaPhi,
                    const volScalarField::Internal& sp,
                    const volScalarField::Internal& su)
    {
        using namespace sc1AlphaBudget;
        const label i=index(phase);
        if (!open_[i] || event_[i]!=mulesChange)
            FatalErrorInFunction << "SC1 alpha audit: unmatched MULES hook"
                << exit(FatalError);
        scalarField div(phase.size(),0);
        fvc::surfaceIntegrate(div,alphaPhi);
        const scalar subDt=mesh_.time().deltaTValue();
        forAll(mesh_.V(),c)
        {
            const auto cell=mulesCell(before_[i][c],old_[i][c],phase[c],subDt,
                div[c],su[c],sp[c],sourceAlpha_[i][c],originalSu_[i][c],
                dgdt_[i][c],requestedDilatation_[i][c]);
            for (label j=0;j<nTerms;++j) add(i,j,c,cell[j]);
        }
        forAll(alphaPhi.boundaryField(),patchi)
        {
            const auto& patch=alphaPhi.boundaryField()[patchi];
            if (!patch.coupled())
                forAll(patch,facei) boundaryVolume_[i]+=subDt*patch[facei];
        }
        ++mulesCalls_[i];
        sourceSeen_[i]=true;
        open_[i]=false;
    }

    void beforeChange(const phaseModel& phase, label event)
    {
        const label i=index(phase);
        if (open_[i])
            FatalErrorInFunction << "SC1 alpha audit: overlapping change hook"
                << exit(FatalError);
        before_[i]=phase.primitiveField();
        event_[i]=event;
        open_[i]=true;
    }

    void afterChange(const phaseModel& phase, label event)
    {
        const label i=index(phase);
        if (!open_[i] || event_[i]!=event)
            FatalErrorInFunction << "SC1 alpha audit: unmatched change hook"
                << exit(FatalError);
        forAll(mesh_.V(),c) add(i,event,c,phase[c]-before_[i][c]);
        open_[i]=false;
    }

    void afterImplicit(const phaseModel& phase, const fvScalarMatrix& equation)
    {
        using namespace sc1AlphaBudget;
        const label i=index(phase);
        afterChange(phase,implicitCorrection);
        const tmp<surfaceScalarField> flux=equation.flux();
        scalarField div(phase.size(),0);
        fvc::surfaceIntegrate(div,flux());
        const scalar subDt=mesh_.time().deltaTValue();
        forAll(mesh_.V(),c)
        {
            const scalar divChange=subDt*div[c];
            const scalar residual=phase[c]-before_[i][c]+divChange;
            add(i,implicitDivergence,c,divChange);
            add(i,implicitEquationResidual,c,residual);
            add(i,implicitEquationResidualL1,c,mag(residual));
        }
    }

    void finish()
    {
        using namespace sc1AlphaBudget;
        fixedMesh();
        Inventories measured{}, unassigned{};
        // One packed reduction per solve, not dozens of collectives per subcycle.
        const label width=nTerms*2+3;
        scalarField packed(3*width,0);
        for (label i=0;i<3;++i)
        {
            if (open_[i] || (i<2 && !sourceSeen_[i]))
                FatalErrorInFunction << "Incomplete SC1 alpha hooks" << exit(FatalError);
            const phaseModel& phase=fluid_.phases()[names_[i]];
            if (i<2)
            {
                // Compare the actual stored averaged flux to the actual fluxes
                // used in every explicit and implicit subcycle update.
                const tmp<surfaceScalarField> flux=phase.alphaPhi();
                scalarField div(phase.size(),0);
                fvc::surfaceIntegrate(div,flux());
                forAll(mesh_.V(),c) add(i,finalTransportDivergence,c,dt_*div[c]);
                for (label k=0;k<2;++k)
                    terms_[i][fluxAveragingDefect][k]=terms_[i][finalTransportDivergence][k]
                        -terms_[i][divergence][k]-terms_[i][implicitDivergence][k];
            }
            for (label j=0;j<nTerms;++j)
                for (label k=0;k<2;++k) packed[i*width+j*2+k]=terms_[i][j][k];
            forAll(mesh_.V(),c)
            {
                const scalar change=(phase[c]-entryAlpha_[i][c])*mesh_.V()[c];
                packed[i*width+nTerms*2]+=change;
                packed[i*width+nTerms*2+1]+=change*entryRho_[i][c];
            }
            packed[i*width+nTerms*2+2]=boundaryVolume_[i];
        }
        Pstream::listCombineGather(packed,plusEqOp<scalar>());
        Pstream::listCombineScatter(packed);
        for (label i=0;i<3;++i)
        {
            for (label j=0;j<nTerms;++j)
                for (label k=0;k<2;++k) terms_[i][j][k]=packed[i*width+j*2+k];
            measured[i][0]=packed[i*width+nTerms*2];
            measured[i][1]=packed[i*width+nTerms*2+1];
            boundaryVolume_[i]=packed[i*width+nTerms*2+2];
            totals_.boundaryVolume[i]+=boundaryVolume_[i];
            for (label k=0;k<2;++k)
                unassigned[i][k]=measured[i][k]-terms_[i][mulesChange][k]
                    -terms_[i][implicitCorrection][k]-terms_[i][referenceSubcycle][k]
                    -terms_[i][clamp][k]-terms_[i][referenceFinal][k];
        }
        aggregateTable(terms_); aggregate(measured); aggregate(unassigned);
        Table& outer=totals_.byCallInStep[callInStep_];
        for (label i=0;i<5;++i)
        {
            for (label k=0;k<2;++k)
            {
                totals_.measured[i][k]+=measured[i][k];
                totals_.unassigned[i][k]+=unassigned[i][k];
            }
            for (label j=0;j<nTerms;++j)
                for (label k=0;k<2;++k)
                {
                    totals_.all[i][j][k]+=terms_[i][j][k];
                    outer[i][j][k]+=terms_[i][j][k];
                }
        }
        if (!mesh_.time().writeTime() && mesh_.time().running()) return;
        static const char* labels[nTerms]={"mulesChange","oldTimeReset","explicitSu",
            "implicitSp","divergence","equationResidual","equationResidualL1",
            "originalSu","rawDgdt","pressureDilatation","couplingDifference",
            "splitAtSourceAlpha","implicitSourceChange","sourceDecompositionRemainder",
            "implicitCorrection","referenceSubcycle","clamp","referenceFinal",
            "implicitDivergence","implicitEquationResidual","implicitEquationResidualL1",
            "finalTransportDivergence","fluxAveragingDefect"};
        const unsigned oldPrecision=Info().precision(15);
        Info<< "SC1_ALPHA_METHOD version=1 fixedMeshGlobalEuler=1"
            << " rhoWeight=frozenAtAlphaSolveEntry pressureDilatation=storedInput"
            << " ratesIntegratedWithActualSubDt=1 physicalSourceReevaluation=0"
            << " sourceTermRows=explicitlySolvedPhasesOnly"
            << " allCallsIncludeOuterIterationResets=1" << endl;
        for (label i=0;i<5;++i)
        {
            Info<< "SC1_ALPHA_SOLVE phase=" << names_[i] << " t_s=" << time_
                << " dt_s=" << dt_ << " callInStep=" << callInStep_
                << " solveCalls=" << totals_.solveCalls
                << " measuredVolume_m3=" << measured[i][0]
                << " measuredRhoBaseMass_kg=" << measured[i][1]
                << " unassignedVolume_m3=" << unassigned[i][0]
                << " unassignedRhoBaseMass_kg=" << unassigned[i][1]
                << " cumulativeMeasuredVolume_m3=" << totals_.measured[i][0]
                << " cumulativeMeasuredRhoBaseMass_kg=" << totals_.measured[i][1]
                << " cumulativeUnassignedVolume_m3=" << totals_.unassigned[i][0]
                << " cumulativeUnassignedRhoBaseMass_kg=" << totals_.unassigned[i][1];
            if (i<3)
                Info<< " mulesCalls=" << mulesCalls_[i]
                    << " mulesBoundaryVolume_m3=" << boundaryVolume_[i]
                    << " mulesDivBoundaryDifference_m3="
                    << terms_[i][divergence][0]-boundaryVolume_[i]
                    << " cumulativeMulesBoundaryVolume_m3=" << totals_.boundaryVolume[i]
                    << " cumulativeMulesDivBoundaryDifference_m3="
                    << totals_.all[i][divergence][0]-totals_.boundaryVolume[i];
            Info<< endl;
            for (label j=0;j<nTerms;++j)
            {
                // Headspace is reconstructed; aggregate source rows would be
                // partial sums over solved phases, not a gas/total source budget.
                if (i>=2 && j!=mulesChange && j!=implicitCorrection
                    && j!=referenceSubcycle && j!=clamp && j!=referenceFinal) continue;
                Info<< "SC1_ALPHA_TERM phase=" << names_[i] << " t_s=" << time_
                    << " callInStep=" << callInStep_ << " term=" << labels[j]
                    << " volume_m3=" << terms_[i][j][0]
                    << " rhoBaseMass_kg=" << terms_[i][j][1]
                    << " cumulativeVolume_m3=" << totals_.all[i][j][0]
                    << " cumulativeRhoBaseMass_kg=" << totals_.all[i][j][1]
                    << " sameCallIndexCumulativeVolume_m3=" << outer[i][j][0]
                    << " sameCallIndexCumulativeRhoBaseMass_kg=" << outer[i][j][1]
                    << endl;
            }
        }
        Info().precision(oldPrecision);
    }
};
}
}
#endif
