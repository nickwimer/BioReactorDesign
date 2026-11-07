// Exact measured inventory changes across paired, non-nested solver blocks.
// Observational only: no sources, corrections, solves or normalization invoked.
// alpha/PB blocks do NOT isolate MULES sources, subcycles or alpha clamping.
#ifndef sc1MassStageAudit_H
#define sc1MassStageAudit_H
#include "fvCFD.H"
#include "phaseSystem.H"
#include "sc1MassStageLedger.H"

namespace Foam
{
class sc1MassStageAudit
{
    phaseSystem& fluid_;
    const fvMesh& mesh_;
    wordList names_;
    autoPtr<sc1MassStage::Ledger> ledger_;
    std::array<scalarField, 3> alphaBefore_, rhoBefore_;

    void fixedMesh() const
    {
        if (mesh_.moving() || mesh_.changing())
            FatalErrorInFunction << "SC1 mass stage audit requires a fixed mesh"
                << exit(FatalError);
    }

    static void aggregate(sc1MassStage::Mass& mass)
    {
        mass[3] = mass[1] + mass[2];
        mass[4] = mass[0] + mass[3];
    }

    sc1MassStage::Mass inventory() const
    {
        fixedMesh();
        sc1MassStage::Mass result{};
        for (label i = 0; i < 3; ++i)
        {
            const phaseModel& phase = fluid_.phases()[names_[i]];
            const tmp<volScalarField> rho = phase.rho();
            result[i] = fvc::domainIntegrate(phase*rho()).value();
        }
        aggregate(result);
        return result;
    }

public:
    sc1MassStageAudit(phaseSystem& fluid, bool LTS)
    : fluid_(fluid), mesh_(fluid.mesh()), names_(5)
    {
        fixedMesh();
        if (LTS || fluid_.phases().size() != 3)
            FatalErrorInFunction
                << "SC1 mass stage audit requires global timestep and three phases"
                << exit(FatalError);

        // Read only configuration, without altering the live discretization.
        // Fail closed on overrides rather than guessing which ddt scheme is used.
        const IOdictionary schemes
        (
            IOobject
            (
                "fvSchemes", mesh_.time().system(), mesh_,
                IOobject::MUST_READ, IOobject::NO_WRITE, false
            )
        );
        const dictionary& ddt = schemes.subDict("ddtSchemes");
        if (ddt.size() != 1 || ddt.lookupOrDefault<word>("default", "none") != "Euler")
            FatalErrorInFunction
                << "SC1 mass stage audit requires only default Euler in ddtSchemes"
                << exit(FatalError);

        names_[0] = "liquid"; names_[1] = "bubbleGas";
        names_[2] = "headspaceGas"; names_[3] = "gas"; names_[4] = "total";
        const auto mass = inventory();
        ledger_.reset(new sc1MassStage::Ledger(mesh_.time().value(), mass));
        const unsigned oldPrecision = Info().precision(15);
        Info<< "SC1_MASS_STAGE_METHOD version=1 inventoryChanges=exactMeasured"
            << " allocation=deltaAlphaTimesRhoBefore"
            << " density=alphaAfterTimesDeltaRho"
            << " allOuterIterations=1 sourcesReevaluated=0"
            << " alphaSourcesAndClamps=unresolved physicalConservationCertificate=0"
            << endl;
        for (label i = 0; i < 5; ++i)
            Info<< "SC1_MASS_STAGE_START phase=" << names_[i]
                << " t_s=" << mesh_.time().value() << " mass_kg=" << mass[i] << endl;
        Info().precision(oldPrecision);
    }

    // After scheduled updates, before runTime++; changes since the prior sample
    // remain visible as unassigned changes instead of silently resetting a baseline.
    void beforeAdvance()
    { ledger_->beforeAdvance(mesh_.time().value(), inventory()); }

    // Blocks must not overlap. Repeated calls with one name (e.g. both outer
    // iterations) accumulate together. Time rates are never multiplied by dt.
    void beforeBlock(const word& stage)
    {
        const auto mass = inventory();
        ledger_->beforeBlock(stage.c_str(), mass);
        for (label i = 0; i < 3; ++i)
        {
            const phaseModel& phase = fluid_.phases()[names_[i]];
            const tmp<volScalarField> rho = phase.rho();
            alphaBefore_[i] = phase.primitiveField();
            rhoBefore_[i] = rho().primitiveField();
        }
    }

    void afterBlock(const word& stage)
    {
        fixedMesh();
        sc1MassStage::Mass allocation{}, density{};
        for (label i = 0; i < 3; ++i)
        {
            const phaseModel& phase = fluid_.phases()[names_[i]];
            const tmp<volScalarField> rho = phase.rho();
            if (alphaBefore_[i].size() != phase.size()
                || rhoBefore_[i].size() != phase.size())
                FatalErrorInFunction << "SC1 mass stage: missing or invalid block snapshot"
                    << exit(FatalError);
            forAll(mesh_.V(), c)
            {
                const auto change = sc1MassStage::cellChange
                (
                    alphaBefore_[i][c], rhoBefore_[i][c],
                    phase[c], rho()[c], mesh_.V()[c]
                );
                allocation[i] += change[0];
                density[i] += change[1];
            }
            reduce(allocation[i], sumOp<scalar>());
            reduce(density[i], sumOp<scalar>());
        }
        aggregate(allocation); aggregate(density);
        ledger_->afterBlock(stage.c_str(), inventory(), allocation, density);
    }

    // Call once EVERY timestep; report controls output only, not accumulation.
    void sample(bool report)
    {
        const auto mass = inventory();
        const scalar time = mesh_.time().value();
        ledger_->finishStep(time, mesh_.time().deltaTValue(), mass);
        if (!report) return;
        const unsigned oldPrecision = Info().precision(15);
        for (const auto& item : ledger_->stages)
        {
            const auto& step = item.second.step;
            const auto& sum = item.second.cumulative;
            for (label i = 0; i < 5; ++i)
                Info<< "SC1_MASS_STAGE phase=" << names_[i]
                    << " stage=" << item.first.c_str() << " t_s=" << time
                    << " stepCalls=" << step.calls << " cumulativeCalls=" << sum.calls
                    << " stepMassDelta_kg=" << step.mass[i]
                    << " stepAllocation_kg=" << step.allocation[i]
                    << " stepDensity_kg=" << step.density[i]
                    << " stepDecompositionRemainder_kg=" << step.remainder[i]
                    << " cumulativeMassDelta_kg=" << sum.mass[i]
                    << " cumulativeAllocation_kg=" << sum.allocation[i]
                    << " cumulativeDensity_kg=" << sum.density[i]
                    << " cumulativeDecompositionRemainder_kg=" << sum.remainder[i]
                    << endl;
        }
        for (label i = 0; i < 5; ++i)
            Info<< "SC1_MASS_STAGE_TOTAL phase=" << names_[i]
                << " t_s=" << time << " steps=" << ledger_->steps
                << " mass_kg=" << mass[i]
                << " beforeAdvanceChange_kg=" << ledger_->beforeAdvanceChange[i]
                << " stepActualDelta_kg=" << ledger_->stepActual[i]
                << " stepAssignedDelta_kg=" << ledger_->stepAssigned[i]
                << " stepUnassignedDelta_kg=" << ledger_->stepUnassigned[i]
                << " cumulativeActualDelta_kg=" << ledger_->cumulativeActual[i]
                << " cumulativeAssignedDelta_kg=" << ledger_->cumulativeAssigned[i]
                << " cumulativeUnassignedDelta_kg=" << ledger_->cumulativeUnassigned[i]
                << endl;
        Info().precision(oldPrecision);
    }
};
}
#endif
