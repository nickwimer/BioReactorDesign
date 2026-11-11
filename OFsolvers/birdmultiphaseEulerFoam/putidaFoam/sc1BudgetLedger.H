// Small, OpenFOAM-independent accumulator used by the scratch SC1 audit.
// Rates are right-endpoint samples, not a claim of time-discrete closure.
#ifndef sc1BudgetLedger_H
#define sc1BudgetLedger_H
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace sc1Budget
{
enum Term
{
    massExternal, massTransfer, massIn, massOut,
    oxygenExternal, oxygenTransfer, oxygenReaction,
    oxygenAdvection, oxygenDiffusion, oxygenInProxy, oxygenOutProxy,
    equationResidual, equationResidualL1,
    oxygenDdt, oxygenStabilization,
    oxygenComponentMismatch, oxygenComponentMismatchL1, nTerms
};
typedef std::array<double, nTerms> Rates;
struct Ledger
{
    double start, previous, windowStart, initialMass, initialOxygen;
    double windowMass, windowOxygen;
    double oxygenJump = 0, windowJump = 0;
    double normalization = 0, windowNormalization = 0;
    unsigned long steps = 0, windowSteps = 0;
    Rates integral{}, window{};
    Ledger(double t, double m, double o)
    : start(t), previous(t), windowStart(t), initialMass(m),
      initialOxygen(o), windowMass(m), windowOxygen(o) {}

    void advance(double t, double dt, const Rates& rates,
                 double jump, double norm)
    {
        if (!std::isfinite(t) || !std::isfinite(dt) || dt <= 0
            || !std::isfinite(jump) || !std::isfinite(norm)
            || std::abs((t-previous)-dt) > 1e-9*std::max(1.0, std::abs(t)))
            throw std::runtime_error("SC1 audit: skipped/duplicate timestep");
        for (unsigned i=0; i<nTerms; ++i)
        {
            if (!std::isfinite(rates[i]))
                throw std::runtime_error("SC1 audit: non-finite rate");
            integral[i] += dt*rates[i];
            window[i] += dt*rates[i];
        }
        oxygenJump += jump; windowJump += jump;
        normalization += norm; windowNormalization += norm;
        previous=t; ++steps; ++windowSteps;
    }
    static double massExpected(const Rates& r)
    { return r[massExternal]+r[massTransfer]+r[massIn]-r[massOut]; }
    static double oxygenExpected(const Rates& r)
    {
        return r[oxygenExternal]+r[oxygenTransfer]+r[oxygenReaction]
             - r[oxygenAdvection]-r[oxygenDiffusion];
    }
    void resetWindow(double t, double m, double o)
    {
        windowStart=t; windowMass=m; windowOxygen=o;
        window={}; windowSteps=0; windowJump=0; windowNormalization=0;
    }
};

// Positive correction predicted by Y <- max(Y,0), using the actual alpha/rho
// weights. A negative phase weight remains visible instead of being hidden.
struct NegativeOxygenStats
{
    double correctionMass = 0;
    double minimumY = std::numeric_limits<double>::infinity();
    unsigned long cells = 0;

    void add(double alpha, double rho, double y, double volume)
    {
        minimumY = std::min(minimumY,y);
        if (y < 0)
        {
            correctionMass -= alpha*rho*y*volume;
            ++cells;
        }
    }
};

// These are stage impulses in kg, not rates: do not multiply them by dt.
// Keep the v1 Ledger::oxygenExpected definition unchanged so old and new
// diagnostics can be compared without silently changing residual meanings.
struct OxygenStageLedger
{
    double oldTimeOffset = 0, oldTimeOffsetSum = 0;
    double ddtInventoryDefect = 0, ddtInventoryDefectSum = 0;
    double postSpeciesChange = 0, postSpeciesChangeSum = 0;
    double negativeMassSum = 0;
    double correctionDefect = 0, correctionDefectSum = 0;
    unsigned long negativePhaseCellSamples = 0;

    void advance
    (
        double startOxygen, double oldOxygen, double beforeCorrectionOxygen,
        double afterCorrectionOxygen, double endOxygen,
        double dt, double ddtRate, double negativeMass,
        unsigned long negativeCells, double correction
    )
    {
        const double values[] = {startOxygen,oldOxygen,beforeCorrectionOxygen,
            afterCorrectionOxygen,endOxygen,dt,ddtRate,negativeMass,correction};
        for (double value : values)
            if (!std::isfinite(value))
                throw std::runtime_error("SC1 audit: non-finite oxygen stage");
        if (dt <= 0)
            throw std::runtime_error("SC1 audit: non-positive oxygen-stage dt");
        oldTimeOffset = oldOxygen-startOxygen;
        ddtInventoryDefect = beforeCorrectionOxygen-oldOxygen-dt*ddtRate;
        postSpeciesChange = endOxygen-afterCorrectionOxygen;
        correctionDefect = correction-negativeMass;
        oldTimeOffsetSum += oldTimeOffset;
        ddtInventoryDefectSum += ddtInventoryDefect;
        postSpeciesChangeSum += postSpeciesChange;
        negativeMassSum += negativeMass;
        correctionDefectSum += correctionDefect;
        negativePhaseCellSamples += negativeCells;
    }

    double explainedResidual(const Rates& integral) const
    {
        return integral[oxygenStabilization]+integral[equationResidual]
             + integral[oxygenComponentMismatch]+oldTimeOffsetSum
             + ddtInventoryDefectSum+postSpeciesChangeSum;
    }
};
}
#endif
