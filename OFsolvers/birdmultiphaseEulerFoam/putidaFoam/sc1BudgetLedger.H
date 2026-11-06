// Small, OpenFOAM-independent accumulator used by the scratch SC1 audit.
// Rates are right-endpoint samples, not a claim of time-discrete closure.
#ifndef sc1BudgetLedger_H
#define sc1BudgetLedger_H
#include <array>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace sc1Budget
{
enum Term
{
    massExternal, massTransfer, massIn, massOut,
    oxygenExternal, oxygenTransfer, oxygenReaction,
    oxygenAdvection, oxygenDiffusion, oxygenInProxy, oxygenOutProxy,
    equationResidual, equationResidualL1, nTerms
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
}
#endif
