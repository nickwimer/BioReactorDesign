// OpenFOAM-independent arithmetic for observational SC1 stage attribution.
// These are measured inventory changes, not integrated physical mass sources.
#ifndef sc1MassStageLedger_H
#define sc1MassStageLedger_H
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace sc1MassStage
{
using Mass = std::array<double, 5>;

// Include the product cross-term in the density contribution by using alphaAfter.
inline std::array<double, 2> cellChange(double alphaBefore, double rhoBefore,
                                      double alphaAfter, double rhoAfter,
                                      double volume)
{
    return {{(alphaAfter - alphaBefore)*rhoBefore*volume,
             alphaAfter*(rhoAfter - rhoBefore)*volume}};
}

inline void finite(const Mass& values)
{
    for (double value : values)
        if (!std::isfinite(value))
            throw std::runtime_error("SC1 mass stage: non-finite inventory change");
}

struct Change
{
    Mass mass{}, allocation{}, density{}, remainder{};
    unsigned long calls = 0;
};

struct Stage
{
    Change step, cumulative;
};

class Ledger
{
    Mass initial_, previous_, blockStart_{};
    double previousTime_;
    bool stepOpen_ = false, blockOpen_ = false;
    std::string block_;

public:
    std::map<std::string, Stage> stages;
    Mass stepActual{}, cumulativeActual{}, stepAssigned{}, cumulativeAssigned{};
    Mass stepUnassigned{}, cumulativeUnassigned{}, beforeAdvanceChange{};
    unsigned long steps = 0;

    Ledger(double time, const Mass& mass)
    : initial_(mass), previous_(mass), previousTime_(time)
    {
        finite(mass);
        if (!std::isfinite(time))
            throw std::runtime_error("SC1 mass stage: non-finite start time");
    }

    void beforeAdvance(double time, const Mass& mass)
    {
        if (stepOpen_ || blockOpen_ || time != previousTime_)
            throw std::runtime_error("SC1 mass stage: duplicate/skipped beforeAdvance");
        finite(mass);
        stepOpen_ = true;
        stepAssigned = {};
        for (auto& stage : stages) stage.second.step = {};
        for (unsigned i = 0; i < 5; ++i)
            beforeAdvanceChange[i] = mass[i] - previous_[i];
    }

    void beforeBlock(const std::string& name, const Mass& mass)
    {
        if (!stepOpen_ || blockOpen_ || name.empty())
            throw std::runtime_error("SC1 mass stage: missing step or nested/empty block");
        finite(mass);
        block_ = name;
        blockStart_ = mass;
        blockOpen_ = true;
    }

    void afterBlock(const std::string& name, const Mass& mass,
                    const Mass& allocation, const Mass& density)
    {
        if (!stepOpen_ || !blockOpen_ || name != block_)
            throw std::runtime_error("SC1 mass stage: unmatched afterBlock");
        finite(mass); finite(allocation); finite(density);
        Stage& stage = stages[name];
        for (unsigned i = 0; i < 5; ++i)
        {
            const double change = mass[i] - blockStart_[i];
            const double remainder = change - allocation[i] - density[i];
            for (Change* sum : {&stage.step, &stage.cumulative})
            {
                sum->mass[i] += change;
                sum->allocation[i] += allocation[i];
                sum->density[i] += density[i];
                sum->remainder[i] += remainder;
            }
            stepAssigned[i] += change;
            cumulativeAssigned[i] += change;
        }
        ++stage.step.calls;
        ++stage.cumulative.calls;
        blockOpen_ = false;
    }

    void finishStep(double time, double dt, const Mass& mass)
    {
        if (!stepOpen_ || blockOpen_ || !std::isfinite(time) || time <= previousTime_
            || !std::isfinite(dt) || dt <= 0
            || std::abs((time - previousTime_) - dt)
                > 1e-9*std::fmax(1.0, std::abs(time)))
            throw std::runtime_error("SC1 mass stage: incomplete/duplicate/skipped timestep");
        finite(mass);
        for (unsigned i = 0; i < 5; ++i)
        {
            stepActual[i] = mass[i] - previous_[i];
            cumulativeActual[i] = mass[i] - initial_[i];
            stepUnassigned[i] = stepActual[i] - stepAssigned[i];
            cumulativeUnassigned[i] = cumulativeActual[i] - cumulativeAssigned[i];
        }
        previous_ = mass;
        previousTime_ = time;
        ++steps;
        stepOpen_ = false;
    }
};
}
#endif
