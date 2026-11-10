// Counterfactual fixed-flux source diagnostics only. Never modify alpha/flux.
#ifndef sc1AlphaSourceProbe_H
#define sc1AlphaSourceProbe_H
#include <algorithm>
#include <array>
#include <cmath>

namespace sc1SourceProbe
{
enum Metric
{
    samples, rawBelow, rawAbove, pressureBelow, pressureAbove,
    rawSourceMass, pressureSourceMass, splitMass, implicitMass,
    rawCorrectionMass, rawExcessMass, pressureExcessMass,
    identityRemainderMass, identityRemainderL1Mass, maxDtSp, nMetrics
};
using Row = std::array<double,nMetrics>;
using Table = std::array<std::array<Row,4>,2>;

inline unsigned bin(double alpha)
{
    if (alpha<0 || alpha>1) return 3;
    if (alpha<1e-4) return 0;
    if (alpha>1-1e-4) return 1;
    return 2;
}

inline Row cell(double old, double after, double dt, double div,
                double su, double sp, double assemblyAlpha,
                double baseSu, double dgdt, double pressure,
                double volume, double rho, double alphaMax)
{
    Row r{};
    const double weight=volume*rho;
    const double raw=old+dt*(baseSu+dgdt-div);
    const double physical=old+dt*(pressure-div);
    // This tolerance classifies diagnostics only; no field/source is clipped.
    const double tolerance=1e-12;
    r[samples]=1;
    r[rawBelow]=(raw < -tolerance);
    r[rawAbove]=(raw > alphaMax+tolerance);
    r[pressureBelow]=(physical < -tolerance);
    r[pressureAbove]=(physical > alphaMax+tolerance);
    r[rawSourceMass]=dt*(baseSu+dgdt)*weight;
    r[pressureSourceMass]=dt*pressure*weight;
    r[splitMass]=dt*(su+sp*assemblyAlpha-baseSu-dgdt)*weight;
    r[implicitMass]=dt*sp*(after-assemblyAlpha)*weight;
    r[rawCorrectionMass]=(raw-after)*weight;
    r[rawExcessMass]=(std::max(-raw,0.0)+std::max(raw-alphaMax,0.0))*weight;
    r[pressureExcessMass]=(std::max(-physical,0.0)
                         +std::max(physical-alphaMax,0.0))*weight;
    const double equationResidual=(after-old-dt*(su+sp*after-div))*weight;
    r[identityRemainderMass]=r[rawCorrectionMass]+r[splitMass]
                           +r[implicitMass]+equationResidual;
    r[identityRemainderL1Mass]=std::abs(r[identityRemainderMass]);
    r[maxDtSp]=std::abs(dt*sp);
    return r;
}

inline void add(Row& target, const Row& value)
{
    for (unsigned j=0;j<maxDtSp;++j) target[j]+=value[j];
    target[maxDtSp]=std::max(target[maxDtSp],value[maxDtSp]);
}

inline const char* label(unsigned metric)
{
    static const char* labels[nMetrics]={"phaseCellSubcycleSamples",
        "rawBelowSamples","rawAboveSamples","pressureBelowSamples",
        "pressureAboveSamples","rawSource_kg","pressureSource_kg",
        "split_kg","implicit_kg","rawCorrection_kg","rawExcess_kg",
        "pressureExcess_kg","identityRemainder_kg","identityRemainderL1_kg",
        "maxDtAbsSp"};
    return labels[metric];
}
}
#endif
