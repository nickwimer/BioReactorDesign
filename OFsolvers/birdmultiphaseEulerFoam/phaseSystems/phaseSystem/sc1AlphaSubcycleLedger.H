// Fixed-mesh, geometricOneField MULES arithmetic; no OpenFOAM dependencies.
#ifndef sc1AlphaSubcycleLedger_H
#define sc1AlphaSubcycleLedger_H
#include <array>
#include <cmath>

namespace sc1AlphaBudget
{
enum Term
{
    mulesChange, oldTimeReset, explicitSu, implicitSp, divergence,
    equationResidual, equationResidualL1, originalSu, rawDgdt,
    pressureDilatation, couplingDifference, splitAtSourceAlpha,
    implicitSourceChange, sourceDecompositionRemainder,
    implicitCorrection, referenceSubcycle, clamp, referenceFinal,
    implicitDivergence, implicitEquationResidual, implicitEquationResidualL1,
    finalTransportDivergence, fluxAveragingDefect, nTerms
};
using CellTerms = std::array<double, nTerms>;

inline CellTerms mulesCell
(
    double before, double old, double after, double dt, double div,
    double su, double sp, double sourceAlpha, double baseSu,
    double dgdt, double requestedDilatation
)
{
    CellTerms v{};
    v[mulesChange] = after-before;
    v[oldTimeReset] = old-before;
    v[explicitSu] = dt*su;
    v[implicitSp] = dt*sp*after;
    v[divergence] = dt*div;
    v[equationResidual] = after-old-dt*(su+sp*after-div);
    v[equationResidualL1] = std::abs(v[equationResidual]);
    v[originalSu] = dt*baseSu;
    v[rawDgdt] = dt*dgdt;
    v[pressureDilatation] = dt*requestedDilatation;
    v[couplingDifference] = dt*(baseSu+dgdt-requestedDilatation);
    v[splitAtSourceAlpha] = dt*(su+sp*sourceAlpha-baseSu-dgdt);
    v[implicitSourceChange] = dt*sp*(after-sourceAlpha);
    v[sourceDecompositionRemainder] = v[explicitSu]+v[implicitSp]
        -v[originalSu]-v[rawDgdt]-v[splitAtSourceAlpha]-v[implicitSourceChange];
    return v;
}
}
#endif
