// Narrow two-phase, nonnegative-rate experiment. No OpenFOAM dependencies.
#ifndef sc1TwoFieldRateAlgebra_H
#define sc1TwoFieldRateAlgebra_H
#include <cmath>
namespace sc1RateSource
{
struct Coefficients { double su, sp; };
inline bool coefficients(double a, double b, double qi, double qj, Coefficients& out)
{
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(qi)
        || !std::isfinite(qj) || !std::isfinite(qi+qj)
        || a<0 || a>1 || b<0 || b>1 || std::abs(a+b-1)>1e-12
        || qi<0 || qj<0) return false;
    out.su=qi;
    out.sp=-(qi+qj);
    return true;
}
}
#endif
