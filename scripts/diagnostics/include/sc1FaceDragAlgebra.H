#ifndef sc1FaceDragAlgebra_H
#define sc1FaceDragAlgebra_H
#include <algorithm>
#include <cmath>

// Solve M phi = rhs, M = [[a+k,-k],[-k,b+k]].
// Return T = inverse(M)*diag(a+k,b+k), acting on diagonal-normalized RHSs.
// Positive a,b and nonnegative k make M positive definite. Scale before
// evaluating ab+k(a+b); never subtract k*k from (a+k)*(b+k).
template<class Scalar>
bool sc1FaceDragTransform
(
    Scalar a, Scalar b, Scalar k,
    Scalar& t00, Scalar& t01, Scalar& t10, Scalar& t11
)
{
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(k)
        || a<=0 || b<=0 || k<0) return false;
    const Scalar scale=std::max(std::max(a,b),k);
    a/=scale; b/=scale; k/=scale;
    const Scalar determinant=a*b+k*(a+b);
    if (!(determinant>0) || !std::isfinite(determinant)) return false;
    t00=((b+k)/determinant)*(a+k);
    t01=(k/determinant)*(b+k);
    t10=(k/determinant)*(a+k);
    t11=((a+k)/determinant)*(b+k);
    return std::isfinite(t00) && std::isfinite(t01)
        && std::isfinite(t10) && std::isfinite(t11);
}
#endif
