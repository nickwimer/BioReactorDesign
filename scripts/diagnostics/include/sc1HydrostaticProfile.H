// Fixed flat-interface reference; SI units. No evolving flow field is reset.
#ifndef sc1HydrostaticProfile_H
#define sc1HydrostaticProfile_H
#include <algorithm>
namespace sc1Hydro
{
inline double alpha(double y, double h, double width)
{
    return std::max(0.0, std::min(1.0, (h+width/2-y)/width));
}
inline double pressure(double y, double h, double width)
{
    const double x=y-h;
    const double integral=x < -width/2 ? x+width/8
        : (x > width/2 ? width/8 : x/2-x*x/(2*width));
    return 100000-9.81*y+999*9.81*(width/8-integral);
}
}
#endif
