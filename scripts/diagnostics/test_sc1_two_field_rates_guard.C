// Minimal API double exercising the actual guard, not an OpenFOAM header build.
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#define forAll(field,i) for (int i=0;i<int((field).size());++i)
#define FatalErrorInFunction std::cerr
namespace Foam
{
using scalar=double; using label=int; using word=std::string;
const int FatalError=1;
inline int exit(int) { throw std::runtime_error("guard rejected"); }
using std::max;using std::endl;
inline double mag(double x) { return std::abs(x); }
struct InfoDouble
{
    std::ostream& operator()() { return std::cout; }
    template<class T>std::ostream& operator<<(const T& value) { return std::cout<<value; }
};
static InfoDouble Info;
struct scalarField:std::vector<double>
{
    using std::vector<double>::vector;
    scalarField()=default;
    scalarField(const scalarField&)=default;
    void operator=(const scalarField& other) { assign(other.begin(),other.end()); }
};
template<class T> struct maxOp {};
template<class T> void reduce(T&,maxOp<T>) {}
inline double gMin(const scalarField& a) { return *std::min_element(a.begin(),a.end()); }
inline double gMax(const scalarField& a) { return *std::max_element(a.begin(),a.end()); }
struct volScalarField:scalarField
{
    using scalarField::scalarField;using Internal=scalarField;
    const scalarField& primitiveField() const { return *this; }
};
template<class T> struct tmp
{
    T value;bool present=true;
    bool valid() const { return present; }
    const T& operator()() const { return value; }
};
struct phaseModel:volScalarField
{
    int id;word label_;double density,rate=0;bool available=true,constant=true;
    phaseModel(int i,word n,double a,double r):volScalarField(1,a),id(i),label_(n),density(r) {}
    int index() const { return id; } const word& name() const { return label_; }
    bool pure() const { return constant; } bool isochoric() const { return constant; }
    bool incompressible() const { return constant; }
    tmp<volScalarField> rho() const { return {volScalarField(1,density),true}; }
    tmp<volScalarField> divU() const { return {volScalarField(1,rate),available}; }
};
struct Time { double value() const { return .001; } };
struct fvMesh
{
    bool dynamic=false;
    bool moving() const { return dynamic; } bool changing() const { return dynamic; }
    scalarField V() const { return scalarField(1,.001); }
    Time time() const { return {}; }
};
struct phaseSystem
{
    std::vector<phaseModel> p{{0,"gas",.00001,1},{1,"liquid",.99999,1000}};
    fvMesh grid; word kind="basicMultiphaseSystem";
    const std::vector<phaseModel>& phases() const { return p; }
    std::vector<int> stationaryPhases() const { return {}; }
    const fvMesh& mesh() const { return grid; }
    template<class T>T lookup(const char*) const { return kind; }
};
}
#include "sc1TwoFieldRateSplit.H"
template<class F> void rejects(F action)
{
    bool rejected=false;
    try { action(); } catch (const std::runtime_error&) { rejected=true; }
    assert(rejected);
}
int main()
{
    using namespace Foam;
    phaseSystem fluid;fluid.p[1].rate=.001;
    sc1TwoFieldRateSplit split(fluid,1,false,1,true);
    scalarField su(1,4),sp(1,0),raw(1,-1e-8);
    split.apply(fluid.p[0],sp,su,raw);
    assert(su[0]==4 && sp[0]==-.001);
    rejects([&]{split.report();});
    rejects([&]{split.apply(fluid.p[0],sp,su,raw);});
    su[0]=4;sp[0]=0;raw[0]=1e-8;
    split.apply(fluid.p[1],sp,su,raw);
    assert(su[0]==4.001 && sp[0]==-.001);
    split.report();
    auto guard=[](const phaseSystem& f){sc1TwoFieldRateSplit s(f,1,false,1,true);};
    auto bad=fluid;bad.p[0].rate=-1e-20;rejects([&]{guard(bad);});
    bad=fluid;bad.p[0].rate=std::numeric_limits<double>::quiet_NaN();rejects([&]{guard(bad);});
    bad=fluid;bad.p[0][0]=.5;rejects([&]{guard(bad);});
    bad=fluid;bad.p[0].available=false;rejects([&]{guard(bad);});
    bad=fluid;bad.p[1].density=999;rejects([&]{guard(bad);});
    bad=fluid;bad.p[0].constant=false;rejects([&]{guard(bad);});
    bad=fluid;bad.grid.dynamic=true;rejects([&]{guard(bad);});
    bad=fluid;bad.kind="other";rejects([&]{guard(bad);});
    bad=fluid;bad.p[0].label_="other";rejects([&]{guard(bad);});
    bad=fluid;bad.p.pop_back();rejects([&]{guard(bad);});
    rejects([&]{sc1TwoFieldRateSplit s(fluid,1,true,1,true);});
    rejects([&]{sc1TwoFieldRateSplit s(fluid,1,false,2,true);});
    rejects([&]{sc1TwoFieldRateSplit s(fluid,1,false,1,false);});
    sc1TwoFieldRateSplit wrong(fluid,1,false,1,true);raw[0]=1;
    rejects([&]{wrong.apply(fluid.p[0],sp,su,raw);});
    volScalarField source(1,.001);
    sc1TwoFieldRateSplit::checkOptionSource(source);
    rejects([&]{sc1TwoFieldRateSplit::checkTransfer(source);});
    source[0]=-1;rejects([&]{sc1TwoFieldRateSplit::checkOptionSource(source);});
    source[0]=0;sc1TwoFieldRateSplit::checkTransfer(source);
    fluid.p[0].available=false;fluid.p[1].available=false;
    sc1TwoFieldRateSplit startup(fluid,1,false,1,true);startup.report();
}
