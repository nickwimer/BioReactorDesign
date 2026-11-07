#!/usr/bin/env python3
"""Compile actual oxygen ledger/assembly expressions; not an OpenFOAM build."""
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT / "OFsolvers/birdmultiphaseEulerFoam/putidaFoam"


class OxygenStageAuditTests(unittest.TestCase):
    def compile_and_run(self, source):
        compiler = shutil.which("g++")
        if not compiler:
            self.skipTest("g++ unavailable; oxygen-stage regression NOT compiled")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "test.cpp").write_text(source)
            subprocess.run(
                [compiler, "-std=c++14", "-Wall", "-Wextra", "-Werror",
                 "-I", str(SOLVER), str(root / "test.cpp"), "-o", str(root / "test")],
                check=True,
            )
            subprocess.run([str(root / "test")], check=True)

    def test_compiled_stage_identity_signs_variable_dt_and_clipping(self):
        self.compile_and_run(r'''
#include "sc1BudgetLedger.H"
#include <cassert>
using namespace sc1Budget;
bool near(double a,double b) { return std::abs(a-b)<1e-12; }
int main() {
    // Known negative cells include an empty phase: count it, but do not
    // invent oxygen mass where alpha=0. Positive and exactly zero Y stay put.
    NegativeOxygenStats negative;
    negative.add(.5,2,-.01,3);  // +.03 kg when clamped
    negative.add(1,1,-.05,1);   // +.05 kg when clamped
    negative.add(0,1,-7,1);    // counted, zero phase mass
    negative.add(1,2,.5,3);
    negative.add(1,2,0,3);
    assert(negative.cells==3 && near(negative.correctionMass,.08));
    assert(near(negative.minimumY,-7));
    // Keep invalid phase weights visible instead of silently clipping alpha.
    NegativeOxygenStats badWeight;
    badWeight.add(-1,2,-.01,3);
    assert(near(badWeight.correctionMass,-.06));

    Ledger ordinary(0,1,10);
    OxygenStageLedger stages;
    Rates rates{};
    rates[oxygenExternal]=-2;
    rates[oxygenStabilization]=.4;
    rates[equationResidual]=.05;
    rates[oxygenComponentMismatch]=.02;
    rates[oxygenDdt]=-1.53;
    // start=old previous end + scheduled jump; old-time offset=.03,
    // ddt/inventory defect=.004, clipping=.08, later stage change=.01.
    ordinary.advance(.2,.2,rates,.1,.08);
    stages.advance(10.1,10.13,9.828,9.908,9.918,.2,-1.53,.08,3,.08);
    const double legacyResidual=9.918-10-Ledger::oxygenExpected(ordinary.integral)
                              -ordinary.oxygenJump-ordinary.normalization;
    assert(near(legacyResidual,.138));
    assert(near(stages.explainedResidual(ordinary.integral),legacyResidual));
    assert(near(stages.correctionDefect,0));
    assert(near(stages.negativeMassSum,.08));
    assert(stages.negativePhaseCellSamples==3);
    // Impulses must not be multiplied by dt, nor accumulated a second time
    // when the output window changes. A second, unequal dt has no stage jump.
    ordinary.resetWindow(.2,1,9.918);
    rates.fill(0); rates[oxygenExternal]=-.1; rates[oxygenDdt]=-.1;
    ordinary.advance(.3,.1,rates,0,0);
    stages.advance(9.918,9.918,9.908,9.908,9.908,.1,-.1,0,0,0);
    assert(near(stages.oldTimeOffsetSum,.03));
    assert(near(stages.ddtInventoryDefectSum,.004));
    assert(near(stages.postSpeciesChangeSum,.01));
    assert(near(stages.explainedResidual(ordinary.integral),.138));
    assert(stages.negativePhaseCellSamples==3);
    // A measured correction can differ from the predicted clip; expose it.
    OxygenStageLedger extraCorrection;
    extraCorrection.advance(1,1,1,1.09,1.09,.1,0,.08,2,.09);
    assert(near(extraCorrection.correctionDefectSum,.01));
    bool rejected=false;
    try { stages.advance(0,0,0,0,0,0,0,0,0,0); }
    catch(const std::runtime_error&) { rejected=true; }
    assert(rejected);
}
''')

    def test_actual_assembly_block_freezes_explicit_term_and_checks_matrix(self):
        source = (SOLVER / "sc1ConservationAudit.H").read_text()
        match = re.search(
            r"(ddt_\.set\(i,fvm::ddt\(alpha,rho,y\)\.ptr\(\)\);.*?)"
            r"\s*speciesSeen_\[i\]=true;", source, re.S,
        )
        self.assertIsNotNone(match)
        # A strict small affine-matrix stand-in exercises the actual assembly
        # block, including a solve that changes Y. It deliberately exposes only
        # the OF9 overload shapes used here; Kestrel must still verify OF ABI.
        self.compile_and_run(r'''
#include <array>
#include <cassert>
#include <cmath>
#include <memory>
template<class T> struct tmp;
struct fvScalarMatrix {
    double a,b; // operator& action a*Y-b, matching the OF matrix convention
    fvScalarMatrix(double aa=0,double bb=0):a(aa),b(bb) {}
    fvScalarMatrix(const tmp<fvScalarMatrix>&);
    double action(double y) const { return a*y-b; }
};
template<class T> struct tmp {
    std::shared_ptr<T> value;
    explicit tmp(const T& x):value(new T(x)) {}
    const T& operator()() const { return *value; }
    T* ptr() const { return new T(*value); }
};
fvScalarMatrix::fvScalarMatrix(const tmp<fvScalarMatrix>& x):a(x().a),b(x().b) {}
fvScalarMatrix operator+(const fvScalarMatrix& x,const fvScalarMatrix& y)
{ return {x.a+y.a,x.b+y.b}; }
fvScalarMatrix operator-(const fvScalarMatrix& x,const fvScalarMatrix& y)
{ return {x.a-y.a,x.b-y.b}; }
tmp<fvScalarMatrix> operator-(double explicitValue,const tmp<fvScalarMatrix>& x)
{ return tmp<fvScalarMatrix>({-x().a,-explicitValue-x().b}); }
struct Field { double value,old; };
struct Weighted { double value,old; };
Weighted operator*(double coefficient,const Field& rho)
{ return {coefficient*rho.value,coefficient*rho.old}; }
const double deltaT=.1;
unsigned explicitCalls=0;
namespace fvc {
double ddt(const Weighted& rho,const Field& y) {
    ++explicitCalls;
    return (rho.value*y.value-rho.old*y.old)/deltaT;
}}
namespace fvm {
tmp<fvScalarMatrix> ddt(const Weighted& rho,const Field& y) {
    return tmp<fvScalarMatrix>({rho.value/deltaT,rho.old*y.old/deltaT});
}
tmp<fvScalarMatrix> ddt(const Field& alpha,const Field& rho,const Field& y) {
    return tmp<fvScalarMatrix>({alpha.value*rho.value/deltaT,
                              alpha.old*rho.old*y.old/deltaT});
}}
struct Matrices {
    std::array<fvScalarMatrix,1> values;
    fvScalarMatrix& operator[](unsigned i) { return values[i]; }
    void set(unsigned i,fvScalarMatrix* p) { values[i]=*p; delete p; }
};
bool near(double a,double b) { return std::abs(a-b)<1e-12; }
int main() {
    const unsigned i=0;
    const Field alpha{.4,.4},rho{2,2};
    Field y{.1,.15};
    const double speciesResidualAlpha_=.01;
    Matrices ddt_,stabilization_,componentMismatch_,advection_,diffusion_,
             reaction_,transfer_,external_,equations_;
    advection_[0]={.4,.01}; diffusion_[0]={.2,.005};
    reaction_[0]={0,-.03}; transfer_[0]={-.02,-.04}; external_[0]={0,.06};
    // Independent known equation: 8.82*Y-1.245 after assembly at Y=.1.
    equations_[0]={8.82,1.245};
    const auto capture=[&]() { ASSEMBLY };
    capture();
    assert(explicitCalls==1);
    y.value=.2; // the solve updates Y after the explicit operand was frozen
    assert(near(stabilization_[0].action(y.value),-.02));
    assert(near(ddt_[0].action(y.value),.4));
    assert(near(componentMismatch_[0].action(y.value),0));
    assert(near(componentMismatch_[0].action(-.1),0));
    assert(explicitCalls==1); // no post-solve fvc reconstruction
    // An unaccounted actual equation source must remain visible with the
    // correct sign: T = physical+stabilization+actualResidual+mismatch.
    y.value=.1; equations_[0].b+=.003;
    capture();
    y.value=.2;
    assert(near(componentMismatch_[0].action(y.value),.003));
    const double physical=external_[0].action(y.value)+transfer_[0].action(y.value)
        +reaction_[0].action(y.value)-advection_[0].action(y.value)
        -diffusion_[0].action(y.value);
    assert(near(ddt_[0].action(y.value),physical+stabilization_[0].action(y.value)
        +equations_[0].action(y.value)+componentMismatch_[0].action(y.value)));
}
'''.replace("ASSEMBLY", match[1]))

    def test_species_coefficient_and_const_existing_old_time_contract(self):
        source = (SOLVER / "sc1ConservationAudit.H").read_text()
        phase = (ROOT / "OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseModel/"
                 "MultiComponentPhaseModel/MultiComponentPhaseModel.C").read_text()
        self.assertRegex(source, r'speciesResidualAlpha_\("residualAlpha",dimless,mesh_\.solverDict\("Yi"\)\)')
        self.assertRegex(phase, r'residualAlpha_\s*\(\s*"residualAlpha",\s*dimless,\s*fluid\.mesh\(\)\.solverDict\("Yi"\)')
        self.assertNotIn("p.residualAlpha()", source)
        self.assertIn("const volScalarField& oldTimeY = y;", source)
        self.assertIn("alpha.oldTime()*rho.oldTime()*oldTimeY.oldTime()", source)
        self.assertNotIn("oldTime().oldTime()", source)
        post_solve = source[source.index("void beforeSpeciesCorrection()"):]
        self.assertNotIn("fvc::ddt", post_solve)
        self.assertIn("SC1_BALANCE_METHOD version=1", source)
        self.assertIn("oxygenDecompositionVersion=2", source)


if __name__ == "__main__":
    unittest.main(verbosity=2)
