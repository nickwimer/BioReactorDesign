#!/usr/bin/env python3
"""Test source-hook placement and compiled exact alpha-accounting arithmetic.

No OpenFOAM library is built and no CFD/physical closure is validated here.
"""
import difflib
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from instrument_sc1_alpha_subcycles import instrument_phase_system_solve

ROOT = Path(__file__).resolve().parents[2]
HEADERS = ROOT / "OFsolvers/birdmultiphaseEulerFoam/phaseSystems/phaseSystem"

CPP = r'''
#include "sc1AlphaSubcycleLedger.H"
#include <algorithm>
#include <cassert>
#include <string>
using namespace sc1AlphaBudget;
void near(double a,double b)
{ assert(std::abs(a-b)<5e-13*std::max(1.0,std::abs(b))); }
void identities(const CellTerms& v)
{
    near(v[mulesChange],v[oldTimeReset]+v[explicitSu]+v[implicitSp]
         -v[divergence]+v[equationResidual]);
    near(v[explicitSu]+v[implicitSp],v[originalSu]+v[rawDgdt]
         +v[splitAtSourceAlpha]+v[implicitSourceChange]
         +v[sourceDecompositionRemainder]);
    near(v[originalSu]+v[rawDgdt]-v[pressureDilatation],v[couplingDifference]);
}
int main(int argc,char** argv)
{
    assert(argc==2); const std::string test(argv[1]);
    if(test=="old_time_reset")
    {
        const double dt=.01,old=.6,before=.7,su=1,sp=-2,div=.3;
        const double after=(old+dt*(su-div))/(1-dt*sp);
        const auto v=mulesCell(before,old,after,dt,div,su,sp,before,.5,-.9,.1);
        identities(v); near(v[oldTimeReset],-.1); near(v[equationResidual],0);
        assert(std::abs(v[mulesChange]-(v[explicitSu]+v[implicitSp]-v[divergence]))>.09);
    }
    else if(test=="positive_denominator_floor")
    {
        const double a=.99999,dgdt=.1,dt=1e-5;
        const double sp=-dgdt/1e-4,su=dgdt/1e-4;
        const double after=(a+dt*su)/(1-dt*sp);
        const auto v=mulesCell(a,a,after,dt,0,su,sp,a,0,dgdt,dgdt);
        identities(v); near(v[equationResidual],0);
        near(v[splitAtSourceAlpha]/dt,-.09);
        assert(v[implicitSourceChange]<0);
    }
    else if(test=="negative_denominator_floor")
    {
        const double a=1e-6,dgdt=-.01,dt=1e-5,sp=dgdt/1e-4;
        const double after=a/(1-dt*sp);
        const auto v=mulesCell(a,a,after,dt,0,0,sp,a,0,dgdt,dgdt);
        identities(v); near(v[equationResidual],0);
        near(v[splitAtSourceAlpha]/dt,.0099);
        assert(v[implicitSourceChange]>0);
    }
    else if(test=="unfloored_split")
    {
        const double a=.4,dgdt=.1,dt=.01,q=.03;
        const double sp=-dgdt/(1-a),su=q+dgdt/(1-a);
        const double after=(a+dt*su)/(1-dt*sp);
        const auto v=mulesCell(a,a,after,dt,0,su,sp,a,q,dgdt,.08);
        identities(v); near(v[equationResidual],0);
        near(v[splitAtSourceAlpha],0);
        near(v[couplingDifference],dt*.05);
        near(v[implicitSourceChange],dt*sp*(after-a));
    }
    else if(test=="subcycle_sum")
    {
        double current=.2; CellTerms sum{}; const double start=current;
        for(int sub=0;sub<6;++sub)
        {
            const double dt=.1/6,old=current,sp=-.2,su=.3,div=.04;
            const double after=(old+dt*(su-div))/(1-dt*sp);
            const auto v=mulesCell(current,old,after,dt,div,su,sp,start,.1,.16,.24);
            identities(v); for(unsigned j=0;j<nTerms;++j) sum[j]+=v[j];
            current=after;
        }
        identities(sum); near(sum[mulesChange],current-start);
        near(sum[divergence],.1*.04); near(sum[oldTimeReset],0);
    }
    else if(test=="weighted_divergence")
    {
        // Internal-face volume flux cancels globally; rho-weighted divergence
        // includes density gradients and is not an external boundary mass flux.
        const double flux=.1,dt=.01,v1=.2,v2=.3,rho1=1000,rho2=1;
        const auto c1=mulesCell(.5,.5,.5-dt*flux/v1,dt,flux/v1,0,0,.5,0,0,0);
        const auto c2=mulesCell(.5,.5,.5+dt*flux/v2,dt,-flux/v2,0,0,.5,0,0,0);
        identities(c1); identities(c2);
        near(c1[divergence]*v1+c2[divergence]*v2,0);
        near(c1[divergence]*v1*rho1+c2[divergence]*v2*rho2,.999);
    }
    else if(test=="equation_defect_visible")
    {
        const auto v=mulesCell(.5,.5,.6,.01,0,0,0,.5,0,0,0);
        identities(v); near(v[equationResidual],.1); near(v[equationResidualL1],.1);
    }
    else return 2;
}
'''


class SourceInstrumentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (HEADERS / "phaseSystemSolve.C").read_text()

    def test_every_original_line_is_retained_in_order(self):
        result = instrument_phase_system_solve(self.source)
        removed = [line for line in difflib.ndiff(self.source.splitlines(), result.splitlines())
                   if line.startswith("- ")]
        self.assertEqual(removed, [])

    def test_mules_capture_uses_actual_limited_flux_and_existing_sources(self):
        result = instrument_phase_system_solve(self.source)
        start = result.index("sc1Alpha->beforeMules(phase)")
        solve = result.index("MULES::explicitSolve", start)
        end = result.index("sc1Alpha->afterMules", solve)
        self.assertLess(start, solve)
        self.assertLess(solve, end)
        self.assertIn("(phase, alphaPhi, Sps[phase.index()], Sus[phase.index()])", result[end:])
        self.assertLess(result.rfind("phase.correctInflowOutflow(alphaPhi);", 0, start), start)
        self.assertEqual(result.count('extern "C" int sc1AlphaAuditLibraryVersion()'), 1)
        self.assertEqual(result.count("sc1Alpha->finish()"), 1)

    def test_implicit_reference_and_clamp_hooks_are_present_once(self):
        result = instrument_phase_system_solve(self.source)
        self.assertEqual(result.count("sc1Alpha->afterImplicit(phase, alphaEqn)"), 1)
        for event in ("referenceSubcycle", "referenceFinal", "clamp"):
            self.assertEqual(result.count("sc1AlphaBudget::" + event), 2)

    def test_missing_duplicate_and_preinstrumented_sources_fail(self):
        with self.assertRaises(ValueError):
            instrument_phase_system_solve(self.source.replace("phase.maxMin(0, 1);", "phase.maxMin(0, 2);"))
        with self.assertRaises(ValueError):
            instrument_phase_system_solve(self.source + "\n            phase.maxMin(0, 1);")
        with self.assertRaises(ValueError):
            instrument_phase_system_solve(instrument_phase_system_solve(self.source))


class CompiledAlphaAccountingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++") or shutil.which("c++")
        if not compiler:
            raise RuntimeError("C++ compiler required; arithmetic tests cannot silently skip")
        cls.temp = tempfile.TemporaryDirectory(prefix="sc1-alpha-accounting-")
        cls.addClassCleanup(cls.temp.cleanup)
        source = Path(cls.temp.name) / "test.C"
        cls.executable = Path(cls.temp.name) / "test"
        source.write_text(CPP)
        subprocess.run([compiler, "-std=c++11", "-Wall", "-Wextra", "-Werror",
                        "-I", str(HEADERS), str(source), "-o", str(cls.executable)],
                       check=True, capture_output=True, text=True)

    def check(self, name):
        subprocess.run([str(self.executable), name], check=True, capture_output=True, text=True)

    def test_old_time_reset_is_not_silently_treated_as_time_advancement(self):
        self.check("old_time_reset")

    def test_positive_source_floor_attenuation(self):
        self.check("positive_denominator_floor")

    def test_negative_source_floor_attenuation(self):
        self.check("negative_denominator_floor")

    def test_unfloored_split_and_implicit_source_evolution(self):
        self.check("unfloored_split")

    def test_six_subcycles_telescope_with_actual_substeps(self):
        self.check("subcycle_sum")

    def test_rho_weighted_internal_flux_is_not_boundary_mass_flux(self):
        self.check("weighted_divergence")

    def test_nonzero_equation_defect_remains_visible(self):
        self.check("equation_defect_visible")


if __name__ == "__main__":
    print("Alpha hook/arithmetic tests; not OpenFOAM compilation or CFD validation.")
    unittest.main(verbosity=2)
