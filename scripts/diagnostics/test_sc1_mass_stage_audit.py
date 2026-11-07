#!/usr/bin/env python3
"""Compile/test the actual independent mass-stage arithmetic helper.

These are analytic inventory-accounting tests. They do not compile the OpenFOAM
wrapper, execute CFD, or validate any hydrodynamic or mass-transfer closure.
"""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
HEADERS = ROOT / "OFsolvers/birdmultiphaseEulerFoam/putidaFoam"

CPP = r'''
#include "sc1MassStageLedger.H"
#include <algorithm>
#include <cassert>
#include <limits>
#include <string>
using namespace sc1MassStage;
Mass m(double l, double b=0, double h=0) { return {{l,b,h,b+h,l+b+h}}; }
void near(double actual, double expected)
{
    assert(std::abs(actual-expected) < 2e-12*std::max(1.0,std::abs(expected)));
}
template<class F> void rejects(F action)
{
    bool rejected=false;
    try { action(); } catch (const std::runtime_error&) { rejected=true; }
    assert(rejected);
}
int main(int argc, char** argv)
{
    assert(argc==2);
    std::string test(argv[1]);
    if (test=="cell_product")
    {
        const auto allocation=cellChange(0.2,1000,0.3,1000,0.01);
        near(allocation[0],1); near(allocation[1],0);
        const auto density=cellChange(0.2,1.2,0.2,1.3,2);
        near(density[0],0); near(density[1],0.04);
        const auto both=cellChange(0.2,1.2,0.3,1.4,2);
        near(both[0],0.24); near(both[1],0.12);
        near(both[0]+both[1],(0.3*1.4-0.2*1.2)*2);
        // Using alphaBefore for both terms would omit the 0.04 kg cross-term.
        near(both[0]+both[1]-((0.3-0.2)*1.2+0.2*(1.4-1.2))*2,0.04);
        const auto second=cellChange(0.8,1000,0.7,998,0.001);
        near(both[0]+both[1]+second[0]+second[1],
             (0.3*1.4-0.2*1.2)*2+(0.7*998-0.8*1000)*0.001);
    }
    else if (test=="all_outer_iterations")
    {
        Ledger ledger(0,m(10)); ledger.beforeAdvance(0,m(10));
        ledger.beforeBlock("alphaPB",m(10));
        ledger.afterBlock("alphaPB",m(12),m(2),m(0));
        ledger.beforeBlock("alphaPB",m(12));
        ledger.afterBlock("alphaPB",m(11),m(-1),m(0));
        ledger.finishStep(0.1,0.1,m(11));
        const auto& stage=ledger.stages.at("alphaPB");
        assert(stage.step.calls==2 && stage.cumulative.calls==2);
        near(stage.step.mass[0],1); near(stage.step.allocation[0],1);
        near(ledger.stepActual[0],1); near(ledger.stepUnassigned[0],0);
    }
    else if (test=="unassigned_gaps")
    {
        Ledger ledger(0,m(10)); ledger.beforeAdvance(0,m(12));
        ledger.beforeBlock("pressure",m(13));
        ledger.afterBlock("pressure",m(15),m(0),m(2));
        ledger.finishStep(0.1,0.1,m(16));
        near(ledger.beforeAdvanceChange[0],2);
        near(ledger.stepActual[0],6); near(ledger.stepAssigned[0],2);
        near(ledger.stepUnassigned[0],4); near(ledger.cumulativeUnassigned[0],4);
    }
    else if (test=="step_reset_and_cumulative")
    {
        Ledger ledger(0,m(10)); ledger.beforeAdvance(0,m(10));
        ledger.beforeBlock("alphaPB",m(10));
        ledger.afterBlock("alphaPB",m(11),m(1),m(0));
        ledger.finishStep(0.1,0.1,m(11));
        ledger.beforeAdvance(0.1,m(11));
        ledger.beforeBlock("pressure",m(11));
        ledger.afterBlock("pressure",m(13),m(0),m(2));
        ledger.finishStep(0.15,0.05,m(13));
        assert(ledger.steps==2);
        assert(ledger.stages.at("alphaPB").step.calls==0);
        assert(ledger.stages.at("alphaPB").cumulative.calls==1);
        near(ledger.stepActual[0],2); near(ledger.cumulativeActual[0],3);
        near(ledger.cumulativeAssigned[0],3); near(ledger.cumulativeUnassigned[0],0);
    }
    else if (test=="phase_aggregation")
    {
        Ledger ledger(0,m(1000,0.01,0.7));
        ledger.beforeAdvance(0,m(1000,0.01,0.7));
        ledger.beforeBlock("alphaPB",m(1000,0.01,0.7));
        ledger.afterBlock("alphaPB",m(1000,0.009,0.701),m(0,-0.001,0.001),m(0));
        ledger.finishStep(0.1,0.1,m(1000,0.009,0.701));
        near(ledger.stepActual[1],-0.001); near(ledger.stepActual[2],0.001);
        near(ledger.stepActual[3],0); near(ledger.stepActual[4],0);
        for (unsigned i=0;i<5;++i)
            near(ledger.stages.at("alphaPB").step.remainder[i],0);
    }
    else if (test=="decomposition_remainder")
    {
        Ledger ledger(0,m(10)); ledger.beforeAdvance(0,m(10));
        ledger.beforeBlock("wrongSplit",m(10));
        ledger.afterBlock("wrongSplit",m(12),m(0.5),m(1));
        ledger.finishStep(0.1,0.1,m(12));
        near(ledger.stages.at("wrongSplit").step.remainder[0],0.5);
        // Complete block coverage and valid product decomposition are distinct checks.
        near(ledger.stepUnassigned[0],0);
    }
    else if (test=="invalid_order")
    {
        Ledger ledger(0,m(10));
        rejects([&]{ledger.beforeBlock("a",m(10));});
        rejects([&]{ledger.finishStep(0.1,0.1,m(10));});
        ledger.beforeAdvance(0,m(10));
        rejects([&]{ledger.beforeAdvance(0,m(10));});
        rejects([&]{ledger.beforeBlock("",m(10));});
        rejects([&]{ledger.afterBlock("a",m(10),m(0),m(0));});
        ledger.beforeBlock("a",m(10));
        rejects([&]{ledger.beforeBlock("nested",m(10));});
        rejects([&]{ledger.afterBlock("wrong",m(10),m(0),m(0));});
        rejects([&]{ledger.finishStep(0.1,0.1,m(10));});
        ledger.afterBlock("a",m(10),m(0),m(0));
        ledger.finishStep(0.1,0.1,m(10));
        rejects([&]{ledger.finishStep(0.1,0.1,m(10));});
    }
    else if (test=="time_and_nonfinite")
    {
        const double nan=std::numeric_limits<double>::quiet_NaN();
        rejects([&]{Ledger bad(nan,m(1));});
        rejects([&]{Ledger bad(0,m(nan));});
        Ledger ledger(0,m(10));
        rejects([&]{ledger.beforeAdvance(0.1,m(10));});
        ledger.beforeAdvance(0,m(10));
        rejects([&]{ledger.finishStep(0.2,0.1,m(10));});
        rejects([&]{ledger.finishStep(0,0,m(10));});
        rejects([&]{ledger.finishStep(0,1e-12,m(10));});
        rejects([&]{ledger.beforeBlock("a",m(nan));});
        ledger.beforeBlock("a",m(10));
        rejects([&]{ledger.afterBlock("a",m(10),m(nan),m(0));});
        ledger.afterBlock("a",m(10),m(0),m(0));
        ledger.finishStep(0.1,0.1,m(10));
    }
    else return 2;
}
'''


class MassStageAccountingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which("g++") or shutil.which("c++")
        if not compiler:
            raise RuntimeError("A C++ compiler is required; these tests must not silently skip")
        cls.temp = tempfile.TemporaryDirectory(prefix="sc1-mass-stage-test-")
        cls.addClassCleanup(cls.temp.cleanup)
        source = Path(cls.temp.name) / "test.C"
        cls.executable = Path(cls.temp.name) / "test"
        source.write_text(CPP)
        subprocess.run([compiler, "-std=c++11", "-Wall", "-Wextra", "-Werror",
                        "-I", str(HEADERS), str(source), "-o", str(cls.executable)],
                       check=True, capture_output=True, text=True)

    def check(self, name):
        subprocess.run([str(self.executable), name], check=True,
                       capture_output=True, text=True)

    def test_cell_product_includes_cross_term_and_unequal_cell_volumes(self):
        self.check("cell_product")

    def test_both_outer_iterations_are_counted(self):
        self.check("all_outer_iterations")

    def test_unmeasured_changes_before_during_and_after_blocks_remain_visible(self):
        self.check("unassigned_gaps")

    def test_current_step_resets_without_losing_cumulative_changes(self):
        self.check("step_reset_and_cumulative")

    def test_phase_transfer_cancels_in_gas_and_total_aggregates(self):
        self.check("phase_aggregation")

    def test_bad_product_split_is_distinct_from_missing_block(self):
        self.check("decomposition_remainder")

    def test_nested_unmatched_or_incomplete_hooks_fail(self):
        self.check("invalid_order")

    def test_invalid_timing_and_nonfinite_values_fail(self):
        self.check("time_and_nonfinite")


if __name__ == "__main__":
    print("Analytic stage-accounting tests; not OpenFOAM compilation or CFD validation.")
    unittest.main(verbosity=2)
