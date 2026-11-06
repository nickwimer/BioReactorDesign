#!/usr/bin/env python3
"""Local arithmetic/provenance/hook regressions; NOT an OpenFOAM/MPI build."""
import copy
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

from prepare_sc1_conservation_solver import ROOT, instrument, audit_make_options
from prepare_sc1_conservation_pair import enable_audit, hashes, validate_pair, CHECKPOINT
from summarize_sc1_conservation import load_log, compare, PHASES

SOLVER=ROOT/'OFsolvers/birdmultiphaseEulerFoam/putidaFoam'


class AuditTests(unittest.TestCase):
    def test_grouped_diffusion_access_contract(self):
        # Compile the actual access block against a strict registry contract.
        # Checks phase selection and fail-closed behavior; NOT an OF9 build.
        source=(SOLVER/'sc1ConservationAudit.H').read_text()
        access=re.search(
            r'const thermophysicalTransportModel& transport\s*=.*?diffusion_\.set\(i,.*?\);',
            source,re.S)
        self.assertIsNotNone(access)
        code=r'''
#include <cassert>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
struct volScalarField { double value; };
struct fvScalarMatrix { double value; };
struct MatrixTmp {
    std::unique_ptr<fvScalarMatrix> value;
    fvScalarMatrix* ptr() { return value.release(); }
};
struct RegisteredObject { virtual ~RegisteredObject() = default; };
struct thermophysicalTransportModel : RegisteredObject {
    static const std::string typeName;
    virtual MatrixTmp divj(volScalarField&) const = 0;
    virtual void correct() = 0;
};
const std::string thermophysicalTransportModel::typeName="thermophysicalTransport";
struct Model : thermophysicalTransportModel {
    double coefficient;
    mutable int calls=0;
    int corrections=0;
    explicit Model(double c) : coefficient(c) {}
    MatrixTmp divj(volScalarField& y) const override {
        ++calls;
        return {std::unique_ptr<fvScalarMatrix>(new fvScalarMatrix{coefficient*y.value})};
    }
    void correct() override { ++corrections; }
};
struct IOobject {
    static std::string groupName(const std::string& base, const std::string& group) {
        return base+"."+group;
    }
};
struct Mesh {
    std::map<std::string,const RegisteredObject*> objects;
    template<class T> const T& lookupObject(const std::string& name) const {
        return dynamic_cast<const T&>(*objects.at(name));
    }
};
struct Flux {
    std::string phase;
    const Flux& operator()() const { return *this; }
    const std::string& group() const { return phase; }
};
// Deliberately no divj on this abstract-phase stand-in.
struct Phase { Flux alphaRhoPhi() const { return {name}; } std::string name; };
struct Matrices {
    std::unique_ptr<fvScalarMatrix> values[3];
    void set(int i,fvScalarMatrix* m) { values[i].reset(m); }
};
struct Audit {
    const Mesh& mesh_;
    Matrices diffusion_;
    explicit Audit(const Mesh& mesh) : mesh_(mesh) {}
    void capture(Phase& p,volScalarField& y,int i) {
        const Flux flux=p.alphaRhoPhi();
        ACCESS
    }
};
int main() {
    Mesh mesh;
    Model liquid(2), bubbles(3), headspace(5);
    Model* models[]={&liquid,&bubbles,&headspace};
    Phase phases[]={{"liquid"},{"bubbleGas"},{"headspaceGas"}};
    for(int i=0;i<3;++i) mesh.objects["thermophysicalTransport."+phases[i].name]=models[i];
    Audit audit(mesh);
    volScalarField y{0.23};
    for(int i=0;i<3;++i) {
        audit.capture(phases[i],y,i);
        assert(audit.diffusion_.values[i]->value==models[i]->coefficient*y.value);
        assert(models[i]->calls==1 && models[i]->corrections==0);
    }
    assert(y.value==0.23);
    mesh.objects.erase("thermophysicalTransport.bubbleGas");
    bool rejected=false;
    try { audit.capture(phases[1],y,1); } catch(const std::out_of_range&) { rejected=true; }
    assert(rejected);
    RegisteredObject wrongType;
    mesh.objects["thermophysicalTransport.bubbleGas"]=&wrongType;
    rejected=false;
    try { audit.capture(phases[1],y,1); } catch(const std::bad_cast&) { rejected=true; }
    assert(rejected);
}
'''.replace('ACCESS',access[0])
        compiler=shutil.which('g++')
        if not compiler:
            self.skipTest('g++ unavailable; transport access contract NOT compiled')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'test.cpp').write_text(code)
            subprocess.run([compiler,'-std=c++14','-Wall','-Wextra','-Werror',str(p/'test.cpp'),'-o',str(p/'test')],check=True)
            subprocess.run([str(p/'test')],check=True)

    def test_audit_build_options_preserve_local_settings(self):
        local='LIB_SRC=/test/OpenFOAM/src\nEXE_INC = -DLOCAL_BIOLOGY -I/local/headers\nEXE_LIBS = -L/local/lib -llocalModel\n'
        modified=audit_make_options(local)
        self.assertTrue(modified.startswith(local))
        old=subprocess.check_output(['git','show','19410f47bbf398da0bbf8f8ccce931449de05c5d:OFsolvers/birdmultiphaseEulerFoam/putidaFoam/Make/options'],cwd=ROOT,text=True)
        self.assertEqual(audit_make_options(old),(SOLVER/'Make/options').read_text())
        make=shutil.which('make')
        if not make:
            self.skipTest('make unavailable; effective appended flags NOT checked')
        recipe="\n.PHONY: sc1-print\nsc1-print:\n\t@printf '%s\\n' '$(EXE_INC)' '$(EXE_LIBS)'\n"
        with tempfile.TemporaryDirectory() as tmp:
            makefile=Path(tmp)/'Makefile'
            makefile.write_text(modified+recipe)
            result=subprocess.check_output([make,'--no-print-directory','-s','-f',str(makefile),'sc1-print'],text=True)
        includes,libraries=result.splitlines()
        self.assertEqual(includes.split(),['-DLOCAL_BIOLOGY','-I/local/headers','-I/test/OpenFOAM/src/ThermophysicalTransportModels/lnInclude'])
        self.assertEqual(libraries.split(),['-L/local/lib','-llocalModel','-lthermophysicalTransportModels'])

    def test_final_report_before_nominal_end_time(self):
        # Compile the actual report predicate from the solver against the
        # authoritative OF9 Time::running() stopping rule. This checks a real
        # lost-final-interval case, not just presence of a source-code string.
        source=(SOLVER/'putidaFoam.C').read_text()
        predicate=re.search(r'sc1Conservation->sample\s*\((.*?)\);',source,re.S)
        self.assertIsNotNone(predicate)
        code=r'''
#include <cassert>
struct Runtime {
    double t, end, dt;
    bool write;
    bool writeTime() const { return write; }
    // OpenFOAM-9 src/OpenFOAM/db/Time/Time.C: no function-object side effects.
    bool running() const { return t < end - 0.5*dt; }
};
bool report(const Runtime& runTime) { return PREDICATE; }
int main() {
    assert(report({0.96, 1.0, 0.1, false})); // terminal, below nominal end
    assert(report({0.95, 1.0, 0.1, false})); // stopping threshold itself
    assert(!report({0.94, 1.0, 0.1, false}));
    assert(report({0.94, 1.0, 0.1, true}));  // ordinary selected write
    assert(report({1.01, 1.0, 0.1, false})); // terminal overshoot
}
'''.replace('PREDICATE',predicate[1].strip())
        compiler=shutil.which('g++')
        if not compiler:
            self.skipTest('g++ unavailable; final-report predicate NOT compiled')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'test.cpp').write_text(code)
            subprocess.run([compiler,'-std=c++11','-Wall','-Wextra','-Werror',str(p/'test.cpp'),'-o',str(p/'test')],check=True)
            subprocess.run([str(p/'test')],check=True)

    def test_compiled_production_ledger(self):
        # Tests the SAME C++ accumulator included by the solver, not a Python
        # reimplementation. An impulsive rate between writes exposes the old
        # sparse trapezoidal error. Also check variable dt, signs and resets.
        code=r'''
#include "sc1BudgetLedger.H"
#include <cassert>
using namespace sc1Budget;
bool near(double a,double b) { return std::abs(a-b)<1e-12; }
int main() {
    Ledger l(31.2,1,.23);
    Rates r{}; r[massExternal]=.01; r[massTransfer]=-.004;
    r[massIn]=1.14; r[massOut]=1.146;
    r[oxygenExternal]=.0023; r[oxygenTransfer]=-.001;
    r[oxygenReaction]=-.0005; r[oxygenAdvection]=.0007;
    r[oxygenDiffusion]=.0001;
    l.advance(31.21,.01,r,0,0);
    assert(near(Ledger::massExpected(l.integral),0));
    assert(near(Ledger::oxygenExpected(l.integral),0));
    r[massTransfer]=.096; // pulse missed by sparse output times
    l.advance(31.24,.03,r,-.002,.0001);
    r[massTransfer]=-.004;
    l.advance(31.3,.06,r,0,0);
    assert(near(Ledger::massExpected(l.integral),.003));
    assert(near(l.integral[massIn],.114));
    assert(near(l.integral[massOut],.1146));
    assert(l.steps==3 && near(l.oxygenJump,-.002));
    assert(near(l.normalization,.0001));
    l.resetWindow(31.3,1.003,.2281);
    l.advance(31.35,.05,r,0,0);
    assert(l.windowSteps==1 && l.steps==4);
    assert(near(Ledger::massExpected(l.window),0));
    assert(near(Ledger::massExpected(l.integral),.003));
    bool rejected=false;
    try { l.advance(31.35,.05,r,0,0); } catch(const std::runtime_error&) { rejected=true; }
    assert(rejected); rejected=false;
    try { l.advance(31.45,.05,r,0,0); } catch(const std::runtime_error&) { rejected=true; }
    assert(rejected);
    // Equal/opposite phase transfer cancels in the combined inventory.
    Rates b{},h{}; b[massTransfer]=.002; h[massTransfer]=-.002;
    b[oxygenTransfer]=.002*.23; h[oxygenTransfer]=-.002*.23;
    assert(near(Ledger::massExpected(b)+Ledger::massExpected(h),0));
    assert(near(Ledger::oxygenExpected(b)+Ledger::oxygenExpected(h),0));
    // Donor lag would be visible as an interphase O2 source imbalance.
    h[oxygenTransfer]=-.002*.229;
    assert(near(Ledger::oxygenExpected(b)+Ledger::oxygenExpected(h),.000002));
}
'''
        compiler=shutil.which('g++')
        if not compiler:
            self.skipTest('g++ unavailable; compiled ledger test NOT run')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'test.cpp').write_text(code)
            subprocess.run([compiler,'-std=c++11','-Wall','-Wextra','-Werror','-I',str(SOLVER),str(p/'test.cpp'),'-o',str(p/'test')],check=True)
            subprocess.run([str(p/'test')],check=True)

    def test_instrumentation_preserves_local_physics(self):
        # Use the parent revision to test installing into the actual old code.
        old=lambda name: subprocess.check_output(['git','show',f'8e10fc70f0048b1e414d7d92701781c5e6122ef0:OFsolvers/birdmultiphaseEulerFoam/putidaFoam/{name}'],cwd=ROOT,text=True)
        main,species=old('putidaFoam.C'),old('YEqns.H')
        modified,ys=instrument(main,species)
        self.assertEqual(modified,(SOLVER/'putidaFoam.C').read_text())
        self.assertEqual(ys,(SOLVER/'YEqns.H').read_text())
        self.assertEqual(modified.count('->sample'),1)
        self.assertEqual(modified.count('->beforePB'),2)
        self.assertIn('lookupOrDefault<Switch>("sc1ConservationAudit", false)',modified)
        self.assertEqual(main[main.index('double pH_highRes_Fix'):main.index('    while (pimple.run')],
                         modified[modified.index('double pH_highRes_Fix'):modified.index('    autoPtr<sc1ConservationAudit>')])
        for marker in ('#include "bioReact.H"','#include "phReact.H"','reaction_time += reaction_update_time.value();'):
            self.assertEqual(main.count(marker),modified.count(marker))
        self.assertEqual(ys.count('fvModels.source(alpha, rho, Y[i])'),1)
        with self.assertRaises(ValueError): instrument(modified,ys)
        with self.assertRaises(ValueError): instrument(main.replace('fluid.solve(rAUs, rAUfs);','fluid.solveChanged();'),species)

    def test_pimple_comments_and_duplicate_guard(self):
        text='// PIMPLE { fake }\nPIMPLE { nOuterCorrectors 2; /* } */ }\n'
        result=enable_audit(text)
        self.assertIn('nOuterCorrectors 2;',result)
        self.assertEqual(result.count('sc1ConservationAudit true;'),1)
        with self.assertRaises(ValueError): enable_audit(result)

    def test_matched_initial_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); c=root/'control'; n=root/'normal'
            for p in (c,n):
                for part in ('constant','system','processors52/constant',f'processors52/{CHECKPOINT}'):
                    (p/part).mkdir(parents=True)
                (p/'processors52/constant/mesh').write_bytes(b'\0mesh')
                for field in ('alpha.liquid','alpha.bubbleGas','alpha.headspaceGas','O2.liquid','O2.bubbleGas','O2.headspaceGas','f.bubbleGas'):
                    (p/'processors52'/CHECKPOINT/field).write_bytes(b'\0same binary field')
                (p/'constant/microbeUpdateProperties').write_text('reaction_update_time 14400;')
                (p/'system/fvSolution').write_text('PIMPLE {}')
            common='maxGlobalEntrainedKgS 0.002; nucleationDiameter 0.002; foamAreaMultiplier 0;'
            (c/'constant/phaseProperties').write_text(common+'entrainmentEfficiency 0; requireDownwardLiquid true;')
            (n/'constant/phaseProperties').write_text(common+'entrainmentEfficiency 0.01; requireDownwardLiquid false; entrainmentVelocityCriterion normalGasIntoLiquid;')
            validate_pair(c,n,CHECKPOINT)
            (n/'processors52'/CHECKPOINT/'O2.headspaceGas').write_bytes(b'different oxygen')
            with self.assertRaises(ValueError): validate_pair(c,n,CHECKPOINT)

    def test_symlink_directory_not_silently_unhashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'data').mkdir(); (root/'data/f').write_text('field')
            (root/'alias').symlink_to(root/'data',target_is_directory=True)
            with self.assertRaises(ValueError): hashes(root)

    def test_log_missing_corrupt_or_restarted(self):
        lines=['SC1_BALANCE_METHOD version=1 exactDiscreteBudget=0']
        for p in PHASES:
            lines.append(f'SC1_BALANCE_START phase={p} t_s=31.2 mass_kg=1 oxygen_kg=0.23')
        for p in PHASES:
            lines.append(f'SC1_BALANCE phase={p} t_s=31.3 start_s=31.2 steps=100 mass_kg=1.001 oxygen_kg=0.23 massResidual_kg=0 oxygenResidual_kg=0 massExternalIntegral_kg=0.001 massTransferIntegral_kg=0 massInIntegral_kg=0.114 massOutIntegral_kg=0.114 oxygenTransferIntegral_kg=0 normalization_kg=0 equationResidualIntegral_kg=0')
        lines.append('End')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'log'; p.write_text('\n'.join(lines))
            valid=load_log(p)
            self.assertEqual(compare(valid,valid)[0],(31.3,0,0))
            changed=copy.deepcopy(valid); changed[0]['gas']['oxygen_kg']=0.24
            with self.assertRaises(ValueError): compare(valid,changed)
            for corrupt in ('\n'.join(lines[:-1]),'BUBBLE_MASS_AUDIT inventory_kg=1\nEnd',
                            '\n'.join(lines).replace('mass_kg=1.001','mass_kg=1.002'),
                            '\n'.join(lines+lines), '\n'.join(lines).replace('oxygenResidual_kg=0','oxygenResidual_kg=nan')):
                p.write_text(corrupt)
                with self.assertRaises(ValueError): load_log(p)


if __name__=='__main__':
    unittest.main(verbosity=2)
