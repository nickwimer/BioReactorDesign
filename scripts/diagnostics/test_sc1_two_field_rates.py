"""Compile the actual source kernel; test private preparation and mocked CFD flow."""
import copy
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import sc1_two_field as base
import sc1_two_field_rates as runner
import sc1_two_field_rates_instrument as instrument
from sc1_two_field_rates_report import summarize
from sc1_two_field_report import summarize as baseline_report
from sc1_two_field_bounds_report import summarize as bound_report
from sc1_two_field_bounds_cases import contract
import test_sc1_two_field as baseline_tests
import test_sc1_two_field_bounds as bound_tests


def log_text(name):
    near=name in runner.BOUND_CASES
    rows=bound_tests.records(name) if near else baseline_tests.records(name)
    if near:
        for kind,row in rows:
            if kind=='SC1_PHASE_TERM':row['volume_m3']=row['frozenRhoMass_kg']=0
    text=baseline_tests.render(rows).replace('SC1_TWO_LIBRARY version=1',
        'SC1_RATE_LIBRARY version=1\nSC1_TWO_LIBRARY version=1')
    steps=contract(name).steps if near else 100
    dt=contract(name).dt if near else .001
    states=[]
    for step in range(steps):
        for call in (1,2):
            valid=int(step>0 or call>1)
            rate=.001 if valid and not name.startswith('closed_') else 0
            states.append(f'SC1_RATE_STATE t_s={round((step+1)*dt,12)} call={call} '
                f'valid={valid} applied={2*valid} qGasMin=0 qGasMax=0 '
                f'qLiquidMin={rate} qLiquidMax={rate} alphaSumMaxError=0 rawIdentityScaled=0')
    return text.replace('\nEnd\n','\n'+'\n'.join(states)+'\nEnd\n')


def original_reports():
    result={}
    for n in runner.CASES:
        if n in runner.BASELINE_CASES:
            r=baseline_report(baseline_tests.render(baseline_tests.records(n)),n)
        else:r=bound_report(baseline_tests.render(bound_tests.records(n)),n)
        r['log_sha256']='original-'+n;result[n]=r
    return result


class Algebra(unittest.TestCase):
    def test_actual_guard_header_with_minimal_api_contract(self):
        compiler=shutil.which('g++') or shutil.which('c++')
        if not compiler:self.skipTest('C++ compiler unavailable')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for name in ('fvCFD.H','phaseSystem.H'):(p/name).write_text('// API double declared in test translation unit.\n')
            done=subprocess.run([compiler,'-std=c++14','-Wall','-Wextra','-Werror',
                '-I'+str(base.HERE/'include'),'-I'+str(p),
                str(base.HERE/'test_sc1_two_field_rates_guard.C'),'-o',str(p/'guard')],text=True,capture_output=True)
            self.assertEqual(done.returncode,0,done.stderr)
            done=subprocess.run([str(p/'guard')],text=True,capture_output=True)
            self.assertEqual(done.returncode,0,done.stderr)
            self.assertEqual(done.stdout.count('SC1_RATE_STATE'),2)

    def test_compiled_actual_kernel_bounds_identity_and_rejections(self):
        compiler=shutil.which('g++') or shutil.which('c++')
        if not compiler:self.skipTest('C++ compiler unavailable')
        code=r'''
#include "sc1TwoFieldRateAlgebra.H"
#include <cassert>
#include <limits>
#include <algorithm>
int main()
{
    using namespace sc1RateSource;
    const double rates[]={0,1e-12,.001,1,100};
    const double dt[]={0,1./12000,1./6000,.001,1,1000};
    for (int k=0;k<=1000;++k)
    for (double qi:rates) for (double qj:rates)
    {
        const double a=k/1000.; Coefficients c;
        assert(coefficients(a,1-a,qi,qj,c));
        assert(std::abs(c.su+c.sp*a-((1-a)*qi-a*qj))<=1e-12*(1+qi+qj));
        for (double h:dt)
        {
            const double x=(a+h*c.su)/(1-h*c.sp);
            const double y=(1-a+h*qj)/(1+h*(qi+qj));
            assert(x>=0 && x<=1 && std::abs(x+y-1)<=1e-14);
        }
    }
    for (double a:{0.,1e-8,1e-5,.99999,1.})
    {
        Coefficients c; assert(coefficients(a,1-a,.001,0,c));
        assert(std::isfinite(c.su) && std::isfinite(c.sp));
        assert(std::abs(c.su+c.sp*a-(1-a)*.001)<1e-18);
    }
    Coefficients c;const double a=.99999,d=(1-a)*.001;
    assert(coefficients(a,1-a,.001,0,c));
    const double legacy=d*(1-a)/std::max(1-a,1e-4);
    assert(std::abs(legacy/d-.1)<1e-10);
    assert(std::abs((c.su+c.sp*a)/d-1)<1e-10);
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double inf=std::numeric_limits<double>::infinity();
    const double big=std::numeric_limits<double>::max();
    for (double bad:{nan,inf,-inf,-1.})
    {
        assert(!coefficients(.5,.5,bad,0,c));
        assert(!coefficients(.5,.5,0,bad,c));
        assert(!coefficients(bad,.5,0,0,c));
    }
    assert(!coefficients(.5,.6,0,0,c));
    assert(!coefficients(1.000001,-.000001,0,0,c));
    assert(!coefficients(.5,.5,big,big,c));
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'check.C').write_text(code)
            done=subprocess.run([compiler,'-std=c++14','-Wall','-Wextra','-Werror',
                '-I'+str(base.HERE/'include'),str(p/'check.C'),'-o',str(p/'check')],text=True,capture_output=True)
            self.assertEqual(done.returncode,0,done.stderr)
            subprocess.run([str(p/'check')],check=True)


class Transformation(unittest.TestCase):
    def test_only_split_replaced_original_observer_and_equations_preserved(self):
        source=(base.SOURCE/'phaseSystems/phaseSystem/phaseSystemSolve.C').read_text()
        result=instrument.phase_source(source)
        restored=result.replace(instrument.MARKER,base.MARKER)
        restored=restored.replace('\n#include "sc1TwoFieldRateSplit.H"','')
        restored=restored.replace('    sc1TwoFieldRateSplit sc1Rates(*this, sc1PhaseTotals.call, LTS, nAlphaCorr, sc1Observer.valid());\n\n','')
        restored=restored.replace('                sc1Rates.apply(phase, Sp, Su, dgdt);',instrument.OLD_SPLIT)
        restored=restored.replace('        sc1Rates.report();\n\n','')
        self.assertEqual(restored,base.phase_source(source))
        with self.assertRaises(ValueError):instrument.phase_source(result)
        with self.assertRaises(ValueError):instrument.phase_source(source.replace('1e-4','1e-5'))

    def test_pressure_source_still_evaluated_once_and_markers_are_distinct(self):
        source=(base.SOURCE/'multiphaseEulerFoam/pEqnComps.H').read_text()
        text=instrument.pressure_source(source)
        self.assertEqual(text.count('fvModels.source(alpha, rho)'),1)
        self.assertLess(text.index('checkOptionSource'),text.index('sc1Stages.optionSource'))
        self.assertIn('pEqnComp -= sc1OptionSource()/rho;',text)
        self.assertIn('checkTransfer(dmdts[phasei])',text)
        main=instrument.solver_source((base.SOURCE/'multiphaseEulerFoam/multiphaseEulerFoam.C').read_text())
        self.assertNotIn(base.MARKER,main)
        self.assertIn(instrument.MARKER+'()!=1',main)
        self.assertIn('SC1_RATE_LIBRARY version=1',main)

    def test_preparation_preserves_originals_and_separates_build_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);installed=p/'installed';installed.mkdir()
            lib=installed/'libbirdphaseSystem.so';lib.write_bytes(b'untouched')
            before=base.selected_hashes(base.SOURCE)
            out=base.prepare_build(p/'rates',installed=installed,recipe=instrument.RECIPE)
            self.assertEqual(before,base.selected_hashes(base.SOURCE));self.assertEqual(lib.read_bytes(),b'untouched')
            m=json.loads((out/'build_manifest.json').read_text())
            self.assertEqual(m['variant'],instrument.TAG)
            self.assertIn('/bin/'+instrument.APP,(out/'src/multiphaseEulerFoam/Make/files').read_text())
            for header in ('sc1TwoFieldRateAlgebra.H','sc1TwoFieldRateSplit.H'):
                self.assertEqual((out/'src/phaseSystems/phaseSystem'/header).read_bytes(),(base.HERE/'include'/header).read_bytes())
            with self.assertRaisesRegex(ValueError,'variant'):base.verify_build(out)
            with self.assertRaisesRegex(ValueError,'Missing/redirected private binary'):
                base.verify_build(out,recipe=instrument.RECIPE)
            m.pop('variant');(out/'build_manifest.json').write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError,'variant'):base.verify_build(out,recipe=instrument.RECIPE)

    def test_all_eight_cases_only_change_executable(self):
        for name in runner.CASES:
            original=runner.baseline_files(name) if name in runner.BASELINE_CASES else runner.bound_files(name)
            files=runner.case_files(name)
            self.assertEqual([p for p in files if files[p]!=original[p]],['system/controlDict'])
            self.assertEqual(files['system/controlDict'].replace(instrument.APP,base.APP),original['system/controlDict'])


class Reporting(unittest.TestCase):
    def test_all_contracts_pass_without_floor_effect_or_relaxed_mass_tolerances(self):
        for name in runner.CASES:
            r=summarize(log_text(name),name)
            self.assertTrue(r['screen_pass'],r['errors']);self.assertTrue(r['diagnostic_pass'])
            self.assertIsNone(r['sourceFloor'])
            if name in runner.BOUND_CASES:
                self.assertFalse(r['sourceSplitObserved'])
                self.assertEqual(r['phases']['gas']['screenTolerance_kg'],1e-15)
                self.assertAlmostEqual(r['phases']['liquid']['screenTolerance_kg'],1.01e-10)

    def test_runtime_marker_missing_duplicate_and_bad_coverage_rejected(self):
        name=runner.CASES[4];text=log_text(name)
        with self.assertRaises(ValueError):summarize(text.replace('SC1_RATE_LIBRARY version=1\n',''),name)
        with self.assertRaises(ValueError):summarize('SC1_RATE_LIBRARY version=1\n'+text,name)
        line=next(l for l in text.splitlines() if l.startswith('SC1_RATE_STATE'))
        for bad in (text.replace(line+'\n','',1),text+line+'\n',text.replace(line,line.replace('t_s=0.001','t_s=0.0011'),1)):
            with self.assertRaises(ValueError):summarize(bad,name)

    def test_rate_guards_and_remaining_floor_defect_fail(self):
        name=runner.CASES[4];text=log_text(name)
        for before,after in (('applied=2','applied=1'),('qGasMax=0','qGasMax=0.001'),
            ('qLiquidMin=0.001','qLiquidMin=-0.001'),('rawIdentityScaled=0','rawIdentityScaled=0.01'),
            ('alphaSumMaxError=0','alphaSumMaxError=0.01')):
            self.assertFalse(summarize(text.replace(before,after,1),name)['diagnostic_pass'])
        old=baseline_tests.render(bound_tests.records(name))
        states='\n'.join(l for l in text.splitlines() if l.startswith('SC1_RATE_STATE'))
        r=summarize('SC1_RATE_LIBRARY version=1\n'+old+states+'\n',name)
        self.assertFalse(r['screen_pass'])
        self.assertIn('liquid: rate split does not preserve the assembled coupled source',r['errors'])

    def test_small_mass_failure_remains_failure_with_valid_diagnostics(self):
        name=runner.CASES[4];text=log_text(name)
        line=next(l for l in text.splitlines() if l.startswith('SC1_TWO_STEP phase=gas'))
        bad=line.replace('boundaryOut_kg=0','boundaryOut_kg=1e-13').replace('physicalResidual_kg=0','physicalResidual_kg=1e-13')
        r=summarize(text.replace(line,bad,1),name)
        self.assertFalse(r['screen_pass']);self.assertTrue(r['diagnostic_pass'],r['errors'])


class Runner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.controls=original_reports()

    @staticmethod
    def fake_execute(command,case,log,env):
        if Path(command[0]).name==instrument.APP:log.write_text(log_text(case.name))
        else:baseline_tests.ClosedContinuation.fake_execute(command,case,log,env)

    def prepared(self,p):
        (p/'build').mkdir();(p/'control').mkdir()
        return runner.prepare(p/'build',p/'control',p/'cases')

    def test_full_regression_and_readonly_comparison_preserve_all_controls(self):
        with tempfile.TemporaryDirectory() as tmp,redirect_stdout(io.StringIO()):
            p=Path(tmp);baseline=p/'baseline'
            m=bound_tests.Isolation().seed_baseline(baseline)
            with patch.object(base,'verify_build',return_value={'assets':m['assets']}):
                near=runner.bounds.prepare(baseline,p/'control')
                with patch.object(base,'execute',side_effect=bound_tests.Isolation.fake_execute):
                    runner.bounds.run(near,'coarse');runner.bounds.run(near,'fine')
                before=baseline_tests.ClosedContinuation.hashes(near)
                old=baseline_tests.ClosedContinuation.hashes(baseline)
                (p/'build').mkdir();root=runner.prepare(p/'build',near,p/'rates')
                with patch.object(base,'execute',side_effect=self.fake_execute) as execute:
                    result=runner.run_all(root)
                self.assertTrue(result['regression_pass']);self.assertEqual(len(execute.call_args_list),32)
                solver_calls=[c for c in execute.call_args_list if Path(c.args[0][0]).name==instrument.APP]
                self.assertEqual([c.args[1].name for c in solver_calls],list(runner.CASES))
                self.assertTrue(all(c.args[0][0]==p/'build/bin'/instrument.APP for c in solver_calls))
                after=baseline_tests.ClosedContinuation.hashes(root)
                manifest,original=runner.verify(root)
                self.assertEqual(result,runner.comparison(root,manifest,original))
                self.assertEqual(after,baseline_tests.ClosedContinuation.hashes(root))
                with patch.object(base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'attempted/modified'):runner.run_all(root)
                    execute.assert_not_called()
                self.assertEqual(before,baseline_tests.ClosedContinuation.hashes(near))
                self.assertEqual(old,baseline_tests.ClosedContinuation.hashes(baseline))

    def test_modified_input_and_changed_control_stop_before_execution(self):
        with tempfile.TemporaryDirectory() as tmp,redirect_stdout(io.StringIO()), \
                patch.object(base,'verify_build',return_value={'assets':{}}), \
                patch.object(runner,'controls',return_value=self.controls) as controls:
            p=Path(tmp);root=self.prepared(p)
            last=root/runner.CASES[-1]/'system/fvSolution';text=last.read_text();last.write_text('changed')
            with patch.object(base,'execute') as execute:
                with self.assertRaisesRegex(ValueError,'inputs changed'):runner.run_all(root)
                execute.assert_not_called()
            last.write_text(text)
            changed=copy.deepcopy(self.controls);changed[runner.CASES[0]]['log_sha256']='changed';controls.return_value=changed
            with patch.object(base,'execute') as execute:
                with self.assertRaisesRegex(ValueError,'control logs changed'):runner.run_all(root)
                execute.assert_not_called()

    def test_failed_closed_screen_stops_before_source_pair(self):
        def fail(command,case,log,env):
            self.fake_execute(command,case,log,env)
            if Path(command[0]).name==instrument.APP:
                text=log.read_text();line=next(l for l in text.splitlines() if l.startswith('SC1_TWO_STEP phase=gas'))
                bad=line.replace('boundaryOut_kg=0','boundaryOut_kg=1e-8').replace('physicalResidual_kg=0','physicalResidual_kg=1e-8')
                log.write_text(text.replace(line,bad,1))
        with tempfile.TemporaryDirectory() as tmp,redirect_stdout(io.StringIO()), \
                patch.object(base,'verify_build',return_value={'assets':{}}), \
                patch.object(runner,'controls',return_value=self.controls):
            root=self.prepared(Path(tmp))
            with patch.object(base,'execute',side_effect=fail) as execute:
                with self.assertRaisesRegex(ValueError,'stage failed screen_pass'):runner.run_all(root)
            self.assertEqual(len(execute.call_args_list),8)
            self.assertFalse((root/runner.CASES[2]/'attempt.json').exists())

    def test_explicit_run_and_allocation_required(self):
        with patch('sys.stderr',new=io.StringIO()),patch.object(runner,'run_all') as run,patch.object(base,'build') as build:
            with self.assertRaises(SystemExit):runner.main(['run-all','--root','/unused'])
            with patch.dict(os.environ,{},clear=True):
                for args in (['run-all','--root','/unused','--run'],['build','--out','/unused']):
                    with self.assertRaisesRegex(ValueError,'active Kestrel'):runner.main(args)
            run.assert_not_called();build.assert_not_called()

    def test_near_mass_failures_allow_fine_but_invalid_diagnostics_stop(self):
        for invalid in (False,True):
            def fail(command,case,log,env):
                self.fake_execute(command,case,log,env)
                if Path(command[0]).name==instrument.APP and case.name.startswith('bound_coarse_'):
                    text=log.read_text()
                    if invalid:text=text.replace('rawIdentityScaled=0','rawIdentityScaled=0.01',1)
                    else:
                        line=next(l for l in text.splitlines() if l.startswith('SC1_TWO_STEP phase=gas'))
                        bad=line.replace('boundaryOut_kg=0','boundaryOut_kg=1e-13').replace('physicalResidual_kg=0','physicalResidual_kg=1e-13')
                        text=text.replace(line,bad,1)
                    log.write_text(text)
            with self.subTest(invalid_diagnostics=invalid),tempfile.TemporaryDirectory() as tmp, \
                    redirect_stdout(io.StringIO()),patch.object(base,'verify_build',return_value={'assets':{}}), \
                    patch.object(runner,'controls',return_value=self.controls):
                root=self.prepared(Path(tmp))
                with patch.object(base,'execute',side_effect=fail) as execute:
                    if invalid:
                        with self.assertRaisesRegex(ValueError,'stage failed diagnostic_pass'):runner.run_all(root)
                        self.assertEqual(len(execute.call_args_list),24)
                        self.assertFalse((root/runner.CASES[-1]/'attempt.json').exists())
                    else:
                        result=runner.run_all(root)
                        self.assertFalse(result['regression_pass'])
                        self.assertTrue(result['cases'][runner.CASES[-1]]['screen_pass'])
                        self.assertEqual(len(execute.call_args_list),32)
        with patch.dict(os.environ,{'SLURM_JOB_ID':'test'}),patch.object(runner,'run_all',return_value={'regression_pass':False}):
            self.assertEqual(runner.main(['run-all','--root','/unused','--run']),2)


if __name__=='__main__':unittest.main()
