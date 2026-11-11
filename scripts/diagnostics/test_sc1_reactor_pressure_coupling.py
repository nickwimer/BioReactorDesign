"""Independent face-matrix algebra, private build/case isolation and acceptance."""
from contextlib import redirect_stdout
from decimal import Decimal, localcontext
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_pressure_coupling as candidate
import sc1_reactor_pressure_coupling_instrument as instrument
import sc1_reactor_failure as failure
import sc1_reactor_watch as watch
from sc1_two_field_rates_instrument import solver_source, RECIPE


def decimal_solve(matrix,rhs):
    """Independent high-precision Gaussian elimination, not the candidate formula."""
    with localcontext() as ctx:
        ctx.prec=100
        a=[[Decimal(str(v)) for v in row]+[Decimal(str(b))] for row,b in zip(matrix,rhs)]
        for i in range(len(a)):
            for j in range(i+1,len(a)):
                multiplier=a[j][i]/a[i][i]
                a[j]=[x-multiplier*y for x,y in zip(a[j],a[i])]
        x=[Decimal(0)]*len(a)
        for i in reversed(range(len(a))):
            x[i]=(a[i][-1]-sum(a[i][j]*x[j] for j in range(i+1,len(a))))/a[i][i]
        return np.array([float(v) for v in x])


class Algebra(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which('g++'):raise unittest.SkipTest('g++ unavailable')
        cls.tmp=tempfile.TemporaryDirectory();p=Path(cls.tmp.name);cls.binary=p/'matrix'
        (p/'test.C').write_text('''#include <iostream>
#include <iomanip>
#include <limits>
#include "sc1FaceDragAlgebra.H"
int main() {
 double a,b,k,t00,t01,t10,t11;
 if (sc1FaceDragTransform(0.,1.,1.,t00,t01,t10,t11)) return 10;
 if (sc1FaceDragTransform(1.,-1.,1.,t00,t01,t10,t11)) return 11;
 if (sc1FaceDragTransform(1.,1.,-1.,t00,t01,t10,t11)) return 12;
 if (sc1FaceDragTransform(1.,1.,std::numeric_limits<double>::infinity(),t00,t01,t10,t11)) return 13;
 if (sc1FaceDragTransform(std::numeric_limits<double>::quiet_NaN(),1.,1.,t00,t01,t10,t11)) return 14;
 while (std::cin>>a>>b>>k) {
   if (!sc1FaceDragTransform(a,b,k,t00,t01,t10,t11)) return 15;
   std::cout<<std::setprecision(17)<<t00<<" "<<t01<<" "<<t10<<" "<<t11<<"\\n";
 }
}
''')
        subprocess.run(['g++','-std=c++14','-Wall','-Wextra','-Werror','-I',str(candidate.HEADERS[2].parent),
            str(p/'test.C'),'-o',str(cls.binary)],check=True,capture_output=True,text=True)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def transforms(self,rows):
        output=subprocess.run([str(self.binary)],input=''.join(f'{a:.17g} {b:.17g} {k:.17g}\n' for a,b,k in rows),
            capture_output=True,text=True,check=True).stdout
        return np.array([[float(v) for v in line.split()] for line in output.splitlines()]).reshape(-1,2,2)

    def test_compiled_inverse_against_independent_elimination(self):
        rows=[(a,b,k) for a,b in ((.02,2e7),(2e4,20),(1,1),(1e-10,1e10))
              for k in (0,.11282816668027722,1,1.4511846946572018e9,1e16,1e25)]
        for (a,b,k),t in zip(rows,self.transforms(rows)):
            with localcontext() as ctx:
                ctx.prec=100
                da,db,dk=map(lambda v:Decimal(str(v)),(a,b,k));d0=da+dk;d1=db+dk
                matrix=[[d0,-dk],[-dk,d1]]
                expected=np.column_stack((decimal_solve(matrix,[d0,0]),decimal_solve(matrix,[0,d1])))
            np.testing.assert_allclose(t,expected,rtol=1e-14,atol=0)
        np.testing.assert_array_equal(self.transforms([(2,3,0)])[0],np.eye(2))

    def test_coupled_pressure_and_phase_flux_match_full_three_equation_system(self):
        # Two momentum equations plus prescribed mixture flux. Covers absent,
        # near-bound and mixed phases; no lagged other-phase flux enters the solve.
        rows=[(max(alpha,1e-6)*2e4,max(1-alpha,1e-6)*2e7,k)
              for alpha in (0,1e-5,.5,1-1e-5,1) for k in (0,1e3,1.451e9,1e20)]
        for index,((a,b,k),t) in enumerate(zip(rows,self.transforms(rows))):
            alpha=(0,1e-5,.5,1-1e-5,1)[index//4]
            weights=np.array([alpha,1-alpha]);response=np.maximum(weights,1e-6)
            rhs=np.array([.1,-.3]);target=1e-8;diag=np.array([a+k,b+k])
            predictor=t@(rhs/diag);pressure=t@(response/diag)
            self.assertGreater(weights@pressure,0)
            q=(target-weights@predictor)/(weights@pressure)
            phi=predictor+pressure*q
            with localcontext() as ctx:
                ctx.prec=100
                da,db,dk=map(lambda v:Decimal(str(v)),(a,b,k))
                expected=decimal_solve([[da+dk,-dk,-response[0]],[-dk,db+dk,-response[1]],
                    [alpha,1-alpha,0]],[.1,-.3,target])
            np.testing.assert_allclose(np.r_[phi,q],expected,rtol=2e-10,atol=1e-17)
            self.assertAlmostEqual(weights@phi,target,places=15)


class Instrumentation(unittest.TestCase):
    def test_inverse_used_in_predictor_and_pressure_without_second_elimination(self):
        source=candidate.base.SOURCE/'multiphaseEulerFoam'
        old=(source/'pUf/pEqn.H').read_text();new=instrument.pressure(old)
        self.assertNotIn('fluid.phiKdPhifs(',new)
        self.assertNotIn('fluid.partialEliminationf(',new)
        self.assertIn('rAUf += alphafs[phasei]*sc1PressureAlpha[phasei];',new)
        self.assertIn('+ sc1PressureAlpha[phase.index()]*mSfGradp;',new)
        self.assertLess(new.index('sc1Coupling.apply(phiHbyAs'),new.index('constrainPhiHbyA'))
        self.assertLess(new.index('sc1Probe.pressure('),new.index('sc1CoupledFluxCheck'))
        self.assertIn('setSnGrad<fixedFluxPressureFvPatchScalarField>',new)
        self.assertNotIn('setMixturePhi',new)
        main=instrument.solver(solver_source((source/'multiphaseEulerFoam.C').read_text()))
        self.assertNotIn('if (sc1Probe.endStep()) return 0;',main)
        self.assertLess(main.index('sc1QuiescentStep(fluid);'),main.index('runTime.write();'))
        self.assertIn('Info<< "End',main)
        with self.assertRaises(ValueError):instrument.pressure(new)


class Build(unittest.TestCase):
    def test_application_only_build_is_pinned_and_preserves_old_sources(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);installed=p/'installed';installed.mkdir()
            parent=candidate.base.prepare_build(p/'rates',installed=installed,recipe=RECIPE)
            before=candidate.base.selected_hashes(parent/'src');out=p/'candidate'
            def execute(command,cwd,log,env):
                self.assertEqual(command,['wmake']);self.assertEqual(cwd,out/'src')
                self.assertEqual(env['FOAM_USER_LIBBIN'],str(parent/'lib'))
                log.write_text('mock');(out/'bin'/candidate.APP).write_text('mock binary')
            identity=dict(binary=str(out/'bin'/candidate.APP),sha256='mock',dependencies={})
            with (patch.object(candidate.startup,'require_allocation'),patch.dict(os.environ,WM_PROJECT_VERSION='9'),
                  patch.object(candidate.shutil,'which',return_value='/mock/wmake'),
                  patch.object(candidate.base,'verify_build',return_value={'assets':{}}),
                  patch.object(candidate.base,'execute',side_effect=execute) as build,
                  patch.object(candidate,'executable_identity',return_value=identity)):
                candidate.build(parent,out);build.assert_called_once()
                self.assertEqual(candidate.verify_build(out),identity)
                self.assertEqual(candidate.base.selected_hashes(parent/'src'),before)
                self.assertFalse((out/'lib').exists())
                self.assertIn('maxCompletedSteps=200',(out/'src/sc1PressureProbe.H').read_text())
                with self.assertRaises(ValueError):candidate.build(parent,out)
                (out/'src/sc1FaceDragAlgebra.H').write_text('changed')
                with self.assertRaisesRegex(ValueError,'coupling recipe'):candidate.verify_build(out)


class Cases(unittest.TestCase):
    def controls(self,p):
        old=p/'coupled';case=old/'quiescent'
        for folder in ('0','constant','system'):(case/folder).mkdir(parents=True)
        for name in ('0/alpha.gas','constant/g','system/controlDict'):(case/name).write_text('original')
        for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/name).write_text('saved')
        (case/'processor0').mkdir();(case/'log.solver').write_text('preserve failed log')
        (p/'candidate').mkdir();(p/'observed').mkdir()
        parent=dict(variant=candidate.startup.COUPLED_TAG,root=str(old),inputs=candidate.base.case_input_hashes(case),
            ranks=8,build=str(p/'rates'),assets={'oldbinary':'oldsha'},initial={},archive=str(p/'archive'),
            mesh_hashes={},regression_root=str(p/'regression'),regression_sha256='proof',
            baseline={'root':str(p/'v1')},failed_control={'root':str(p/'v2')},face_control={'root':str(p/'v3')})
        control=dict(root=str(old),log_sha256=candidate.base.sha(case/'log.solver'),manifest_sha256='manifest')
        (p/'candidate/coupling_build_manifest.json').write_text(json.dumps(dict(base_build=parent['build'],base_assets=parent['assets'])))
        identity=dict(binary=str(p/'candidate/bin'/candidate.APP),sha256='candidate',dependencies={})
        evidence=dict(root=str(p/'observed'),log_sha256=candidate.OBSERVED_LOG,manifest_sha256='proof')
        return old,parent,control,identity,evidence

    def test_fresh_same_inputs_with_received_probe_and_reuse_rejection(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,control,identity,evidence=self.controls(p)
            before=candidate.startup.tree_hashes(old)
            with (patch.object(candidate,'failed_control',return_value=(parent,control)),
                  patch.object(candidate,'observed_probe',return_value=evidence),
                  patch.object(candidate,'verify_build',return_value=identity),patch.object(candidate.startup,'require_allocation')):
                out=candidate.prepare(old,p/'candidate',p/'new',p/'observed');m,_=candidate.verify(out)
                self.assertEqual(m['inputs'],parent['inputs']);self.assertEqual(m['pressure_probe'],evidence)
                self.assertEqual(candidate.startup.tree_hashes(old),before)
                self.assertFalse((out/'quiescent/processor0').exists())
                (out/'quiescent/attempt.json').write_text('attempt')
                with self.assertRaisesRegex(ValueError,'attempted'):candidate.run(out)
                (out/'quiescent/system/controlDict').write_text('changed')
                with self.assertRaisesRegex(ValueError,'inputs changed'):candidate.verify(out)

    def test_numerical_stop_collects_evidence_without_reconstruction(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old,parent,control,identity,evidence=self.controls(p)
            with (patch.object(candidate,'failed_control',return_value=(parent,control)),
                  patch.object(candidate,'observed_probe',return_value=evidence),
                  patch.object(candidate,'verify_build',return_value=identity),patch.object(candidate.startup,'require_allocation')):
                out=candidate.prepare(old,p/'candidate',p/'new',p/'observed')
                def execute(command,cwd,log,env):
                    self.assertEqual(command[-3:],['-parallel','-fileHandler','uncollated'])
                    self.assertIn('--time=00:30:00',command)
                    log.write_text('SC1_PRESSURE_STOP reason=unclippedPressureBelow1000Pa\nFOAM FATAL ERROR\n')
                    raise ValueError('solver failure')
                with (patch.object(candidate.base,'execute') as utility,patch.object(watch,'execute',side_effect=execute),
                      patch.object(failure,'collect') as collect):
                    with self.assertRaisesRegex(ValueError,'solver failure'):candidate.run(out)
                    self.assertEqual(utility.call_count,1)
                    self.assertEqual(utility.call_args.args[0][0],'decomposePar')
                    collect.assert_called_once_with(out,out.with_name(out.name+'_evidence'))
                self.assertIn(candidate.TAG,failure.VARIANTS)

    def test_observed_failure_hash_required(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'quiescent').mkdir();(p/'quiescent/log.solver').write_text('changed')
            with patch.object(candidate.observation,'verify',return_value=({},{})):
                with self.assertRaisesRegex(ValueError,'received, unchanged'):candidate.observed_probe(p,{})


class Acceptance(unittest.TestCase):
    def log(self):
        lines=['SC1_PRESSURE_PROBE version=1 maxCompletedSteps=200',
            'SC1_PRESSURE_DONE completedSteps=200 meaning=requiresStartupAnalysis']
        lines+=['SC1_FACE_COUPLING version=1 order=beforePressure']*400
        lines+=['SC1_PRESSURE_FIELD term=pUnclipped min=99983 nonfinite=0']*800
        lines+=['SC1_COUPLED_FLUX maxAbsError_m3_s=1e-20 maxAbsMixture_m3_s=1e-8']*800
        lines += [f'SC1_QUIESCENT_STEP t_s={i*5e-5:.12g} maxLiquidPresent_m_s=1e-5' for i in range(1,201)]
        return '\n'.join(lines)+'\n'

    def test_requires_old_startup_acceptance_plus_complete_new_path(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);log=p/'log.solver';log.write_text(self.log())
            m=dict(pressure_probe={'log_sha256':candidate.OBSERVED_LOG})
            with patch.object(candidate.startup,'analyze',side_effect=lambda case,meta:dict(startup_pass=False)) as old:
                report=candidate.analyze(p,m)
                self.assertFalse(report['startup_pass'])
                self.assertEqual(report['pressureCoupling']['mixtureFluxChecks'],800)
                self.assertEqual(old.call_args.args[1]['variant'],candidate.startup.FACE_TAG)
                for bad in (self.log().replace('min=99983','min=-100',1),
                            self.log().replace('maxAbsError_m3_s=1e-20','maxAbsError_m3_s=1e-3',1),
                            self.log().replace('maxLiquidPresent_m_s=1e-5','maxLiquidPresent_m_s=0.002',1),
                            self.log()+'Inverting drag system: Min face det = 1\n',
                            self.log().replace('SC1_FACE_COUPLING version=1 order=beforePressure\n','',1)):
                    log.write_text(bad)
                    with self.assertRaises(ValueError):candidate.analyze(p,m)


if __name__=='__main__':unittest.main()
