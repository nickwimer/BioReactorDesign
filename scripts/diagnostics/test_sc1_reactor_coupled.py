"""Test the coupled-drag branch selection, actual scalar solve loops and isolation."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_startup as run
import sc1_reactor_startup_case as recipe
import sc1_reactor_startup_flow as flow
import sc1_reactor_watch as watch
import test_sc1_reactor_drag as drag_tests
import test_sc1_reactor_startup as startup_tests
from test_sc1_reactor_startup import array_field


class CoupledAlgebra(unittest.TestCase):
    def test_actual_two_phase_elimination_and_mixture_restoration_loops(self):
        compiler=shutil.which('g++')
        if not compiler:self.skipTest('C++ compiler unavailable')
        source=(run.base.SOURCE/'phaseSystems/PhaseSystems/MomentumTransferPhaseSystem/MomentumTransferPhaseSystem.C').read_text()
        source=source.split('::partialEliminationf\n',1)[1]
        decompose=source.split('        // Decompose\n',1)[1].split('\n        {\n            surfaceScalarField detPhiKdfs',1)[0]
        solve=source.split('        // Solve for the fluxes\n',1)[1].split('\n    }\n\n    this->setMixturePhi',1)[0]
        restore=(run.base.SOURCE/'phaseSystems/phaseSystem/phaseSystem.C').read_text().split('void Foam::phaseSystem::setMixturePhi\n',1)[1]
        restore=restore.split('\n{\n',1)[1].split('\n}\n',1)[0]
        self.assertIn('-rAUfs[phase1i]*K',source)
        code=r'''
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>
using label=int;
using surfaceScalarField=double;
struct Phase {
    int id; double flux;
    int index()const{return id;}
    bool stationary()const{return false;}
    double phi()const{return flux;}
    double& phiRef(){return flux;}
};
#define forAll(v,i) for(label i=0;i<static_cast<label>(v.size());++i)
void solve(std::array<std::array<double,2>,2>& phiKdfs,std::array<Phase,2>& phases) {
'''+decompose+solve+r'''
}
void restore(std::array<Phase,2>& movingPhaseModels_,const std::array<double,2>& alphafs,double phim0) {
'''+restore+r'''
}
int main() {
    for(double g:{0.,.01,.5,.9,.999,.999999})
    for(double l:{0.,.02,.6,.95,.9999,.999999})
    for(double b:{-3.,0.,.125,2.}) {
        std::array<std::array<double,2>,2> matrix{{{{1.,-g}},{{-l,1.}}}};
        std::array<Phase,2> phases{{{0,b},{1,.7}}};
        long double det=1.L-static_cast<long double>(g)*l;
        long double e0=(static_cast<long double>(b)+static_cast<long double>(g)*.7)/det;
        long double e1=(.7L+static_cast<long double>(l)*b)/det;
        solve(matrix,phases);
        assert(matrix[0][0]*matrix[1][1]>0);
        assert(std::abs(phases[0].flux-e0)<1e-9L*std::max(1.L,std::abs(e0)));
        assert(std::abs(phases[1].flux-e1)<1e-9L*std::max(1.L,std::abs(e1)));
    }
    for(double a:{0.,1e-5,.3,.99,1.}) {
        std::array<Phase,2> phases{{{0,20.},{1,-1.}}};
        restore(phases,{a,1-a},.7);
        assert(std::abs(a*phases[0].flux+(1-a)*phases[1].flux-.7)<1e-12);
        assert(std::abs(phases[0].flux-phases[1].flux-21.)<1e-12);
    }
}
'''
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'test.C').write_text(code)
            subprocess.run([compiler,'-std=c++17',str(p/'test.C'),'-o',str(p/'test')],check=True,capture_output=True,text=True)
            subprocess.run([str(p/'test')],check=True,capture_output=True,text=True)

    def test_one_dictionary_change_with_legacy_recipes_retained(self):
        old=recipe.case_files(8,segregated_drag=True,face_momentum=True)
        new=recipe.case_files(8,segregated_drag=True,face_momentum=True,coupled_drag=True)
        self.assertEqual([k for k in old if old[k]!=new[k]],['system/fvSolution'])
        self.assertEqual(new['system/fvSolution'].replace(' partialElimination true;',''),old['system/fvSolution'])
        with self.assertRaises(ValueError):recipe.case_files(8,segregated_drag=True,coupled_drag=True)


class CoupledWorkflow(unittest.TestCase):
    def test_prepare_and_verify_entire_control_chain_and_detect_tampering(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old=p/'old';old.mkdir();m=drag_tests.SavedControl().fixture(old);case=old/'quiescent'
            c=np.array([[0,y,0] for y in (.999,1.01,1.022,1.033)]);v=np.full(4,.5)
            fields,m['initial']=recipe.initialized_fields(c,v)
            for rel,text in fields.items():(case/rel).write_text(text)
            (case/'0/C').write_text(array_field('C',c));(case/'0/V').write_text(array_field('V',v))
            (case/'log.checkMesh').write_text('Mesh OK.');m['inputs']=run.base.case_input_hashes(case)
            (p/'archive/constant/polyMesh').mkdir(parents=True)
            report={'log_sha256':run.base.sha(case/'log.solver')}
            def attempt(root,text):
                meta=run.verify(root);cp=root/'quiescent'
                (cp/'attempt.json').write_text(json.dumps({k:meta[k] for k in ('build','assets','ranks')}))
                (cp/'log.solver').write_text(text)
                return run.base.sha(cp/'log.solver')
            with (patch.object(run,'saved_startup',return_value=(m,report)),patch.object(run,'require_allocation'),
                  patch.object(run,'CELLS',4),patch.object(run,'regression',return_value=({},'proof')),
                  patch.object(run.base,'verify_build',return_value={'assets':{}})):
                drag=run.prepare_drag(old,p/'drag');dh=attempt(drag,'failed drag\n')
                with patch.object(run,'FAILED_DRAG_LOG',dh):
                    face=run.prepare_face(drag,p/'face');fh=attempt(face,'failed face\n')
                    (face/'quiescent/processor0/0').mkdir(parents=True)
                    before={str(root):run.tree_hashes(root) for root in (old,drag,face)}
                    with patch.object(run,'FAILED_FACE_LOG',fh):
                        coupled=run.prepare_face(face,p/'coupled',coupled=True);cm=run.verify(coupled)
                        fm=json.loads((face/'startup_manifest.json').read_text())
                        self.assertEqual(cm['variant'],run.COUPLED_TAG)
                        self.assertEqual([k for k in fm['inputs'] if fm['inputs'][k]!=cm['inputs'][k]],['system/fvSolution'])
                        self.assertFalse((coupled/'quiescent/processor0').exists())
                        self.assertFalse((coupled/'quiescent/log.solver').exists())
                        for root in (old,drag,face):self.assertEqual(run.tree_hashes(root),before[str(root)])
                        with self.assertRaises(ValueError):run.prepare_face(face,coupled,coupled=True)
                        field=coupled/'quiescent/system/controlDict'
                        field.write_text(field.read_text().replace('deltaT 5e-05;','deltaT 1e-05;'))
                        cm['inputs']['system/controlDict']=run.base.sha(field)
                        (coupled/'startup_manifest.json').write_text(json.dumps(cm))
                        with self.assertRaisesRegex(ValueError,'recipe'):run.verify(coupled)
                    with self.assertRaisesRegex(ValueError,'received evidence'):run.failed_face(face)
                    with self.assertRaisesRegex(ValueError,'failed v3'):run.failed_face(drag)

    def test_report_requires_complete_positive_coupled_solve_evidence(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=startup_tests.Reports().setup_case(p);m['variant']=run.COUPLED_TAG
            for s in range(recipe.FACE_WRITE,recipe.STEPS+1,recipe.FACE_WRITE):
                dest=p/f'{s*recipe.DT:.12g}'
                if not dest.exists():shutil.copytree(p/'0.002',dest)
            log=p/'log.solver'
            base='Selecting drag blending method: none\n'+('Constructing face momentum equations\n'*400)+log.read_text()
            line='Inverting drag system: Min face det = 0.25\n'
            log.write_text(line*800+base);r=run.analyze(p,m)
            self.assertTrue(r['startup_pass']);self.assertEqual(len(r['snapshots']),40)
            self.assertEqual(r['dragCoupling']['calls'],800)
            self.assertEqual(r['dragCoupling']['minimumFaceDeterminant'],.25)
            log.write_text(line*799+base)
            with self.assertRaisesRegex(ValueError,'determinant coverage'):run.analyze(p,m)
            log.write_text(line.replace('0.25','0')+line*799+base)
            with self.assertRaisesRegex(ValueError,'nonpositive'):run.analyze(p,m)

    def test_runner_keeps_monitor_and_records_immediate_face_control(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d);case=root/'quiescent';case.mkdir()
            for n in ('0','constant','system'):(case/n).mkdir()
            for n in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/n).write_text('saved')
            m=dict(variant=run.COUPLED_TAG,ranks=8,build='/private',assets={},face_control={'root':'/face'},
                inputs=run.base.case_input_hashes(case))
            commands=[]
            def execute(cmd,cwd,log,env):commands.append(cmd);log.write_text('mock')
            with (patch.object(run,'verify',return_value=m),patch.object(run,'require_allocation'),
                  patch.object(run.base,'environment',return_value={}),patch.object(run.base,'execute',side_effect=execute),
                  patch.object(watch,'execute',side_effect=execute) as monitored,
                  patch.object(run,'analyze',return_value={'startup_pass':False}),patch.object(flow,'flow',return_value={})):
                run.run(root)
            monitored.assert_called_once()
            self.assertIn('--kill-on-bad-exit=1',commands[1]);self.assertIn('0.00025:0.01',commands[2])
            comparison=json.loads((root/'coupled_comparison.json').read_text())
            self.assertEqual(comparison['failed_control'],m['face_control'])
            self.assertEqual(comparison['numerical_change'],'PIMPLE partialElimination true')

    def test_monitor_rejects_invalid_determinants_without_a_small_positive_cutoff(self):
        self.assertIsNone(watch.Monitor().feed('Inverting drag system: Min face det = 1e-15'))
        for value in ('0','-1','nan','inf','bad'):
            self.assertIn('face drag determinant',watch.Monitor().feed('Inverting drag system: Min face det = '+value)['reason'])
        self.assertIn('malformed',watch.Monitor().feed('Inverting drag system: missing determinant')['reason'])


if __name__=='__main__':unittest.main()
