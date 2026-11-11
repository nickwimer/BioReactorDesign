"""Face-only numerical change, preserved controls, complete reports and live stops."""
from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shutil
import sys
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


class FaceRecipe(unittest.TestCase):
    def test_only_momentum_option_and_write_frequency_change(self):
        old=recipe.case_files(8,segregated_drag=True)
        new=recipe.case_files(8,segregated_drag=True,face_momentum=True)
        self.assertEqual({k for k in old if old[k]!=new[k]},{'system/fvSolution','system/controlDict'})
        self.assertEqual(new['system/fvSolution'].replace(' faceMomentum true;',''),old['system/fvSolution'])
        self.assertEqual(new['system/controlDict'].replace('writeInterval 5;','writeInterval 40;'),old['system/controlDict'])
        self.assertNotIn('partialElimination',new['system/fvSolution'])
        with self.assertRaises(ValueError):recipe.case_files(8,face_momentum=True)
        source=run.base.SOURCE/'multiphaseEulerFoam'
        self.assertIn('#include "pEqnComps.H"',(source/'pUf/pEqn.H').read_text())
        self.assertIn('Constructing face momentum equations',(source/'pUf/UEqns.H').read_text())

    def test_prepare_verify_preserves_controls_and_rejects_extra_numerical_change(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);old=p/'old';old.mkdir();m=drag_tests.SavedControl().fixture(old);case=old/'quiescent'
            c=np.array([[0,y,0] for y in (.999,1.01,1.022,1.033)]);v=np.full(4,.5)
            fields,m['initial']=recipe.initialized_fields(c,v)
            for rel,text in fields.items():(case/rel).write_text(text)
            (case/'0/C').write_text(array_field('C',c));(case/'0/V').write_text(array_field('V',v))
            (case/'log.checkMesh').write_text('Mesh OK.')
            m['inputs']=run.base.case_input_hashes(case)
            (p/'archive/constant/polyMesh').mkdir(parents=True)
            report={'log_sha256':run.base.sha(case/'log.solver')}
            with (patch.object(run,'saved_startup',return_value=(m,report)),patch.object(run,'require_allocation'),
                  patch.object(run,'CELLS',4),patch.object(run,'regression',return_value=({},'proof')),
                  patch.object(run.base,'verify_build',return_value={'assets':{}})):
                failed=run.prepare_drag(old,p/'drag');fm=run.verify(failed);fc=failed/'quiescent'
                (fc/'attempt.json').write_text(json.dumps({k:fm[k] for k in ('build','assets','ranks')}))
                (fc/'log.solver').write_text('saved failed solver log\n')
                (fc/'processor0/0').mkdir(parents=True);(fc/'0.00025').mkdir()
                before=run.tree_hashes(failed);old_before=run.tree_hashes(old)
                with patch.object(run,'FAILED_DRAG_LOG',run.base.sha(fc/'log.solver')):
                    out=run.prepare_face(failed,p/'face');new=run.verify(out)
                    self.assertEqual(new['variant'],run.FACE_TAG)
                    self.assertEqual({k for k in fm['inputs'] if fm['inputs'][k]!=new['inputs'][k]},
                                     {'system/fvSolution','system/controlDict'})
                    self.assertFalse((out/'quiescent/processor0').exists())
                    self.assertFalse((out/'quiescent/0.00025').exists())
                    self.assertFalse((out/'quiescent/log.solver').exists())
                    self.assertEqual(run.tree_hashes(failed),before)
                    self.assertEqual(run.tree_hashes(old),old_before)
                    with self.assertRaises(ValueError):run.prepare_face(failed,out)
                    sol=out/'quiescent/system/fvSolution'
                    sol.write_text(sol.read_text().replace('faceMomentum true;','faceMomentum true; partialElimination true;'))
                    new['inputs']['system/fvSolution']=run.base.sha(sol)
                    (out/'startup_manifest.json').write_text(json.dumps(new))
                    with self.assertRaisesRegex(ValueError,'recipe'):run.verify(out)
                with self.assertRaisesRegex(ValueError,'received evidence'):run.failed_drag(failed)
                with self.assertRaisesRegex(ValueError,'failed v2'):run.failed_drag(out)

    def test_face_report_requires_runtime_path_and_all_early_snapshots(self):
        with tempfile.TemporaryDirectory() as d,patch.object(run,'CELLS',4):
            p=Path(d);m=startup_tests.Reports().setup_case(p);m['variant']=run.FACE_TAG
            log=p/'log.solver'
            text='Selecting drag blending method: none\n'+('Constructing face momentum equations\n'*400)+log.read_text()
            log.write_text(text)
            with self.assertRaises(FileNotFoundError):run.analyze(p,m)
            for s in range(recipe.FACE_WRITE,recipe.STEPS+1,recipe.FACE_WRITE):
                dest=p/f'{s*recipe.DT:.12g}'
                if not dest.exists():shutil.copytree(p/'0.002',dest)
            r=run.analyze(p,m)
            self.assertTrue(r['startup_pass']);self.assertEqual(len(r['snapshots']),40)
            self.assertEqual(r['snapshots'][0]['time_s'],.00025)
            log.write_text(text.replace('Constructing face momentum equations\n','',1))
            with self.assertRaisesRegex(ValueError,'runtime coverage'):run.analyze(p,m)
            log.write_text(text+'limitPressure: p min -73364\n')
            with self.assertRaisesRegex(ValueError,'pressure limiter'):run.analyze(p,m)

    def test_runner_uses_monitor_and_reconstructs_early_fields(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d);case=root/'quiescent';case.mkdir()
            for n in ('0','constant','system'):(case/n).mkdir()
            for n in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/n).write_text('saved')
            m=dict(variant=run.FACE_TAG,ranks=8,build='/private',assets={},failed_control={},inputs=run.base.case_input_hashes(case))
            commands=[]
            def execute(cmd,cwd,log,env):commands.append(cmd);log.write_text('mock')
            with (patch.object(run,'verify',return_value=m),patch.object(run,'require_allocation'),
                  patch.object(run.base,'environment',return_value={}),patch.object(run.base,'execute',side_effect=execute),
                  patch.object(watch,'execute',side_effect=execute) as monitored,
                  patch.object(run,'analyze',return_value={'startup_pass':False}),patch.object(flow,'flow',return_value={})):
                r=run.run(root);self.assertFalse(r['startup_pass'])
                self.assertEqual(monitored.call_count,1)
                self.assertIn('--kill-on-bad-exit=1',commands[1]);self.assertIn('--time=00:30:00',commands[1])
                self.assertEqual(commands[2],['reconstructPar','-time','0.00025:0.01','-fileHandler','uncollated'])
                self.assertTrue((root/'face_comparison.json').is_file())
                with self.assertRaisesRegex(ValueError,'attempted'):run.run(root)

    def test_guard_failure_collects_evidence_and_never_reconstructs_or_analyzes(self):
        import sc1_reactor_failure
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d)/'face';case=root/'quiescent';case.mkdir(parents=True)
            for n in ('0','constant','system'):(case/n).mkdir()
            for n in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/n).write_text('saved')
            m=dict(variant=run.FACE_TAG,ranks=8,build='/private',assets={},inputs=run.base.case_input_hashes(case))
            def stop(cmd,cwd,log,env):
                log.write_text('Time = 0.0003\nlimitPressure: p min -1\n')
                raise ValueError('pressure limiter')
            with (patch.object(run,'verify',return_value=m),patch.object(run,'require_allocation'),
                  patch.object(run.base,'environment',return_value={}),patch.object(run.base,'execute') as execute,
                  patch.object(watch,'execute',side_effect=stop),patch.object(run,'analyze') as analyze,
                  patch.object(sc1_reactor_failure,'collect') as collect):
                with self.assertRaisesRegex(ValueError,'pressure limiter'):run.run(root)
                self.assertEqual(execute.call_count,1)
                self.assertEqual(execute.call_args.args[0][0],'decomposePar')
                analyze.assert_not_called();collect.assert_called_once_with(root,root.with_name('face_failure'))
                self.assertTrue((case/'log.solver').is_file())


class MonitorTests(unittest.TestCase):
    def test_boundaries_banner_and_fatal_conditions(self):
        normal=['sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).',
                'Selecting finite volume constraint type limitPressure','Courant Number mean: 0 max: 0.2',
                'liquid fraction, min, max = .6 -1e-10 1.0000000001']
        m=watch.Monitor()
        for line in normal:self.assertIsNone(m.feed(line))
        for line in ('limitPressure: p min -1','Courant Number mean: 0 max: .21',
                     'gas fraction, min, max = .1 -1e-6 .2','liquid fraction, min, max = nan -inf inf',
                     'SC1_TWO_STEP mass_kg=nan','srun: error: Floating point exception',
                     'FOAM FATAL IO ERROR: missing dictionary'):
            m=watch.Monitor();m.feed('Time = 0.0003')
            self.assertIsNotNone(m.feed(line),line);self.assertEqual(m.reason['lastPrintedTime'],'0.0003')

    def test_live_abort_keeps_log_and_stops_before_hanging_child(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            p=Path(d);code="import time; print('Time = 0.0003',flush=True); print('limitPressure: p min -1',flush=True); time.sleep(30)"
            with self.assertRaisesRegex(ValueError,'pressure limiter'):
                watch.execute([sys.executable,'-c',code],p,p/'log.solver',dict(os.environ),poll_seconds=.01)
            r=json.loads((p/'guard_stop.json').read_text())
            self.assertEqual(r['lastPrintedTime'],'0.0003');self.assertEqual(r['line'],2)
            self.assertNotEqual(r['returncode'],0)
            self.assertIn('limitPressure:',(p/'log.solver').read_text())

    def test_fast_exit_fragment_nonzero_and_timeout(self):
        with tempfile.TemporaryDirectory() as d,redirect_stdout(io.StringIO()):
            root=Path(d)
            for name,code,expect in (
                ('ok',"print('End')",None),
                ('fragment',"import sys; sys.stdout.write('Courant Number mean: 0 max: 1')",'Courant'),
                ('badexit',"raise SystemExit(7)",'Command failed'),
                ('timeout',"import time; time.sleep(30)",'wall-clock')):
                p=root/name;p.mkdir()
                if expect:
                    with self.assertRaisesRegex(ValueError,expect):
                        watch.execute([sys.executable,'-c',code],p,p/'log.solver',dict(os.environ),wall_seconds=.2,poll_seconds=.01)
                else:watch.execute([sys.executable,'-c',code],p,p/'log.solver',dict(os.environ),poll_seconds=.01)
            with self.assertRaises(FileExistsError):
                watch.execute([sys.executable,'-c',"print('overwritten')"],root/'ok',root/'ok/log.solver',dict(os.environ))
            self.assertEqual((root/'ok/log.solver').read_text(),'End\n')


if __name__=='__main__':unittest.main()
