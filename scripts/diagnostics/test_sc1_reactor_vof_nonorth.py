"""Pressure-iteration coverage, unchanged budgets, and isolated fresh-case runs."""
from contextlib import ExitStack, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof_nonorth as candidate
import sc1_reactor_vof_nonorth_report as reporting
import sc1_reactor_watch as watch
from test_sc1_reactor_vof import good_log, INITIAL
import test_sc1_reactor_vof_gradients as fixtures


def corrected_log():
    lines=[]
    for line in good_log().splitlines():
        if line.startswith('SC1_VOF_PRESSURE '):
            for i in range(10):
                lines.append('DICPCG:  Solving for p_rgh, Initial residual = '+('1e-4' if i==0 else '5e-13')+
                    ', Final residual = 5e-13, No Iterations '+('100' if i==0 else '0'))
                lines.append(line)
        else:lines.append(line)
    return '\n'.join(lines)+'\n'


class Pressure(unittest.TestCase):
    def test_every_pressure_pass_validated_and_original_mass_budget_preserved(self):
        text=corrected_log().replace('min_Pa=99983','min_Pa=90000',1)
        r,_=reporting.budget(text,INITIAL)
        self.assertEqual(r['pressureIteration']['pressureRecords'],8000)
        self.assertEqual(r['pressureIteration']['convergedCorrections'],800)
        self.assertEqual(r['minimumPressure_Pa'],90000)
        self.assertFalse(r['errors'])
        self.assertAlmostEqual(r['phases']['liquid']['boundaryOut_kg'],.0001)

    def test_incomplete_misordered_nonfinite_or_unconverged_records_rejected(self):
        text=corrected_log();first=next(x for x in text.splitlines() if x.startswith('SC1_VOF_PRESSURE '))
        variants=[text.replace(first+'\n','',1),text.replace(first,first+'\n'+first,1),
            text.replace('min_Pa=99983','min_Pa=nan',1),text.replace('min_Pa=99983','min_Pa=999',1),
            text.replace('nonfinite=0','nonfinite=1',1),text.replace('Initial residual = 1e-4','Initial residual = nan',1),
            text.replace('Final residual = 5e-13','Final residual = 2e-12',1),
            text.replace(first,first+' t_s=5e-5',1),text.replace('End\n',''),
            text.replace('Time = 5e-05','Time = 0.0001',1)]
        for i,t in enumerate(variants):
            with self.subTest(i=i),self.assertRaises(ValueError):reporting.budget(t,INITIAL)
        with self.assertRaises(ValueError):reporting.budget(good_log(),INITIAL)

    def test_final_linear_solve_is_not_enough_for_explicit_correction_convergence(self):
        text=corrected_log().replace('Initial residual = 5e-13','Initial residual = 2e-12')
        r,_=reporting.budget(text,INITIAL)
        self.assertFalse(r['pressureIteration']['pressure_iteration_pass'])
        self.assertEqual(r['pressureIteration']['convergedCorrections'],0)
        text=corrected_log().replace('No Iterations 0','No Iterations 1')
        r,_=reporting.budget(text,INITIAL)
        self.assertFalse(r['pressureIteration']['pressure_iteration_pass'])

    def test_failed_prefix_can_report_convergence_without_accepting_completion(self):
        text=corrected_log().split('Time = 0.0001\n')[0]+'FOAM FATAL ERROR:\n'
        _,r=reporting.pressure_coverage(text,require_complete=False)
        self.assertEqual(r['corrections'],4);self.assertEqual(r['convergedCorrections'],4)
        self.assertFalse(r['complete']);self.assertFalse(r['pressure_iteration_pass'])
        with self.assertRaises(ValueError):reporting.budget(text,INITIAL)

    def test_native_snapshot_checks_both_phase_speeds(self):
        with tempfile.TemporaryDirectory() as d:
            case=Path(d);(case/'log.solver').write_text(corrected_log())
            def read(p,n,vector=False):
                t=float(p.parent.name)
                return {'C':np.array([[0,0,0],[0,1.5,0]]),'V':np.array([1.,.642]),
                    'alpha.liquid':np.array([1-t*1e-5,0.]),'U':np.tile([.0001,0,0],(2,1)),
                    'p_rgh':np.ones(2)*110000}[p.name]
            masks={'walls':np.array([True,True]),'propeller':np.array([False,False])}
            with patch.object(reporting,'internal',side_effect=read),patch.object(reporting,'patch_cells',return_value=masks):
                r=reporting.analyze(case,INITIAL,2)
                self.assertTrue(r['reference_pass']);self.assertEqual(len(r['snapshots']),40)
                (case/'log.solver').write_text(corrected_log().replace('gasMax_m_s=.0001','gasMax_m_s=.0002'))
                with self.assertRaisesRegex(ValueError,'phase speed'):reporting.analyze(case,INITIAL,2)


class Cases(unittest.TestCase):
    def fixture(self,path):
        root,spatial,m=fixtures.Gradients().fixture(path)
        case=root/'quiescent'
        scheme=case/'system/fvSchemes';scheme.write_text(candidate.gradients.schemes(scheme.read_text()))
        m.update(inputs=candidate.base.case_input_hashes(case),variant=candidate.gradients.TAG,
            gradient_control=dict(root=str(path/'original_vof'),spatial_root=str(spatial)),gradient_recipe=candidate.gradients.recipe())
        (root/'startup_manifest.json').write_text(json.dumps(m))
        log=good_log().split('Time = 0.00845\n')[0].replace('liquidMax_m_s=.0001','liquidMax_m_s=.001003',168)
        # Only the final of the 168 records is the rejected liquid-speed state.
        log=log.replace('liquidMax_m_s=.001003','liquidMax_m_s=.0001',167)
        (case/'log.solver').write_text(log+'FOAM FATAL ERROR:\n')
        return root,spatial,m

    def patches(self,root,spatial,m):
        stack=ExitStack()
        stack.enter_context(patch.object(candidate.gradients,'verify',return_value=m))
        stack.enter_context(patch.object(candidate,'CONTROL_LOG',candidate.base.sha(root/'quiescent/log.solver')))
        stack.enter_context(patch.object(candidate.gradients,'SPATIAL_SHA',candidate.base.sha(spatial/'report.json')))
        stack.enter_context(patch.object(candidate.startup,'require_allocation'))
        stack.enter_context(redirect_stdout(io.StringIO()))
        return stack

    def test_only_fvSolution_changes_and_original_case_is_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,spatial,m=self.fixture(p)
            with self.patches(root,spatial,m):
                before=candidate.startup.tree_hashes(root)
                out=candidate.prepare(root,p/'candidate');n=candidate.verify(out)
                self.assertEqual(candidate.startup.tree_hashes(root),before)
                self.assertEqual([k for k in m['inputs'] if n['inputs'][k]!=m['inputs'][k]],['system/fvSolution'])
                old=(root/'quiescent/system/fvSolution').read_text();new=(out/'quiescent/system/fvSolution').read_text()
                self.assertEqual(new.replace('nNonOrthogonalCorrectors 9;','nNonOrthogonalCorrectors 0;'),old)
                self.assertFalse((out/'quiescent/processor0').exists())
                self.assertEqual(n['reference_identity'],m['reference_identity'])
                with self.assertRaises(ValueError):candidate.prepare(root,out)
                with self.assertRaises(ValueError):candidate.prepare(root,root/'nested')
                (out/'quiescent/system/fvSchemes').write_text('changed')
                with self.assertRaises(ValueError):candidate.verify(out)
        with self.assertRaises(ValueError):candidate.solution('nNonOrthogonalCorrectors 1;')

    def test_control_log_and_attempt_metadata_are_required(self):
        for name in ('log.solver','attempt.json'):
            with tempfile.TemporaryDirectory() as d:
                root,spatial,m=self.fixture(Path(d))
                with self.patches(root,spatial,m):
                    (root/'quiescent'/name).write_text('{}')
                    with self.assertRaises(ValueError):candidate.control(root)

    def test_pressure_convergence_gas_and_original_screens_all_required(self):
        for key in ('reference_pass','gas','pressure'):
            r,steps,control=fixtures.Gradients().reports()
            r['pressureIteration']={'pressure_iteration_pass':True}
            if key=='reference_pass':r[key]=False
            elif key=='gas':r['maxGas_m_s']=.002
            else:r['pressureIteration']['pressure_iteration_pass']=False
            self.assertFalse(candidate.finish_report(r,steps,control)['nonorth_pass'])

    def test_run_and_failure_preserve_control_and_collect_evidence(self):
        for mode in ('complete','failed','missing'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,spatial,m=self.fixture(p)
                with self.patches(root,spatial,m):
                    out=candidate.prepare(root,p/'candidate');before=candidate.startup.tree_hashes(root)
                    def execute(command,case,log,env):
                        self.assertIn(m['reference_identity']['binary'],command)
                        log.write_text(corrected_log())
                        if mode=='failed':raise ValueError('controlled stop')
                        for rank in range(8):
                            for i in range(5,201,5):
                                for name in ('U','alpha.liquid','p_rgh'):
                                    if mode=='missing' and rank==7 and i==200 and name=='U':continue
                                    f=case/f'processor{rank}'/f'{i*5e-5:.12g}'/name
                                    f.parent.mkdir(parents=True,exist_ok=True);f.write_text('saved')
                    r,_,_=fixtures.Gradients().reports();r['pressureIteration']={'pressure_iteration_pass':True}
                    with patch.object(watch,'execute',side_effect=execute),patch.object(candidate.base,'execute') as shell,patch.object(reporting,'analyze',return_value=r) as analyze:
                        if mode=='complete':
                            self.assertTrue(candidate.run(out)['nonorth_pass'])
                            self.assertEqual(shell.call_count,2)
                        else:
                            with self.assertRaises(ValueError):candidate.run(out)
                            self.assertEqual(shell.call_count,1);analyze.assert_not_called()
                    self.assertEqual(candidate.startup.tree_hashes(root),before)
                    self.assertTrue((p/'candidate_evidence/reactor_failure_evidence.tar.gz').is_file())
                    with self.assertRaisesRegex(ValueError,'attempted'):candidate.run(out)


if __name__=='__main__':unittest.main()
