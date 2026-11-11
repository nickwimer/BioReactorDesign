"""One-setting convergence extension; physical limits and control evidence stay fixed."""
from contextlib import ExitStack,redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import sc1_reactor_vof_hydro_convergence as candidate
import sc1_reactor_vof_hydro_convergence_report as report
import sc1_reactor_watch as watch
import test_sc1_reactor_vof_hydro as old
from test_sc1_reactor_vof import good_log,INITIAL


def corrected_log(passes=24,failed=False):
    lines=[old.HEADER.strip()];group=0
    for line in good_log().splitlines():
        if line.startswith('SC1_VOF_PRESSURE '):
            group+=1
            for i in range(passes):
                bad=failed and i==passes-1 and (group%4==1 or group in (2,6,10))
                initial='8e-10' if bad else ('1e-4' if i==0 else '5e-13')
                its=100 if i==0 else (2 if bad else 0)
                lines.append(f'DICPCG:  Solving for pPert, Initial residual = {initial}, Final residual = 5e-13, No Iterations {its}')
                lines.append(line)
        else:lines.append(line)
    return '\n'.join(lines)+'\n'


def snapshots():
    return [dict(time_s=i*5e-5,liquidVolume_L=1000.,phases={p:{'present':{'max_m_s':.0001}} for p in ('gas','liquid')}) for i in range(5,201,5)]


class Reports(unittest.TestCase):
    def test_requires_all_19200_passes_and_all_original_budgets(self):
        text=corrected_log();r,steps=report.budget(text,INITIAL)
        self.assertEqual(r['pressureIteration']['pressureRecords'],19200);self.assertEqual(r['pressureIteration']['convergedCorrections'],800)
        self.assertTrue(r['complete']);self.assertFalse(r['errors']);self.assertEqual(len(steps),200)
        for damaged in (text.replace('End\n',''),corrected_log(10),text.replace(old.HEADER,''),
            text.replace('rhoError_kg_m3=0','rhoError_kg_m3=nan'),text.replace('min_Pa=99983','min_Pa=999',1),
            text.replace('Final residual = 5e-13','Final residual = 2e-12',1)):
            with self.assertRaises(ValueError):report.budget(damaged,INITIAL)
        p=next(x for x in text.splitlines() if x.startswith('SC1_VOF_PRESSURE '))
        for damaged in (text.replace(p+'\n','',1),text.replace(p,p+'\n'+p,1)):
            with self.assertRaises(ValueError):report.budget(damaged,INITIAL)

    def test_final_solved_residual_does_not_replace_unchanged_initial_residual_gate(self):
        r,_=report.budget(corrected_log(failed=True),INITIAL)
        self.assertFalse(r['pressureIteration']['pressure_iteration_pass']);self.assertEqual(r['pressureIteration']['convergedCorrections'],597)
        text=corrected_log().replace('No Iterations 0','No Iterations 1')
        self.assertFalse(report.budget(text,INITIAL)[0]['pressureIteration']['pressure_iteration_pass'])

    def test_native_pressure_identity_and_both_phase_velocity_limits(self):
        # Reuse the independent hydrostatic snapshot test against the new parser.
        with patch.object(old,'report',report),patch.object(old,'log',corrected_log):
            old.Reports().test_snapshots_require_pressure_identity_and_both_phase_speeds()

    def test_matched_time_comparison_does_not_change_acceptance(self):
        a={'snapshots':snapshots(),'hydrostatic_pass':False};b={'snapshots':snapshots()}
        r=candidate.compare(a,b);self.assertFalse(r['hydrostatic_pass']);self.assertEqual(len(r['comparison']),40)
        self.assertEqual(r['comparison'][-1]['gasPeakDifference_m_s'],0)
        b['snapshots'][0]['time_s']+=.0001
        with self.assertRaises(ValueError):candidate.compare(a,b)


class Cases(unittest.TestCase):
    def fixture(self,p):
        source,m,fr,identity,br,binary=old.Cases().fixture(p)
        with old.Cases().patches(m,identity,binary):root=old.candidate.prepare(fr,br,p/'hydro')
        m=json.loads((root/'startup_manifest.json').read_text());case=root/'quiescent';text=corrected_log(10,True)
        (case/'log.solver').write_text(text);(case/'attempt.json').write_text(json.dumps({k:m[k] for k in ('build','assets','ranks')}))
        r,_=old.report.budget(text,m['initial']);r.update(hydrostatic_pass=False,snapshots=snapshots(),errors=['Nonorthogonal pressure corrections did not converge'])
        (case/'result.json').write_text(json.dumps(r))
        return root,m,r

    def patches(self,root,m):
        s=ExitStack();s.enter_context(patch.object(candidate.hydro,'verify',return_value=m))
        s.enter_context(patch.multiple(candidate,CONTROL_LOG=candidate.base.sha(root/'quiescent/log.solver'),
            CONTROL_MANIFEST=candidate.base.sha(root/'startup_manifest.json'),CONTROL_RESULT=candidate.base.sha(root/'quiescent/result.json')))
        s.enter_context(patch.object(candidate.startup,'require_allocation'));s.enter_context(redirect_stdout(io.StringIO()))
        return s

    def test_only_pressure_pass_count_changes_no_build_or_control_mutation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,m,_=self.fixture(p)
            with self.patches(root,m):
                before=candidate.startup.tree_hashes(root);out=candidate.prepare(root,p/'candidate');n=candidate.verify(out)
                self.assertEqual(n['hydro_identity'],m['hydro_identity']);self.assertEqual(before,candidate.startup.tree_hashes(root))
                self.assertEqual([k for k in m['inputs'] if m['inputs'][k]!=n['inputs'][k]],['system/fvSolution'])
                new=(out/'quiescent/system/fvSolution').read_text();old_text=(root/'quiescent/system/fvSolution').read_text()
                self.assertEqual(new.replace('nNonOrthogonalCorrectors 23;','nNonOrthogonalCorrectors 9;'),old_text)
                self.assertFalse((out/'quiescent/processor0').exists())
                for target in (out,root/'nested',Path(m['build'])/'nested'):
                    with self.assertRaises(ValueError):candidate.prepare(root,target)
                (out/'quiescent/0/pPert').write_text('tamper')
                with self.assertRaises(ValueError):candidate.verify(out)

    def test_changed_control_log_result_attempt_or_manifest_is_rejected(self):
        for name in ('startup_manifest.json','quiescent/log.solver','quiescent/result.json','quiescent/attempt.json'):
            with self.subTest(name=name),tempfile.TemporaryDirectory() as d:
                root,m,_=self.fixture(Path(d))
                with self.patches(root,m):
                    (root/name).write_text('{}')
                    with self.assertRaises(ValueError):candidate.control(root)

    def test_complete_missing_output_and_failed_runs_preserve_control_and_collect(self):
        for mode in ('complete','missing','failed'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,m,r=self.fixture(p)
                with self.patches(root,m):
                    out=candidate.prepare(root,p/'candidate');before=candidate.startup.tree_hashes(root)
                    def execute(command,case,logfile,env):
                        self.assertIn(m['hydro_identity']['binary'],command);logfile.write_text(corrected_log())
                        if mode=='failed':raise ValueError('controlled stop')
                        for rank in range(8):
                            for i in range(5,201,5):
                                for name in ('U','alpha.liquid','p_rgh','p','pPert'):
                                    if mode=='missing' and rank==7 and i==200 and name=='pPert':continue
                                    f=case/f'processor{rank}'/f'{i*5e-5:.12g}'/name;f.parent.mkdir(parents=True,exist_ok=True);f.write_text('saved')
                    final=dict(hydrostatic_pass=True,snapshots=snapshots())
                    with patch.object(watch,'execute',side_effect=execute),patch.object(candidate.base,'execute') as shell,patch.object(report,'analyze',return_value=final) as analyze:
                        if mode=='complete':
                            self.assertTrue(candidate.run(out)['hydrostatic_pass']);self.assertEqual(shell.call_count,2)
                        else:
                            with self.assertRaises(ValueError):candidate.run(out)
                            self.assertEqual(shell.call_count,1);analyze.assert_not_called()
                    self.assertTrue((p/'candidate_evidence/reactor_failure_evidence.tar.gz').is_file())
                    self.assertEqual(before,candidate.startup.tree_hashes(root))
                    with self.assertRaisesRegex(ValueError,'attempted'):candidate.run(out)


if __name__=='__main__':unittest.main()
