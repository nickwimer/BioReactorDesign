"""Exercise incomplete logs and evidence preservation without OpenFOAM or CFD."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import sc1_reactor_failure as failure

LOG='''Selecting drag blending method: none
Courant Number mean: 0 max: 0
Time = 5e-05
PIMPLE: Iteration 1
liquid fraction, min, max = 0.6 0 1
DICPCG: Solving for p_rgh, Initial residual = 1, Final residual = 1e-13, No Iterations 23
SC1_TWO_STEP phase=gas t_s=0.00005 alphaMin=0 alphaMax=1
SC1_TWO_STEP phase=liquid t_s=0.0000500000000000001 alphaMin=0 alphaMax=1
Courant Number mean: 0.01 max: 0.3
Time = 0.0001
PIMPLE: Iteration 1
liquid fraction, min, max = 0.6 -0.0002 1.02
limitPressure: p min 20
PIMPLE: Iteration 2
SC1_PHASE_TERM phase=liquid t_s=0.0001 call=2 term=clamp volume_m3=3e-5 frozenRhoMass_kg=0.03
SC1_TWO_STEP phase=gas t_s=0.0001 alphaMin=0 alphaMax=1
SC1_TWO_STEP phase=liquid t_s=0.0001 alphaMin=0 alphaMax=1
Time = 0.0009
PIMPLE: Iteration 2
liquid fraction, min, max = -4e120 -1.4e142 9.9e141
[2] #1 Foam::sigFpe::sigHandler(int) at ??:?
srun: error: task 2: Floating point exception (core dumped)
*** STEP CANCELLED DUE TO TIME LIMIT ***
'''


def fixture(root):
    case=root/'quiescent';case.mkdir(parents=True)
    files={'0/alpha.liquid':'field\n','constant/phaseProperties':'drag config\n',
           'constant/polyMesh/owner':'mesh\n','system/fvSolution':'PIMPLE config\n'}
    inputs={}
    for rel,text in files.items():
        p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text);inputs[rel]=failure.sha(p)
    (case/'log.solver').write_text(LOG)
    manifest=dict(variant=failure.VARIANTS[1],root=str(root),build=str(root.parent/'build'),
        assets={'solver':'saved_hash'},ranks=8,inputs=inputs,source_commit='0876843')
    (root/'startup_manifest.json').write_text(json.dumps(manifest))
    (case/'attempt.json').write_text(json.dumps({k:manifest[k] for k in ('build','assets','ranks')}))
    return manifest


class FailureEvidence(unittest.TestCase):
    def test_trapping_banner_is_not_an_exception_and_does_not_hide_later_crash(self):
        banner='sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).\n'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'log';path.write_text(banner+'Time = 5e-05\nEnd\n')
            r=failure.scan(path)
            self.assertNotIn('floating_point_exception',r['firstEvents'])
            self.assertNotIn('floating_point_exception',r['eventCounts'])
            path.write_text(banner+LOG);r=failure.scan(path)
        event=r['firstEvents']['floating_point_exception']
        self.assertEqual(event['text'],'[2] #1 Foam::sigFpe::sigHandler(int) at ??:?')
        self.assertEqual(event['lastPrintedTime_s'],.0009)
        self.assertEqual(event['outer'],2)
        self.assertEqual(r['eventCounts']['floating_point_exception'],2)

    def test_earliest_events_keep_courant_context_and_unclamped_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'log';path.write_text(LOG);r=failure.scan(path)
        self.assertFalse(r['endMarkerSeen'])
        self.assertEqual(r['pairedStepRecordCount'],2)
        self.assertEqual(r['lastPairedStepRecordTime_s'],.0001)
        self.assertEqual(r['firstEvents']['courant_above_0_2']['lastPrintedTime_s'],.00005)
        self.assertEqual(r['firstEvents']['alpha_out_of_bounds']['lastPrintedTime_s'],.0001)
        self.assertEqual(r['firstEvents']['alpha_magnitude_above_10']['lastPrintedTime_s'],.0009)
        self.assertEqual(r['eventCounts']['alpha_out_of_bounds'],2)
        self.assertEqual(r['firstEvents']['clamp_above_1e_10_m3']['outer'],2)
        self.assertEqual(r['timeline'][0]['pressureIterationsMax'],23)
        self.assertEqual(r['timeline'][1]['pressureLimitCalls'],1)
        self.assertIn('allocation_time_limit',r['firstEvents'])
        self.assertIn('floating_point_exception',r['firstEvents'])
        self.assertNotIn('startup_pass',r)

    def test_nonfinite_and_duplicate_records_are_evidence_not_json_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'log';path.write_text(LOG+'''Time = nan
Courant Number mean: nan max: inf
liquid fraction, min, max = nan -inf inf
SC1_TWO_STEP phase=gas t_s=0.0001 alphaMin=0 alphaMax=1
SC1_TWO_STEP phase=gas t_s=nan mass_kg=nan
End
''');r=failure.scan(path)
        json.dumps(r,allow_nan=False)
        self.assertTrue(r['endMarkerSeen'])
        self.assertEqual(r['duplicateStepRecords'],1)
        self.assertEqual(r['malformedStepRecords'],1)
        for kind in ('invalid_time','nonfinite_courant','nonfinite_alpha','nonfinite_audit'):
            self.assertIn(kind,r['firstEvents'])
        self.assertNotIn('screen_pass',r)

    def test_bundle_is_complete_log_but_no_mesh_or_field_copy_and_case_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'case';fixture(root);out=Path(tmp)/'report'
            before={str(p):failure.sha(p) for p in root.rglob('*') if p.is_file()}
            with redirect_stdout(io.StringIO()):r=failure.collect(root,out)
            self.assertEqual(before,{str(p):failure.sha(p) for p in root.rglob('*') if p.is_file()})
            self.assertEqual(r['inputIdentity'],dict(checked=4,changed=[],missing=[],extra=[]))
            with tarfile.open(out/'reactor_failure_evidence.tar.gz') as tar:
                self.assertEqual(tar.extractfile('case_evidence/quiescent/log.solver').read().decode(),LOG)
                self.assertIn('case_evidence/quiescent/system/fvSolution',tar.getnames())
                self.assertNotIn('case_evidence/quiescent/constant/polyMesh/owner',tar.getnames())
                self.assertNotIn('case_evidence/quiescent/0/alpha.liquid',tar.getnames())
            self.assertIn('liquid fraction, min, max = 0.6 -0.0002 1.02',(out/'first_events.txt').read_text())
            with self.assertRaisesRegex(ValueError,'Destination exists'):failure.collect(root,out)

    def test_missing_changed_extra_inputs_are_reported_without_destroying_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'case';fixture(root);case=root/'quiescent'
            (case/'0/alpha.liquid').unlink()
            (case/'system/fvSolution').write_text('changed')
            (case/'system/extra').write_text('extra')
            with redirect_stdout(io.StringIO()):r=failure.collect(root,Path(tmp)/'report')
            self.assertEqual(r['inputIdentity']['changed'],['system/fvSolution'])
            self.assertEqual(r['inputIdentity']['missing'],['0/alpha.liquid'])
            self.assertEqual(r['inputIdentity']['extra'],['system/extra'])

    def test_reject_redirected_evidence_path_escape_overlap_and_wrong_attempt(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'case';m=fixture(root);case=root/'quiescent'
            with self.assertRaisesRegex(ValueError,'overlaps'):failure.collect(root,case/'report')
            m['baseline']={'root':str(Path(tmp)/'control')}
            (root/'startup_manifest.json').write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError,'overlaps'):
                failure.collect(root,Path(m['baseline']['root'])/'report')
            real=case/'log.solver';real.rename(case/'log.saved');real.symlink_to('log.saved')
            with self.assertRaisesRegex(ValueError,'Redirected'):failure.collect(root,Path(tmp)/'report')
            real.unlink();(case/'log.saved').rename(real)
            m['inputs']['system/../../../secret']='ignored'
            (root/'startup_manifest.json').write_text(json.dumps(m))
            with self.assertRaisesRegex(ValueError,'Unsafe'):failure.collect(root,Path(tmp)/'report')
            m['inputs'].pop('system/../../../secret');(root/'startup_manifest.json').write_text(json.dumps(m))
            (case/'attempt.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'identity mismatch'):failure.collect(root,Path(tmp)/'report')

    def test_log_mutation_during_archive_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'case';fixture(root);original=tarfile.TarFile.add
            def mutate(tar,path,*args,**kwargs):
                original(tar,path,*args,**kwargs)
                if Path(path).name=='log.solver':Path(path).write_text(LOG+'changed\n')
            with patch.object(tarfile.TarFile,'add',mutate),redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError,'changed during collection'):
                    failure.collect(root,Path(tmp)/'report')

    def test_empty_log_is_not_successful_evidence_collection(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'log';path.write_text('')
            with self.assertRaisesRegex(ValueError,'No recognizable'):failure.scan(path)


if __name__=='__main__':unittest.main()
