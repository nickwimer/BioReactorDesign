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


def vof_step(**changes):
    values=dict(t_s=.00005,dt_s=.00005,liquid_m3=1,gas_m3=.64,
        liquidOut_m3_s=0,gasOut_m3_s=0,postPressureOut_m3_s=0,
        liquidResidual_kg=0,gasResidual_kg=0,liquidL1_kg=0,gasL1_kg=0,
        alphaMin=0,alphaMax=1,complement=0,liquidMax_m_s=.0001,gasMax_m_s=.1,
        alphaChangeL1_m3=0,Co=.001,alphaCalls=2)
    values.update(changes)
    return 'SC1_VOF_STEP '+' '.join(f'{k}={v}' for k,v in values.items())+'\n'


VOF_INITIAL='SC1_VOF_INITIAL liquid_m3=1 gas_m3=.64 version=1\n'
VOF_PRESSURE='SC1_VOF_PRESSURE t_s=.00005 min_Pa=99983 max_Pa=111800 nonfinite=0\n'


class FailureEvidence(unittest.TestCase):
    def scan_text(self,text):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'log';path.write_text('Time = 5e-05\n'+text)
            return failure.scan(path)

    def test_vof_records_are_distinct_from_bird_pairs_and_not_accepted_steps(self):
        r=self.scan_text(VOF_INITIAL+VOF_PRESSURE+vof_step()+vof_step(t_s=.0001,liquidMax_m_s=.00101))
        v=r['vof']
        self.assertEqual(r['pairedStepRecordCount'],0)
        self.assertEqual(v['recordCounts'],dict(INITIAL=1,PRESSURE=1,STEP=2))
        self.assertEqual(v['uniqueStepTimeCount'],2)
        self.assertEqual(v['maximumGasSpeed_m_s'],.1)
        self.assertEqual(v['minimumFinitePressure_Pa'],99983)
        self.assertEqual(v['lastValidStep']['record']['liquidL1_kg'],0)
        self.assertEqual(v['firstLimitEvents']['liquidMax_m_s']['t_s'],.0001)
        self.assertEqual(set(v['firstLimitEvents']),{'liquidMax_m_s'})
        self.assertEqual(r['largestFiniteAlphaMagnitude'],1)
        self.assertFalse(r['endMarkerSeen'])
        self.assertNotIn('startup_pass',v)
        self.assertNotIn('vof',self.scan_text(LOG))

    def test_vof_reports_each_existing_step_guard_and_initial_mass_limits(self):
        cases=[('alphaMin',-2e-10),('alphaMax',1+2e-10),('complement',2e-10),
               ('Co',.21),('liquidMax_m_s',.0011),('alphaChangeL1_m3',1.1e-5),
               ('alphaCalls',1),('liquidL1_kg',1.1e-5),('gasL1_kg',7e-9)]
        for key,value in cases:
            with self.subTest(key=key):
                v=self.scan_text(VOF_INITIAL+vof_step(**{key:value}))['vof']
                self.assertEqual(set(v['firstLimitEvents']),{key})
        at_limits=vof_step(alphaMin=-1e-10,alphaMax=1+1e-10,complement=1e-10,
            Co=.2,liquidMax_m_s=.001,alphaChangeL1_m3=1e-5,
            liquidL1_kg=1e-10+1e-8*1000,gasL1_kg=1e-10+1e-8*.64)
        self.assertEqual(self.scan_text(VOF_INITIAL+at_limits)['vof']['firstLimitEvents'],{})

    def test_vof_pressure_guard_does_not_require_a_completed_step(self):
        text=VOF_INITIAL+VOF_PRESSURE.replace('min_Pa=99983','min_Pa=-20').replace('nonfinite=0','nonfinite=3')
        r=self.scan_text(text+'[1] --> FOAM FATAL ERROR:\nInvalid VOF unclipped pressure\n')
        v=r['vof']
        self.assertEqual(v['validStepRecordCount'],0)
        self.assertIsNone(v['lastValidStep'])
        self.assertEqual(set(v['firstLimitEvents']),{'nonfinite_pressure','pressure_below_1000_Pa'})
        self.assertIn('foam_fatal',r['firstEvents'])

    def test_vof_malformed_nonfinite_and_duplicate_records_are_json_safe(self):
        text=(VOF_INITIAL+vof_step()+vof_step(t_s=.0000500000000000001)
              +vof_step(Co='nan')+vof_step().replace('alphaCalls=2','')
              +vof_step().replace('Co=0.001','Co=0.001 Co=0.002')
              +VOF_PRESSURE.replace('99983','inf')+vof_step(dt_s=0))
        v=self.scan_text(text)['vof'];json.dumps(v,allow_nan=False)
        self.assertEqual(v['duplicateStepRecords'],1)
        self.assertEqual(v['uniqueStepTimeCount'],1)
        self.assertEqual(v['validStepRecordCount'],2)
        self.assertEqual(v['malformedRecordCount'],5)
        self.assertEqual(v['firstMalformedRecord']['record']['Co'],'nan')
        self.assertIsNone(v['minimumFinitePressure_Pa'])

    def test_vof_ambiguous_or_missing_initial_record_does_not_invent_mass_limits(self):
        for initial in ('',VOF_INITIAL*2,VOF_INITIAL.replace('version=1','version=2')):
            v=self.scan_text(initial+vof_step(liquidL1_kg=1))['vof']
            self.assertIsNone(v['initial'])
            self.assertNotIn('liquidL1_kg',v['firstLimitEvents'])

    def test_vof_stock_alpha_output_and_bundle_keep_guard_context(self):
        stock='Phase-1 volume fraction = 0.61  Min(alpha.liquid) = -0.01  Max(alpha.liquid) = 1.02\n'
        r=self.scan_text(stock)
        self.assertEqual(r['largestFiniteAlphaMagnitude'],1.02)
        self.assertEqual(r['timeline'][0]['alphaMin'],-.01)
        self.assertIn('alpha_out_of_bounds',r['firstEvents'])
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'case';m=fixture(root);m['variant']='sc1_meml_vof_reference_v1'
            (root/'startup_manifest.json').write_text(json.dumps(m))
            text='Time = 5e-05\n'+VOF_INITIAL+vof_step(liquidMax_m_s=.002)
            (root/'quiescent/log.solver').write_text(text)
            before={str(p):failure.sha(p) for p in root.rglob('*') if p.is_file()}
            out=Path(tmp)/'report'
            with redirect_stdout(io.StringIO()):r=failure.collect(root,out)
            self.assertEqual(before,{str(p):failure.sha(p) for p in root.rglob('*') if p.is_file()})
            self.assertIn('SC1_VOF_STEP',(out/'first_events.txt').read_text())
            self.assertEqual(r['vof']['limitEventCounts'],{'liquidMax_m_s':1})

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
