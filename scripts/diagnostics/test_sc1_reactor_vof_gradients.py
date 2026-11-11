"""Scalar-gradient-only preparation, immutable controls and strict run/report gates."""
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import sc1_reactor_vof_gradients as candidate
import sc1_reactor_vof_case as cases
import sc1_reactor_watch as watch


class Gradients(unittest.TestCase):
    def fixture(self,path):
        root=path/'control';case=root/'quiescent'
        files=dict(cases.case_files(),**{'0/C':'centres','0/V':'volumes','0/alpha.liquid':'initial alpha',
            '0/p_rgh':'initial pressure','0/p':'static pressure','constant/polyMesh/points':'mesh',
            'constant/g':'gravity','system/decomposeParDict':'eight ranks'})
        for rel,text in files.items():
            p=case/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
        for name in ('log.checkMesh','log.writeCellCentres','log.writeCellVolumes'):(case/name).write_text('control log')
        m=dict(root=str(root),variant='sc1_meml_vof_reference_v1',ranks=8,assets={'solver':'hash'},
            build=str(path/'build'),archive=str(path/'archive'),regression_root=str(path/'rates'),
            baseline={'root':str(path/'baseline')},control={'root':str(path/'original_control')},
            initial={'liquid_volume_m3':1,'gas_volume_m3':.642},mesh_hashes={'points':'mesh'},
            reference_identity={'binary':str(path/'build/bin/sc1InterFoamReference')},
            recipe=candidate.reference.recipe(),build_manifest_sha256='build manifest',
            inputs=candidate.base.case_input_hashes(case),source_commit='control commit')
        (root/'startup_manifest.json').write_text(json.dumps(m));(case/'log.solver').write_text('Time = 0.0084\n')
        (case/'attempt.json').write_text(json.dumps({k:m[k] for k in ('build','assets','ranks')}))
        (case/'processor0').mkdir();(case/'processor0/saved').write_text('native saved field')
        report=dict(suite='sc1_vof_saved_spatial_v1',sourceRoot=str(root),
            sourceHashes={rel:candidate.base.sha(root/rel) for rel in ('startup_manifest.json','quiescent/log.solver','quiescent/attempt.json')},
            fieldHashes=dict(m['inputs'],**{'processor0/saved':candidate.base.sha(case/'processor0/saved')}),
            snapshots=self.reports()[2]['snapshots'])
        spatial=path/'spatial';spatial.mkdir();(spatial/'report.json').write_text(json.dumps(report))
        return root,spatial,m

    def patches(self,root,spatial,m):
        from contextlib import ExitStack
        stack=ExitStack()
        stack.enter_context(patch.multiple(candidate.saved,MANIFEST_SHA=candidate.base.sha(root/'startup_manifest.json'),
            LOG_SHA=candidate.base.sha(root/'quiescent/log.solver')))
        stack.enter_context(patch.object(candidate,'SPATIAL_SHA',candidate.base.sha(spatial/'report.json')))
        stack.enter_context(patch.object(candidate.reference,'verify',return_value=m))
        stack.enter_context(patch.object(candidate.startup,'require_allocation'))
        stack.enter_context(redirect_stdout(io.StringIO()))
        return stack

    def test_exact_four_scalar_overrides_preserve_all_other_schemes(self):
        original=cases.case_files()['system/fvSchemes'];text=candidate.schemes(original)
        expected='''gradSchemes {
    default Gauss linear;
    grad(p_rgh) leastSquares;
    grad(rho) leastSquares;
    grad(alpha.liquid) leastSquares;
    nHat leastSquares;
}'''
        self.assertEqual(text.replace(expected,'gradSchemes { default Gauss linear; }'),original)
        self.assertNotIn('grad(U) leastSquares',text)
        self.assertIn('nHat leastSquares;',text)  # OF9 curvature uses this named lookup.
        for changed in ('gradSchemes {}',original+original):
            with self.assertRaises(ValueError):candidate.schemes(changed)

    def test_fresh_t0_copy_changes_only_fvSchemes_and_keeps_control_untouched(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,spatial,parent=self.fixture(p)
            before=candidate.startup.tree_hashes(root)
            with self.patches(root,spatial,parent):
                out=candidate.prepare(root,spatial,p/'candidate');m=candidate.verify(out)
                self.assertEqual(before,candidate.startup.tree_hashes(root))
                differences=[k for k in parent['inputs'] if parent['inputs'][k]!=m['inputs'][k]]
                self.assertEqual(differences,['system/fvSchemes'])
                self.assertEqual(m['reference_identity'],parent['reference_identity'])
                self.assertFalse((out/'quiescent/processor0').exists())
                self.assertFalse((out/'quiescent/log.solver').exists())
                with self.assertRaises(ValueError):candidate.prepare(root,spatial,out)
                with self.assertRaises(ValueError):candidate.prepare(root,spatial,root/'nested')
                (out/'quiescent/0/p_rgh').write_text('altered')
                with self.assertRaisesRegex(ValueError,'identity changed'):candidate.verify(out)

    def test_control_report_native_fields_and_candidate_metadata_are_revalidated(self):
        for change in ('spatial','native','metadata'):
            with self.subTest(change=change),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,spatial,parent=self.fixture(p)
                with self.patches(root,spatial,parent):
                    out=candidate.prepare(root,spatial,p/'candidate')
                    if change=='spatial':(spatial/'report.json').write_text('{}')
                    elif change=='native':(root/'quiescent/processor0/saved').write_text('altered')
                    else:
                        f=out/'startup_manifest.json';m=json.loads(f.read_text());m['ranks']=16;f.write_text(json.dumps(m))
                    with self.assertRaises(ValueError):candidate.verify(out)

    def test_failed_solver_collected_without_reconstructing_or_changing_control(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);root,spatial,parent=self.fixture(p)
            with self.patches(root,spatial,parent):
                out=candidate.prepare(root,spatial,p/'candidate');before=candidate.startup.tree_hashes(root)
                def fail(command,case,log,env):
                    self.assertIn(parent['reference_identity']['binary'],command)
                    self.assertEqual(command[-3:],['-parallel','-fileHandler','uncollated'])
                    log.write_text('Time = 0.00005\nFOAM FATAL ERROR:\n')
                    raise ValueError('controlled solver stop')
                with patch.object(watch,'execute',side_effect=fail),patch.object(candidate.base,'execute') as execute:
                    with self.assertRaisesRegex(ValueError,'controlled solver stop'):candidate.run(out)
                self.assertEqual(execute.call_count,1)
                self.assertEqual(execute.call_args.args[0][0],'decomposePar')
                self.assertTrue((p/'candidate_evidence/reactor_failure_evidence.tar.gz').is_file())
                self.assertEqual(before,candidate.startup.tree_hashes(root))
                with self.assertRaisesRegex(ValueError,'attempted'):candidate.run(out)

    def reports(self):
        times=(.00025,.001,.004,.00825)
        def snap(t,speed):
            return dict(time_s=t,phases={p:{'present':dict(max_m_s=speed,rms_m_s=speed/2)} for p in ('gas','liquid')})
        r=dict(reference_pass=True,maxGas_m_s=.0001,snapshots=[snap(t,.0001) for t in times])
        steps=[dict(t_s=t,gasMax_m_s=.0001) for t in times]
        control=dict(snapshots=[snap(t,.9) for t in times])
        return r,steps,control

    def test_gas_rest_is_required_even_when_original_liquid_screen_passes(self):
        r,steps,control=self.reports();r['maxGas_m_s']=.002
        result=candidate.finish_report(r,steps,control)
        self.assertTrue(result['reference_pass']);self.assertFalse(result['gradient_pass'])
        self.assertFalse(result['gas_rest_pass']);self.assertEqual(result['gasRestLimit_m_s'],.001)
        self.assertEqual(len(result['comparison']),4)
        r,steps,control=self.reports();r['maxGas_m_s']=.001
        self.assertTrue(candidate.finish_report(r,steps,control)['gradient_pass'])
        r,steps,control=self.reports();r['reference_pass']=False
        self.assertFalse(candidate.finish_report(r,steps,control)['gradient_pass'])

    def test_native_gas_disagreement_is_rejected_before_acceptance(self):
        r,steps,control=self.reports();steps[-1]['gasMax_m_s']=.0005
        with self.assertRaisesRegex(ValueError,'Native gas speed'):candidate.finish_report(r,steps,control)

    def test_normal_run_requires_all_outputs_then_reports_comparison_and_gas_rest(self):
        from test_sc1_reactor_vof import good_log
        for missing in (False,True):
            with self.subTest(missing=missing),tempfile.TemporaryDirectory() as d:
                p=Path(d);root,spatial,parent=self.fixture(p)
                with self.patches(root,spatial,parent):
                    out=candidate.prepare(root,spatial,p/'candidate')
                    def complete(command,case,log,env):
                        log.write_text(good_log())
                        for rank in range(8):
                            for i in range(5,201,5):
                                for name in ('U','alpha.liquid','p_rgh'):
                                    if missing and rank==7 and i==200 and name=='U':continue
                                    f=case/f'processor{rank}'/f'{i*5e-5:.12g}'/name
                                    f.parent.mkdir(parents=True,exist_ok=True);f.write_text('saved field')
                    report,_,_=self.reports()
                    with (patch.object(watch,'execute',side_effect=complete),
                          patch.object(candidate.base,'execute') as execute,
                          patch.object(candidate.reporting,'analyze',return_value=report) as analyze):
                        if missing:
                            with self.assertRaisesRegex(ValueError,'Missing evidence'):candidate.run(out)
                            self.assertEqual(execute.call_count,1);analyze.assert_not_called()
                        else:
                            result=candidate.run(out)
                            self.assertTrue(result['gradient_pass']);self.assertEqual(len(result['comparison']),4)
                            self.assertEqual([c.args[0][0] for c in execute.call_args_list],['decomposePar','reconstructPar'])
                            self.assertTrue((out/'quiescent/result.json').is_file())
                        self.assertTrue((p/'candidate_evidence/reactor_failure_evidence.tar.gz').is_file())


if __name__=='__main__':unittest.main()
