#!/usr/bin/env python3
"""One-variable numerical controls, immutable state, and launch guards."""
import contextlib
import io
import json
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import prepare_sc1_conservation_pair as pair
import sc1_alpha_convergence as convergence
import test_sc1_alpha_normal as fixtures


class NumericalControlsTests(unittest.TestCase):
    def setUp(self):
        fixture=fixtures.PreparedNormalTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.baseline_pair=fixture.root
        self.baseline=fixture.control
        self.source=fixture.source
        self.out=self.baseline_pair.parent/'numerical_checks'
        control=self.baseline/'system/controlDict'
        control.write_text(control.read_text()+'\ndeltaT 5e-05;\nmaxDeltaT 5e-05;\nadjustTimeStep yes;\n')
        solution=self.baseline/'system/fvSolution'
        text=solution.read_text().replace('PIMPLE {','PIMPLE {\n nOuterCorrectors 2;')
        solution.write_text(text)
        (self.baseline/'scaler_test.txt').write_text('unchanged scaler')

    def prepare(self):
        return convergence.prepare(self.baseline_pair,self.out)

    def test_only_reviewed_numerical_entries_change_and_all_state_is_preserved(self):
        before=pair.hashes(self.baseline_pair)
        original=pair.hashes(self.source)
        manifest=self.prepare()
        bcontrol=(self.baseline/'system/controlDict').read_text()
        bsolution=(self.baseline/'system/fvSolution').read_text()
        for name in convergence.VARIANTS:
            case=self.out/name
            self.assertEqual(pair.hashes(case/'constant'),pair.hashes(self.baseline/'constant'))
            for section,digests in manifest['initial_state_sha256'].items():
                self.assertEqual(pair.hashes(case/section),digests)
            self.assertEqual((case/'scaler_test.txt').read_bytes(),b'unchanged scaler')
            self.assertEqual(convergence.preflight(self.out,name),manifest)
        four=(self.out/'four_outer/system/fvSolution').read_text()
        start,end=convergence.pimple_bounds(four)
        self.assertEqual(convergence.value(four[start+1:end],'nOuterCorrectors'),'4')
        self.assertEqual((self.out/'four_outer/system/controlDict').read_text(),bcontrol)
        half=(self.out/'half_dt/system/controlDict').read_text()
        self.assertEqual(convergence.value(half,'deltaT'),'2.5e-05')
        self.assertEqual(convergence.value(half,'maxDeltaT'),'2.5e-05')
        self.assertEqual((self.out/'half_dt/system/fvSolution').read_text(),bsolution)
        self.assertEqual(pair.hashes(self.baseline_pair),before)
        self.assertEqual(pair.hashes(self.source),original)

    def test_original_initial_state_is_used_even_when_completed_control_co_was_rewritten(self):
        path=self.baseline/f'processors52/{pair.CHECKPOINT}/Co.liquid'
        path.write_bytes(path.read_bytes().replace(str(self.source).encode(),str(self.baseline).encode()))
        manifest=self.prepare()
        for name in convergence.VARIANTS:
            self.assertEqual(pair.hashes(self.out/name/f'processors52/{pair.CHECKPOINT}'),
                             manifest['initial_state_sha256'][f'processors52/{pair.CHECKPOINT}'])

    def test_changed_original_checkpoint_blocks_creation(self):
        path=self.source/f'processors52/{pair.CHECKPOINT}/O2.liquid'
        path.write_bytes(path.read_bytes()+b'changed')
        with self.assertRaisesRegex(ValueError,'Original common checkpoint/mesh changed'):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_active_entrainment_or_changed_baseline_controls_block_creation(self):
        path=self.baseline/'constant/phaseProperties'
        original=path.read_text()
        path.write_text(original.replace('entrainmentEfficiency 0;','entrainmentEfficiency 0.01;'))
        with self.assertRaisesRegex(ValueError,'zero-entrainment'):
            self.prepare()
        path.write_text(original)
        control=self.baseline/'system/controlDict'
        control.write_text(control.read_text().replace('maxDeltaT 5e-05;','maxDeltaT 1e-4;'))
        with self.assertRaisesRegex(ValueError,'timestep baseline'):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_existing_and_overlapping_destinations_are_never_reused(self):
        self.prepare()
        before=pair.hashes(self.out)
        with self.assertRaisesRegex(ValueError,'Destination exists'):
            self.prepare()
        for root in (self.source,self.baseline):
            with self.assertRaisesRegex(ValueError,'overlaps'):
                convergence.prepare(self.baseline_pair,root/'must_not_create')
            self.assertFalse((root/'must_not_create').exists())
        self.assertEqual(pair.hashes(self.out),before)

    def test_changed_case_input_or_model_asset_is_rejected(self):
        self.prepare()
        case=self.out/'four_outer'
        for rel in ('constant/microbeUpdateProperties','system/fvSolution','scaler_test.txt'):
            path=case/rel
            before=path.read_bytes()
            path.write_bytes(before+b'changed')
            with self.assertRaises(ValueError):
                convergence.preflight(self.out,'four_outer')
            path.write_bytes(before)

    def test_attempted_case_or_extra_time_is_not_restarted(self):
        self.prepare()
        case=self.out/'four_outer'
        time=case/'processors52/31.3'
        time.mkdir()
        with self.assertRaisesRegex(ValueError,'Unexpected case times'):
            convergence.preflight(self.out,'four_outer')
        time.rmdir()
        (case/'log.failed').write_text('preserve')
        with self.assertRaisesRegex(ValueError,'already attempted'):
            convergence.preflight(self.out,'four_outer')

    def test_redirected_checkpoint_is_rejected(self):
        self.prepare()
        field=self.out/'four_outer'/f'processors52/{pair.CHECKPOINT}/O2.liquid'
        field.unlink()
        field.symlink_to(self.source/f'processors52/{pair.CHECKPOINT}/O2.liquid')
        with self.assertRaisesRegex(ValueError,'redirected input'):
            convergence.preflight(self.out,'four_outer')

    def test_run_is_opt_in_and_checks_fail_before_launch(self):
        self.prepare()
        with (patch.object(convergence,'scratch',lambda p:Path(p).resolve()),
              patch.object(convergence,'check_build'),
              patch.object(convergence,'run_case',side_effect=AssertionError('Must not launch')),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(convergence.main(['run','--root',str(self.out),'--case','four_outer']),0)
            (self.out/'four_outer/system/controlDict').write_text('changed')
            self.assertEqual(convergence.main(['run','--root',str(self.out),'--case','four_outer','--run']),1)
        self.assertFalse((self.out/'four_outer/log.conservation').exists())

    def test_launch_is_one_case_with_private_process_environment_and_exclusive_log(self):
        manifest=self.prepare()
        baseline_before=pair.hashes(self.baseline_pair)
        half_before=pair.hashes(self.out/'half_dt')
        runner=Mock(return_value=SimpleNamespace(returncode=17))
        with patch.dict(os.environ,{'SLURM_JOB_ID':'test','LD_LIBRARY_PATH':'/original'}):
            with contextlib.redirect_stdout(io.StringIO()):
                status=convergence.run_case(self.out,'four_outer',self.out/'build',manifest,runner)
            self.assertEqual(status,17)
            self.assertEqual(os.environ['LD_LIBRARY_PATH'],'/original')
            self.assertEqual(runner.call_args.kwargs['cwd'],self.out/'four_outer')
            self.assertEqual(runner.call_args.args[0][:3],['srun','-n','52'])
            with self.assertRaises(FileExistsError),contextlib.redirect_stdout(io.StringIO()):
                convergence.run_case(self.out,'four_outer',self.out/'build',manifest,runner)
        self.assertEqual(runner.call_count,1)
        self.assertEqual(pair.hashes(self.baseline_pair),baseline_before)
        self.assertEqual(pair.hashes(self.out/'half_dt'),half_before)

    def test_no_allocation_leaves_case_unrun(self):
        manifest=self.prepare()
        with patch.dict(os.environ,{},clear=True),self.assertRaisesRegex(ValueError,'No active allocation'):
            convergence.run_case(self.out,'four_outer',self.out/'build',manifest,Mock())
        self.assertFalse((self.out/'four_outer/log.conservation').exists())

    def test_comments_and_nested_outer_count_do_not_override_pimple_setting(self):
        control=(self.baseline/'system/controlDict').read_text()
        solution=(self.baseline/'system/fvSolution').read_text()
        with_comments='// PIMPLE { nOuterCorrectors 200; }\n'+solution
        result=convergence.variants(control,with_comments)['four_outer'][1]
        self.assertIn('// PIMPLE { nOuterCorrectors 200; }',result)
        nested=solution.replace('nOuterCorrectors 2;','nested {\n nOuterCorrectors 2;\n }')
        with self.assertRaisesRegex(ValueError,'one explicit top-level nOuterCorrectors'):
            convergence.variants(control,nested)


if __name__=='__main__':
    unittest.main(verbosity=2)
