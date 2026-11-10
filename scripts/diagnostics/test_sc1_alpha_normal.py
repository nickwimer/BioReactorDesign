#!/usr/bin/env python3
"""Prepared-pair lifecycle and one-shot launcher regressions; no CFD launches."""
import contextlib
import io
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import prepare_sc1_conservation_pair as pair
import run_sc1_alpha_normal as normal_run
import test_sc1_checkpoint_metadata as fixtures
from test_sc1_conservation_reporting import audit_text


class PreparedNormalTests(unittest.TestCase):
    def setUp(self):
        fixture=fixtures.MatchedCheckpointCopyTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.source=fixture.control
        for case in (fixture.control,fixture.normal):
            (case/'system/fvSolution').write_text(
                'solvers { "alpha.*" { nAlphaCorr 1; nAlphaSubCycles 6; } }\nPIMPLE {}\n')
        self.root=fixture.root/'prepared_pair'
        argv=['prepare','--control',str(fixture.control),'--active',str(fixture.normal),
              '--allow-courant-location-metadata','--duration','0.1','--alpha-audit',
              '--out',str(self.root)]
        with (patch.object(sys,'argv',argv),patch.object(pair,'scratch',lambda p:Path(p).resolve()),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(pair.main(),0)
        self.control,self.normal=self.root/'control',self.root/'normal'
        self.section=f'processors52/{pair.CHECKPOINT}'
        self.write_control_log()

    def write_control_log(self):
        start=float(pair.CHECKPOINT)
        (self.control/'log.conservation').write_text(
            'SC1_ALPHA_LIBRARY version=1\n'+audit_text([start+0.1],start=start))

    def rewrite_location(self,name='Co.liquid'):
        path=self.control/self.section/name
        old=str((self.source/self.section/name).resolve()).encode()
        path.write_bytes(path.read_bytes().replace(old,str(path.resolve()).encode(),1))
        return path

    def test_pristine_pair_passes_without_mutation(self):
        before=pair.hashes(self.root)
        report=normal_run.preflight_pair(self.root)
        self.assertEqual(report['normal_initial_files'],60)  # 59 fields, one fixture mesh file
        self.assertEqual(report['control_location_proofs'],{})
        self.assertEqual(pair.hashes(self.root),before)

    def test_completed_control_location_changes_pass_full_payload_proof(self):
        for name in fixtures.CO_FIELDS:
            self.rewrite_location(name)
        # Reproduce the old preflight failure without changing the normal case.
        with self.assertRaisesRegex(ValueError,'Inputs are not identical'):
            pair.validate_pair(self.control,self.normal,pair.CHECKPOINT)
        before=pair.hashes(self.root)
        report=normal_run.preflight_pair(self.root)
        self.assertEqual(set(report['control_location_proofs']),set(fixtures.CO_FIELDS))
        for name,proof in report['control_location_proofs'].items():
            self.assertEqual(proof['normal_location'],str((self.source/self.section/name).resolve()))
            self.assertEqual(proof['control_location'],str((self.control/self.section/name).resolve()))
            self.assertNotEqual(proof['control_sha256'],proof['normal_sha256'])
        self.assertEqual(pair.hashes(self.root),before)

    def test_normal_changes_never_get_a_metadata_exception(self):
        for name in ('Co.liquid','O2.liquid'):
            with self.subTest(name=name):
                path=self.normal/self.section/name
                before=path.read_bytes()
                path.write_bytes(before+b'changed')
                with self.assertRaisesRegex(ValueError,'Normal initial state changed'):
                    normal_run.preflight_pair(self.root)
                path.write_bytes(before)

    def test_non_courant_and_mesh_control_changes_are_rejected(self):
        for relative in (self.section+'/O2.liquid','processors52/constant/mesh'):
            with self.subTest(relative=relative):
                path=self.control/relative
                before=path.read_bytes()
                path.write_bytes(before+b'changed')
                with self.assertRaisesRegex(ValueError,'non-Courant input changed'):
                    normal_run.preflight_pair(self.root)
                path.write_bytes(before)

    def test_courant_payload_change_is_not_mistaken_for_a_location_change(self):
        path=self.rewrite_location()
        path.write_bytes(path.read_bytes()+b'changed payload')
        with self.assertRaisesRegex(ValueError,'NOT proved location-only.*after the outer header'):
            normal_run.preflight_pair(self.root)

    def test_courant_other_header_change_is_rejected(self):
        path=self.rewrite_location()
        path.write_bytes(path.read_bytes().replace(b'version 2.0;',b'version 2.1;'))
        with self.assertRaisesRegex(ValueError,'NOT proved location-only.*beyond its location'):
            normal_run.preflight_pair(self.root)

    def test_completed_header_must_identify_its_actual_file(self):
        path=self.rewrite_location()
        path.write_bytes(path.read_bytes().replace(str(path).encode(),b'/unexpected/Co.liquid'))
        with self.assertRaisesRegex(ValueError,'NOT proved location-only.*source file'):
            normal_run.preflight_pair(self.root)

    def test_manifest_source_must_match_pristine_copied_header(self):
        self.rewrite_location()
        path=self.root/'comparison_manifest.json'
        data=json.loads(path.read_text())
        data['initial_state_source']=data['control_source']='/another/original/source'
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError,'NOT proved location-only.*source file'):
            normal_run.preflight_pair(self.root)

    def test_missing_or_extra_control_files_are_rejected(self):
        field=self.control/self.section/'state00'
        content=field.read_bytes()
        field.unlink()
        with self.assertRaisesRegex(ValueError,'file set changed'):
            normal_run.preflight_pair(self.root)
        field.write_bytes(content)
        (self.control/self.section/'extra').write_text('extra')
        with self.assertRaisesRegex(ValueError,'file set changed'):
            normal_run.preflight_pair(self.root)

    def test_extra_normal_time_or_log_prevents_restart(self):
        time=self.normal/'processors52/31.3'
        time.mkdir()
        with self.assertRaisesRegex(ValueError,'unexpected time'):
            normal_run.preflight_pair(self.root)
        time.rmdir()
        (self.normal/'log.failed').write_text('prior attempt')
        with self.assertRaisesRegex(ValueError,'Normal log already exists'):
            normal_run.preflight_pair(self.root)

    def test_settings_and_control_completion_remain_required(self):
        path=self.normal/'constant/microbeUpdateProperties'
        before=path.read_bytes()
        path.write_bytes(before+b'changed schedule')
        with self.assertRaisesRegex(ValueError,'Constants/biology differ'):
            normal_run.preflight_pair(self.root)
        path.write_bytes(before)
        log=self.control/'log.conservation'
        log.write_text(log.read_text().replace('\nEnd\n','\n'))
        with self.assertRaisesRegex(ValueError,'no solver End'):
            normal_run.preflight_pair(self.root)

    def test_cli_default_is_read_only_and_failure_cannot_launch(self):
        with (patch.object(normal_run,'scratch',lambda p:Path(p).resolve()),
              patch.object(normal_run,'check_build'),
              patch.object(normal_run,'launch',side_effect=AssertionError('Must not launch')),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(normal_run.main(['--pair-root',str(self.root)]),0)
            (self.normal/self.section/'O2.liquid').write_text('changed')
            self.assertEqual(normal_run.main(['--pair-root',str(self.root),'--run']),1)
        self.assertFalse((self.normal/'log.conservation').exists())

    def test_launch_scopes_library_path_and_never_overwrites_a_log(self):
        build=self.root/'build'
        control_before=pair.hashes(self.control)
        runner=Mock(return_value=SimpleNamespace(returncode=17))
        with patch.dict(os.environ,{'SLURM_JOB_ID':'test','LD_LIBRARY_PATH':'/original/lib'}):
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(normal_run.launch(self.root,build,runner),17)
            self.assertEqual(os.environ['LD_LIBRARY_PATH'],'/original/lib')
            args,kwargs=runner.call_args
            self.assertEqual(args[0][:3],['srun','-n','52'])
            self.assertEqual(kwargs['cwd'],self.normal)
            self.assertEqual(kwargs['env']['LD_LIBRARY_PATH'],str(build/'lib')+':/original/lib')
            with self.assertRaises(FileExistsError),contextlib.redirect_stdout(io.StringIO()):
                normal_run.launch(self.root,build,runner)
        self.assertEqual(runner.call_count,1)
        self.assertEqual(pair.hashes(self.control),control_before)

    def test_no_allocation_cannot_create_log(self):
        with patch.dict(os.environ,{},clear=True),self.assertRaisesRegex(ValueError,'No active'):
            normal_run.launch(self.root,self.root/'build',Mock())
        self.assertFalse((self.normal/'log.conservation').exists())


if __name__=='__main__':
    unittest.main(verbosity=2)
