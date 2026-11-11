#!/usr/bin/env python3
"""Probe provenance, immutable checkpoint, one-shot launch and regression gates."""
import contextlib
import io
import os
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import prepare_sc1_conservation_pair as pair
import run_sc1_alpha_source_probe as probe
import test_sc1_alpha_convergence as fixtures


class SourceProbeRunTests(unittest.TestCase):
    def setUp(self):
        fixture = fixtures.NumericalControlsTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.baseline_pair, self.baseline = fixture.baseline_pair, fixture.baseline
        self.source = fixture.source
        self.out = self.baseline_pair.parent / 'probe_case'
        self.build = self.baseline_pair.parent / 'probe_build'
        self.build.mkdir()
        for name in ('controlDict', 'fvSolution'):
            (self.baseline_pair / 'normal/system' / name).write_bytes((self.baseline / 'system' / name).read_bytes())
        self.log = self.baseline / 'log.conservation'
        self.log.write_text(self.log.read_text().replace('steps=10000', 'steps=2000'))
        self.build_check = patch.object(probe, 'check_build').start()
        self.addCleanup(patch.stopall)

    def prepare(self):
        return probe.prepare(self.baseline_pair, self.out, self.build)

    def observed_log(self):
        endpoint = float(pair.CHECKPOINT) + 0.1
        lines = ['SC1_SOURCE_PROBE_LIBRARY version=1']
        for call in (1, 2):
            lines.append('SC1_SOURCE_PROBE_METHOD version=1 observerOnly=1')
            for phase in ('liquid', 'bubbleGas'):
                for bin_name in ('lower', 'upper', 'interior', 'outside'):
                    lines.append(f'SC1_SOURCE_PROBE phase={phase} t_s={endpoint:.17g} '
                                 f'callInStep={call} bin={bin_name} phaseCellSubcycleSamples=10 '
                                 'rawBelowSamples=0 rawAboveSamples=0 pressureBelowSamples=0 '
                                 'pressureAboveSamples=0 rawSource_kg=1 pressureSource_kg=1 '
                                 'split_kg=0 implicit_kg=0 rawCorrection_kg=0 rawExcess_kg=0 '
                                 'pressureExcess_kg=0 identityRemainder_kg=0 identityRemainderL1_kg=0 '
                                 'maxDtAbsSp=0')
        return self.log.read_text() + '\n'.join(lines) + '\n'

    def test_only_observer_switch_changes_and_original_inputs_survive(self):
        before, source_before = pair.hashes(self.baseline_pair), pair.hashes(self.source)
        manifest = self.prepare()
        self.assertEqual(probe.preflight(self.out), manifest)
        for section, expected in manifest['initial_state_sha256'].items():
            self.assertEqual(pair.hashes(self.out / section), expected)
        self.assertEqual(pair.hashes(self.out / 'constant'), pair.hashes(self.baseline / 'constant'))
        self.assertEqual((self.out / 'system/controlDict').read_bytes(), (self.baseline / 'system/controlDict').read_bytes())
        observed = (self.out / 'system/fvSolution').read_text()
        self.assertEqual(observed.replace('\n        sc1AlphaSourceProbe true;\n    ', ''),
                         (self.baseline / 'system/fvSolution').read_text())
        self.assertEqual(probe.assets(self.out), probe.assets(self.baseline))
        self.assertEqual(pair.hashes(self.baseline_pair), before)
        self.assertEqual(pair.hashes(self.source), source_before)

    def test_completed_control_co_is_not_used_as_pristine_state(self):
        path = self.baseline / f'processors52/{pair.CHECKPOINT}/Co.liquid'
        path.write_bytes(path.read_bytes().replace(str(self.source).encode(), str(self.baseline).encode()))
        manifest = self.prepare()
        self.assertEqual(pair.hashes(self.out / f'processors52/{pair.CHECKPOINT}'),
                         manifest['initial_state_sha256'][f'processors52/{pair.CHECKPOINT}'])

    def test_changed_original_state_prevents_creation(self):
        (self.source / f'processors52/{pair.CHECKPOINT}/O2.liquid').write_text('changed')
        with self.assertRaisesRegex(ValueError, 'Original common checkpoint/mesh changed'):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_changed_biology_and_numerics_prevent_creation(self):
        for name in ('constant/microbeUpdateProperties', 'system/fvSolution', 'system/controlDict'):
            path = self.baseline / name
            before = path.read_bytes()
            path.write_bytes(before + b'changed')
            with self.assertRaises(ValueError):
                self.prepare()
            path.write_bytes(before)
        self.assertFalse(self.out.exists())

    def test_build_failure_prevents_creation(self):
        self.build_check.side_effect = ValueError('binary changed')
        with self.assertRaisesRegex(ValueError, 'binary changed'):
            self.prepare()
        self.assertFalse(self.out.exists())

    def test_existing_symlink_and_overlapping_outputs_are_rejected(self):
        for parent in (self.baseline_pair, self.source, self.build):
            with self.assertRaisesRegex(ValueError, 'overlaps'):
                probe.prepare(self.baseline_pair, parent / 'probe', self.build)
        self.out.symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'Destination exists'):
            self.prepare()
        self.out.unlink()
        self.prepare()
        before = pair.hashes(self.out)
        with self.assertRaisesRegex(ValueError, 'Destination exists'):
            self.prepare()
        self.assertEqual(pair.hashes(self.out), before)

    def test_inputs_assets_and_baseline_log_remain_checked(self):
        self.prepare()
        for path in (self.out / 'constant/microbeUpdateProperties', self.out / 'system/fvSolution',
                     self.out / 'scaler_test.txt', self.log):
            before = path.read_bytes()
            path.write_bytes(before + b'changed')
            with self.assertRaises(ValueError):
                probe.preflight(self.out)
            path.write_bytes(before)

    def test_any_attempt_extra_time_or_redirect_blocks_launch(self):
        self.prepare()
        time = self.out / 'processors52/31.3'
        time.mkdir()
        with self.assertRaisesRegex(ValueError, 'Unexpected case times'):
            probe.preflight(self.out)
        time.rmdir()
        log = self.out / 'log.failed'
        log.write_text('preserve')
        with self.assertRaisesRegex(ValueError, 'already attempted'):
            probe.preflight(self.out)
        log.unlink()
        field = self.out / f'processors52/{pair.CHECKPOINT}/O2.liquid'
        field.unlink()
        field.symlink_to(self.source / f'processors52/{pair.CHECKPOINT}/O2.liquid')
        with self.assertRaisesRegex(ValueError, 'redirected input'):
            probe.preflight(self.out)

    def test_no_allocation_cannot_create_log(self):
        manifest = self.prepare()
        with patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(ValueError, 'No active allocation'):
            probe.launch(self.out, manifest, Mock())
        self.assertFalse((self.out / 'log.conservation').exists())

    def test_launch_is_single_case_scoped_environment_and_exclusive_log(self):
        manifest = self.prepare()
        runner = Mock(return_value=SimpleNamespace(returncode=17))
        before = pair.hashes(self.baseline_pair)
        with patch.dict(os.environ, {'SLURM_JOB_ID': 'test', 'LD_LIBRARY_PATH': '/original'}), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(probe.launch(self.out, manifest, runner), 17)
            self.assertEqual(os.environ['LD_LIBRARY_PATH'], '/original')
            self.assertEqual(runner.call_args.args[0][:3], ['srun', '-n', '52'])
            self.assertEqual(runner.call_args.kwargs['cwd'], self.out)
            self.assertEqual(runner.call_args.kwargs['env']['LD_LIBRARY_PATH'], str(self.build / 'lib') + ':/original')
            with self.assertRaises(FileExistsError):
                probe.launch(self.out, manifest, runner)
        self.assertEqual(runner.call_count, 1)
        self.assertEqual(pair.hashes(self.baseline_pair), before)

    def test_run_default_is_read_only_and_mutation_blocks_explicit_run(self):
        self.prepare()
        with (patch.object(probe, 'scratch', lambda p: Path(p).resolve()),
              patch.object(probe, 'launch', side_effect=AssertionError('Must not launch')),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(probe.main(['run', '--case', str(self.out)]), 0)
            (self.out / 'system/fvSolution').write_text('changed')
            self.assertEqual(probe.main(['run', '--case', str(self.out), '--run']), 1)
        self.assertFalse((self.out / 'log.conservation').exists())

    def test_regression_requires_complete_finite_unique_records_and_marker(self):
        manifest = self.prepare()
        text = self.observed_log()
        path = self.out / 'log.conservation'
        path.write_text(text)
        with contextlib.redirect_stdout(io.StringIO()):
            probe.regression(self.out, manifest)
        for changed in (text.replace('SC1_SOURCE_PROBE_LIBRARY version=1\n', ''),
                        text.replace('bin=lower ', 'bin=interior ', 1),
                        text.replace('rawSource_kg=1 ', 'rawSource_kg=nan ', 1),
                        text.replace('phaseCellSubcycleSamples=10 ', '', 1),
                        text.replace('SC1_SOURCE_PROBE phase=liquid', 'LOST phase=liquid', 1)):
            path.write_text(changed)
            with self.assertRaises((ValueError, KeyError)):
                probe.regression(self.out, manifest)

    def test_legacy_comparison_retains_multiplicity_and_detects_solution_changes(self):
        self.log.write_text(self.log.read_text() + 'SC1_ALPHA_TERM phase=liquid term=effectiveSource rhoBaseMass_kg=1\n')
        manifest = self.prepare()
        text = self.observed_log()
        for changed in (text + 'SC1_ALPHA_LIBRARY version=1\n', text.replace('rhoBaseMass_kg=1', 'rhoBaseMass_kg=2', 1)):
            (self.out / 'log.conservation').write_text(changed)
            with self.assertRaisesRegex(ValueError, 'Printed legacy regression differs'):
                probe.regression(self.out, manifest)

    def test_probe_flag_parser_ignores_comments_and_rejects_ambiguity(self):
        text = (self.baseline / 'system/fvSolution').read_text()
        self.assertIn('// sc1AlphaSourceProbe false;', probe.enable_probe('// sc1AlphaSourceProbe false;\n' + text))
        for changed in (probe.enable_probe(text), text.replace('nAlphaSubCycles 6;', 'nAlphaSubCycles 12;'),
                        text.replace('nAlphaCorr 1;', 'nested { nAlphaCorr 1; }'),
                        text.replace('nAlphaCorr 1;', '$inherited;'),
                        text.replace('"alpha.*"', '"alpha.other"'),
                        text + '\nsolvers { alpha {} }'):
            with self.assertRaises(ValueError):
                probe.enable_probe(changed)

    def test_wrong_executable_and_library_hashes_fail_before_runtime_check(self):
        # Restore the actual build checker; fake files never invoke nm/ldd.
        patch.stopall()
        (self.build / 'bin').mkdir()
        (self.build / 'lib').mkdir()
        exe = self.build / 'bin/putidaFoamSC1AlphaAudit'
        lib = self.build / 'lib/libbirdphaseSystem.so'
        exe.write_text('wrong executable')
        lib.write_text('wrong library')
        with self.assertRaisesRegex(ValueError, 'Binary differs.*bin/'):
            probe.check_build(self.build)
        with patch.object(probe, 'EXE_SHA', probe.sha256(exe)), self.assertRaisesRegex(ValueError, 'Binary differs.*lib/'):
            probe.check_build(self.build)


if __name__ == '__main__':
    unittest.main(verbosity=2)
