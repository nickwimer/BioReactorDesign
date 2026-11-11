#!/usr/bin/env python3
"""Strict byte-proof/copy regressions; synthetic containers, no field decoding."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import prepare_sc1_conservation_pair as pair
from sc1_checkpoint_metadata import verify_courant_location_metadata


CO_FIELDS = ("Co.bubbleGas", "Co.headspaceGas", "Co.liquid", "Co.mixture")


def header(path):
    # Include semicolons in a quoted arch value and comments containing braces
    # so a simplistic split/brace regex cannot stand in for header tokenization.
    return (
        '/* Synthetic byte fixture: { not a dictionary } */\n'
        'FoamFile\n{\n'
        '    version 2.0;\n'
        '    format binary;\n'
        '    arch "LSB;label=32;scalar=64";\n'
        '    class decomposedBlockData; // outer container\n'
        f'    location "{path.resolve()}";\n'
        f'    object {path.name};\n'
        '}'
    ).encode("ascii")


def write_container(path, payload=b"\n// synthetic tail\n(\x00\xff\x01\x02)\n"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header(path) + payload)
    return path


def digest(content):
    return hashlib.sha256(content).hexdigest()


class CourantLocationProofTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.control = self.root / "run_control_05"
        self.normal = self.root / "run_normalactive_05"
        self.control.mkdir()
        self.normal.mkdir()

    def files(self, name="Co.liquid", payload=b"\n\x00\xffidentical tail\n"):
        return (
            write_container(self.control / name, payload),
            write_container(self.normal / name, payload),
        )

    def test_all_four_fields_prove_only_five_byte_location_length_difference(self):
        payload = b"\n// payload, not parsed\n" + bytes(range(256))
        for name in CO_FIELDS:
            with self.subTest(name=name):
                control, normal = self.files(name, payload)
                before = (control.read_bytes(), normal.read_bytes())
                self.assertEqual(len(before[1]) - len(before[0]), 5)
                proof = verify_courant_location_metadata(control, normal)
                self.assertEqual(proof["control_sha256"], digest(before[0]))
                self.assertEqual(proof["normal_sha256"], digest(before[1]))
                self.assertNotEqual(proof["control_sha256"], proof["normal_sha256"])
                self.assertEqual(proof["control_location"], str(control.resolve()))
                self.assertEqual(proof["normal_location"], str(normal.resolve()))
                self.assertEqual(proof["payload_sha256"], digest(payload))
                self.assertEqual(proof["payload_bytes"], len(payload))
                self.assertRegex(proof["normalized_header_sha256"], r"^[0-9a-f]{64}$")
                self.assertEqual((control.read_bytes(), normal.read_bytes()), before)

    def test_late_binary_difference_is_rejected(self):
        payload = b"\n" + b"\x00" * (4 * 1024 * 1024 + 31)
        control, normal = self.files(payload=payload)
        normal.write_bytes(header(normal) + payload[:-1] + b"\x01")
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)

    def test_non_location_header_change_is_rejected(self):
        control, normal = self.files()
        normal.write_bytes(normal.read_bytes().replace(b"version 2.0;", b"version 2.1;"))
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)

    def test_comment_or_whitespace_change_is_not_silently_normalized(self):
        control, normal = self.files()
        normal.write_bytes(normal.read_bytes().replace(b"    version", b"     version", 1))
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)

    def test_missing_duplicate_unquoted_and_wrong_location_are_rejected(self):
        for mode in ("missing", "duplicate", "unquoted", "directory", "other_file"):
            with self.subTest(mode=mode):
                control, normal = self.files()
                entry = f'    location "{normal.resolve()}";\n'.encode("ascii")
                replacement = {
                    "missing": b"",
                    "duplicate": entry + entry,
                    "unquoted": f"    location {normal.resolve()};\n".encode("ascii"),
                    "directory": f'    location "{normal.parent.resolve()}";\n'.encode("ascii"),
                    "other_file": f'    location "{control.resolve()}";\n'.encode("ascii"),
                }[mode]
                normal.write_bytes(normal.read_bytes().replace(entry, replacement, 1))
                with self.assertRaises(ValueError):
                    verify_courant_location_metadata(control, normal)

    def test_wrong_type_format_object_and_malformed_headers_are_rejected(self):
        variants = (
            (b"class decomposedBlockData;", b"class volScalarField;"),
            (b"format binary;", b"format ascii;"),
            (b"object Co.liquid;", b"object Co.mixture;"),
            (b"version 2.0;", b"version 2.0; version 2.0;"),
            (b"version 2.0;", b"nested { version 2.0; }"),
            (b"FoamFile", b"unexpectedToken FoamFile"),
            (b"// outer container\n", b"/* unterminated comment\n"),
        )
        for old, new in variants:
            with self.subTest(replacement=new):
                control, normal = self.files()
                normal.write_bytes(normal.read_bytes().replace(old, new, 1))
                with self.assertRaises(ValueError):
                    verify_courant_location_metadata(control, normal)

    def test_oversized_header_is_rejected(self):
        control, normal = self.files()
        for path in (control, normal):
            path.write_bytes(b"/*" + b"x" * 8192 + b"*/\n" + path.read_bytes())
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)

    def test_forbidden_field_and_compressed_filename_are_rejected(self):
        for name in ("O2.liquid", "alpha.bubbleGas", "Co.other", "Co.liquid.gz"):
            with self.subTest(name=name):
                control, normal = self.files(name)
                with self.assertRaises(ValueError):
                    verify_courant_location_metadata(control, normal)

    def test_payload_location_text_is_compared_without_replacement(self):
        control, normal = self.files()
        for path in (control, normal):
            path.write_bytes(header(path) + f'\nlocation "{path.resolve()}";\n'.encode("ascii"))
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)

    def test_extra_tail_bytes_are_rejected(self):
        control, normal = self.files()
        normal.write_bytes(normal.read_bytes() + b"\x00")
        with self.assertRaises(ValueError):
            verify_courant_location_metadata(control, normal)


class MatchedCheckpointCopyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.control = self.root / "run_control_05"
        self.normal = self.root / "run_normalactive_05"
        self.time_section = f"processors52/{pair.CHECKPOINT}"
        required = (
            "alpha.liquid", "alpha.bubbleGas", "alpha.headspaceGas", "O2.liquid",
            "O2.bubbleGas", "O2.headspaceGas", "f.bubbleGas",
        )
        for case, efficiency, gate in (
            (self.control, "0", "requireDownwardLiquid true;"),
            (self.normal, "0.01", "requireDownwardLiquid false; entrainmentVelocityCriterion normalGasIntoLiquid;"),
        ):
            for section in ("constant", "system", "processors52/constant", self.time_section):
                (case / section).mkdir(parents=True)
            (case / "processors52/constant/mesh").write_bytes(b"same mesh bytes\x00")
            for field in required:
                (case / self.time_section / field).write_bytes(b"same physical-state bytes\x00")
            # Preserve the observed 59-file checkpoint size, including all
            # four derived diagnostics and every otherwise identical file.
            for index in range(48):
                (case / self.time_section / f"state{index:02d}").write_bytes(bytes([index]))
            for field in CO_FIELDS:
                write_container(case / self.time_section / field)
            (case / "constant/phaseProperties").write_text(
                f"entrainmentEfficiency {efficiency}; {gate} "
                "maxGlobalEntrainedKgS 0.002; nucleationDiameter 0.002; foamAreaMultiplier 0;"
            )
            (case / "constant/microbeUpdateProperties").write_text(
                "reaction_update_time 4;\nfluid_update_time 200;\nph_update_time 0.1;\n"
            )
            (case / "system/fvSolution").write_text("PIMPLE {}\n")
            (case / "system/fvSchemes").write_text("ddtSchemes { default Euler; }\n")
            (case / "system/controlDict").write_text(
                "startFrom startTime;\nstartTime 31.199982372783058;\nendTime 31.7;\n"
                "writeControl runTime;\nwriteInterval 0.1;\n"
            )

    def validate(self, proofs=None):
        return pair.validate_pair(
            self.control, self.normal, pair.CHECKPOINT,
            allow_courant_location_metadata=True, metadata_proofs=proofs,
        )

    def test_default_remains_strict_and_opt_in_records_four_proofs(self):
        with self.assertRaises(ValueError):
            pair.validate_pair(self.control, self.normal, pair.CHECKPOINT)
        proofs = {}
        evidence = self.validate(proofs)
        self.assertEqual(set(proofs), set(CO_FIELDS))
        self.assertEqual(len(evidence[self.time_section]), 59)
        self.assertEqual(evidence[self.time_section], pair.hashes(self.control / self.time_section))

    def test_opt_in_cannot_hide_physical_field_or_mesh_difference(self):
        for relative in (f"{self.time_section}/O2.liquid", "processors52/constant/mesh"):
            with self.subTest(file=relative):
                path = self.normal / relative
                original = path.read_bytes()
                path.write_bytes(original + b"changed")
                with self.assertRaises(ValueError):
                    self.validate({})
                path.write_bytes(original)

    def test_opt_in_cannot_hide_missing_files(self):
        (self.normal / self.time_section / "state00").unlink()
        with self.assertRaises(ValueError):
            self.validate({})

    def test_main_copies_all_control_bytes_to_both_without_changing_sources(self):
        before = {str(case): pair.hashes(case) for case in (self.control, self.normal)}
        destination = self.root / "new_pair"
        argv = [
            "prepare", "--control", str(self.control), "--active", str(self.normal),
            "--out", str(destination), "--allow-courant-location-metadata",
        ]
        output = io.StringIO()
        with (
            patch.object(sys, "argv", argv),
            patch.object(pair, "scratch", side_effect=lambda p: Path(p).resolve()),
            contextlib.redirect_stdout(output),
        ):
            status = pair.main()
        self.assertEqual(status, 0, output.getvalue())
        expected = pair.hashes(self.control / self.time_section)
        self.assertEqual(len(expected), 59)
        for case_name in ("control", "normal"):
            copied = destination / case_name
            self.assertEqual(pair.hashes(copied / self.time_section), expected)
            self.assertEqual(
                pair.hashes(copied / "processors52/constant"),
                pair.hashes(self.control / "processors52/constant"),
            )
            self.assertEqual(
                (copied / "constant/microbeUpdateProperties").read_bytes(),
                (self.control / "constant/microbeUpdateProperties").read_bytes(),
            )
        self.assertEqual({str(case): pair.hashes(case) for case in (self.control, self.normal)}, before)
        manifest = json.loads((destination / "comparison_manifest.json").read_text())
        self.assertEqual(manifest["initial_state_source"], str(self.control.resolve()))
        self.assertEqual(set(manifest["courant_location_metadata_proofs"]), set(CO_FIELDS))
        self.assertEqual(manifest["initial_state_sha256"][self.time_section], expected)

    def test_short_diagnostic_changes_only_new_duration_and_preserves_biology(self):
        before = {str(case): pair.hashes(case) for case in (self.control, self.normal)}
        destination = self.root / "short_pair"
        argv = ["prepare", "--control", str(self.control), "--active", str(self.normal),
                "--out", str(destination), "--allow-courant-location-metadata", "--duration", "0.1"]
        with (patch.object(sys, "argv", argv),
              patch.object(pair, "scratch", side_effect=lambda p: Path(p).resolve()),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(pair.main(), 0)
        manifest = json.loads((destination / "comparison_manifest.json").read_text())
        self.assertEqual(manifest["duration_s"], 0.1)
        for name, original in (("control", self.control), ("normal", self.normal)):
            copied = destination / name
            compact = ' '.join((copied / "system/controlDict").read_text().split())
            self.assertIn(f"endTime {float(pair.CHECKPOINT)+0.1:.14g};", compact)
            self.assertEqual(pair.hashes(copied / "constant"), pair.hashes(original / "constant"))
            self.assertEqual(pair.hashes(copied / self.time_section),
                             pair.hashes(self.control / self.time_section))
        self.assertEqual({str(case): pair.hashes(case) for case in (self.control, self.normal)}, before)

    def test_alpha_audit_is_explicit_and_only_changes_new_solver_controls(self):
        settings = 'solvers { "alpha.*" { nAlphaCorr 1; nAlphaSubCycles 6; meanFluxReference true; } }\nPIMPLE {}\n'
        for case in (self.control, self.normal):
            (case / 'system/fvSolution').write_text(settings)
        before = {str(case): pair.hashes(case) for case in (self.control, self.normal)}
        destination = self.root / 'alpha_pair'
        argv = ['prepare', '--control', str(self.control), '--active', str(self.normal),
                '--out', str(destination), '--allow-courant-location-metadata',
                '--duration', '0.1', '--alpha-audit']
        with (patch.object(sys, 'argv', argv),
              patch.object(pair, 'scratch', side_effect=lambda p: Path(p).resolve()),
              contextlib.redirect_stdout(io.StringIO())):
            self.assertEqual(pair.main(), 0)
        manifest = json.loads((destination / 'comparison_manifest.json').read_text())
        self.assertTrue(manifest['alpha_audit'])
        for name, source in (('control', self.control), ('normal', self.normal)):
            copied = destination / name
            self.assertEqual(pair.hashes(copied/'constant'), pair.hashes(source/'constant'))
            self.assertEqual(pair.hashes(copied/self.time_section), pair.hashes(self.control/self.time_section))
            self.assertEqual((copied/'system/fvSolution').read_text(),
                             pair.enable_alpha_audit(pair.enable_audit(settings)))
        self.assertEqual({str(case): pair.hashes(case) for case in (self.control, self.normal)}, before)

    def test_unreviewed_duration_is_rejected_before_preparation(self):
        destination = self.root / "must_not_exist"
        with (patch.object(sys, "argv", ["prepare", "--out", str(destination), "--duration", "200"]),
              contextlib.redirect_stderr(io.StringIO())):
            with self.assertRaises(SystemExit) as raised:
                pair.main()
        self.assertEqual(raised.exception.code, 2)
        self.assertFalse(destination.exists())

    def test_failed_proof_does_not_create_destination(self):
        path = self.normal / self.time_section / "Co.liquid"
        path.write_bytes(path.read_bytes() + b"payload change")
        destination = self.root / "must_not_exist"
        argv = [
            "prepare", "--control", str(self.control), "--active", str(self.normal),
            "--out", str(destination), "--allow-courant-location-metadata",
        ]
        with (
            patch.object(sys, "argv", argv),
            patch.object(pair, "scratch", side_effect=lambda p: Path(p).resolve()),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            self.assertEqual(pair.main(), 1)
        self.assertFalse(destination.exists())

    def test_raw_inspection_cannot_be_combined_with_metadata_acceptance(self):
        with (
            patch.object(sys, "argv", ["prepare", "--inspect-inputs", "--allow-courant-location-metadata"]),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            with self.assertRaises(SystemExit) as raised:
                pair.main()
        self.assertEqual(raised.exception.code, 2)

    def test_raw_inspection_rejects_alpha_audit(self):
        with (patch.object(sys, 'argv', ['prepare', '--inspect-inputs', '--alpha-audit']),
              contextlib.redirect_stderr(io.StringIO())):
            with self.assertRaises(SystemExit) as raised:
                pair.main()
        self.assertEqual(raised.exception.code, 2)


class AlphaControlsTests(unittest.TestCase):
    def test_actual_quoted_alpha_pattern_and_comment_braces(self):
        original = '// solvers { alpha { } }\nsolvers { "alpha.*" { nAlphaCorr 1; /* } */ } p { solver GAMG; } }\nPIMPLE {}'
        result = pair.enable_alpha_audit(original)
        self.assertIn('nAlphaCorr 1; /* } */ \n        sc1AlphaAudit true;', result)
        self.assertEqual(result.replace('\n        sc1AlphaAudit true;\n    ', ''), original)

    def test_plain_alpha_is_supported(self):
        self.assertIn('sc1AlphaAudit true;', pair.enable_alpha_audit('solvers { alpha {} }'))

    def test_missing_duplicate_ambiguous_and_inherited_controls_rejected(self):
        for original in ('PIMPLE {}', 'solvers {}', 'solvers { alpha {} alpha {} }',
                         'solvers { alpha {} "alpha.*" {} }', 'solvers {} solvers {}',
                         'solvers { "alpha.+" {} }', 'solvers { nested { alpha {} } }',
                         'solvers { alpha { sc1AlphaAudit false; } }',
                         'solvers { alpha { "sc1AlphaAudit" false; } }',
                         'solvers { alpha { $base; } }',
                         'solvers { #include "alpha"\n alpha {} }'):
            with self.subTest(original=original), self.assertRaises(ValueError):
                pair.enable_alpha_audit(original)


if __name__ == "__main__":
    unittest.main(verbosity=2)
