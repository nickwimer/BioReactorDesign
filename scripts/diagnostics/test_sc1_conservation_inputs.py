#!/usr/bin/env python3
"""Read-only input-inspection regressions; no OpenFOAM or numeric field decoding."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import prepare_sc1_conservation_pair as pair


class InputInspectionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.control = self.root / "control"
        self.normal = self.root / "normal"
        self.checkpoint = pair.CHECKPOINT
        self.sections = (
            "processors52/constant",
            f"processors52/{self.checkpoint}",
            "constant",
            "system",
        )
        self.time_section = self.sections[1]
        for case in (self.control, self.normal):
            for section in self.sections:
                folder = case / section
                folder.mkdir(parents=True)
                (folder / "same").write_bytes(b"identical container bytes\n")

    def write(self, case, relative, content):
        path = case / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return path

    def inspect(self):
        return pair.inspect_pair_inputs(
            self.control, self.normal, self.checkpoint
        )

    def assert_metadata(self, metadata, path):
        content = path.read_bytes()
        self.assertEqual(metadata["size_bytes"], len(content))
        self.assertEqual(metadata["sha256"], hashlib.sha256(content).hexdigest())
        self.assertEqual(metadata["resolved_path"], str(path.resolve()))

    def assert_excerpt(self, excerpt, content):
        # The evidence is a bounded byte excerpt, never interpreted as scalars.
        offset = excerpt["offset"]
        decoded = bytes.fromhex(excerpt["hex"])
        self.assertGreaterEqual(offset, 0)
        self.assertLessEqual(offset, len(content))
        self.assertLessEqual(len(decoded), 96)
        self.assertEqual(decoded, content[offset:offset + len(decoded)])

    @staticmethod
    def snapshot(root):
        # Exclude atime: reading a file can legitimately update access time.
        return {
            str(path.relative_to(root)): (
                "directory" if path.is_dir() else "file",
                path.stat().st_size,
                path.stat().st_mtime_ns,
                None if path.is_dir() else hashlib.sha256(path.read_bytes()).hexdigest(),
            )
            for path in sorted(root.rglob("*"))
        }

    def test_sets_nested_missing_files_sizes_and_hashes(self):
        control_path = self.write(
            self.control, f"{self.time_section}/O2.liquid", b"control\x00field"
        )
        normal_path = self.write(
            self.normal, f"{self.time_section}/O2.liquid", b"normal\x00field"
        )
        only_control = self.write(
            self.control, f"{self.time_section}/uniform/time", b"timeIndex 123;\n"
        )
        only_normal = self.write(
            self.normal, f"{self.time_section}/uniform/extra", b"extra state\n"
        )
        report = self.inspect()
        self.assertTrue(report["inspection_only"])
        self.assertEqual(report["checkpoint"], self.checkpoint)
        self.assertEqual(report["control_source"], str(self.control.resolve()))
        self.assertEqual(report["normal_source"], str(self.normal.resolve()))
        self.assertEqual(report["read_errors"], 0)
        self.assertEqual(set(report["sections"]), set(self.sections))
        section = report["sections"][self.time_section]
        self.assertEqual(section["files_control"], 3)
        self.assertEqual(section["files_normal"], 3)
        self.assertEqual(section["identical_files"], 1)
        self.assertEqual([x["file"] for x in section["changed"]], ["O2.liquid"])
        changed = section["changed"][0]
        self.assert_metadata(changed["control"], control_path)
        self.assert_metadata(changed["normal"], normal_path)
        self.assertEqual(changed["classification"], "undetermined")
        self.assertEqual([x["file"] for x in section["only_control"]], ["uniform/time"])
        self.assertEqual([x["file"] for x in section["only_normal"]], ["uniform/extra"])
        self.assert_metadata(section["only_control"][0], only_control)
        self.assert_metadata(section["only_normal"][0], only_normal)

    def test_first_difference_beyond_hash_chunk_is_bounded_and_exact(self):
        offset = 4 * 1024 * 1024 + 17
        original = b"\x00" * (offset + 80)
        modified = original[:offset] + b"\xff" + original[offset + 1:]
        self.write(self.control, f"{self.time_section}/binary", original)
        self.write(self.normal, f"{self.time_section}/binary", modified)
        changed = self.inspect()["sections"][self.time_section]["changed"][0]
        first = changed["first_difference"]
        self.assertEqual(first["offset"], offset)
        self.assert_excerpt(first["control"], original)
        self.assert_excerpt(first["normal"], modified)
        self.assertEqual(changed["classification"], "undetermined")

    def test_early_difference_never_claims_header_only(self):
        original = b"FoamFile header_A\n" + b"\x00" * 65536 + b"payload_A"
        modified = b"FoamFile header_B\n" + b"\x00" * 65536 + b"payload_B"
        self.write(self.control, f"{self.time_section}/binary", original)
        self.write(self.normal, f"{self.time_section}/binary", modified)
        changed = self.inspect()["sections"][self.time_section]["changed"][0]
        self.assertEqual(changed["first_difference"]["offset"], 16)
        self.assertEqual(changed["classification"], "undetermined")
        self.assertNotEqual(changed["control"]["sha256"], changed["normal"]["sha256"])

    def test_length_difference_reports_end_of_common_prefix(self):
        self.write(self.control, f"{self.time_section}/short", b"abc")
        self.write(self.normal, f"{self.time_section}/short", b"abcdef")
        changed = self.inspect()["sections"][self.time_section]["changed"][0]
        first = changed["first_difference"]
        self.assertEqual(first["offset"], 3)
        self.assert_excerpt(first["control"], b"abc")
        self.assert_excerpt(first["normal"], b"abcdef")
        self.assertEqual(changed["classification"], "undetermined")

    def test_missing_section_reports_error_without_hiding_other_sections(self):
        shutil.rmtree(self.control / "system")
        report = self.inspect()
        self.assertEqual(report["read_errors"], 1)
        self.assertIn("error", report["sections"]["system"])
        self.assertTrue(report["sections"]["system"]["error"])
        self.assertEqual(report["sections"][self.time_section]["identical_files"], 1)

    def test_inspect_cli_with_differences_is_read_only_and_needs_no_out(self):
        self.write(self.normal, f"{self.time_section}/same", b"different bytes\x1b\x00")
        before = self.snapshot(self.root)
        output = io.StringIO()
        argv = [
            "prepare_sc1_conservation_pair.py", "--inspect-inputs",
            "--control", str(self.control), "--active", str(self.normal),
        ]
        forbidden = AssertionError("Inspection attempted a preparation/write operation")
        with (
            patch.object(sys, "argv", argv),
            patch.object(pair, "scratch", side_effect=lambda p: Path(p).resolve()),
            patch.object(pair, "validate_pair", side_effect=forbidden),
            patch.object(pair, "enable_audit", side_effect=forbidden),
            patch.object(shutil, "copytree", side_effect=forbidden),
            patch.object(shutil, "copy2", side_effect=forbidden),
            patch.object(Path, "mkdir", side_effect=forbidden),
            patch.object(Path, "write_text", side_effect=forbidden),
            patch.object(Path, "write_bytes", side_effect=forbidden),
            contextlib.redirect_stdout(output),
        ):
            status = pair.main()
        self.assertEqual(status, 0)
        report = json.loads(output.getvalue())
        self.assertTrue(report["inspection_only"])
        self.assertEqual(len(report["sections"][self.time_section]["changed"]), 1)
        self.assertNotIn("\x1b", output.getvalue())
        self.assertNotIn("\x00", output.getvalue())
        self.assertEqual(self.snapshot(self.root), before)

    def test_cli_rejects_out_in_inspection_mode(self):
        with (
            patch.object(sys, "argv", ["prepare", "--inspect-inputs", "--out", str(self.root / "unused")]),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            with self.assertRaises(SystemExit) as raised:
                pair.main()
        self.assertEqual(raised.exception.code, 2)
        self.assertFalse((self.root / "unused").exists())

    def test_preparation_still_requires_out(self):
        with (
            patch.object(sys, "argv", ["prepare"]),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            with self.assertRaises(SystemExit) as raised:
                pair.main()
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
