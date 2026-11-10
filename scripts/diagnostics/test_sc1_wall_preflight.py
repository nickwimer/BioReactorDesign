"""Adversarial mesh-format and surface-measurement tests; no CFD required."""
import gzip
import contextlib
import io
import json
from pathlib import Path
import tempfile
import tarfile
import unittest
from unittest.mock import patch

import numpy as np

import collect_sc1_wall_preflight as collect
from sc1_wall_mesh import Reader, zone_adjacency


def header(name, binary=False, bits=32, endian='LSB'):
    cls = {'owner': 'labelList', 'neighbour': 'labelList',
           'cellZones': 'regIOobject', 'boundary': 'polyBoundaryMesh'}[name]
    return (f'/* comment {{ignored}} */\nFoamFile {{ version 2.0; '
            f'format {"binary" if binary else "ascii"}; '
            f'arch "{endian};label={bits};scalar=64"; '
            f'class {cls}; object {name}; }}\n').encode()


def labels(values, binary=False, bits=32, endian='LSB'):
    if binary:
        if not values:
            return b'0\n'
        dtype = ('<' if endian == 'LSB' else '>') + 'i' + str(bits // 8)
        return str(len(values)).encode() + b'\n(' + np.array(values, dtype=dtype).tobytes() + b')'
    return str(len(values)).encode() + b'(' + ' '.join(map(str, values)).encode() + b')'


def mesh_fixture(root, binary=False, bits=32, endian='LSB'):
    opts = dict(binary=binary, bits=bits, endian=endian)
    data = {
        'owner': labels([0, 1, 0, 1, 2, 3], **opts),
        'neighbour': labels([4, 4], **opts),
        'cellZones': b'2( empty { type cellZone; cellLabels List<label> '
                     + labels([], **opts) + b'; } rotor { type cellZone; cellLabels List<label> '
                     + labels([0, 2, 4], **opts) + b'; })',
        'boundary': b'3( wall { type wall; inGroups 1(wall); nFaces 2; startFace 2; } '
                    b'propeller { type wall; nFaces 1; startFace 4; } '
                    b'outlet { type patch; nFaces 1; startFace 5; })',
    }
    for name, payload in data.items():
        (root / name).write_bytes(header(name, **opts) + payload + b'\n// end\n')


class MeshTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_known_adjacency_ascii_binary_both_widths_and_endianness(self):
        for binary, bits, endian in [(False, 32, 'LSB'), (True, 32, 'LSB'),
                                     (True, 64, 'LSB'), (True, 32, 'MSB'), (True, 64, 'MSB')]:
            with self.subTest(binary=binary, bits=bits, endian=endian):
                mesh_fixture(self.root, binary, bits, endian)
                before = {p.name: p.read_bytes() for p in self.root.iterdir()}
                result = zone_adjacency(self.root)
                self.assertEqual(result['zone_cells'], 3)
                self.assertEqual([p['zone_adjacent_faces'] for p in result['patches']], [1, 1, 0])
                self.assertEqual(before, {p.name: p.read_bytes() for p in self.root.iterdir()})

    def test_binary_first_byte_whitespace_is_not_skipped(self):
        for first in (9, 10, 13, 32):
            reader = Reader(header('owner', True) + labels([first, 400], True))
            np.testing.assert_array_equal(reader.label_list(), [first, 400])
            reader.finish()

    def test_ascii_uniform_and_comments(self):
        r = Reader(header('owner') + b'3{17}')
        np.testing.assert_array_equal(r.label_list(), [17, 17, 17])
        r = Reader(header('owner') + b'2(3 /* 99 */ 4)')
        np.testing.assert_array_equal(r.label_list(), [3, 4])

    def test_gzipped_mesh(self):
        mesh_fixture(self.root, True)
        for p in list(self.root.iterdir()):
            p.with_name(p.name + '.gz').write_bytes(gzip.compress(p.read_bytes()))
            p.unlink()
        self.assertEqual(zone_adjacency(self.root)['patches'][0]['zone_adjacent_faces'], 1)

    def test_truncated_wrong_count_and_unknown_arch_fail(self):
        bad = [header('owner', True) + b'3(\x00\x00)',
               header('owner') + b'3(0 1)',
               header('owner') + b'2(0 -1)',
               header('owner', True).replace(b'label=32', b'label=16') + labels([0], True)]
        for data in bad:
            with self.subTest(data=data), self.assertRaises(ValueError):
                Reader(data).label_list()

    def test_invalid_patch_range_missing_zone_and_directive_fail(self):
        mesh_fixture(self.root)
        with self.assertRaises(ValueError):
            zone_adjacency(self.root, 'absent')
        p = self.root / 'boundary'
        original = p.read_bytes()
        for replacement in (original.replace(b'startFace 2', b'startFace 3'),
                            original.replace(b'nFaces 2', b'nFaces 50'),
                            original.replace(b'type wall;', b'#include "something"')):
            p.write_bytes(replacement)
            with self.assertRaises(ValueError):
                zone_adjacency(self.root)


class SurfaceTests(unittest.TestCase):
    def test_inventory_is_volume_weighted_and_never_clipped(self):
        c = np.array([[0.515, 1.03, 0], [0.515, 1.04, 0], [0, 0.5, 0]])
        data = collect.surface_profiles(c, [1, 3, 2], [0.2, 0.8, 1.1])
        self.assertAlmostEqual(data['total_liquid_L'], 4800)
        self.assertAlmostEqual(data['liquid_above_1_02_L'], 2600)
        self.assertEqual(data['out_of_bounds_cells'], 1)
        row = [r for r in data['sectors'] if r['radial_band_m'] == [0.50, 0.53] and r['sector'] == 0][0]
        self.assertAlmostEqual(next(x for x in row['mean_alpha'] if x is not None), 0.65)

    def test_all_crossings_retained_and_missing_bins_not_bridged(self):
        profile = np.full(len(collect.HEIGHT_CENTERS), np.nan)
        profile[:6] = [1, 0, np.nan, 1, 0, 1]
        rows = collect.all_crossings(profile, 0.5)
        self.assertEqual([r['direction'] for r in rows], ['down', 'down', 'up'])
        self.assertAlmostEqual(rows[0]['height_m'], 0.325)

    def test_azimuthal_variation_survives_sector_measurement(self):
        # Opposite sectors have distinct surfaces; a radial mean would hide them.
        c, v, a = [], [], []
        for angle, height in [(0.1, 1.0), (3.2, 1.2)]:
            for y in collect.HEIGHT_CENTERS:
                c.append([0.515 * np.cos(angle), y, 0.515 * np.sin(angle)])
                v.append(0.001)
                a.append(float(y < height))
        data = collect.surface_profiles(c, v, a)
        rows = [r for r in data['sectors'] if r['radial_band_m'] == [0.5, 0.53]]
        self.assertAlmostEqual(rows[0]['crossings']['0.5'][0]['height_m'], 1.0)
        self.assertAlmostEqual(rows[6]['crossings']['0.5'][0]['height_m'], 1.2)
        self.assertEqual(rows[1]['occupied_height_bins'], 0)
        self.assertEqual(rows[1]['crossings']['0.5'], [])

    def test_invalid_cells_rejected_instead_of_silently_dropped(self):
        for v, a in [([0], [1]), ([1], [float('nan')]), ([float('inf')], [1])]:
            with self.assertRaises(ValueError):
                collect.surface_profiles([[0, 1, 0]], v, a)


class PreservationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_output_refuses_reuse_overlap_and_symlink(self):
        case = self.root / 'case'
        case.mkdir()
        out = self.root / 'new'
        collect.prepare_output(out, [case], scratch_roots=[self.root])
        for target in (out, case / 'report', case):
            with self.assertRaises(ValueError):
                collect.prepare_output(target, [case], scratch_roots=[self.root])
        link = self.root / 'link'
        link.symlink_to(case / 'missing')
        with self.assertRaises(ValueError):
            collect.prepare_output(link, [case], scratch_roots=[self.root])

    def test_vtk_selection_rejects_patch_only_and_ambiguous_volumes(self):
        root = self.root / 'EMM/VTK'
        (root / 'boundary').mkdir(parents=True)
        (root / 'boundary/wall_9.9.vtk').touch()
        with self.assertRaises(ValueError):
            collect.existing_vtk(self.root, 'EMM', '9.9')
        volume = root / 'EMM_9.9.vtk'
        volume.touch()
        self.assertEqual(collect.existing_vtk(self.root, 'EMM', '9.9'), volume)
        (root / 'other_9.9.vtk').touch()
        with self.assertRaises(ValueError):
            collect.existing_vtk(self.root, 'EMM', '9.9')

    def test_boundary_query_has_no_expand_or_mutation_and_inputs_unchanged(self):
        case, out = self.root / 'case', self.root / 'report'
        (case / '9.9').mkdir(parents=True)
        out.mkdir()
        field = case / '9.9/U.liquid'
        field.write_text('fixture')
        before = field.read_bytes()
        result = type('Result', (), dict(stdout='boundaryField {}', stderr='', returncode=0))()
        with patch.object(collect.subprocess, 'run', return_value=result) as run:
            collect.collect_configuration(case, 'EMM', out, '/bin/foamDictionary', [])
            self.assertEqual(run.call_args.args[0],
                             ['/bin/foamDictionary', str(field), '-entry', 'boundaryField'])
        self.assertEqual(before, field.read_bytes())

    def test_missing_vtks_still_leave_inspectable_evidence_bundle(self):
        case, out = self.root / 'case', self.root / 'report'
        mesh = case / 'constant/polyMesh'
        mesh.mkdir(parents=True)
        mesh_fixture(mesh, binary=True)
        before = {p.name: p.read_bytes() for p in mesh.iterdir()}
        original_prepare = collect.prepare_output
        with patch.object(collect, 'prepare_output',
                          side_effect=lambda p, sources: original_prepare(p, sources, [self.root])), \
                contextlib.redirect_stdout(io.StringIO()):
            status = collect.main(['--out', str(out), '--emm', str(case), '--meml', str(case),
                                   '--vtk-root', str(self.root / 'missing_vtks')])
        self.assertEqual(status, 0)
        report = json.loads((out / 'report.json').read_text())
        self.assertEqual(report['cases']['EMM']['zone_adjacency']['zone_cells'], 3)
        self.assertEqual(report['cases']['EMM']['surfaces'], {})
        self.assertTrue(any('surface analysis unavailable' in item for item in report['review_items']))
        with tarfile.open(out / 'evidence.tar.gz') as archive:
            self.assertEqual(set(archive.getnames()),
                             {'summary.txt', 'report.json', 'EMM_configuration.txt', 'MEML_BASE_configuration.txt'})
        self.assertEqual(before, {p.name: p.read_bytes() for p in mesh.iterdir()})


if __name__ == '__main__':
    unittest.main()
