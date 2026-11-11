"""Selection ordering, center uncertainty and provenance checks."""
import tempfile
from pathlib import Path
import unittest

import numpy as np

from estimate_sc1_prescribed_fill import estimate, parse_fill, verify_vtk


HEADER = b'FoamFile {format ascii; class dictionary; object setFieldsDict;}\n'
FILL = (b'defaultFieldValues(volScalarFieldValue alpha.liquid 0); regions('
        b'boxToCell {box(-1 -1 -1)(1 1 1);fieldValues(volScalarFieldValue alpha.liquid 0.99);}'
        b'boxToCell {box(-1 0.9 -1)(1 1.1 1);fieldValues(volScalarFieldValue alpha.liquid 0.5);}'
        b');')


class PrescribedFillTests(unittest.TestCase):
    def test_later_region_overrides_and_uses_all_three_coordinates(self):
        centers = [[0, 0, 0], [0, 0.95, 0], [2, 0, 0], [0, 2, 0]]
        c = np.array(centers)
        bounds = np.column_stack([c[:, 0]-.01,c[:, 0]+.01,c[:, 1]-.01,c[:, 1]+.01,c[:, 2]-.01,c[:, 2]+.01])
        result = estimate(c, bounds, [1, 2, 1, 1], parse_fill(HEADER + FILL))
        self.assertAlmostEqual(result['prescribed_liquid_center_estimate_L'], 1990)
        self.assertEqual(result['prescribed_liquid_box_selection_bounds_L'], [1990, 1990])

    def test_uncertain_centroid_selection_is_bounded(self):
        result = estimate([[0, 0.85, 0]], [[-.1, .1, .8, 1.2, -.1, .1]], [1], parse_fill(HEADER + FILL))
        self.assertAlmostEqual(result['prescribed_liquid_center_estimate_L'], 990)
        self.assertEqual(result['prescribed_liquid_box_selection_bounds_L'], [0, 990])

    def test_comments_and_exponents(self):
        x = parse_fill(HEADER + FILL.replace(b'0.99', b'9.9e-1 /* ignore boxToCell */'))
        self.assertEqual(x['regions'][0]['values']['alpha.liquid'], .99)

    def test_unsupported_or_dynamic_configuration_rejected(self):
        for data in [FILL.replace(b'boxToCell',b'cylinderToCell'),
                     FILL.replace(b'0.99',b'$alpha'),
                     FILL.replace(b'0.99',b'nan'),
                     FILL.replace(b'(1 1 1)',b'(-2 -2 -2)'),
                     FILL + b'#include "later"']:
            with self.subTest(data=data), self.assertRaises(ValueError):
                parse_fill(HEADER + data)

    def test_invalid_geometry_or_alpha_rejected(self):
        with self.assertRaises(ValueError):
            estimate([[0, 0, 0]], [[-1, 1, -1, 1, -1, 1]], [0], parse_fill(HEADER + FILL))
        with self.assertRaises(ValueError):
            estimate([[0, 0, 0]], [[-1, 1, -1, 1, -1, 1]], [1],
                     parse_fill(HEADER + FILL.replace(b'0.99', b'1.1')))

    def test_changed_vtk_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'mesh.vtk'
            p.write_bytes(b'changed')
            with self.assertRaises(ValueError):
                verify_vtk(p, dict(sha256='0'*64, bytes=7))


if __name__ == '__main__':
    unittest.main()
