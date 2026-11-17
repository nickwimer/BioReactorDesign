"""The active MRF zone must take whole mesh layers.

topoSet selects the zone by cell centre in a box. A box edge that lies on a
layer of cell centres splits that layer by rounding, so rotating and
stationary cells sit side by side. These tests use the NLR pilot-plant
yaml, whose nominal 3-cell buffer (1.5 in) is 3.5 real layers of the
3/7 in buffer block and so ends exactly on a layer of centres.
"""

import re
from pathlib import Path

import numpy as np
import pytest

from bird.meshing.mrf_zone_from_yaml import (
    active_mrf_bounds,
    impeller_edges,
    layer_mrf_bounds,
    mesh_layer_planes,
    write_toposet,
)
from bird.meshing.stirred_tank_mesh import get_reactor_geom
from bird.utilities.parser import parse_yaml

PILOT_YAML = (
    Path(__file__).parents[2]
    / "bird"
    / "meshing"
    / "stirred_tank_mesh_templates"
    / "nlr_pilot_plant"
    / "tank_par.yaml"
)
INCH = 0.0254


def pilot_planes():
    """Heights (m) of the pilot tank's horizontal mesh planes.

    Built the way the mesh writer builds them: block heights from the yaml
    and int(nz * height) + 1 equal layers per block.
    """
    reactor = get_reactor_geom(str(PILOT_YAML))
    heights = [0.0]
    for z0, z1, n in zip(
        reactor.reacthts[:-1], reactor.reacthts[1:], reactor.meshz
    ):
        heights.extend(np.linspace(z0, z1, int(n) + 1)[1:])
    return INCH * np.array(heights)


def centres(planes):
    return 0.5 * (planes[:-1] + planes[1:])


def test_nominal_bounds_end_on_cell_centres():
    # Documents the defect the layer-based bounds remove.
    cfg = parse_yaml(str(PILOT_YAML))
    bounds, nbuf, _, _, _ = active_mrf_bounds(cfg)
    mid = centres(pilot_planes())
    lower = {name: (lo, hi) for name, lo, hi in bounds}["lower"]
    assert nbuf == 3
    assert np.min(np.abs(mid - lower[0])) < 1e-9
    assert np.min(np.abs(mid - lower[1])) < 1e-9


def test_layer_bounds_take_whole_layers():
    cfg = parse_yaml(str(PILOT_YAML))
    planes = pilot_planes()
    layer = np.min(np.diff(planes))
    edges = impeller_edges(cfg)
    bounds = layer_mrf_bounds(edges, planes, 3)

    for (name, lo_edge, hi_edge), (_, lo, hi) in zip(edges, bounds):
        # Ends are mesh planes, three layers beyond the blade edges.
        i_lo = int(np.argmin(np.abs(planes - lo)))
        i_hi = int(np.argmin(np.abs(planes - hi)))
        assert planes[i_lo] == pytest.approx(lo, abs=1e-12)
        assert planes[i_hi] == pytest.approx(hi, abs=1e-12)
        assert planes[i_lo + 3] == pytest.approx(lo_edge, abs=1e-9)
        assert planes[i_hi - 3] == pytest.approx(hi_edge, abs=1e-9)
        # No cell centre is near either end.
        mid = centres(planes)
        assert np.min(np.abs(mid - lo)) > 0.4 * layer, name
        assert np.min(np.abs(mid - hi)) > 0.4 * layer, name


def test_upper_zone_unchanged_lower_zone_whole_layers():
    # The upper nominal span already rounded to whole layers in practice;
    # the layer-based span reproduces it and gives the lower one the same
    # three layers.
    cfg = parse_yaml(str(PILOT_YAML))
    planes = pilot_planes()
    bounds = {
        name: (lo, hi)
        for name, lo, hi in layer_mrf_bounds(impeller_edges(cfg), planes, 3)
    }
    assert bounds["upper"][0] == pytest.approx(0.713471, abs=2e-6)
    assert bounds["upper"][1] == pytest.approx(0.931180, abs=2e-6)
    assert bounds["lower"][0] == pytest.approx(0.185148, abs=2e-6)
    assert bounds["lower"][1] == pytest.approx(0.341903, abs=2e-6)


def test_tie_goes_outward():
    planes = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    (_, lo, hi) = layer_mrf_bounds([("imp", 1.5, 2.5)], planes, 0)[0]
    assert (lo, hi) == (1.0, 3.0)


def test_buffer_past_mesh_raises():
    planes = np.array([0.0, 1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        layer_mrf_bounds([("imp", 1.0, 2.0)], planes, 2)


def write_points(case, points):
    mesh_dir = case / "constant" / "polyMesh"
    mesh_dir.mkdir(parents=True)
    lines = "\n".join(f"({x:.9g} {y:.9g} {z:.9g})" for x, y, z in points)
    (mesh_dir / "points").write_text(
        "FoamFile\n{\n    format      ascii;\n    class       vectorField;\n"
        '    location    "constant/polyMesh";\n    object      points;\n}\n'
        f"// * * //\n\n{len(points)}\n(\n{lines}\n)\n\n// ** //\n",
        encoding="utf-8",
    )


def test_mesh_layer_planes_reads_points(tmp_path):
    # Three full planes of 16 points and a curved-bottom ring of 4.
    ring = [(np.cos(a), np.sin(a)) for a in np.linspace(0, 6, 16)]
    points = [(x, y, z) for y in (0.1, 0.2, 0.35) for x, z in ring]
    points += [(0.5 * x, -0.05, 0.5 * z) for x, z in ring[:4]]
    write_points(tmp_path, points)
    planes = mesh_layer_planes(tmp_path)
    np.testing.assert_allclose(planes, [0.1, 0.2, 0.35])


def test_mesh_layer_planes_without_mesh(tmp_path):
    assert mesh_layer_planes(tmp_path) is None


def test_write_toposet_uses_mesh_planes(tmp_path):
    cfg = parse_yaml(str(PILOT_YAML))
    planes = pilot_planes()
    ring = [(0.3 * np.cos(a), 0.3 * np.sin(a)) for a in np.arange(8)]
    write_points(tmp_path, [(x, y, z) for y in planes for x, z in ring])
    write_toposet(cfg, tmp_path)
    text = (tmp_path / "system" / "topoSetDict.mrfFromYaml").read_text()
    lows = [float(v) for v in re.findall(r"box \(-1 ([-\d.]+) -1\)", text)]
    mid = centres(planes)
    assert len(lows) == 2
    for lo in lows:
        assert np.min(np.abs(mid - lo)) > 0.004
