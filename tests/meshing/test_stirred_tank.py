import argparse
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

from bird.meshing.stirred_tank_mesh import (
    get_reactor_geom,
    write_blocks,
    write_edges,
    write_ofoam_preamble,
    write_patches,
    write_vertices,
)


def test_stirred_tank():
    BIRD_STIRRED_TANK_MESH_TEMP_DIR = os.path.join(
        Path(__file__).parent,
        "..",
        "..",
        "bird",
        "meshing",
        "stirred_tank_mesh_templates",
    )
    inp = os.path.join(
        BIRD_STIRRED_TANK_MESH_TEMP_DIR, "base_tank", "tank_par.yaml"
    )
    # Output to temporary directory and delete when done
    with tempfile.TemporaryDirectory() as tmpdirname:
        output_folder = tmpdirname
        outfilename = os.path.join(tmpdirname, "blockMeshDict")
        with open(outfilename, "w") as outfile:
            react = get_reactor_geom(inp)
            write_ofoam_preamble(outfile, react)
            write_vertices(outfile, react)
            write_edges(outfile, react)
            write_blocks(outfile, react)
            write_patches(outfile, react)


def _nlr_reactor(tmp_path, **mesh):
    """NLR pilot tank with extra mesh settings, written to a temp yaml."""
    from ruamel.yaml import YAML

    src = (
        Path(__file__).parent
        / ".."
        / ".."
        / "bird"
        / "meshing"
        / "stirred_tank_mesh_templates"
        / "nlr_pilot_plant"
        / "tank_par.yaml"
    )
    yaml = YAML()
    with open(src) as stream:
        data = yaml.load(stream)
    data["mesh"].update(mesh)
    out = tmp_path / "tank_par.yaml"
    with open(out, "w") as stream:
        yaml.dump(data, stream)
    return get_reactor_geom(str(out))


def _block_mesh_dict(react, path):
    with open(path, "w") as outfile:
        write_ofoam_preamble(outfile, react)
        write_vertices(outfile, react)
        write_edges(outfile, react)
        write_blocks(outfile, react)
        write_patches(outfile, react)
    return path.read_text()


def test_square_core_topology(tmp_path):
    """The square core has no degenerate blocks and closes up."""
    import re
    from collections import Counter

    react = _nlr_reactor(tmp_path, square_core=True)
    text = _block_mesh_dict(react, tmp_path / "blockMeshDict")

    start = text.index("vertices\n(")
    vertices = np.array(
        re.findall(
            r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)",
            text[start : text.index(");", start)],
        ),
        dtype=float,
    )
    start = text.index("blocks\n(")
    blocks = [
        [int(v) for v in ids.split()]
        for ids in re.findall(
            r"hex \(([\d ]+)\)", text[start : text.index("\n);", start)]
        )
    ]
    patch_faces = Counter(
        frozenset(int(v) for v in face.split())
        for face in re.findall(
            r"\(\s*([\d ]+)\)", text[text.index("patches\n(") :]
        )
    )

    hex_faces = [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 7, 6, 2),
        (0, 3, 2, 1),
        (4, 5, 6, 7),
    ]
    faces = Counter()
    for block in blocks:
        assert len(set(block)) == 8, "degenerate block"
        assert max(block) < len(vertices)
        p = vertices[block]
        volume = np.dot(np.cross(p[1] - p[0], p[3] - p[0]), p[4] - p[0])
        assert volume > 0, "left-handed block"
        for face in hex_faces:
            faces[frozenset(block[i] for i in face)] += 1

    assert max(faces.values()) == 2
    boundary = {face for face, count in faces.items() if count == 1}
    assert boundary <= set(patch_faces)
    assert set(patch_faces) <= set(faces)


def test_square_core_fraction_limits(tmp_path):
    """The square's corners must stay inside the shaft-radius polygon."""
    import pytest

    with pytest.raises(ValueError):
        _nlr_reactor(tmp_path, square_core=True, square_core_fraction=0.75)
