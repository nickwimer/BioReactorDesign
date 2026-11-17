#!/usr/bin/env python3
"""Active MRF cellZone and MRFProperties from a stirred-tank yaml.

The active zone is every cell of the mesh's rotor envelope whose centre
lies in one axial span per impeller: the blade span widened by
``mrf.vertical_buffer_cells`` mesh layers above and below.

When the case already has a mesh (``constant/polyMesh/points``), each span
starts and ends on a horizontal mesh plane, so topoSet's centre-in-box test
takes or leaves whole layers. Without a mesh the span falls back to the
nominal ``vertical_buffer_cells / mesh.nz`` geometry units, which can put a
box edge exactly on a layer of cell centres; the mesh writer gives each
block ``int(nz * height) + 1`` layers, so its layers are thinner than
``1 / nz``. topoSet then splits that layer cell by cell, by rounding, and
rotating and stationary cells end up side by side.
"""

import argparse
import gzip
import math
import re
from pathlib import Path

import numpy as np
from ruamel.yaml import YAML


UNIT_TO_M = {
    "m": 1.0,
    "cm": 0.01,
    "mm": 0.001,
    "in": 0.0254,
    "ft": 0.3048,
}


def _vector(values, name):
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values")
    return tuple(float(v) for v in values)


def _format_vector(values):
    return "(" + " ".join(f"{v:g}" for v in values) + ")"


def _buffer_cells(cfg):
    """Active-zone buffer, in cells beyond the blade edges, and its key."""
    mesh = cfg["mesh"]
    mrf = cfg.get("mrf", {})

    # Preferred split schema:
    #
    # mesh.mrf_mesh_vertical_buffer_cells
    #     controls broad blockMesh topology only.
    #
    # mrf.vertical_buffer_cells
    #     controls the active rotating cellZone.
    #
    # The old mesh.mrf_vertical_buffer_cells key remains a compatibility
    # fallback for cases created before the schema was separated.
    if "vertical_buffer_cells" in mrf:
        nbuf = int(mrf["vertical_buffer_cells"])
        source = "mrf.vertical_buffer_cells"
    elif "mrf_vertical_buffer_cells" in mesh:
        nbuf = int(mesh["mrf_vertical_buffer_cells"])
        source = "mesh.mrf_vertical_buffer_cells (legacy)"
    elif "mrf_buffer_cells" in mesh:
        nbuf = int(mesh["mrf_buffer_cells"])
        source = "mesh.mrf_buffer_cells (legacy)"
    else:
        nbuf = 0
        source = "default"

    if nbuf < 0:
        raise ValueError("active MRF vertical_buffer_cells must be >= 0")

    return nbuf, source


def _scale(cfg):
    units = cfg["geometry"]["units"]
    if units not in UNIT_TO_M:
        raise ValueError(f"Unsupported geometry.units: {units}")
    return UNIT_TO_M[units], units


def impeller_edges(cfg):
    """Lower and upper edge of each impeller's blade span, in metres.

    Returns a list of (name, lo, hi). The span is the blades' axial extent,
    or 1.05 hub heights if that is larger, as in the mesh writer.
    """
    scale, _ = _scale(cfg)
    edges = []

    for imp in cfg["geometry"]["impellers"]:
        blades = imp["blades"]

        name = imp["name"]
        zc = float(imp["center_z"])
        hub_height = float(imp["hub_height"])
        width = float(blades["axial_width"])
        pitch = math.radians(float(blades["pitch_deg"]))

        dz = width * abs(math.cos(pitch))
        dz = max(dz, 1.05 * hub_height)

        edges.append((name, (zc - 0.5 * dz) * scale, (zc + 0.5 * dz) * scale))

    return edges


def active_mrf_bounds(cfg):
    """Nominal active spans: blade edges -/+ vertical_buffer_cells / nz.

    This ignores where the mesh layers are and is used only when the case
    has no mesh to read. Returns (bounds, nbuf, buffer_height, units,
    source), with bounds a list of (name, lo, hi) in metres.
    """
    nz = float(cfg["mesh"]["nz"])

    if nz <= 0:
        raise ValueError("mesh.nz must be positive")

    scale, units = _scale(cfg)
    nbuf, source = _buffer_cells(cfg)
    buffer_height = nbuf / nz

    bounds = [
        (name, lo - buffer_height * scale, hi + buffer_height * scale)
        for name, lo, hi in impeller_edges(cfg)
    ]

    return bounds, nbuf, buffer_height, units, source


def mesh_layer_planes(case, axis=1, resolution=1e-7):
    """Heights of the horizontal mesh planes, from the case's points.

    The stirred-tank mesh is built from stacked blocks, so in the
    cylindrical part of the tank every point lies on one of a set of planes
    of constant height; the twist in the impeller region is a rotation
    about the axis and keeps heights. A plane is a height (to within
    ``resolution`` m) shared by at least half as many points as the most
    populated one, which leaves out the single rings of a curved bottom.

    Returns the sorted plane heights in metres, or None when
    ``constant/polyMesh/points`` is missing or not in ASCII format.
    """
    mesh_dir = Path(case) / "constant" / "polyMesh"

    if (mesh_dir / "points").is_file():
        text = (mesh_dir / "points").read_text(encoding="utf-8")
    elif (mesh_dir / "points.gz").is_file():
        with gzip.open(mesh_dir / "points.gz", "rt", encoding="utf-8") as f:
            text = f.read()
    else:
        return None

    if re.search(r"format\s+binary\s*;", text[:2000]):
        return None

    head = re.search(r"^(\d+)\s*\(", text, re.M)
    if not head:
        raise ValueError(f"cannot find the point list in {mesh_dir}/points")

    count = int(head.group(1))
    body = text[head.end() : text.rindex(")")]
    values = np.array(
        body.replace("(", " ").replace(")", " ").split(), dtype=float
    )

    if values.size != 3 * count:
        raise ValueError(
            f"{mesh_dir}/points: expected {count} points, "
            f"read {values.size / 3:g}"
        )

    heights = values.reshape(count, 3)[:, axis]
    keys = np.round(heights / resolution).astype(np.int64)
    _, inverse, counts = np.unique(
        keys, return_inverse=True, return_counts=True
    )
    sums = np.bincount(inverse.ravel(), weights=heights)
    planes = sums / counts

    return np.sort(planes[counts >= 0.5 * counts.max()])


def _nearest_plane(planes, height, outward):
    """Index of the plane nearest height; ties go outward (-1 down, +1 up)."""
    dist = np.abs(planes - height)
    ties = np.flatnonzero(dist <= dist.min() + 1e-9)
    return int(ties[0] if outward < 0 else ties[-1])


def layer_mrf_bounds(edges, planes, nbuf):
    """Active spans that start and end on mesh planes.

    Each blade edge is moved to its nearest plane (blade edges are block
    boundaries, so normally they already lie on one), and the span is
    widened by nbuf whole layers above and below. Every cell centre is then
    half a layer or more from a span's ends.

    edges: list of (name, lo, hi) in metres, from impeller_edges.
    planes: sorted plane heights in metres, from mesh_layer_planes.
    Returns a list of (name, lo, hi) in metres.
    """
    planes = np.asarray(planes, dtype=float)
    bounds = []

    for name, lo_edge, hi_edge in edges:
        i_lo = _nearest_plane(planes, lo_edge, -1) - nbuf
        i_hi = _nearest_plane(planes, hi_edge, +1) + nbuf

        if i_lo < 0 or i_hi >= len(planes):
            raise ValueError(
                f"{name} impeller: a buffer of {nbuf} layers runs past the "
                "mesh planes"
            )

        bounds.append((name, float(planes[i_lo]), float(planes[i_hi])))

    return bounds


def write_toposet(cfg, case):
    mrf = cfg.get("mrf", {})

    parent_zone = str(mrf.get("parent_zone", "rotorEnvelope"))
    active_zone = str(mrf.get("cell_zone", "rotor"))

    planes = mesh_layer_planes(case)

    if planes is not None and len(planes) > 1:
        nbuf, source = _buffer_cells(cfg)
        edges = impeller_edges(cfg)
        bounds = layer_mrf_bounds(edges, planes, nbuf)
        buffer_text = (
            f"{nbuf} mesh layers beyond the blade edges, ends on mesh "
            f"planes; source={source}"
        )
    else:
        bounds, nbuf, buffer_height, units, source = active_mrf_bounds(cfg)
        buffer_text = (
            f"{nbuf} nominal cells ({buffer_height:g} {units}); "
            f"source={source}"
        )
        print(
            "WARNING: no ASCII mesh points in "
            f"{Path(case) / 'constant' / 'polyMesh'}; the active MRF spans "
            "use nominal heights and may end on a layer of cell centres"
        )

    out = case / "system" / "topoSetDict.mrfFromYaml"
    out.parent.mkdir(parents=True, exist_ok=True)

    parent_set = "__mrfEnvelopeCells"
    active_set = "__activeMRFCells"

    text = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}}

actions
(
    // Preserve the broad radial/axial MRF envelope created by blockMesh.
    {{
        name    {parent_set};
        type    cellSet;
        action  new;
        source  zoneToCell;
        sourceInfo
        {{
            name {parent_zone};
        }}
    }}

"""

    for i, (name, lo, hi) in enumerate(bounds):
        action = "new" if i == 0 else "add"

        text += f"""    // {name} impeller active axial span
    {{
        name    {active_set};
        type    cellSet;
        action  {action};
        source  boxToCell;
        sourceInfo
        {{
            box (-1 {lo:.8f} -1) (1 {hi:.8f} 1);
        }}
    }}

"""

    text += f"""    // Preserve the parent envelope's radial definition exactly.
    {{
        name    {active_set};
        type    cellSet;
        action  subset;
        source  cellToCell;
        sourceInfo
        {{
            set {parent_set};
        }}
    }}

    // Production active rotating cellZone.
    {{
        name    {active_zone};
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        sourceInfo
        {{
            set {active_set};
        }}
    }}
);
"""

    out.write_text(text, encoding="utf-8")

    print(f"Wrote: {out}")
    print(f"parent zone: {parent_zone}")
    print(f"active zone: {active_zone}")
    print(f"active vertical buffer = {buffer_text}")

    for name, lo, hi in bounds:
        print(f"{name:>8s}: {lo:.8f} <= y <= {hi:.8f} m")


def write_mrf_properties(cfg, case):
    mrf = cfg.get("mrf")

    if not mrf or "omega_rad_s" not in mrf:
        print("MRFProperties not generated: mrf.omega_rad_s is not configured")
        return

    cell_zone = str(mrf.get("cell_zone", "rotor"))

    origin = _vector(mrf.get("origin", [0, 0, 0]), "mrf.origin")
    axis = _vector(mrf.get("axis", [0, 1, 0]), "mrf.axis")

    omega = float(mrf["omega_rad_s"])
    ramp_time = float(mrf.get("ramp_time_s", 0.0))
    hold_until = float(mrf.get("hold_until_s", 10.0))

    if ramp_time < 0:
        raise ValueError("mrf.ramp_time_s must be >= 0")

    if hold_until <= ramp_time:
        raise ValueError("mrf.hold_until_s must be greater than ramp_time_s")

    non_rotating = mrf.get("nonRotatingPatches", [])
    non_rotating_text = " ".join(str(v) for v in non_rotating)

    if ramp_time > 0:
        omega_entry = f"""omega table
    (
        (0 {0.0:g})
        ({ramp_time:g} {omega:g})
        ({hold_until:g} {omega:g})
    );"""
    else:
        omega_entry = f"omega {omega:g};"

    out = case / "constant" / "MRFProperties"
    out.parent.mkdir(parents=True, exist_ok=True)

    text = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      MRFProperties;
}}

MRF1
{{
    cellZone            {cell_zone};
    active              yes;
    nonRotatingPatches  ({non_rotating_text});

    origin              {_format_vector(origin)};
    axis                {_format_vector(axis)};

    {omega_entry}
}}
"""

    out.write_text(text, encoding="utf-8")

    print(f"Wrote: {out}")
    print(
        f"MRF: zone={cell_zone}, omega={omega:g} rad/s, "
        f"ramp={ramp_time:g} s"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml", required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument(
        "--no-mrf-properties",
        action="store_true",
        help="Generate only topoSetDict.mrfFromYaml",
    )
    args = parser.parse_args()

    yaml_path = Path(args.yaml).resolve()
    case = Path(args.case).resolve()

    yaml = YAML(typ="safe")

    with yaml_path.open("r", encoding="utf-8") as f:
        cfg = yaml.load(f)

    write_toposet(cfg, case)

    if not args.no_mrf_properties:
        write_mrf_properties(cfg, case)


if __name__ == "__main__":
    main()
