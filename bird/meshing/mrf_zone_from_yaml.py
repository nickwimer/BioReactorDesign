#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

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


def active_mrf_bounds(cfg):
    geom = cfg["geometry"]
    mesh = cfg["mesh"]
    mrf = cfg.get("mrf", {})

    units = geom["units"]
    if units not in UNIT_TO_M:
        raise ValueError(f"Unsupported geometry.units: {units}")

    scale = UNIT_TO_M[units]
    nz = float(mesh["nz"])

    if nz <= 0:
        raise ValueError("mesh.nz must be positive")

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

    buffer_height = nbuf / nz

    bounds = []

    for imp in geom["impellers"]:
        blades = imp["blades"]

        name = imp["name"]
        zc = float(imp["center_z"])
        hub_height = float(imp["hub_height"])
        width = float(blades["axial_width"])
        pitch = math.radians(float(blades["pitch_deg"]))

        dz = width * abs(math.cos(pitch))
        dz = max(dz, 1.05 * hub_height)

        lo = (zc - 0.5 * dz - buffer_height) * scale
        hi = (zc + 0.5 * dz + buffer_height) * scale

        bounds.append((name, lo, hi))

    return bounds, nbuf, buffer_height, units, source


def write_toposet(cfg, case):
    mrf = cfg.get("mrf", {})

    parent_zone = str(mrf.get("parent_zone", "rotorEnvelope"))
    active_zone = str(mrf.get("cell_zone", "rotor"))

    bounds, nbuf, buffer_height, units, source = active_mrf_bounds(cfg)

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
    print(
        f"active vertical buffer = {nbuf} nominal cells "
        f"({buffer_height:g} {units}); source={source}"
    )

    for name, lo, hi in bounds:
        print(
            f"{name:>8s}: "
            f"{lo:.8f} <= y <= {hi:.8f} m"
        )


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
