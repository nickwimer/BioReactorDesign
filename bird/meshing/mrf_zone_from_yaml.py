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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml", required=True)
    parser.add_argument("--case", required=True)
    args = parser.parse_args()

    yaml_path = Path(args.yaml).resolve()
    case = Path(args.case).resolve()

    yaml = YAML(typ="safe")
    with yaml_path.open("r", encoding="utf-8") as f:
        cfg = yaml.load(f)

    geom = cfg["geometry"]
    mesh = cfg["mesh"]

    units = geom["units"]
    scale = UNIT_TO_M[units]

    nz = float(mesh["nz"])
    nbuf = int(mesh["mrf_vertical_buffer_cells"])

    # Same convention used in BiRD:
    # buffer height = requested nominal cells / nz
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

    out = case / "system" / "topoSetDict.mrfFromYaml"

    text = """FoamFile
{
    format      ascii;
    class       dictionary;
    object      topoSetDict;
}

actions
(
    // Save the full validated radial MRF parent zone.
    {
        name    __rotorOriginal;
        type    cellSet;
        action  new;
        source  zoneToCell;
        sourceInfo
        {
            name rotor;
        }
    }

"""

    for i, (name, lo, hi) in enumerate(bounds):
        action = "new" if i == 0 else "add"

        text += f"""    // {name} impeller
    {{
        name    __rotorV3Cells;
        type    cellSet;
        action  {action};
        source  boxToCell;
        sourceInfo
        {{
            box (-1 {lo:.8f} -1) (1 {hi:.8f} 1);
        }}
    }}

"""

    text += """    // Keep the parent rotor's radial definition exactly.
    {
        name    __rotorV3Cells;
        type    cellSet;
        action  subset;
        source  cellToCell;
        sourceInfo
        {
            set __rotorOriginal;
        }
    }

    {
        name    rotorV3;
        type    cellZoneSet;
        action  new;
        source  setToCellZone;
        sourceInfo
        {
            set __rotorV3Cells;
        }
    }
);
"""

    out.write_text(text, encoding="utf-8")

    print(f"Wrote: {out}")
    print(
        f"mrf_vertical_buffer_cells = {nbuf} "
        f"({buffer_height:g} {units})"
    )

    for name, lo, hi in bounds:
        print(
            f"{name:>8s}: "
            f"{lo:.8f} <= y <= {hi:.8f} m"
        )


if __name__ == "__main__":
    main()
