#!/usr/bin/env python3

"""Build a complete BiRD stirred-tank OpenFOAM mesh.

Workflow
--------
1. Write blockMeshDict from the reactor YAML.
2. Run blockMesh.
3. Perfect-stitch the duplicated conformal hub interfaces.
4. Rotate the native z-axis mesh onto the production +y vertical axis.
5. Generate the active MRF cellZone from YAML.
6. Run topoSet.
7. Run checkMesh.

The broad blockMesh MRF envelope is rotorEnvelope.  The physically active
MRF zone is generated separately (normally rotor).
"""

import argparse
import shutil
import subprocess
from pathlib import Path

from bird.meshing._stirred_tank_reactor import StirredTankReactor
from bird.meshing.mrf_zone_from_yaml import write_mrf_properties, write_toposet
from bird.meshing.stirred_tank_mesh import (
    write_blocks,
    write_edges,
    write_ofoam_preamble,
    write_patches,
    write_vertices,
)
from bird.utilities.parser import parse_yaml


STITCH_PAIRS = (
    ("inside_to_hub", "inside_to_hub_copy"),
    ("hub_to_rotor", "hub_to_rotor_copy"),
)


def run_command(command, log_path):
    """Run a command, stream output to stdout, and preserve a case log."""

    print()
    print("+ " + " ".join(str(x) for x in command))

    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(x) for x in command],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert process.stdout is not None

        for line in process.stdout:
            print(line, end="")
            log.write(line)

        return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            f"Command failed with return code {return_code}: "
            + " ".join(str(x) for x in command)
        )


def require_case_dictionaries(case):
    """stitchMesh constructs an fvMesh and requires these dictionaries."""

    required = (
        case / "system" / "controlDict",
        case / "system" / "fvSchemes",
        case / "system" / "fvSolution",
    )

    missing = [str(path) for path in required if not path.is_file()]

    if missing:
        joined = "\n  ".join(missing)
        raise FileNotFoundError(
            "The stirred-tank mesh workflow requires a complete OpenFOAM "
            "case before stitching. Missing:\n  "
            + joined
        )


def write_block_mesh_dict(yaml_path, case):
    reactor = StirredTankReactor.from_file(str(yaml_path))

    path = case / "system" / "blockMeshDict"

    with path.open("w", encoding="utf-8") as stream:
        write_ofoam_preamble(stream, reactor)
        write_vertices(stream, reactor)
        write_edges(stream, reactor)
        write_blocks(stream, reactor)
        write_patches(stream, reactor)

    print(f"Wrote: {path}")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--yaml",
        required=True,
        help="Stirred-tank YAML geometry/mesh configuration",
    )
    parser.add_argument(
        "--case",
        required=True,
        help="Existing OpenFOAM case directory",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete an existing constant/polyMesh before rebuilding",
    )
    parser.add_argument(
        "--skip-check",
        action="store_true",
        help="Do not run the final checkMesh",
    )

    args = parser.parse_args()

    yaml_path = Path(args.yaml).resolve()
    case = Path(args.case).resolve()

    if not yaml_path.is_file():
        raise FileNotFoundError(yaml_path)

    if not case.is_dir():
        raise FileNotFoundError(case)

    require_case_dictionaries(case)

    poly_mesh = case / "constant" / "polyMesh"

    if poly_mesh.exists():
        if not args.force:
            raise RuntimeError(
                f"{poly_mesh} already exists. Refusing to rebuild/rotate an "
                "existing mesh. Use --force only when a fresh rebuild is intended."
            )

        print(f"Removing existing mesh: {poly_mesh}")
        shutil.rmtree(poly_mesh)

    write_block_mesh_dict(yaml_path, case)

    run_command(
        ("blockMesh", "-case", case),
        case / "log.blockMesh",
    )

    for master, slave in STITCH_PAIRS:
        run_command(
            (
                "stitchMesh",
                "-case",
                case,
                "-overwrite",
                "-perfect",
                master,
                slave,
            ),
            case / f"log.stitchMesh.{master}",
        )

    # The Python mesh writer is natively z-up. Production BiRD convention
    # is +y vertical. Apply this exactly once to a freshly generated mesh.
    run_command(
        (
            "transformPoints",
            "-case",
            case,
            "rotate=((0 0 1) (0 1 0))",
        ),
        case / "log.transformPoints",
    )

    cfg = parse_yaml(str(yaml_path))

    write_toposet(cfg, case)
    write_mrf_properties(cfg, case)

    run_command(
        (
            "topoSet",
            "-case",
            case,
            "-dict",
            "system/topoSetDict.mrfFromYaml",
        ),
        case / "log.topoSet.mrf",
    )

    if not args.skip_check:
        run_command(
            (
                "checkMesh",
                "-case",
                case,
                "-allTopology",
                "-allGeometry",
            ),
            case / "log.checkMesh",
        )

    print()
    print("Stirred-tank mesh generation complete.")
    print(f"Case: {case}")


if __name__ == "__main__":
    main()
