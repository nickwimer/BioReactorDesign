#!/usr/bin/env python3
"""Identify headspace vs bubble gas at the SC1 upper propeller wall.

Reads the *original MEML* 9.9 s scratch VTK mirror only.
With --convert, runs foamToVTK on the existing scratch-only case mirror,
adding alpha.headspaceGas and alpha.bubbleGas to the saved VTK export.
It never runs the solver or writes in the production case.

Usage:
    python sc1_meml_propeller_gas_identity.py --convert
    python sc1_meml_propeller_gas_identity.py
"""

import argparse
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pyvista as pv


DEFAULT_CASE = Path("/scratch/nwimer/sc1_surface_evolution_20261103/MEML_BASE")
FIELDS = ("alpha.liquid", "alpha.headspaceGas", "alpha.bubbleGas")
LOWER = (0.179705, 0.347345)
UPPER = (0.70802753, 0.93662247)


def check_scratch(case: Path, time: str) -> bool:
    """Disallow conversion outside the user's scratch workspace."""
    scratch = Path("/scratch/nwimer").resolve()
    real = case.resolve()
    if not real.is_relative_to(scratch) or real == scratch:
        print(f"Refusing operation outside scratch: {real}")
        return False
    for name in FIELDS:
        if not (case / time / name).is_file():
            print(f"Required source field unavailable: {case / time / name}")
            return False
    if not (case / "constant/polyMesh/points").is_file():
        print(f"Scratch mirror lacks mesh: {case}")
        return False
    return True


def convert(case: Path, time: str) -> bool:
    if not shutil.which("foamToVTK"):
        print("foamToVTK is unavailable; initialize the RHEL9 OpenFOAM environment.")
        return False
    command = [
        "foamToVTK", "-case", str(case), "-time", time,
        "-fields", "(" + " ".join(FIELDS) + ")", "-useTimeName"
    ]
    print("Converting three alpha fields in scratch only:", case, flush=True)
    result = subprocess.run(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=False,
    )
    log = case / f"foamToVTK.threePhase.{time}.log"
    log.write_text(result.stdout)
    print(f"Conversion return code={result.returncode}; log={log}")
    if result.returncode != 0:
        print("\n".join(result.stdout.splitlines()[-35:]))
        return False
    return True


def get_cell_field(mesh: pv.DataSet, name: str) -> np.ndarray:
    if name in mesh.cell_data:
        return np.asarray(mesh.cell_data[name], dtype=float)
    if name in mesh.point_data:
        cells = mesh.point_data_to_cell_data(pass_point_data=False)
        return np.asarray(cells.cell_data[name], dtype=float)
    raise KeyError(f"{name} missing: cell={list(mesh.cell_data)}; point={list(mesh.point_data)}")


def evaluate(patch_file: Path) -> None:
    if not patch_file.is_file():
        print("Expected propeller patch VTK missing:", patch_file)
        return

    original = pv.read(patch_file)
    mesh = original if isinstance(original, pv.PolyData) else original.extract_surface()
    try:
        al, ah, ab = [get_cell_field(mesh, name) for name in FIELDS]
    except KeyError as exc:
        print("Phase-resolved patch fields unavailable:", exc)
        return

    face_area = np.asarray(
        mesh.compute_cell_sizes(length=False, area=True, volume=False)
        .cell_data["Area"], dtype=float,
    )
    center = mesh.cell_centers().points
    x, y, z = center[:, 0], center[:, 1], center[:, 2]
    radius = np.hypot(x, z)
    good = (
        (face_area > 0) & np.isfinite(face_area)
        & np.isfinite(al) & np.isfinite(ah) & np.isfinite(ab)
    )

    # For radial blades, broad vertical blade faces are approximately
    # circumferential-normal. A geometric orientation proxy is NOT
    # an authoritative separation of hub, shaft and blade surfaces.
    candidate = np.zeros(len(face_area), dtype=bool)
    try:
        normal_mesh = mesh.compute_normals(
            cell_normals=True, point_normals=False,
            auto_orient_normals=False, consistent_normals=False,
        )
        normal = np.asarray(normal_mesh.cell_data["Normals"], dtype=float)
        r_safe = np.maximum(radius, 1e-8)
        abs_normal_theta = np.abs((-z*normal[:, 0] + x*normal[:, 2]) / r_safe)
        candidate = abs_normal_theta >= 0.8
        print("Geometric blade-face proxy: |normal dot e_theta| >= 0.8")
    except (KeyError, AttributeError, ValueError, TypeError) as exc:
        print("Face-normal proxy unavailable:", exc)

    all_area = float(np.sum(face_area[good]))
    sum_error = np.max(np.abs(al[good] + ah[good] + ab[good] - 1.0))
    print(f"Patch VTK: {patch_file}")
    print(f"Patch faces: {mesh.n_cells}; valid area: {all_area:.5f} m2")
    print(f"Max per-face |alphaL+alphaH+alphaB-1|: {sum_error:.6g}")

    print(
        "Region / surface type               Area m2    alphaL"
        "    alphaH    alphaB   H>0.5%   B>0.5%   L>0.5%"
    )

    def report(title: str, y0: float, y1: float,
               r0: float, r1: float, oriented: bool = False) -> None:
        chosen = good & (y >= y0) & (y < y1) & (radius >= r0) & (radius < r1)
        if oriented:
            chosen &= candidate
        a = float(np.sum(face_area[chosen]))
        if a == 0:
            print(f"{title:35s}  no faces")
            return
        means = [
            float(np.sum(face_area[chosen] * field[chosen]) / a)
            for field in (al, ah, ab)
        ]
        fractions = [
            float(100*np.sum(face_area[chosen & (field >= 0.5)]) / a)
            for field in (ah, ab, al)
        ]
        print(
            f"{title:35s} {a:8.5f}"
            f"   {means[0]:7.3f}  {means[1]:7.3f}  {means[2]:7.3f}"
            f"    {fractions[0]:6.1f}   {fractions[1]:6.1f}   {fractions[2]:6.1f}"
        )

    report("Lower window, all r<0.30 m", *LOWER, 0, 0.30)
    report("Upper window, all r<0.30 m", *UPPER, 0, 0.30)
    report("Upper 0.708-0.780, all r<0.30", UPPER[0], .780, 0, .30)
    report("Upper 0.780-0.850, all r<0.30", .780, .850, 0, .30)
    report("Upper 0.850-0.920, all r<0.30", .850, .920, 0, .30)
    report("Upper r=0.10-0.175, y=.78-.85", .780, .850, .10, .175)
    report("Upper r=0.175-0.30, y=.78-.85", .780, .850, .175, .30)
    report("Upper r=0.10-0.175, y=.85-.92", .850, .920, .10, .175)
    report("Upper r=0.175-0.30, y=.85-.92", .850, .920, .175, .30)
    report("Upper tangential-normal faces", *UPPER, .10, .30, True)
    report("Lower tangential-normal faces", *LOWER, .10, .30, True)
    print(
        "All classifications use patch-face phase fractions. "
        "The normal-orientation subset is only a blade-face proxy; "
        "confirm actual mesh surface topology before interpreting wetting."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=DEFAULT_CASE)
    parser.add_argument("--time", default="9.9")
    parser.add_argument("--convert", action="store_true")
    args = parser.parse_args()
    if not check_scratch(args.case, args.time):
        return
    if args.convert and not convert(args.case, args.time):
        return
    patch = args.case / "VTK/propeller" / f"propeller_{args.time}.vtk"
    evaluate(patch)


if __name__ == "__main__":
    main()
