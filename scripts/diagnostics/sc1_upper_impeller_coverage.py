#!/usr/bin/env python3
"""Assess upper-impeller liquid coverage from existing SC1 scratch VTK exports.

Reads only scratch volume VTKs and optional propeller boundary VTKs.
Reports cell-volume-weighted liquid fractions in specified axial/radial bands
and liquid-fraction-weighted propeller surface areas (if available).
Does not run OpenFOAM or write/change any file.

Run:
  /scratch/nwimer/.conda-envs/of-post/bin/python \
    scripts/diagnostics/sc1_upper_impeller_coverage.py
"""

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv


DEFAULT_ROOT = Path("/scratch/nwimer/sc1_surface_evolution_20261103")
# SC1 MEML's active upper MRF cellZone axial extent from topoSetDict.mrfFromYaml.
UPPER_MRF_LO = 0.70802753
UPPER_MRF_HI = 0.93662247

# These are analysis sub-bands, not confirmed locations of the physical blades.
HEIGHT_BANDS = (
    (0.70, 0.78),
    (0.78, 0.85),
    (0.85, 0.92),
    (UPPER_MRF_LO, UPPER_MRF_HI),
)
RADIAL_BANDS = (
    (0.00, 0.05),
    (0.05, 0.10),
    (0.10, 0.175),
    (0.175, 0.30),
)


def liquid_cell_field(mesh: pv.DataSet) -> np.ndarray | None:
    if "alpha.liquid" in mesh.cell_data:
        return np.asarray(mesh.cell_data["alpha.liquid"], dtype=float)
    if "alpha.liquid" in mesh.point_data:
        converted = mesh.point_data_to_cell_data(pass_point_data=False)
        return np.asarray(converted.cell_data["alpha.liquid"], dtype=float)
    return None


def report_volume(vtk: Path) -> None:
    if not vtk.is_file():
        print(f"Volume VTK not found: {vtk}")
        return

    mesh = pv.read(vtk)
    al = liquid_cell_field(mesh)
    if al is None:
        print(f"Volume VTK missing alpha.liquid: {vtk}")
        return

    centers = mesh.cell_centers().points
    y = centers[:, 1]
    r = np.hypot(centers[:, 0], centers[:, 2])
    vol = np.asarray(
        mesh.compute_cell_sizes(length=False, area=False, volume=True)
        .cell_data["Volume"], dtype=float
    )
    good = np.isfinite(al) & np.isfinite(vol) & (vol > 0)

    print("Volume-weighted liquid fraction within upper-impeller region:")
    print("Height m       Radius m      Volume L    mean alphaL   alphaL>=0.5 vol %")
    for y0, y1 in HEIGHT_BANDS:
        for r0, r1 in RADIAL_BANDS:
            mask = (
                good & (y >= y0) & (y < y1) & (r >= r0) & (r < r1)
            )
            v = np.sum(vol[mask])
            if v <= 0:
                print(f"{y0:5.3f}-{y1:5.3f}   {r0:5.3f}-{r1:5.3f}   empty")
                continue
            mean_al = float(np.sum(al[mask] * vol[mask]) / v)
            liquid_dominant = (
                100.0 * float(np.sum(vol[mask & (al >= 0.5)])) / v
            )
            print(
                f"{y0:5.3f}-{y1:5.3f}   {r0:5.3f}-{r1:5.3f}"
                f"   {1000*v:9.3f}   {mean_al:11.3f}"
                f"   {liquid_dominant:15.1f}"
            )


def report_patch(vtk: Path) -> None:
    if not vtk.is_file():
        print(f"Propeller patch VTK not found: {vtk}")
        return

    mesh = pv.read(vtk)
    print(
        f"Propeller patch: {mesh.n_cells} faces; "
        f"y bounds [{mesh.bounds[2]:.3f}, {mesh.bounds[3]:.3f}] m"
    )
    al = liquid_cell_field(mesh)
    if al is None:
        print(
            "Propeller patch has no alpha.liquid field; cell_data:",
            list(mesh.cell_data.keys()),
            "; point_data:", list(mesh.point_data.keys()),
        )
        return

    face_area = np.asarray(
        mesh.compute_cell_sizes(length=False, area=True, volume=False)
        .cell_data["Area"], dtype=float
    )
    centers = mesh.cell_centers().points
    y = centers[:, 1]
    r = np.hypot(centers[:, 0], centers[:, 2])
    valid = np.isfinite(al) & np.isfinite(face_area) & (face_area > 0)

    print(
        "Upper MRF height band: "
        f"{UPPER_MRF_LO:.3f}-{UPPER_MRF_HI:.3f} m"
    )
    print("Propeller face radius m   Area m2   area-weighted alphaL")
    for r0, r1 in RADIAL_BANDS:
        mask = (
            valid
            & (y >= UPPER_MRF_LO) & (y < UPPER_MRF_HI)
            & (r >= r0) & (r < r1)
        )
        area = float(np.sum(face_area[mask]))
        if area <= 0:
            print(f"{r0:5.3f}-{r1:5.3f}                no faces")
            continue
        mean_al = float(np.sum(al[mask] * face_area[mask]) / area)
        print(f"{r0:5.3f}-{r1:5.3f}            {area:8.5f}   {mean_al:10.3f}")

    print(
        "Patch alpha is a liquid-area proxy, not a resolved "
        "wetting/contact-line measurement."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--times", nargs="+", default=["9.9", "20"])
    parser.add_argument("--cases", nargs="+", default=["MEML_BASE", "EMM"])
    args = parser.parse_args()

    for name in args.cases:
        vtk_dir = args.root / name / "VTK"
        for time in args.times:
            print(f"\n=== {name} at t={time} s ===")
            report_volume(vtk_dir / f"{name}_{time}.vtk")
            report_patch(
                vtk_dir / "propeller" / f"propeller_{time}.vtk"
            )


if __name__ == "__main__":
    main()
