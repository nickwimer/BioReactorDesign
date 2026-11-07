#!/usr/bin/env python3
"""Read-only SC1 9.9/20 s liquid-interface evolution from root time fields.

Builds a scratch-only OpenFOAM mirror of original EMM and MEML cases, converts
only alpha.liquid with foamToVTK, and computes volume-weighted radial profiles.
Source case data (including checkpoints) are never written or modified.
Requires foamToVTK in PATH and numpy, pyvista in the Python environment.

Examples:
  python3 scripts/diagnostics/sc1_surface_evolution.py
  python3 scripts/diagnostics/sc1_surface_evolution.py --analyze-only
"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pyvista as pv


BASE = Path(
    "/projects/microenviro/final_runs/reacting_Julia_08282026/"
    "meml/putida_sc1_final"
)
EMM = Path(
    "/projects/microenviro/final_runs/reacting_Julia_08282026/"
    "emm/putida_sc1"
)
DEFAULT_OUT = Path("/scratch/nwimer/sc1_surface_evolution_20261103")

RADIAL_EDGES = (0.0, 0.05, 0.10, 0.20, 0.35, 0.45, 0.53)
HEIGHT_EDGES = np.arange(0.30, 1.551, 0.025)
HEIGHT_CENTERS = 0.5 * (HEIGHT_EDGES[:-1] + HEIGHT_EDGES[1:])


def link_children(source: Path, destination: Path) -> None:
    """Make real scratch dirs whose entries link to input files/dirs."""
    destination.mkdir(parents=True, exist_ok=True)
    for source_entry in source.iterdir():
        dest_entry = destination / source_entry.name
        target = source_entry.resolve()
        if dest_entry.is_symlink():
            if dest_entry.resolve() != target:
                raise RuntimeError(
                    f"Existing link points elsewhere: {dest_entry}"
                )
        elif dest_entry.exists():
            raise RuntimeError(
                f"Refusing to replace existing scratch entry: {dest_entry}"
            )
        else:
            dest_entry.symlink_to(target)


def mirror_case(source: Path, destination: Path, times: list[str]) -> None:
    if not (source / "constant/polyMesh/points").is_file():
        raise FileNotFoundError(f"Source mesh not found: {source}")
    if not (source / "system/controlDict").is_file():
        raise FileNotFoundError(f"Source controlDict not found: {source}")

    src_real = source.resolve()
    out_real = destination.resolve()
    if src_real == out_real or src_real in out_real.parents or out_real in src_real.parents:
        raise RuntimeError("Input and scratch paths must not overlap.")

    link_children(source / "constant", destination / "constant")
    link_children(source / "system", destination / "system")
    for time in times:
        source_time = source / time
        if not (source_time / "alpha.liquid").is_file():
            raise FileNotFoundError(
                f"Missing reconstructed alpha.liquid: {source_time}"
            )
        link_children(source_time, destination / time)


def vtk_for_time(vtk_dir: Path, time: str) -> Path | None:
    """Select the main volume mesh, never patch meshes in VTK subdirectories."""
    expected = vtk_dir / f"{vtk_dir.parent.name}_{time}.vtk"
    if expected.is_file():
        return expected

    candidates = sorted(vtk_dir.glob(f"*_{time}.vtk"))
    if len(candidates) == 1:
        return candidates[0]

    print(f"  Cannot uniquely identify volume VTK for time {time!r}.")
    print("  Root-level VTK candidates:", [str(f) for f in candidates])
    return None


def crossing(profile: np.ndarray, threshold: float = 0.5) -> str:
    ok = np.isfinite(profile)
    indices = np.flatnonzero(
        ok[:-1] & ok[1:]
        & (profile[:-1] >= threshold)
        & (profile[1:] < threshold)
    )
    if not indices.size:
        return "not found"
    i = indices[-1]
    denom = profile[i + 1] - profile[i]
    if abs(denom) < 1e-14:
        return "not found"
    portion = (threshold - profile[i]) / denom
    height = HEIGHT_CENTERS[i] + portion * (
        HEIGHT_CENTERS[i + 1] - HEIGHT_CENTERS[i]
    )
    return f"{height:.3f}{'*' if indices.size > 1 else ''}"


def evaluate_vtk(label: str, time: str, vtk: Path) -> None:
    mesh = pv.read(str(vtk))
    if "alpha.liquid" not in mesh.cell_data:
        raise KeyError(f"alpha.liquid not found in {vtk}")

    al = np.asarray(mesh.cell_data["alpha.liquid"], dtype=float)
    vol = np.asarray(
        mesh.compute_cell_sizes(
            length=False, area=False, volume=True
        ).cell_data["Volume"],
        dtype=float,
    )
    centers = mesh.cell_centers().points
    y = centers[:, 1]
    radius = np.hypot(centers[:, 0], centers[:, 2])
    good = np.isfinite(al) & np.isfinite(vol) & (vol > 0)

    total = 1000 * np.sum(al[good] * vol[good])
    above = good & (y >= 1.02)
    above_vol = 1000 * np.sum(al[above] * vol[above])

    print(f"\n=== {label} at {time} s ===")
    print(f"Liquid inventory: {total:.3f} L")
    print(f"Liquid above y=1.02 m: {above_vol:.3f} L")
    print("Radial band (m)   alphaL=0.5 (m)   Liquid above 1.02 m (L)")
    for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:]):
        annulus = good & (radius >= lo) & (radius < hi)
        upper = annulus & (y >= 1.02)
        elevated = 1000 * np.sum(al[upper] * vol[upper])

        within = annulus & (y >= HEIGHT_EDGES[0]) & (y < HEIGHT_EDGES[-1])
        bins = np.searchsorted(
            HEIGHT_EDGES, y[within], side="right"
        ) - 1
        volumes = np.bincount(
            bins, weights=vol[within], minlength=len(HEIGHT_CENTERS)
        )
        weighted = np.bincount(
            bins, weights=al[within] * vol[within],
            minlength=len(HEIGHT_CENTERS)
        )
        profile = np.full(len(HEIGHT_CENTERS), np.nan)
        np.divide(weighted, volumes, out=profile, where=volumes > 0)
        print(
            f"{lo:.2f}-{hi:.2f}            "
            f"{crossing(profile):>10s}       {elevated:10.3f}"
        )
    print("* indicates multiple annularly averaged crossings")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meml", type=Path, default=BASE)
    parser.add_argument("--emm", type=Path, default=EMM)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--times", nargs="+", default=["9.9", "20"],
        help="Reconstructed root checkpoint times to compare",
    )
    parser.add_argument(
        "--analyze-only", action="store_true",
        help="Read previously converted scratch VTK files; skip conversion",
    )
    args = parser.parse_args()

    if not args.analyze_only and not shutil.which("foamToVTK"):
        print("foamToVTK not found: load the RHEL9 OpenFOAM environment.")
        return

    out = args.out.resolve()
    # Production lives under /projects; the output must live in user scratch.
    if not str(out).startswith(("/scratch/", "/kfs3/scratch/")):
        print(f"Refusing to write outside scratch: {out}")
        return

    for name, source in (("MEML_BASE", args.meml), ("EMM", args.emm)):
        destination = out / name
        if args.analyze_only:
            print(f"\nAnalyzing existing {name} VTK: {destination / 'VTK'}", flush=True)
        else:
            print(f"\nPreparing {name}: {source} -> {destination}", flush=True)
            try:
                mirror_case(source, destination, args.times)
            except (FileNotFoundError, RuntimeError, OSError) as exc:
                print(f"  Could not prepare {name}: {exc}")
                continue

            cmd = [
                "foamToVTK",
                "-case", str(destination),
                "-time", ",".join(args.times),
                "-fields", "(alpha.liquid)",
                "-useTimeName",
            ]
            result = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, check=False,
            )
            log = destination / "foamToVTK.surface.log"
            log.write_text(result.stdout)
            print(f"  foamToVTK return code: {result.returncode}; log: {log}")
            if result.returncode != 0:
                print("\n".join(result.stdout.splitlines()[-25:]))
                continue

        for time in args.times:
            vtk = vtk_for_time(destination / "VTK", time)
            if vtk is None:
                continue
            try:
                evaluate_vtk(name, time, vtk)
            except (KeyError, ValueError, OSError) as exc:
                print(f"  Could not analyze {vtk}: {exc}")

    print(f"\nAll VTK and converter logs, if produced, are under: {out}")
    print("All original production case files were read-only.")


if __name__ == "__main__":
    main()
