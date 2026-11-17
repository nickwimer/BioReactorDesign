#!/usr/bin/env python3
"""Shape of the free surface in a two-phase tank run.

The tank is divided into vertical columns (radial rings times azimuthal
sectors). In every column the script finds

* the surface: the highest point where the liquid fraction crosses 0.5;
* the wetted height: the highest layer with a liquid fraction above
  ``--wet`` (default 0.05), which is where a thin film or spray shows up;

and summarises them over the tank and over the ring next to the wall. A
surface that is higher at the wall and lower at the shaft is a vortex; a
wetted height far above the surface is a film climbing the wall.

The fields must be ASCII and reconstructed::

    postProcess -funcs '(writeCellCentres writeCellVolumes)' -time 0
    reconstructPar -fields '(alpha.liquid)'
    surface_metrics.py CASE 0 2 4 6 8 10 --surface 1.0186 --radius 0.5334

Heights are printed in centimetres relative to ``--surface``. A column
with no liquid above ``--ymin`` is counted as dry and reported as lying
at ``--ymin`` or below.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

AXES = {"x": 0, "y": 1, "z": 2}


class FieldError(Exception):
    """A field file is missing or cannot be read."""


def read_field(path, n_cells=None):
    """Return the internal values of an ASCII OpenFOAM field.

    Scalars come back with shape (n,), vectors with shape (n, 3). A uniform
    field needs ``n_cells`` to be expanded.
    """
    path = Path(path)
    if not path.is_file():
        raise FieldError(f"field not found: {path}")
    raw = path.read_bytes()
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        raise FieldError(f"{path} is not an ASCII field") from exc
    if re.search(r"format\s+binary\s*;", text[:1000]):
        raise FieldError(f"{path} is written in binary format")

    match = re.search(
        r"\b(internalField|value)\s+(uniform|nonuniform)\s+", text
    )
    if not match:
        raise FieldError(f"no internalField in {path}")
    rest = text[match.end() :]

    if match.group(2) == "uniform":
        entry = rest[: rest.index(";")].replace("(", " ").replace(")", " ")
        values = np.array(entry.split(), dtype=float)
        if n_cells is None:
            raise FieldError(
                f"{path} is uniform and the cell count is unknown"
            )
        if values.size == 1:
            return np.full(n_cells, values[0])
        return np.tile(values, (n_cells, 1))

    head = re.match(r"List<(\w+)>\s*(\d+)\s*\(", rest)
    if not head:
        raise FieldError(f"cannot read the list header in {path}")
    kind, count = head.group(1), int(head.group(2))
    body = rest[head.end() :]
    stop = body.find("boundaryField")
    if stop >= 0:
        body = body[:stop]
    body = body[: body.rindex(")")]
    width = {"scalar": 1, "vector": 3}.get(kind)
    if width is None:
        raise FieldError(f"{path} holds {kind} values, not scalars or vectors")
    if width == 3:
        body = body.replace("(", " ").replace(")", " ")
    values = np.array(body.split(), dtype=float)
    if values.size != count * width:
        raise FieldError(
            f"{path}: expected {count * width} numbers, found {values.size}"
        )
    return values if width == 1 else values.reshape(count, 3)


def column_metrics(
    centres,
    alpha,
    radius,
    volumes=None,
    axis=1,
    ymin=None,
    wet=0.05,
    rings=10,
    sectors=24,
    layer=None,
):
    """Surface and wetted height in columns of rings times sectors.

    Cells are averaged, weighted by volume, into bins of ring, sector and
    horizontal layer, so the result does not depend on how the mesh is
    built. Layers follow the cell-centre heights; give ``layer`` (a
    thickness) to use uniform slabs instead.

    Returns a dict of arrays with one entry per column that holds cells:
    ``r`` (mid radius of the ring), ``theta`` (mid angle of the sector),
    ``surface`` (NaN where the column has no layer with a fraction of 0.5
    or more), ``wetted`` (NaN where the column is dry) and ``bottom``, the
    lowest layer of the column.
    """
    centres = np.asarray(centres, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    weight = (
        np.ones(alpha.size)
        if volumes is None
        else np.asarray(volumes, dtype=float)
    )
    height = centres[:, axis]
    keep = np.ones(height.size, dtype=bool) if ymin is None else height >= ymin
    height, alpha, weight = height[keep], alpha[keep], weight[keep]
    horizontal = np.delete(centres[keep], axis, axis=1)
    r = np.hypot(horizontal[:, 0], horizontal[:, 1])
    theta = np.arctan2(horizontal[:, 1], horizontal[:, 0]) % (2.0 * np.pi)

    # Horizontal layers: the distinct cell-centre heights where the mesh is
    # layered, uniform slabs of thickness ``layer`` otherwise.
    level, iy = np.unique(np.round(height, 6), return_inverse=True)
    iy = iy.ravel()
    if layer is not None or level.size > 2000:
        if layer is None:
            layer = (height.max() - height.min()) / 200.0
        base = height.min() - 0.5 * layer
        iy = np.floor((height - base) / layer).astype(int)
        level = base + layer * (np.arange(int(iy.max()) + 1) + 0.5)
    ir = np.minimum((r / radius * rings).astype(int), rings - 1)
    it = np.minimum((theta / (2.0 * np.pi) * sectors).astype(int), sectors - 1)
    n_layers = level.size

    shape = (rings, sectors, n_layers)
    flat = np.ravel_multi_index((ir, it, iy), shape)
    size = rings * sectors * n_layers
    total = np.bincount(flat, weights=weight, minlength=size).reshape(shape)
    liquid = np.bincount(flat, weights=weight * alpha, minlength=size)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = liquid.reshape(shape) / total

    out = {"r": [], "theta": [], "surface": [], "wetted": [], "bottom": []}
    for i in range(rings):
        for j in range(sectors):
            filled = total[i, j] > 0
            if not filled.any():
                continue
            y = level[filled]
            a = mean[i, j][filled]
            surface = wetted = np.nan
            full = np.nonzero(a >= 0.5)[0]
            if full.size:
                k = full[-1]
                surface = y[k]
                if k < y.size - 1:
                    surface += (
                        (a[k] - 0.5) / (a[k] - a[k + 1]) * (y[k + 1] - y[k])
                    )
            damp = np.nonzero(a >= wet)[0]
            if damp.size:
                wetted = y[damp[-1]]
            out["r"].append((i + 0.5) * radius / rings)
            out["theta"].append((j + 0.5) * 2.0 * np.pi / sectors)
            out["surface"].append(surface)
            out["wetted"].append(wetted)
            out["bottom"].append(y[0])
    return {key: np.array(value) for key, value in out.items()}


def summarise(columns, radius, rings=10, floor=None):
    """Reduce the column values to the numbers that describe the surface.

    The wall ring is the outermost ring. Dry columns count as lying at
    ``floor`` (default: the lowest layer of the column) in the means and
    minima, so a vortex that reaches below the analysed range is not lost.
    """
    r, wetted = columns["r"], columns["wetted"]
    dry = np.isnan(columns["surface"])
    low = columns["bottom"] if floor is None else np.full(r.size, floor)
    surface = np.where(dry, low, columns["surface"])
    wall = r > (1.0 - 1.0 / rings) * radius
    core = r < radius / rings * 2.0

    def stat(func, values):
        values = values[~np.isnan(values)]
        return float(func(values)) if values.size else float("nan")

    return {
        "columns": int(r.size),
        "dry_columns": int(dry.sum()),
        "surface_min": stat(np.min, surface),
        "surface_max": stat(np.max, surface),
        "core_surface_mean": stat(np.mean, surface[core]),
        "wall_surface_mean": stat(np.mean, surface[wall]),
        "wall_surface_max": stat(np.max, surface[wall]),
        "wall_wetted_mean": stat(np.mean, wetted[wall]),
        "wall_wetted_max": stat(np.max, wetted[wall]),
        "wetted_max": stat(np.max, wetted),
    }


def ring_table(columns, floor=None):
    """Surface and wetted height by ring, as a list of rows."""
    rows = []
    for radius in np.unique(columns["r"]):
        inside = columns["r"] == radius
        dry = np.isnan(columns["surface"][inside])
        low = columns["bottom"][inside] if floor is None else floor
        surface = np.where(dry, low, columns["surface"][inside])
        wetted = columns["wetted"][inside]
        wetted = wetted[~np.isnan(wetted)]
        rows.append(
            {
                "r": float(radius),
                "columns": int(inside.sum()),
                "dry": int(dry.sum()),
                "surface_mean": float(surface.mean()),
                "surface_min": float(surface.min()),
                "surface_max": float(surface.max()),
                "wetted_max": float(wetted.max()) if wetted.size else np.nan,
            }
        )
    return rows


def find_geometry(case, name, times):
    """Locate the cell-centre or cell-volume field written by postProcess."""
    for folder in ["0", "constant", *times]:
        path = Path(case) / folder / name
        if path.is_file():
            return path
    return None


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("case", help="case directory, reconstructed")
    parser.add_argument("times", nargs="+", help="time directories to read")
    parser.add_argument("--liquid", default="liquid", help="liquid phase name")
    parser.add_argument("--axis", default="y", choices=sorted(AXES))
    parser.add_argument(
        "--surface",
        type=float,
        required=True,
        help="height of the initial surface, m",
    )
    parser.add_argument(
        "--radius", type=float, required=True, help="tank radius, m"
    )
    parser.add_argument(
        "--ymin",
        type=float,
        default=None,
        help="ignore cells below this height (default: surface - 0.4 m)",
    )
    parser.add_argument("--wet", type=float, default=0.05)
    parser.add_argument("--rings", type=int, default=10)
    parser.add_argument("--sectors", type=int, default=24)
    parser.add_argument(
        "--csv", help="write the columns of the last time here"
    )
    args = parser.parse_args(argv)

    axis = AXES[args.axis]
    ymin = args.surface - 0.4 if args.ymin is None else args.ymin
    case = Path(args.case)
    cm = 100.0

    try:
        centre_file = find_geometry(case, "C", args.times)
        if centre_file is None:
            raise FieldError(
                "cell centres not found; run "
                "postProcess -funcs '(writeCellCentres writeCellVolumes)' "
                "-time 0"
            )
        centres = read_field(centre_file)
        volume_file = find_geometry(case, "V", args.times)
        volumes = read_field(volume_file) if volume_file else None
        if volumes is None:
            print("note: cell volumes not found, cells are weighted equally")

        print(
            f"Heights in cm above the initial surface ({args.surface:g} m); "
            f"columns are {args.rings} rings x {args.sectors} sectors, "
            f"cells below {ymin:g} m ignored."
        )
        print(
            "time[s]  liquid[L] | surface: lowest  shaft  wall mean  wall max"
            " | wetted: wall mean  wall max  anywhere | dry columns"
        )
        columns = None
        for time in args.times:
            alpha = read_field(
                case / time / f"alpha.{args.liquid}", n_cells=centres.shape[0]
            )
            if alpha.shape[0] != centres.shape[0]:
                raise FieldError(
                    f"time {time}: {alpha.shape[0]} values for "
                    f"{centres.shape[0]} cells"
                )
            columns = column_metrics(
                centres,
                alpha,
                args.radius,
                volumes,
                axis,
                ymin,
                args.wet,
                args.rings,
                args.sectors,
            )
            s = summarise(columns, args.radius, args.rings)
            litres = (
                f"{1000.0 * float(np.sum(alpha * volumes)):9.3f}"
                if volumes is not None
                else "        -"
            )

            def rel(key, summary=s):
                return (summary[key] - args.surface) * cm

            print(
                f"{time:>7}  {litres} | "
                f"{rel('surface_min'):15.1f} {rel('core_surface_mean'):6.1f} "
                f"{rel('wall_surface_mean'):10.1f} "
                f"{rel('wall_surface_max'):9.1f} | "
                f"{rel('wall_wetted_mean'):17.1f} "
                f"{rel('wall_wetted_max'):9.1f} {rel('wetted_max'):9.1f} | "
                f"{s['dry_columns']:5d} of {s['columns']}"
            )
    except FieldError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print()
    print(f"By ring at time {args.times[-1]} (cm above the initial surface):")
    print("  r/R    surface mean  lowest  highest  wetted max  dry sectors")
    for row in ring_table(columns):
        print(
            f"  {row['r'] / args.radius:4.2f}"
            f"  {(row['surface_mean'] - args.surface) * cm:12.1f}"
            f" {(row['surface_min'] - args.surface) * cm:7.1f}"
            f" {(row['surface_max'] - args.surface) * cm:8.1f}"
            f" {(row['wetted_max'] - args.surface) * cm:11.1f}"
            f"  {row['dry']:5d} of {row['columns']}"
        )

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["r", "theta", "surface", "wetted", "bottom"])
            for row in zip(
                columns["r"],
                columns["theta"],
                columns["surface"],
                columns["wetted"],
                columns["bottom"],
            ):
                writer.writerow([f"{value:.6g}" for value in row])
        print(f"\nColumns written to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
