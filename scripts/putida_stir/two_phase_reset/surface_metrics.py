#!/usr/bin/env python3
"""Shape of the free surface in a two-phase tank run.

For every vertical column of cells the script finds

* the surface: the highest point where the liquid fraction crosses 0.5;
* the wetted height: the highest cell with a liquid fraction above
  ``--wet`` (default 0.05), which is where a thin film or spray shows up;

and summarises them over the tank and over the band next to the wall. A
surface that sits a few centimetres higher at the wall is a normal stirred
tank; a wetted height far above the surface is liquid climbing the wall.

The fields must be ASCII and reconstructed::

    postProcess -func writeCellCentres -func writeCellVolumes -time 0
    reconstructPar -fields '(alpha.liquid)'
    surface_metrics.py CASE 0 2 4 6 8 10 --surface 1.0185 --radius 0.5334

Heights are printed in centimetres relative to ``--surface``.
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


def column_metrics(centres, alpha, axis=1, ymin=None, wet=0.05, decimals=5):
    """Surface and wetted height of every vertical column of cells.

    Returns a dict of arrays, one entry per column: the two horizontal
    coordinates ``a`` and ``b``, the radius ``r``, ``surface`` (NaN where
    the column holds no cell with a fraction of 0.5 or more), ``wetted``
    (NaN where the column is dry) and ``cells``.
    """
    centres = np.asarray(centres, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    height = centres[:, axis]
    keep = np.ones(height.size, dtype=bool) if ymin is None else height >= ymin
    height, alpha = height[keep], alpha[keep]
    horizontal = np.delete(centres[keep], axis, axis=1)

    keys, column = np.unique(
        np.round(horizontal, decimals), axis=0, return_inverse=True
    )
    column = column.ravel()
    order = np.lexsort((height, column))
    height, alpha, column = height[order], alpha[order], column[order]
    starts = np.searchsorted(column, np.arange(keys.shape[0] + 1))

    surface = np.full(keys.shape[0], np.nan)
    wetted = np.full(keys.shape[0], np.nan)
    for i in range(keys.shape[0]):
        y = height[starts[i] : starts[i + 1]]
        a = alpha[starts[i] : starts[i + 1]]
        full = np.nonzero(a >= 0.5)[0]
        if full.size:
            j = full[-1]
            if j == y.size - 1:
                surface[i] = y[j]
            else:
                surface[i] = y[j] + (a[j] - 0.5) / (a[j] - a[j + 1]) * (
                    y[j + 1] - y[j]
                )
        damp = np.nonzero(a >= wet)[0]
        if damp.size:
            wetted[i] = y[damp[-1]]

    return {
        "a": keys[:, 0],
        "b": keys[:, 1],
        "r": np.hypot(keys[:, 0], keys[:, 1]),
        "surface": surface,
        "wetted": wetted,
        "cells": np.diff(starts),
    }


def summarise(columns, radius, wall_band=0.1):
    """Reduce the column values to the numbers that describe the surface."""
    r, surface, wetted = columns["r"], columns["surface"], columns["wetted"]
    wall = r >= (1.0 - wall_band) * radius
    core = ~wall

    def stat(func, values):
        values = values[~np.isnan(values)]
        return float(func(values)) if values.size else float("nan")

    return {
        "columns": int(r.size),
        "dry_columns": int(np.isnan(surface).sum()),
        "surface_mean": stat(np.mean, surface),
        "surface_min": stat(np.min, surface),
        "surface_max": stat(np.max, surface),
        "core_surface_mean": stat(np.mean, surface[core]),
        "wall_surface_mean": stat(np.mean, surface[wall]),
        "wall_surface_max": stat(np.max, surface[wall]),
        "wall_wetted_mean": stat(np.mean, wetted[wall]),
        "wall_wetted_max": stat(np.max, wetted[wall]),
        "wetted_max": stat(np.max, wetted),
    }


def band_table(columns, radius, edges=(0.0, 0.25, 0.5, 0.75, 0.9, 1.0)):
    """Surface and wetted height by radial band, as a list of rows."""
    r = columns["r"] / radius
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        last = hi == edges[-1]
        inside = (r >= lo) & ((r <= hi + 1e-9) if last else (r < hi))
        surface = columns["surface"][inside]
        wetted = columns["wetted"][inside]
        surface = surface[~np.isnan(surface)]
        wetted = wetted[~np.isnan(wetted)]
        rows.append(
            {
                "lo": lo,
                "hi": hi,
                "columns": int(inside.sum()),
                "surface_mean": surface.mean() if surface.size else np.nan,
                "surface_min": surface.min() if surface.size else np.nan,
                "surface_max": surface.max() if surface.size else np.nan,
                "wetted_max": wetted.max() if wetted.size else np.nan,
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
    parser.add_argument(
        "--wall-band",
        type=float,
        default=0.1,
        help="fraction of the radius counted as the wall band",
    )
    parser.add_argument(
        "--csv", help="write the columns of the last time here"
    )
    args = parser.parse_args(argv)

    axis = AXES[args.axis]
    ymin = args.surface - 0.4 if args.ymin is None else args.ymin
    case = Path(args.case)

    try:
        centre_file = find_geometry(case, "C", args.times)
        if centre_file is None:
            raise FieldError(
                "cell centres not found; run "
                "postProcess -func writeCellCentres -time 0"
            )
        centres = read_field(centre_file)
        volume_file = find_geometry(case, "V", args.times)
        volumes = read_field(volume_file) if volume_file else None

        to_cm = 100.0
        print(
            "time[s]  liquid[L]   surface: mean    min    max |"
            "  wall: surface mean    max  wetted mean    max   dry cols"
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
            columns = column_metrics(centres, alpha, axis, ymin, args.wet)
            s = summarise(columns, args.radius, args.wall_band)
            litres = (
                f"{1000.0 * float(np.sum(alpha * volumes)):9.3f}"
                if volumes is not None
                else "        -"
            )

            def rel(key):
                return (s[key] - args.surface) * to_cm

            print(
                f"{time:>7}  {litres}   "
                f"{rel('surface_mean'):13.2f} {rel('surface_min'):6.2f} "
                f"{rel('surface_max'):6.2f} |"
                f"  {rel('wall_surface_mean'):18.2f} "
                f"{rel('wall_surface_max'):6.2f} "
                f"{rel('wall_wetted_mean'):12.2f} {rel('wall_wetted_max'):6.2f}"
                f"   {s['dry_columns']:8d}"
            )
    except FieldError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print()
    print(
        f"Radial bands at time {args.times[-1]} (cm above the initial surface):"
    )
    print("  r/R          columns  surface mean     min     max  wetted max")
    for row in band_table(columns, args.radius):
        print(
            f"  {row['lo']:.2f}-{row['hi']:.2f}  {row['columns']:9d}"
            f"  {(row['surface_mean'] - args.surface) * to_cm:12.2f}"
            f" {(row['surface_min'] - args.surface) * to_cm:7.2f}"
            f" {(row['surface_max'] - args.surface) * to_cm:7.2f}"
            f" {(row['wetted_max'] - args.surface) * to_cm:11.2f}"
        )

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["a", "b", "r", "surface", "wetted", "cells"])
            for row in zip(
                columns["a"],
                columns["b"],
                columns["r"],
                columns["surface"],
                columns["wetted"],
                columns["cells"],
            ):
                writer.writerow([f"{value:.6g}" for value in row])
        print(f"\nColumns written to {args.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
