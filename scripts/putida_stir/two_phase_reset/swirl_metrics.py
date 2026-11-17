#!/usr/bin/env python3
"""Swirl of the liquid in a stirred-tank run.

Prints the mean tangential velocity of the liquid by radius and height, in
m/s and as a fraction of the impeller's solid-body speed ``omega * r``, and
the rotation rate of each height band (the rate of a rigid body with the
same angular momentum).

In a baffled tank the liquid outside the impeller region should turn at a
small fraction of the impeller speed. Liquid that turns as a body digs a
vortex: a body turning at rate W has a surface ``W^2 R^2 / 2g`` higher at
the wall than on the axis, and that number is printed for comparison with
the surface shape from ``surface_metrics.py``.

The fields must be ASCII and reconstructed::

    postProcess -funcs '(writeCellCentres writeCellVolumes)' -time 0
    reconstructPar -fields '(alpha.liquid U.liquid)'
    swirl_metrics.py CASE 2 4 --omega 17.38 --radius 0.5334 \\
        --heights 0.1 0.4 0.7 1.0
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from surface_metrics import AXES, FieldError, find_geometry, read_field

GRAVITY = 9.81


def tangential_velocity(centres, velocity, axis=1):
    """Radius and tangential velocity about the tank axis.

    The tangential velocity is positive in the direction of a positive
    rotation about the axis (right-hand rule).
    """
    centres = np.asarray(centres, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    unit = np.zeros(3)
    unit[axis] = 1.0
    radial = centres - np.outer(centres @ unit, unit)
    r = np.linalg.norm(radial, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        direction = np.cross(unit, radial) / r[:, None]
    direction[r == 0] = 0.0
    return r, np.einsum("ij,ij->i", velocity, direction)


def swirl_table(r, height, u_t, weight, radius, heights, rings=5):
    """Weighted mean tangential velocity by height band and ring.

    Returns one dict per height band with ``lo``, ``hi``, the arrays
    ``velocity`` and ``mean_r`` (one value per ring, NaN where the band
    holds no liquid) and ``rate``, the rigid-body rotation rate with the
    same angular momentum.
    """
    ring = np.minimum((r / radius * rings).astype(int), rings - 1)
    rows = []
    for lo, hi in zip(heights[:-1], heights[1:]):
        band = (height >= lo) & (height < hi)
        velocity = np.full(rings, np.nan)
        mean_r = np.full(rings, np.nan)
        for i in range(rings):
            pick = band & (ring == i)
            total = weight[pick].sum()
            if total > 0:
                velocity[i] = (weight[pick] * u_t[pick]).sum() / total
                mean_r[i] = (weight[pick] * r[pick]).sum() / total
        inertia = (weight[band] * r[band] ** 2).sum()
        rate = (
            (weight[band] * r[band] * u_t[band]).sum() / inertia
            if inertia > 0
            else np.nan
        )
        rows.append(
            {
                "lo": lo,
                "hi": hi,
                "velocity": velocity,
                "mean_r": mean_r,
                "rate": float(rate),
            }
        )
    return rows


def rise(rate, radius):
    """Wall-to-axis surface rise of a rigid body turning at ``rate``."""
    return rate**2 * radius**2 / (2.0 * GRAVITY)


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
        "--omega", type=float, required=True, help="impeller speed, rad/s"
    )
    parser.add_argument(
        "--radius", type=float, required=True, help="tank radius, m"
    )
    parser.add_argument(
        "--heights",
        type=float,
        nargs="+",
        required=True,
        help="edges of the height bands, m, lowest first",
    )
    parser.add_argument("--rings", type=int, default=5)
    args = parser.parse_args(argv)

    if len(args.heights) < 2 or sorted(args.heights) != args.heights:
        parser.error("--heights needs at least two increasing values")
    axis = AXES[args.axis]
    case = Path(args.case)
    name = args.liquid

    try:
        centre_file = find_geometry(case, "C", args.times)
        volume_file = find_geometry(case, "V", args.times)
        if centre_file is None or volume_file is None:
            raise FieldError(
                "cell centres or volumes not found; run "
                "postProcess -funcs '(writeCellCentres writeCellVolumes)' "
                "-time 0"
            )
        centres = read_field(centre_file)
        volumes = read_field(volume_file)
        height = centres[:, axis]
        n = centres.shape[0]

        edges = np.arange(args.rings + 1) / args.rings
        header = "  ".join(
            f"r/R {lo:.1f}-{hi:.1f}" for lo, hi in zip(edges[:-1], edges[1:])
        )
        for time in args.times:
            alpha = read_field(case / time / f"alpha.{name}", n_cells=n)
            velocity = read_field(case / time / f"U.{name}", n_cells=n)
            r, u_t = tangential_velocity(centres, velocity, axis)
            weight = np.clip(alpha, 0.0, 1.0) * volumes
            rows = swirl_table(
                r, height, u_t, weight, args.radius, args.heights, args.rings
            )
            print(
                f"time {time}: tangential velocity of the liquid, m/s "
                "(fraction of impeller solid-body speed)"
            )
            print(f"  height band, m   {header}   band rate/impeller")
            for row in reversed(rows):
                cells = []
                for u, mean_r in zip(row["velocity"], row["mean_r"]):
                    if np.isnan(u):
                        cells.append("     -      ")
                    else:
                        cells.append(
                            f"{u:5.2f} ({u / (args.omega * mean_r):4.2f})"
                        )
                print(
                    f"  {row['lo']:5.2f} - {row['hi']:5.2f}    "
                    + "  ".join(cells)
                    + f"   {row['rate'] / args.omega:6.2f}"
                )
            inertia = (weight * r**2).sum()
            rate = (weight * r * u_t).sum() / inertia
            print(
                f"  whole liquid: {rate:.2f} rad/s, "
                f"{rate / args.omega:.2f} of the impeller; a rigid body at "
                f"that rate has a surface {100 * rise(rate, args.radius):.0f}"
                " cm higher at the wall than on the axis"
            )
            print()
    except FieldError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
