#!/usr/bin/env python3
"""Read a putidaFoam log and report whether the liquid inventory is conserved.

    check_run.py LOG [--volume M3] [--tol 1e-3] [--csv FILE]

The solver prints, for every transported phase, its volume-averaged fraction
after each phase-fraction solve. Averaged liquid fraction times the tank
volume is the liquid inventory, so the log alone shows any gain or loss.
Without a reference phase the solver also prints the phase sum before it
rescales the fractions; a sum away from one is volume being created or
removed by that rescaling.

Exit status is 0 when the relative change in liquid inventory over the log is
within the tolerance, 1 when it is not, and 2 when the log has no usable data.
"""

import argparse
import csv
import re
import sys

NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)"
TIME = re.compile(rf"^Time = {NUM}s?\s*$")
FRACTION = re.compile(rf"^(\S+) fraction, min, max = {NUM} {NUM} {NUM}")
PHASE_SUM = re.compile(
    rf"^Phase-sum volume fraction, min, max = {NUM} {NUM} {NUM}"
)
DELTA_T = re.compile(rf"^deltaT = {NUM}")


def read_log(lines, liquid="liquid"):
    """One record per time step, holding the last values printed in it."""
    records, time, step = [], None, {}

    def flush():
        # The time step size is printed before its time line, so a step
        # without phase-fraction data yet is left open.
        if "fraction" in step:
            if time is not None:
                records.append({"time": time, **step})
            step.clear()

    for line in lines:
        m = TIME.match(line)
        if m:
            # putidaFoam prints the time twice per step; only a step that
            # already holds data is closed by the next time line.
            flush()
            time = float(m.group(1))
            continue
        m = FRACTION.match(line)
        if m and m.group(1) == liquid:
            step["fraction"] = float(m.group(2))
            step["alpha_min"] = float(m.group(3))
            step["alpha_max"] = float(m.group(4))
            continue
        m = PHASE_SUM.match(line)
        if m:
            low, high = float(m.group(2)), float(m.group(3))
            step["sum_min"] = min(low, step.get("sum_min", low))
            step["sum_max"] = max(high, step.get("sum_max", high))
            continue
        m = DELTA_T.match(line)
        if m:
            step["deltaT"] = float(m.group(1))
    flush()
    return records


def summarise(records, volume=None):
    first, last = records[0], records[-1]
    f0 = first["fraction"]
    out = {
        "steps": len(records),
        "t_start": first["time"],
        "t_end": last["time"],
        "fraction_start": f0,
        "fraction_end": last["fraction"],
        "relative_change": (last["fraction"] - f0) / f0,
        "largest_relative_excursion": max(
            abs(r["fraction"] - f0) / f0 for r in records
        ),
        "alpha_min": min(r["alpha_min"] for r in records),
        "alpha_max": max(r["alpha_max"] for r in records),
        "rescaled": any("sum_min" in r for r in records),
    }
    if out["rescaled"]:
        out["sum_min"] = min(r["sum_min"] for r in records if "sum_min" in r)
        out["sum_max"] = max(r["sum_max"] for r in records if "sum_max" in r)
    steps = [r["deltaT"] for r in records if "deltaT" in r]
    if steps:
        out["deltaT_min"], out["deltaT_max"] = min(steps), max(steps)
    if volume is not None:
        litres = (last["fraction"] - f0) * volume * 1000
        out["change_L"] = litres
        span = last["time"] - first["time"]
        out["rate_L_per_s"] = litres / span if span > 0 else float("nan")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("log")
    p.add_argument("--liquid", default="liquid", help="liquid phase name")
    p.add_argument(
        "--volume",
        type=float,
        help="total mesh volume in m3, to report litres",
    )
    p.add_argument(
        "--tol",
        type=float,
        default=1e-3,
        help="allowed relative inventory change",
    )
    p.add_argument("--csv", help="write the per-step history to this file")
    opts = p.parse_args(argv)

    with open(opts.log, errors="replace") as handle:
        records = read_log(handle, opts.liquid)
    if not records:
        print(
            f"check_run: no '{opts.liquid} fraction' lines in {opts.log}",
            file=sys.stderr,
        )
        return 2

    s = summarise(records, opts.volume)
    print(
        f"steps                       {s['steps']}  ({s['t_start']:g} s to {s['t_end']:g} s)"
    )
    print(
        f"mean liquid fraction        {s['fraction_start']:.6g} -> {s['fraction_end']:.6g}"
    )
    print(
        f"relative inventory change   {s['relative_change']:+.3e}  (largest {s['largest_relative_excursion']:.3e})"
    )
    if "change_L" in s:
        print(
            f"liquid volume change        {s['change_L']:+.4g} L  ({s['rate_L_per_s']:+.4g} L/s)"
        )
    print(
        f"liquid fraction bounds      {s['alpha_min']:g} to {s['alpha_max']:g}"
    )
    if s["rescaled"]:
        print(
            f"phase sum before rescaling  {s['sum_min']:g} to {s['sum_max']:g}  (no reference phase)"
        )
    else:
        print(
            "phase sum                   not printed: a reference phase is in use"
        )
    if "deltaT_min" in s:
        print(
            f"time step                   {s['deltaT_min']:g} to {s['deltaT_max']:g} s"
        )
    print("note: the log prints fractions to about six digits, so changes")
    print("      below about 2e-6 of the inventory are not resolved.")

    if opts.csv:
        names = [
            "time",
            "fraction",
            "alpha_min",
            "alpha_max",
            "sum_min",
            "sum_max",
            "deltaT",
        ]
        with open(opts.csv, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=names, restval="")
            writer.writeheader()
            writer.writerows(records)

    ok = abs(s["relative_change"]) <= opts.tol
    print(
        f"liquid inventory: {'PASS' if ok else 'FAIL'} at tolerance {opts.tol:g}"
    )
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
