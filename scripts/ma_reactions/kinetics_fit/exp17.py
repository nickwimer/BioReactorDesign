"""Exp 17 (DO-controlled P. putida batch with off-gas analysis).

Two parts:

- ``extract`` reads the three Exp 17 workbooks and writes small CSV files.
  It needs pandas and openpyxl. Run it once::

      python exp17.py extract --xlsx-dir "<Exp17 folder>" --out-dir local/exp17

- everything else reads those CSVs with numpy only.

Sign convention: uptake and production rates are positive.
"""

import argparse
from pathlib import Path

import numpy as np

REACTORS = (49, 50, 52, 54)
T0 = "2024-11-06 08:30:00"  # first sample, from the gas workbook
MW_GLUCOSE = 180.156
MW_MUCONIC_ACID = 142.11

# Column positions in the "<reactor> Gas Analysis" sheets.
COL_HOURS = 97
COL_CPR = 98
COL_CPR_NOAC = 100
COL_OUR = 103
COL_OUR_NOAC = 105
COL_O2_IN = 74
COL_O2_OUT = 76


def _trapz(y, x):
    f = getattr(np, "trapezoid", None) or np.trapz
    return float(f(y, x))


# --------------------------------------------------------------- extraction
def extract(xlsx_dir, out_dir):
    import pandas as pd

    xlsx_dir = Path(xlsx_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    books = [p for p in xlsx_dir.glob("*.xlsx") if not p.name.startswith("~")]
    gas = next(p for p in books if "Gas" in p.name)
    par = next(p for p in books if "Parameters" in p.name)
    raw = next(p for p in books if p not in (gas, par))
    t0 = pd.Timestamp(T0)

    s = pd.read_excel(raw, "Raw Data Exp 17")
    s.columns = [str(c).strip() for c in s.columns]
    s = s.dropna(subset=["time"])
    samples = pd.DataFrame(
        {
            "reactor": s["name"].astype(float).astype(int),
            "t_h": s["time"].astype(float),
            "X_kg_m3": s["CDW"].astype(float),
            "G_mol_m3": s["glucose"].astype(float) / MW_GLUCOSE * 1000.0,
            "M_mol_m3": s["Muconic Acid"].astype(float)
            / MW_MUCONIC_ACID
            * 1000.0,
            "V_L": s["Reactor Volume"].astype(float) / 1000.0,
            "base_mL": s["Base Addition"].astype(float),
        }
    )
    samples.to_csv(out_dir / "samples.csv", index=False)
    print("samples:", len(samples), "rows")

    g = pd.ExcelFile(gas)
    for r in REACTORS:
        d = g.parse(f"{r} Gas Analysis", header=None).iloc[2:]
        cols = {
            "t_h": COL_HOURS,
            "our_mmol_h": COL_OUR,
            "our_noac_mmol_h": COL_OUR_NOAC,
            "cpr_mmol_h": COL_CPR,
            "cpr_noac_mmol_h": COL_CPR_NOAC,
            "o2_in_frac": COL_O2_IN,
            "o2_out_frac": COL_O2_OUT,
        }
        out = pd.DataFrame(
            {k: pd.to_numeric(d[c], errors="coerce") for k, c in cols.items()}
        ).dropna(subset=["t_h", "our_noac_mmol_h"])
        out.to_csv(out_dir / f"offgas_{r}.csv", index=False)
        print(f"offgas {r}: {len(out)} rows to {out.t_h.max():.2f} h")

    p = pd.ExcelFile(par)
    for sheet in p.sheet_names:
        digits = "".join(ch for ch in sheet[-3:] if ch.isdigit())
        if not digits or int(digits) not in REACTORS:
            continue
        d = p.parse(sheet, header=0).iloc[:, :7]
        d.columns = ["stamp", "age", "air", "pH", "pO2", "stir", "temp"]
        stamp = pd.to_datetime(
            d["stamp"], format="%m/%d/%Y %I:%M:%S %p", errors="coerce"
        )
        d["t_h"] = (stamp - t0).dt.total_seconds() / 3600.0
        d = d.dropna(subset=["t_h"]).sort_values("t_h")
        num = d[["air", "pH", "pO2", "stir", "temp"]].apply(
            pd.to_numeric, errors="coerce"
        )
        d[num.columns] = num.ffill()
        d = d[(d.t_h >= -0.5) & (d.t_h <= 37.0)]
        binned = (
            d.assign(bin=(d.t_h * 12).round() / 12)
            .groupby("bin")[["air", "pH", "pO2", "stir", "temp"]]
            .median()
            .reset_index()
            .rename(
                columns={
                    "bin": "t_h",
                    "air": "air_ccm",
                    "pO2": "do_pct",
                    "stir": "stir_rpm",
                    "temp": "temp_C",
                }
            )
        )
        binned.to_csv(out_dir / f"online_{int(digits)}.csv", index=False)
        print(f"online {digits}: {len(binned)} rows (5 min bins)")


# ------------------------------------------------------------------ loading
def _read(path):
    return np.genfromtxt(path, delimiter=",", names=True)


def samples(local_dir, reactor):
    s = _read(Path(local_dir) / "samples.csv")
    return s[s["reactor"] == reactor]


def offgas(local_dir, reactor, gap=0.05):
    """Off-gas rates, one row per measurement burst.

    The mass spectrometer cycles between reactors and takes a short burst
    on each. Bursts are reduced to their median.
    """
    d = _read(Path(local_dir) / f"offgas_{reactor}.csv")
    start = np.flatnonzero(np.diff(d["t_h"], prepend=-1.0) > gap)
    out = {k: [] for k in d.dtype.names}
    for chunk in np.split(np.arange(len(d)), start[1:]):
        for k in out:
            out[k].append(np.nanmedian(d[k][chunk]))
    return {k: np.array(v) for k, v in out.items()}


def online(local_dir, reactor):
    d = _read(Path(local_dir) / f"online_{reactor}.csv")
    return {k: d[k] for k in d.dtype.names}


def volumetric_our(local_dir, reactor, which="our_noac_mmol_h", zero=False):
    """(t, OUR in mol/m3/h). ``zero`` shifts the trace so it starts at 0.

    The raw trace starts slightly negative (inlet and outlet O2 readings do
    not agree exactly), so the true uptake is a little higher than
    reported. Shifting to zero at the start is a lower bound on that
    correction.
    """
    s = samples(local_dir, reactor)
    g = offgas(local_dir, reactor)
    v = np.interp(g["t_h"], s["t_h"], s["V_L"])
    rate = g[which].copy()
    if zero:
        rate = rate - np.nanmedian(rate[g["t_h"] < 0.5])
    return g["t_h"], rate / v


def hourly(t, y, t_end, width=1.0):
    """Mean of y in consecutive windows; returns (centres, means)."""
    edges = np.arange(0.0, t_end + 1e-9, width)
    tc, yc = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (t >= a) & (t < b) & np.isfinite(y)
        if m.sum() >= 2:
            tc.append(0.5 * (a + b))
            yc.append(float(np.mean(y[m])))
    return np.array(tc), np.array(yc)


def load_reactor(local_dir, reactor, o2sat):
    """Dataset for the fit: samples, prescribed DO and measured OUR."""
    from kinetics0d import Dataset

    s = samples(local_dir, reactor)
    o = online(local_dir, reactor)
    t_end = s["t_h"][-1]
    # Hourly means: the 5 min series only makes the integrator work harder.
    td, do = hourly(o["t_h"], np.clip(o["do_pct"], 0.0, None), t_end + 1.0)
    t_our, our = volumetric_our(local_dir, reactor)
    tc, oc = hourly(t_our, our, t_end)
    return Dataset(
        name=f"exp17_{reactor}",
        group="exp17",
        t=s["t_h"],
        X=s["X_kg_m3"],
        G=s["G_mol_m3"],
        M=s["M_mol_m3"],
        o2sat=o2sat,
        do_t=td,
        do_frac=do / 100.0,
        our_t=tc,
        our=oc,
    )


# ----------------------------------------------------------------- analysis
def yields(local_dir, reactor, intervals=None, zero=False):
    """Oxygen and CO2 per glucose from off-gas, by sampling interval.

    Everything is per unit liquid volume, so sample withdrawal does not
    enter. Returns a list of dicts.
    """
    s = samples(local_dir, reactor)
    g = offgas(local_dir, reactor)
    t_o, our = volumetric_our(local_dir, reactor, zero=zero)
    _, cpr = volumetric_our(local_dir, reactor, which="cpr_noac_mmol_h")
    n = len(s)
    if intervals is None:
        intervals = [(0, 4), (4, 5), (5, 6), (6, 8), (0, n - 2), (0, n - 1)]
    rows = []
    for i, j in intervals:
        ta, tb = s["t_h"][i], s["t_h"][j]
        m = (t_o >= ta) & (t_o <= tb) & np.isfinite(our)
        o2 = _trapz(our[m], t_o[m])
        co2 = _trapz(cpr[m], g["t_h"][m])
        dG = s["G_mol_m3"][i] - s["G_mol_m3"][j]
        dM = s["M_mol_m3"][j] - s["M_mol_m3"][i]
        dX = s["X_kg_m3"][j] - s["X_kg_m3"][i]
        x_mean = _trapz(
            np.interp(t_o[m], s["t_h"], s["X_kg_m3"]),
            t_o[m],
        ) / (t_o[m][-1] - t_o[m][0])
        rows.append(
            {
                "reactor": reactor,
                "t_start": float(ta),
                "t_end": float(tb),
                "glucose_used": float(dG),
                "o2_used": o2,
                "co2_made": co2,
                "y_os": o2 / dG,
                "rq": co2 / o2,
                "y_ms": float(dM / dG),
                "y_x": float(dX / dG),
                "qs": float(dG / (tb - ta) / x_mean),
                "our_spec": o2 / (tb - ta) / x_mean,
                "mu": float(
                    np.log(s["X_kg_m3"][j] / s["X_kg_m3"][i]) / (tb - ta)
                ),
            }
        )
    return rows


def kla_series(local_dir, reactor, o2sat, min_our=5.0):
    """kLa = OUR / (C* - C) during the run, with the stir rate.

    C* follows the outlet gas (well-mixed gas phase): o2sat scaled by the
    outlet-to-inlet O2 ratio. Only points with OUR above ``min_our``
    (mol/m3/h) are returned, where the driving force is well resolved.
    """
    g = offgas(local_dir, reactor)
    o = online(local_dir, reactor)
    t, our = volumetric_our(local_dir, reactor)
    do = np.interp(t, o["t_h"], o["do_pct"]) / 100.0
    rpm = np.interp(t, o["t_h"], o["stir_rpm"])
    cstar = o2sat * g["o2_out_frac"] / g["o2_in_frac"]
    drive = cstar - o2sat * do
    m = np.isfinite(our) & (our > min_our) & (drive > 0.02 * o2sat)
    return {
        "t": t[m],
        "rpm": rpm[m],
        "kla": our[m] / drive[m],
        "our": our[m],
        "do": do[m],
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    ex = sub.add_parser("extract", help="workbooks -> CSV")
    ex.add_argument("--xlsx-dir", required=True)
    ex.add_argument("--out-dir", required=True)
    yl = sub.add_parser("yields", help="print off-gas yields")
    yl.add_argument("--local-dir", required=True)
    args = ap.parse_args(argv)
    if args.cmd == "extract":
        extract(args.xlsx_dir, args.out_dir)
    else:
        for r in REACTORS:
            for row in yields(args.local_dir, r):
                print(
                    "{reactor} {t_start:5.1f}-{t_end:5.1f} h  y_os {y_os:.2f}"
                    "  RQ {rq:.2f}  qs {qs:.2f}  mu {mu:.3f}".format(**row)
                )


if __name__ == "__main__":
    main()
