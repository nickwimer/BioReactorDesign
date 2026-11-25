"""Run the fit matrix and write the findings report.

Everything goes to ``out/`` (git-ignored): REPORT.md, figures, tables and one
``microbeUpdateProperties`` block per candidate parameter set. Needs
local.yaml (surrogate files and the Exp 17 extract).

    python make_report.py            # full run, several minutes
    python make_report.py --quick    # fewer starts, no profiles
"""

import argparse
import datetime
import hashlib
import json
import subprocess
from dataclasses import replace
from pathlib import Path

import exp17
import fit_kinetics as fk
import kinetics0d as k
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "0dWillComp"
REF_STATE = (7.33, 0.2528, 228.8, 0.0, 0.0)  # first CFD bio update
K_O_CFD = 0.003125  # mol/m3 = 0.1 mg/L, the O2-sink floor used in the CFD
K_S_PAPER = 0.4008
RPM_BIN = 200  # rpm, width of the stir-rate bins of the kLa table
CAP_GRID = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 30.0)  # cap_n profile
MODEL_FILE = "scripted_model_biomass_h_muconate.pt"
# sha256 of MODEL_FILE in the 1000 L reacting cases on Kestrel (runs/rx_meml)
# on 2026-11-25. It is not the file of the April sc1_test case.
KESTREL_MODEL_SHA256 = (
    "c46adcaa7565fbeb1834c7732e5dff8aa39b9470b67ba5687d220ab893a10a51"
)
C = {
    "blue": "#2a78d6",
    "orange": "#eb6834",
    "aqua": "#1baf7a",
    "yellow": "#eda100",
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "surface": "#fcfcfb",
}
SERIES = [C["blue"], C["orange"], C["aqua"], C["yellow"]]


def style():
    plt.rcParams.update(
        {
            "figure.facecolor": C["surface"],
            "axes.facecolor": C["surface"],
            "savefig.facecolor": C["surface"],
            "axes.edgecolor": C["axis"],
            "axes.labelcolor": C["ink2"],
            "text.color": C["ink"],
            "xtick.color": C["muted"],
            "ytick.color": C["muted"],
            "axes.grid": True,
            "grid.color": C["grid"],
            "grid.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "lines.linewidth": 2.0,
            "lines.markersize": 6,
            "font.size": 10,
            "axes.titlesize": 10,
            "legend.frameon": False,
            "figure.dpi": 110,
        }
    )


def say(*a):
    print(*a, flush=True)


def f(x, n=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "-"
    return f"{x:.{n}g}"


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def truncate(ds, t_max):
    m = ds.t <= t_max
    kw = dict(t=ds.t[m], X=ds.X[m], G=ds.G[m], M=ds.M[m])
    if ds.do is not None:
        kw["do"] = ds.do[m]
    if ds.our_t is not None:
        mo = ds.our_t <= t_max
        kw["our_t"], kw["our"] = ds.our_t[mo], ds.our[mo]
    return replace(ds, **kw)


def group_of(name):
    return "fixed-stir" if name.startswith("fixedstir") else "Exp 17"


def score(values, datasets, sigma, model="emm", sur=None):
    """chi2 and RMSE / noise level by (experiment, quantity)."""
    pb = fk.Problem(datasets, (), dict(values), model, sur, sigma)
    blocks = pb.blocks(values)
    by = {}
    for (name, key), b in blocks.items():
        if key != "hull":
            by.setdefault((group_of(name), key), []).append(b)
    chi2 = sum(float(np.sum(np.concatenate(v) ** 2)) for v in by.values())
    rm = {
        kk: float(np.sqrt(np.mean(np.concatenate(v) ** 2)))
        for kk, v in by.items()
    }
    return chi2, rm


def rm_cells(rm):
    keys = [
        ("fixed-stir", "X"),
        ("fixed-stir", "G"),
        ("fixed-stir", "M"),
        ("fixed-stir", "DO"),
        ("Exp 17", "X"),
        ("Exp 17", "G"),
        ("Exp 17", "M"),
        ("Exp 17", "OUR"),
    ]
    return [f(rm.get(kk), 2) for kk in keys]


RM_HEAD = [
    "fs X",
    "fs G",
    "fs M",
    "fs DO",
    "E17 X",
    "E17 G",
    "E17 M",
    "E17 OUR",
]


def interval(prof, dof):
    """Parameter range with chi2 within 3.84 reduced-chi2 units of the best.

    The crossing is interpolated between grid points in log-parameter. An
    end that never crosses is reported as open.
    """
    v = np.log(np.array([p[0] for p in prof]))
    c = np.array([p[1] for p in prof])
    limit = c.min() + 3.84 * c.min() / max(dof, 1)
    best = int(np.argmin(c))

    def edge(step):
        j = best
        while 0 <= j + step < len(c) and c[j + step] <= limit:
            j += step
        if not 0 <= j + step < len(c):
            return ("<= " if step < 0 else ">= ") + f(np.exp(v[j]))
        n = j + step
        x = v[j] + (limit - c[j]) * (v[n] - v[j]) / (c[n] - c[j])
        return f(np.exp(x))

    return f(np.exp(v[best])), edge(-1), edge(1)


def vals_row(v, names):
    out = []
    for n in names:
        at_floor = n in fk.BOUNDS and v[n] <= 1.05 * fk.BOUNDS[n][0]
        out.append(f"{f(v[n], 4)} (bound)" if at_floor else f(v[n], 4))
    return out


# ------------------------------------------------------------------ figures
def fig_profiles(profiles, dof, path):
    names = list(profiles)
    fig, axes = plt.subplots(
        1, len(names), figsize=(3.0 * len(names), 2.8), squeeze=False
    )
    for ax, name in zip(axes[0], names):
        v = np.array([p[0] for p in profiles[name]])
        c = np.array([p[1] for p in profiles[name]])
        ax.plot(v, c - c.min(), color=C["blue"], marker="o")
        ax.axhline(
            3.84 * c.min() / max(dof, 1), color=C["muted"], lw=1, ls="--"
        )
        if v.max() / v.min() > 20:
            ax.set_xscale("log")
        ax.set_xlabel(name)
        ax.set_ylim(bottom=0)
    axes[0][0].set_ylabel("chi2 - best")
    fig.suptitle(
        "Misfit against each parameter, the others refitted "
        "(dashed: 95 % level)",
        fontsize=10,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_fits(datasets, models, sur, path, quantities):
    """Small multiples: one column per quantity, one row per run.

    ``models`` holds (label, values, model) or (label, values, model,
    colour); without a colour the series colours are taken in order.
    """
    n = len(datasets)
    fig, axes = plt.subplots(
        n,
        len(quantities),
        figsize=(3.1 * len(quantities), 2.1 * n),
        squeeze=False,
        sharex=True,
    )
    for i, ds in enumerate(datasets):
        t = np.linspace(ds.t[0], ds.t[-1], 200)
        for j, (key, label) in enumerate(quantities):
            ax = axes[i][j]
            for entry, col in zip(models, SERIES):
                name, values, model = entry[:3]
                col = entry[3] if len(entry) > 3 else col
                p = fk.params_of(values)
                try:
                    s = k.simulate(
                        p, ds, model, sur, values["theta"], t_dense=t
                    )
                except RuntimeError:
                    continue
                y = s["O2"] / ds.o2sat if key == "DO" else s[key]
                ax.plot(t, y, color=col, label=name if i + j == 0 else None)
            if key == "DO":
                if ds.do is not None:
                    ax.plot(ds.t, ds.do, "o", color=C["ink"])
                else:
                    ax.plot(ds.do_t, ds.do_frac, ".", color=C["ink"], ms=3)
            elif key == "OUR":
                if ds.our is not None:
                    ax.plot(ds.our_t, ds.our, ".", color=C["ink"], ms=4)
            else:
                ax.plot(
                    ds.t,
                    getattr(ds, key),
                    "o",
                    color=C["ink"],
                    label="measured" if i + j == 0 else None,
                )
            if i == 0:
                ax.set_title(label, loc="left")
            if j == 0:
                ax.set_ylabel(ds.name)
            if i == n - 1:
                ax.set_xlabel("time (h)")
            ax.set_ylim(bottom=0)
    fig.legend(loc="upper right", ncol=len(models) + 1, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path)
    plt.close(fig)


def fig_offgas(local, o2sat, path_our, path_kla):
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for r, col in zip(exp17.REACTORS, SERIES):
        t, our = exp17.volumetric_our(local, r)
        ax.plot(t, our, color=col, lw=1.5, label=f"reactor {r}")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("oxygen uptake (mol/m3/h)")
    ax.set_title("Exp 17: oxygen uptake from off-gas", loc="left")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path_our)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for r, col in zip(exp17.REACTORS, SERIES):
        s = exp17.kla_series(local, r, o2sat)
        ax.plot(s["rpm"], s["kla"], "o", color=col, ms=4, label=f"reactor {r}")
    ax.set_xlabel("stir rate (rpm)")
    ax.set_ylabel("kLa (1/h)")
    ax.set_ylim(bottom=0)
    ax.set_title(
        "Exp 17: kLa = OUR / (C* - C) during the fermentation", loc="left"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(path_kla)
    plt.close(fig)


def fig_plane(sur, rows, rays, path):
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    ax.plot(
        -sur.points[:, 0],
        -sur.points[:, 1],
        ".",
        color=C["grid"],
        ms=3,
        label="surrogate training points",
    )
    hull = np.vstack([sur.hull_vertices, sur.hull_vertices[:1]])
    ax.plot(-hull[:, 0], -hull[:, 1], color=C["axis"], lw=1)
    qs = np.linspace(0, 9.5, 50)
    for (label, y_os), col in zip(rays, SERIES):
        ax.plot(qs, y_os * qs, color=col, label=label)
    ax.plot(
        [r["qs"] for r in rows],
        [r["our_spec"] for r in rows],
        "o",
        color=C["ink"],
        label="Exp 17 measured, 9 h onwards",
    )
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 21)
    ax.set_xlabel("specific glucose uptake (mol/kg/h)")
    ax.set_ylabel("specific O2 uptake (mol/kg/h)")
    ax.set_title(
        "Where each parameter set evaluates the surrogate", loc="left"
    )
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_qs_o2(cands, path, o2_set):
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    o2 = np.linspace(0.0, 0.26, 200)
    for (name, v, _), col in cands:
        p = fk.params_of(v)
        ax.plot(o2, k.uptake(228.8, o2, p), color=col, label=name)
    ax.axvspan(0.08, 0.18, color=C["grid"], alpha=0.5, lw=0)
    ax.text(
        0.13,
        0.3,
        "fixed-stir DO range",
        ha="center",
        color=C["ink2"],
        fontsize=8,
    )
    ax.axvline(o2_set, color=C["muted"], lw=1, ls="--")
    ax.text(
        o2_set + 0.002, 0.9, "Exp 17 set point", color=C["ink2"], fontsize=8
    )
    ax.set_xlabel("dissolved O2 (mol/m3)")
    ax.set_ylabel("glucose uptake qs (mol/kg/h)")
    ax.set_ylim(bottom=0)
    ax.set_title(
        "Uptake against dissolved O2 at 229 mol/m3 glucose", loc="left"
    )
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ------------------------------------------------------------------- matrix
def run(quick=False):
    style()
    out = HERE / "out"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "candidates").mkdir(exist_ok=True)
    cfg = fk.local_config()
    local = cfg["exp17_dir"]
    allds = k.load_datasets(HERE / "datasets.yaml", local_dir=local)
    fs = [d for d in allds if d.group == "fixedstir"]
    ex = [d for d in allds if d.group == "exp17"]
    o2sat = fs[0].o2sat
    # DO set point of Exp 17, read off its DO log: fraction of
    # saturation, mol/m3 and per cent for the text.
    do_all = np.concatenate([d.do_frac for d in ex])
    do_set = round(float(np.median(do_all)), 2)
    o2_set = do_set * o2sat
    do_pct = f"{100 * do_set:.0f}"
    sigma = fk.noise_levels(allds)
    sur = k.Surrogate(cfg["surrogate_dir"])
    model_sha = hashlib.sha256(
        (Path(cfg["surrogate_dir"]) / MODEL_FILE).read_bytes()
    ).hexdigest()
    ns = 4 if quick else 10
    md = []
    res = {}
    names = fk.NAMES

    def fit(data, free, fixed, model="emm", n=None):
        fixed = {
            **{q: fk.START[q] for q in names if q not in free},
            **fixed,
        }
        pb = fk.Problem(data, free, fixed, model, sur, sigma)
        return pb, fk.fit(pb, n_start=n or ns)

    prod = dict(k.PRODUCTION.as_dict(), theta=1.195 / 336.0)

    # ---------------------------------------------------------------- R0
    say("R0 legacy")
    ds1 = fs[0]
    stored = []
    for label, folder, kLa, x_max in (
        ("git (julia/of9_working_branch)", "gDCW_conv", 540.0, 7.65),
        ("production", "production_kla336", 336.0, 7.6),
    ):
        v = np.loadtxt(OLD / folder / "fixedstir1" / "fitted_param.dat")
        p = k.Params(**dict(zip(fk.LEGACY_ORDER, v)), x_max=x_max)
        stored.append(
            [label, f(kLa), f(x_max)]
            + [f(getattr(p, q), 4) for q in fk.LEGACY_ORDER]
            + [f(fk.legacy_loss(p, ds1, kLa), 3)]
        )
    start = np.loadtxt(
        OLD / "production_kla336" / "fixedstir1" / "paper_param.dat"
    )
    legacy = []
    for kLa in (254.0, 336.0, 428.0, 540.0):
        p, loss = fk.legacy_fit(ds1, start, kLa, 7.6)
        legacy.append((kLa, p, loss))
    _, new1 = fit(
        fs[:1],
        ("q_max", "K_o", "y_xs", "y_ms", "theta"),
        {"K_s": K_S_PAPER, "x_max": 7.6, "y_os": 1.195},
    )
    new1_loss = fk.legacy_loss(new1.params, ds1, 1.195 / new1.values["theta"])
    res["legacy"] = legacy
    res["new1"] = new1

    # ---------------------------------------------------------------- R1
    say("R1 EMM on fixed-stir")
    free_fs = ("q_max", "K_o", "K_s", "y_xs", "y_ms", "x_max", "theta")
    pb_fs, r_fs = fit(fs, free_fs, {"y_os": 2.4})
    rows_r1 = []
    for d in fs:
        if d.notes.get("do_suspect"):
            free = tuple(q for q in free_fs if q != "theta")
            _, r = fit([d], free, {"y_os": 2.4, "theta": r_fs.values["theta"]})
        else:
            _, r = fit([d], free_fs, {"y_os": 2.4})
        rows_r1.append((d.name, r))
    rows_r1.append(("all three", r_fs))
    prof_fs = {}
    if not quick:
        grids = {
            "q_max": np.geomspace(2.5, 20, 10),
            "K_o": np.geomspace(1e-4, 0.3, 10),
            "K_s": np.geomspace(0.01, 30, 10),
            "x_max": np.linspace(7.2, 9.0, 10),
            "theta": np.geomspace(0.0025, 0.0065, 9),
        }
        for q, g in grids.items():
            say("  profile", q)
            prof_fs[q] = fk.profile(pb_fs, r_fs, q, g)
        fig_profiles(
            prof_fs,
            r_fs.n_residuals - len(free_fs),
            out / "figures" / "profiles_fixedstir.png",
        )

    # ---------------------------------------------------------------- R2
    say("R2 Exp 17 off-gas")
    phases = [
        ("0-23 h (growth)", (0, 4)),
        ("23-30 h", (4, 6)),
        ("30-34.5 h", (6, 8)),
        ("whole run", (0, 9)),
    ]
    yl = {
        ph: [exp17.yields(local, r, [iv])[0] for r in exp17.REACTORS]
        for ph, iv in phases
    }
    yz = [
        exp17.yields(local, r, [(0, 9)], zero=True)[0]["y_os"]
        for r in exp17.REACTORS
    ]
    whole = yl["whole run"]
    mean = {
        q: float(np.mean([w[q] for w in whole]))
        for q in ("y_os", "rq", "y_ms", "y_x")
    }
    bal = k.electron_balance(mean["y_x"], mean["y_ms"], mean["y_os"])
    carbon = 6 * mean["y_ms"] + mean["y_x"] * 1000 / k.MW_BIOMASS_CMOL
    carbon = (carbon + mean["y_os"] * mean["rq"]) / 6.0
    kla_bins = []

    def rpm_bin(i):
        lo = kla_bins[i][0]
        return f"{lo}-{lo + RPM_BIN}"

    series = [exp17.kla_series(local, r, o2sat) for r in exp17.REACTORS]
    rpm = np.concatenate([s["rpm"] for s in series])
    kla = np.concatenate([s["kla"] for s in series])
    # Bins of RPM_BIN from the low end of the stir rates upwards.
    rpm_lo = 100 * int(np.percentile(rpm, 5) // 100)
    for lo in range(rpm_lo, int(rpm.max()), RPM_BIN):
        m = (rpm >= lo) & (rpm < lo + RPM_BIN)
        if m.sum() >= 5:
            kla_bins.append(
                (
                    lo,
                    int(m.sum()),
                    float(np.median(kla[m])),
                    float(np.percentile(kla[m], 10)),
                    float(np.percentile(kla[m], 90)),
                )
            )
    fig_offgas(
        local,
        o2sat,
        out / "figures" / "exp17_our.png",
        out / "figures" / "exp17_kla.png",
    )
    res["y_os_meas"] = mean["y_os"]

    # ------------------------------------------------------------ R3, R4
    say("R3 EMM on Exp 17 and on everything")
    free17 = ("q_max", "K_o", "K_s", "y_xs", "y_ms", "x_max", "y_os")
    _, r_17 = fit(ex, free17, {"theta": r_fs.values["theta"]})
    free_all = free17 + ("theta",)
    pb_all, r_all = fit(allds, free_all, {})
    say("R4 EMM with K_o and K_s fixed")
    free_ko = tuple(q for q in free_all if q != "K_o")
    _, r_ko = fit(allds, free_ko, {"K_o": K_O_CFD})
    free_kk = tuple(q for q in free_ko if q != "K_s")
    _, r_kk = fit(allds, free_kk, {"K_o": K_O_CFD, "K_s": K_S_PAPER})
    prof_all = {}
    if not quick:
        grids = {
            "q_max": np.geomspace(3.0, 12, 9),
            "K_o": np.geomspace(1e-4, 0.15, 9),
            "K_s": np.geomspace(0.05, 30, 9),
            "y_os": np.linspace(1.2, 3.4, 9),
        }
        for q, g in grids.items():
            say("  profile", q)
            prof_all[q] = fk.profile(pb_all, r_all, q, g)
        fig_profiles(
            prof_all,
            r_all.n_residuals - len(free_all),
            out / "figures" / "profiles_all.png",
        )

    # ---------------------------------------------------------------- R5
    say("R5 surrogate against Exp 17")
    r5 = []
    plane_rows = []
    for label, iv in (
        ("9-23 h", (2, 4)),
        ("23-26.5 h", (4, 5)),
        ("26.5-30 h", (5, 6)),
        ("30-34.5 h", (6, 8)),
    ):
        rows = [exp17.yields(local, r, [iv])[0] for r in exp17.REACTORS]
        plane_rows += rows
        pred = np.array([sur(w["qs"], w["our_spec"]) for w in rows])
        qs = np.mean([w["qs"] for w in rows])
        r5.append(
            [
                label,
                f(qs),
                f(np.mean([w["our_spec"] for w in rows])),
                f(np.mean([w["y_os"] for w in rows])),
                f(np.mean([w["mu"] for w in rows]), 2),
                f(np.mean(pred[:, 0]), 2),
                f(np.mean([w["y_ms"] for w in rows]), 2),
                f(np.mean(pred[:, 2]) / qs, 2),
            ]
        )

    # ---------------------------------------------------------------- R6
    say("R6 MEML")
    grow = [truncate(d, 23.0) for d in allds]
    grow_fs = [d for d in grow if d.group == "fixedstir"]
    base = {
        "K_o": K_O_CFD,
        "K_s": r_ko.values["K_s"],
        "y_xs": r_ko.values["y_xs"],
        "y_ms": r_ko.values["y_ms"],
        "x_max": r_ko.values["x_max"],
    }
    free_m = ("q_max", "y_os", "theta")
    _, m_fs = fit(grow_fs, free_m, base, "meml")
    _, m_grow = fit(grow, free_m, base, "meml")
    _, m_full = fit(allds, free_m + ("K_s",), base, "meml")
    _, e_grow = fit(
        grow,
        ("q_max", "y_xs", "y_ms", "y_os", "theta"),
        {**base, "x_max": 12.0},
    )

    # ---------------------------------------------------------------- R7
    say("R7 MEML with a biomass cap")
    base_cap = {q: r_ko.values[q] for q in ("K_o", "y_xs", "y_ms")}
    free_c = ("q_max", "K_s", "y_os", "x_max", "theta")
    _, c_one = fit(allds, free_c, base_cap, "meml_cap")
    pb_cap, c_fit = fit(allds, free_c + ("cap_n",), base_cap, "meml_cap")
    _, e_cap = fit(allds, free_ko + ("cap_n",), {"K_o": K_O_CFD})
    _, e_cap_fs = fit(fs, free_fs + ("cap_n",), {"y_os": 2.4})
    _, e_cap_17 = fit(ex, free17 + ("cap_n",), {"theta": r_fs.values["theta"]})
    prof_cap = []
    if not quick:
        say("  profile cap_n")
        prof_cap = fk.profile(pb_cap, c_fit, "cap_n", CAP_GRID)

    # ------------------------------------------------------- candidates
    say("candidates")
    cands = [
        ("A production", prod, "emm"),
        ("B EMM, everything fitted", r_all.values, "emm"),
        ("C EMM, K_o at the CFD floor", r_ko.values, "emm"),
        ("D EMM, K_o and K_s fixed", r_kk.values, "emm"),
    ]
    meml_c = ("E MEML, fitted to growth phase", m_grow.values, "meml")
    cand_rows = []
    cand_blocks = []
    for name, v, model in cands + [meml_c]:
        p = fk.params_of(v)
        chi_e, rm_e = score(v, allds, sigma, "emm")
        chi_m, _ = score(v, allds, sigma, "meml", sur)
        chi_g, _ = score(v, grow, sigma, "meml", sur)
        st_e, our_e, _ = k.advance_solver(REF_STATE, p)
        st_m, our_m, _ = k.advance_solver(REF_STATE, p, model="meml", sur=sur)
        sim = k.simulate(p, fs[0], "emm", theta=v["theta"])
        y_x = (sim["X"][-1] - sim["X"][0]) / (sim["G"][0] - sim["G"][-1])
        eb = k.electron_balance(float(y_x), p.y_ms, p.y_os)
        lo, hi = sur.qs_range(p.y_os)
        cand_rows.append(
            {
                "name": name,
                "v": v,
                "chi_emm": chi_e,
                "chi_meml": chi_m,
                "chi_meml_growth": chi_g,
                "rm": rm_e,
                "kla_fs": p.y_os / v["theta"],
                "emm": (st_e[0], st_e[2], our_e),
                "meml": (st_m[0], st_m[2], our_m),
                "recovery": eb["electron_recovery"],
                "qs_set": float(k.uptake(228.8, o2_set, p)),
                "qs_sat": float(k.uptake(228.8, 0.2528, p)),
                "hull": (lo, hi),
            }
        )
        block = p.to_foam(kLa=p.y_os / v["theta"], use_ml=model == "meml")
        tag = name.split()[0]
        (out / "candidates" / f"{tag}.txt").write_text(block + "\n")
        cand_blocks.append((name, block))
    q_fs = [
        ("X", "biomass (kg/m3)"),
        ("G", "glucose (mol/m3)"),
        ("M", "muconate (mol/m3)"),
        ("DO", "DO (fraction of sat.)"),
    ]
    q_ex = q_fs[:3] + [("OUR", "O2 uptake (mol/m3/h)")]
    # Capped MEML is not production, so it takes the fourth colour.
    cap_models = [
        ("EMM refit (C)", r_ko.values, "emm", C["orange"]),
        ("MEML, no cap (E)", m_grow.values, "meml", C["aqua"]),
        ("MEML with cap", c_fit.values, "meml_cap", C["yellow"]),
    ]
    fig_fits(fs, cap_models, sur, out / "figures" / "cap_fixedstir.png", q_fs)
    fig_fits(ex, cap_models, sur, out / "figures" / "cap_exp17.png", q_ex)
    fig_fits(
        fs,
        [
            ("production", prod, "emm"),
            ("EMM refit (C)", r_ko.values, "emm"),
            ("MEML, growth fit (E)", m_grow.values, "meml"),
        ],
        sur,
        out / "figures" / "fits_fixedstir.png",
        q_fs,
    )
    fig_fits(
        ex,
        [
            ("production", prod, "emm"),
            ("EMM refit (C)", r_ko.values, "emm"),
            ("MEML, growth fit (E)", m_grow.values, "meml"),
        ],
        sur,
        out / "figures" / "fits_exp17.png",
        q_ex,
    )
    fig_plane(
        sur,
        plane_rows,
        [
            ("production, y_os 1.195", 1.195),
            (
                f"refit (C), y_os {r_ko.values['y_os']:.2f}",
                r_ko.values["y_os"],
            ),
        ],
        out / "figures" / "surrogate_plane.png",
    )
    # Same colour per parameter set as in the fit figures; aqua stays MEML.
    fig_qs_o2(
        [
            (cands[0], C["blue"]),
            (cands[2], C["orange"]),
            (cands[1], C["yellow"]),
        ],
        out / "figures" / "qs_vs_o2.png",
        o2_set,
    )

    # =============================================================== text
    say("writing report")
    try:
        rev = subprocess.run(
            ["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
        ).stdout.strip()
    except OSError:
        rev = "?"
    A, B, Cc, D, E = cand_rows
    th_fs = r_fs.values["theta"]
    qs_e17 = float(np.mean([w["qs"] for w in yl[phases[1][0]]]))
    kla_all = [mean["y_os"] / th_fs] + [c["kla_fs"] for c in (B, Cc, D)]
    kla_lo, kla_hi = min(kla_all), max(kla_all)
    y_grow = float(np.mean([w["y_os"] for w in yl[phases[0][0]]]))
    y_end = float(np.mean([w["y_os"] for w in yl[phases[2][0]]]))
    # Fastest growth of the surrogate along the production O2/glucose
    # ratio.
    qs_line = np.linspace(0.0, 12.0, 481)
    mu_line = sur(qs_line, k.PRODUCTION.y_os * qs_line)[:, 0]
    mu_prod_max = float(np.max(mu_line))
    md.append(
        "# EMM / MEML kinetics: how the parameters were fit, and refits\n"
    )
    md.append(
        f"Generated {datetime.date.today()} by `make_report.py` at commit "
        f"`{rev}`{' (quick mode)' if quick else ''}. Units: biomass kg/m3, "
        "glucose, muconate and O2 mol/m3, time h, specific rates mol/kg/h.\n"
    )
    md.append("## Summary\n")
    md.append(
        "\n".join(
            [
                "- **The production parameters are a fit to one run with "
                "an assumed kLa of 336 1/h.** Case files label it 336 or "
                "428; the solver does not use that entry.",
                "- **The oxygen yield y_os was never determined by the "
                "fit.** Dissolved-O2 data only fix y_os / kLa "
                f"({f(th_fs)} h for the fixed-stir runs), so y_os came out "
                "as whatever kLa was typed in.",
                "- **Off-gas analysis measures it: "
                f"{mean['y_os']:.2f} mol O2 per mol glucose** over the whole "
                "run (four Exp 17 reactors), against 1.195 in production.",
                "- **The production fit stopped at its starting guess.** "
                f"A real search gives q_max {r_all.values['q_max']:.2f} "
                "instead of 8.26 and no measurable O2 limitation down to "
                f"{do_pct} % DO. K_o = 0.1134 mol/m3 (3.6 mg/L) is not "
                "supported: it predicts an uptake of "
                f"{A['qs_set']:.1f} mol/kg/h at {do_pct} % DO where Exp 17 "
                f"shows {qs_e17:.1f}.",
                "- **With y_os at its measured value the fixed-stir runs "
                f"imply a kLa of {kla_lo:.0f}-{kla_hi:.0f} 1/h** "
                "(depending on the fit). In Exp 17 the same reactor "
                f"reaches {kla_bins[-2][2]:.0f} 1/h at {rpm_bin(-2)} rpm "
                f"and {kla_bins[-1][2]:.0f} 1/h at {rpm_bin(-1)} rpm.",
                "- **MEML with the production numbers cannot reproduce the "
                f"lab runs** (it grows at most {mu_prod_max:.2f} 1/h; the "
                f"cultures grew at {r5[0][4]} 1/h). Refitted, MEML matches "
                "the growth phase "
                f"with y_os {m_grow.values['y_os']:.2f}, but it has no "
                "mechanism for the biomass plateau and overshoots after "
                "23 h.",
                "- **A biomass cap fixes that, if it is sharper than "
                "EMM's.** MEML with growth "
                "multiplied by 1 - (X / x_max)^n fits both experiments "
                f"over the whole run at n {c_fit.values['cap_n']:.1f} "
                f"(chi2 {c_fit.chi2:.0f}, EMM set C {r_ko.chi2:.0f}). "
                f"With EMM's own form (n = 1) it does not ({c_one.chi2:.0f}"
                "). EMM itself also fits better with the sharper cap "
                f"({e_cap.chi2:.0f}). putidaFoam takes n as an optional "
                "`cap_n` entry (sets F and G in section 6).",
                "- Candidate sets B-E below fit both experiments with the "
                "solver as it was; A is the current production set for "
                "comparison. F and G need a build that reads `cap_n`.",
            ]
        )
        + "\n"
    )

    md.append("## 1. Where the production parameters come from\n")
    md.append(
        "Both stored parameter sets regenerate their own `soln.dat` and "
        "loss. Loss is the original metric: sum over biomass, glucose, "
        "muconate and dissolved O2 of MSE / mean(obs^2), fixedstir1 only.\n"
    )
    md.append(
        table(
            ["set", "kLa", "x_max"] + list(fk.LEGACY_ORDER) + ["loss"], stored
        )
        + "\n"
    )
    md.append(
        "Repeating the original procedure (same script logic, same "
        "starting guess) at four assumed kLa values:\n"
    )
    md.append(
        table(
            ["kLa assumed"] + list(fk.LEGACY_ORDER) + ["loss", "y_os / kLa"],
            [
                [f(kLa)]
                + [f(getattr(p, q), 4) for q in fk.LEGACY_ORDER]
                + [f(loss, 3), f(p.y_os / kLa, 3)]
                for kLa, p, loss in legacy
            ],
        )
        + "\n"
    )
    md.append(
        f"The loss does not change with kLa, y_os scales with it, and q_max "
        f"and K_s stay at their starting values ({start[0]:.3f} and "
        f"{start[5]:.3f}, after Jamshidzadeh et al. 2025). A multi-start "
        "search in log-parameters on the same run, with K_s and x_max held "
        f"at the production values, reaches loss {f(new1_loss, 3)} with "
        f"q_max {new1.values['q_max']:.2f} and K_o at "
        f"{f(new1.values['K_o'], 2)}"
        + (
            ", its lower bound (no O2 limitation).\n"
            if "K_o" in new1.at_bound
            else ".\n"
        )
    )

    md.append("## 2. What the fixed-stir runs determine\n")
    md.append(
        "EMM fitted to the three fixed-stir runs with quasi-steady "
        "dissolved O2, theta = y_os / kLa free. Residuals are scaled by the "
        "replicate scatter (noise levels in the appendix), so chi2 near the "
        "number of residuals means a fit as good as the replicates agree.\n"
    )
    head = ["data", "chi2 / n"] + list(free_fs) + ["at bound"]
    md.append(
        table(
            head,
            [
                [name, f"{r.chi2:.1f} / {r.n_residuals}"]
                + vals_row(r.values, free_fs)
                + [", ".join(r.at_bound) or "-"]
                for name, r in rows_r1
            ],
        )
        + "\n"
    )
    if prof_fs:
        dof = r_fs.n_residuals - len(free_fs)
        md.append(
            table(
                ["parameter", "best", "lower", "upper"],
                [[q, *interval(prof_fs[q], dof)] for q in prof_fs],
            )
            + "\n"
        )
        md.append("![profiles](figures/profiles_fixedstir.png)\n")
    md.append(
        "- q_max, the yields, x_max and theta are determined. K_o only has "
        "an upper limit: the data show no O2 dependence of uptake over the "
        "measured DO range (45-95 % of saturation).\n"
        "- K_s differs between replicates and is driven by how the late "
        "slow-down of uptake is shared between K_s and x_max; fixedstir2 "
        "never runs out of glucose and says nothing about it.\n"
        "- Data problems found on the way:\n"
        "  - fixedstir1 was stored with 100 % DO = 0.188 mol/m3, fixedstir2 "
        "and 3 with 0.214. `datasets.yaml` now carries DO as a fraction.\n"
        "  - fixedstir2 (BR2) has the same DO % as BR1 at every sample in "
        "the source spreadsheet. Its DO is left out of the fits.\n"
        "  - Stir rate and air flow of the fixed-stir runs are not "
        "recorded in the project files.\n"
    )

    md.append("## 3. Oxygen yield and kLa measured in Exp 17\n")
    md.append(
        "Exp 17 is a batch with DO held at a set point by the stir rate "
        "and off-gas mass spectrometry on all four reactors. Oxygen and "
        "CO2 are integrated between samples and divided by the glucose "
        "used.\n"
    )
    md.append(
        table(
            [
                "phase",
                "y_os (O2 / glucose)",
                "range",
                "RQ",
                "qs",
                "growth 1/h",
            ],
            [
                [
                    ph,
                    f(np.mean([w["y_os"] for w in yl[ph]])),
                    f"{min(w['y_os'] for w in yl[ph]):.2f}-"
                    f"{max(w['y_os'] for w in yl[ph]):.2f}",
                    f(np.mean([w["rq"] for w in yl[ph]])),
                    f(np.mean([w["qs"] for w in yl[ph]])),
                    f(np.mean([w["mu"] for w in yl[ph]]), 2),
                ]
                for ph, _ in phases
            ],
        )
        + "\n"
    )
    md.append(
        f"- Whole-run y_os is {mean['y_os']:.2f} (reactors "
        f"{min(w['y_os'] for w in whole):.2f}-"
        f"{max(w['y_os'] for w in whole):.2f}). The off-gas trace starts "
        "slightly negative; shifting it to zero at the start raises the "
        f"value to {np.mean(yz):.2f}, so the true number is if anything "
        "higher.\n"
        "- y_os is not constant: it rises as growth stops. A single y_os "
        "is a simplification in both EMM and MEML.\n"
        f"- Balances on the whole-run yields (muconate {mean['y_ms']:.3f} "
        f"mol/mol, biomass {mean['y_x'] * 1000:.1f} g/mol): "
        f"{carbon * 100:.0f} % of the carbon and "
        f"{bal['electron_recovery'] * 100:.0f} % of the electrons are "
        "accounted for. A closed electron balance would need y_os "
        f"{bal['y_os_closed']:.2f}.\n"
    )
    md.append(
        "kLa follows from the same data as OUR / (C* - C), with C* taken "
        "from the outlet gas:\n"
    )
    md.append(
        table(
            ["stir rate (rpm)", "points", "kLa median (1/h)", "10-90 %"],
            [
                [f"{lo}-{lo + RPM_BIN}", n, f"{med:.0f}", f"{a:.0f}-{b:.0f}"]
                for lo, n, med, a, b in kla_bins
            ],
        )
        + "\n"
    )
    md.append(
        "![our](figures/exp17_our.png)\n\n![kla](figures/exp17_kla.png)\n"
    )
    md.append(
        "What each kLa assumption means for the fixed-stir runs, "
        f"with theta = {f(th_fs)} h from the fixed-stir fit of section 2 "
        "(theta moves between "
        f"{min(r.values['theta'] for _, r in rows_r1):.4f} and "
        f"{max(r.values['theta'] for _, r in rows_r1):.4f} h across the "
        "replicates):\n"
    )
    md.append(
        table(
            ["assumption", "kLa (1/h)", "y_os"],
            [
                [f"kLa {kk:.0f}", f"{kk:.0f}", f(th_fs * kk)]
                for kk in (93, 336, 428, 540)
            ]
            + [
                [
                    "y_os measured in Exp 17",
                    f"{mean['y_os'] / th_fs:.0f}",
                    f(mean["y_os"]),
                ],
                [
                    "y_os of the growth phase",
                    f"{np.mean([w['y_os'] for w in yl[phases[0][0]]]) / th_fs:.0f}",
                    f(np.mean([w["y_os"] for w in yl[phases[0][0]]])),
                ],
            ],
        )
        + "\n"
    )

    md.append("## 4. EMM fitted to both experiments\n")
    md.append(
        "In Exp 17 the dissolved O2 is prescribed from the log and the "
        "measured oxygen uptake enters the residual, so no kLa is involved "
        "and y_os is determined directly. The last columns are RMSE divided "
        "by the noise level, per experiment (fs = fixed-stir, E17 = "
        "Exp 17).\n"
    )
    rows = []
    for label, v, r in (
        ("production (A)", prod, None),
        ("Exp 17 only", r_17.values, r_17),
        ("both, all free (B)", r_all.values, r_all),
        ("both, K_o fixed (C)", r_ko.values, r_ko),
        ("both, K_o and K_s fixed (D)", r_kk.values, r_kk),
    ):
        chi, rm = score(v, allds, sigma)
        rows.append(
            [label, f"{chi:.0f}"]
            + vals_row(v, fk.KINETIC)
            + [f"{v['y_os'] / v['theta']:.0f}"]
            + rm_cells(rm)
        )
    md.append(
        table(
            ["fit", "chi2 (all data)"]
            + list(fk.KINETIC)
            + ["kLa fs"]
            + RM_HEAD,
            rows,
        )
        + "\n"
    )
    if prof_all:
        dof = r_all.n_residuals - len(free_all)
        md.append(
            table(
                ["parameter", "best", "lower", "upper"],
                [[q, *interval(prof_all[q], dof)] for q in prof_all],
            )
            + "\n"
        )
        md.append(
            "Ranges are where chi2 stays within the 95 % level of the best "
            "fit, treating residuals as independent. Hourly oxygen-uptake "
            "points are correlated in time and the model misses more than "
            "the replicates scatter (chi2 is about "
            f"{r_all.chi2 / dof:.1f} per degree of freedom), so the real "
            "ranges are wider.\n"
        )
        md.append("![profiles](figures/profiles_all.png)\n")
    md.append(
        "![fixed-stir](figures/fits_fixedstir.png)\n\n"
        "![exp17](figures/fits_exp17.png)\n"
    )
    md.append(
        "- Fitted separately, the two experiments give similar kinetics: "
        f"q_max {r_fs.values['q_max']:.2f} (fixed-stir) and "
        f"{r_17.values['q_max']:.2f} (Exp 17), y_xs "
        f"{r_fs.values['y_xs']:.3f} and {r_17.values['y_xs']:.3f}, y_ms "
        f"{r_fs.values['y_ms']:.2f} and {r_17.values['y_ms']:.2f}.\n"
        f"- Uptake is the same at {do_pct} % DO (Exp 17) and at 45-95 % DO "
        "(fixed-stir), which is what bounds K_o. The production K_o would "
        f"cut uptake at {do_pct} % DO to {A['qs_set']:.1f} mol/kg/h; "
        f"{qs_e17:.1f} was measured (23-30 h).\n"
        f"- x_max differs between experiments ({r_fs.values['x_max']:.1f} "
        f"fixed-stir, {r_17.values['x_max']:.1f} Exp 17). It is a property "
        "of the batch (whatever stops growth is not in the model), not a "
        "kinetic constant, and should be chosen per reactor. Most of this "
        "gap comes from the form of the cap (section 6).\n"
        "- A large K_s improves the fit because specific uptake falls late "
        "in every run while glucose is still present. That is a stand-in "
        "for whatever slows the culture, not a glucose affinity. Set D "
        "keeps the literature K_s and accepts a worse late-phase fit.\n"
    )

    md.append("## 5. The surrogate against Exp 17\n")
    md.append(
        "The measured specific glucose and oxygen uptake of Exp 17 are fed "
        "straight to the surrogate (no kinetic parameters involved), "
        "averaged over the four reactors. The first 9 h are left out: the "
        "off-gas signal there is of the size of its offset.\n"
    )
    md.append(
        table(
            [
                "interval",
                "qs",
                "O2 uptake / X",
                "ratio",
                "growth measured",
                "growth surrogate",
                "muconate yield measured",
                "surrogate",
            ],
            r5,
        )
        + "\n"
    )
    md.append("![plane](figures/surrogate_plane.png)\n")
    md.append(
        "- In the growth phase the surrogate comes close to the measured "
        f"growth rate ({r5[0][5]} against {r5[0][4]} 1/h) when given the "
        f"measured oxygen uptake, and its fixed muconate yield ({r5[0][7]}) "
        "matches the data.\n"
        "- After 23 h the cultures stop growing while still taking up "
        "glucose and oxygen at a ratio the surrogate associates with fast "
        "growth. Nothing in the (glucose, O2) inputs tells it to stop.\n"
        "- Along the production ratio (1.195) the surrogate is evaluated "
        "at the low-oxygen edge of its training data; the measured ratios "
        "sit in the middle of it.\n"
    )

    md.append("## 6. MEML\n")
    rows = []
    for label, v, data, r in (
        ("production numbers, EMM", prod, allds, None),
        ("production numbers, MEML", prod, allds, None),
        ("MEML fit, fixed-stir growth phase only", m_fs.values, grow_fs, m_fs),
        ("MEML fit, both, growth phase (E)", m_grow.values, grow, m_grow),
        ("EMM fit, both, growth phase", e_grow.values, grow, e_grow),
        ("MEML fit, both, whole run", m_full.values, allds, m_full),
        ("EMM fit, both, whole run (C)", r_ko.values, allds, r_ko),
    ):
        model = "meml" if "MEML" in label else "emm"
        chi, _ = score(v, data, sigma, model, sur)
        p = fk.params_of(v)
        x_end = [
            k.simulate(p, d, model, sur, v["theta"])["X"][-1]
            for d in (fs[0], ex[2])
        ]
        n = fk.Problem(data, (), dict(v), model, sur, sigma).n_residuals
        rows.append(
            [
                label,
                f"{chi:.0f} / {n}",
                f(v["q_max"]),
                f(v["K_s"]),
                f(v["y_os"]),
                f"{v['y_os'] / v['theta']:.0f}",
                f"{x_end[0]:.1f} ({fs[0].X[-1]:.1f})",
                f"{x_end[1]:.1f} ({ex[2].X[-1]:.1f})",
            ]
        )
    md.append(
        table(
            [
                "case",
                "chi2 / n",
                "q_max",
                "K_s",
                "y_os",
                "kLa fs",
                "final X fixedstir1 (measured)",
                "final X Exp 17 r52 (measured)",
            ],
            rows,
        )
        + "\n"
    )
    md.append(
        "- MEML run with the EMM production numbers (current practice) "
        "misses both experiments by a wide margin.\n"
        "- Fitted to the growth phase, MEML does as well as EMM. Fitted to "
        f"the fixed-stir runs alone, with no off-gas data, it returns y_os "
        f"{m_fs.values['y_os']:.2f}: the growth rate itself tells MEML how "
        "much oxygen per glucose the cells use. Off-gas gives "
        f"{y_grow:.2f} in the growth phase and {mean['y_os']:.2f} over "
        "the whole run.\n"
        "- Over the whole run MEML is much worse than EMM because nothing "
        "caps the biomass. In putidaFoam a MEML culture near x_max keeps "
        "growing (first production update: 7.33 to 9.25 kg/m3).\n"
        "- EMM's y_xs, y_ms and x_max are unused by MEML; set E carries "
        "the values of set C so the file stays complete.\n"
    )
    if model_sha == KESTREL_MODEL_SHA256:
        md.append(
            "- **Model file.** These fits use the same surrogate model "
            f"file as the 1000 L cases on Kestrel (sha256 {model_sha[:12]}"
            ").\n"
        )
    else:
        md.append(
            "- **Model file.** These fits use the surrogate model file "
            f"with sha256 {model_sha[:12]}. The 1000 L cases on Kestrel "
            "carried a different training of the same network on "
            f"2026-11-25 ({KESTREL_MODEL_SHA256[:12]}; scalers and hull "
            "files are the same). Its growth rates differ by 0.6 % at the "
            "median over the training inputs, so the MEML sets should be "
            "refitted with that file before they are used there.\n"
        )

    md.append("### A biomass cap for MEML, and a sharper one for EMM\n")
    md.append(
        "The cap multiplies growth by 1 - (X / x_max)^n and leaves glucose "
        "uptake, oxygen uptake and muconate as they are. n = 1 is the cap "
        "EMM always had; a larger n leaves early growth alone and stops it "
        "more abruptly near x_max. putidaFoam reads n as the optional "
        "entry `cap_n`. All rows are scored on both experiments over the "
        "whole run with K_o fixed as in set C. The MEML refits free q_max, "
        "K_s, y_os, x_max and theta; EMM frees y_xs and y_ms as well.\n"
    )
    cap_cases = (
        ("EMM (C), n = 1", r_ko.values, "emm"),
        ("EMM refitted, n fitted (G)", e_cap.values, "emm"),
        ("MEML, no cap, refitted", m_full.values, "meml"),
        ("MEML set E, cap switched on, n = 1", m_grow.values, "meml_cap"),
        ("MEML refitted, n = 1", c_one.values, "meml_cap"),
        ("MEML refitted, n fitted (F)", c_fit.values, "meml_cap"),
    )
    rows, rows_rm = [], []
    cap_ref = {}
    for label, v, model in cap_cases:
        chi, rm = score(v, allds, sigma, model, sur)
        p = fk.params_of(v)
        x_end = [
            k.simulate(p, d, model, sur, v["theta"])["X"][-1]
            for d in (fs[0], ex[2])
        ]
        st, our, _ = k.advance_solver(REF_STATE, p, model=model, sur=sur)
        cap_ref[label] = (chi, st[0], our)
        rows.append(
            [
                label,
                f"{chi:.0f}",
                f(v["q_max"]),
                f(v["K_s"]),
                f(v["y_os"]),
                "-" if model == "meml" else f(v["x_max"]),
                "-" if model == "meml" else f(p.cap_n),
                f"{v['y_os'] / v['theta']:.0f}",
                f"{x_end[0]:.1f} ({fs[0].X[-1]:.1f})",
                f"{x_end[1]:.1f} ({ex[2].X[-1]:.1f})",
                f"{st[0]:.2f}, {our:.1f}",
            ]
        )
        rows_rm.append([label] + rm_cells(rm))
    md.append(
        table(
            [
                "case",
                "chi2 (all data)",
                "q_max",
                "K_s",
                "y_os",
                "x_max",
                "n",
                "kLa fs",
                "final X fixedstir1 (measured)",
                "final X Exp 17 r52 (measured)",
                "CFD ref. update: X, OUR",
            ],
            rows,
        )
        + "\n"
    )
    md.append("RMSE divided by the noise level, per experiment:\n")
    md.append(table(["case"] + RM_HEAD, rows_rm) + "\n")
    if prof_cap:
        md.append(
            "chi2 against the sharpness n, the other parameters refitted:\n"
        )
        md.append(
            table(
                ["n", "chi2", "q_max", "K_s", "y_os", "x_max"],
                [
                    [f(n, 3), f"{chi:.0f}"]
                    + [f(v[q]) for q in ("q_max", "K_s", "y_os", "x_max")]
                    for n, chi, v in prof_cap
                ],
            )
            + "\n"
        )
        md.append(
            "The profile walks outwards from the optimum with one start "
            "per point, so its n = 1 row can sit a little above the "
            "ten-start fit in the first table "
            f"({cap_ref['MEML refitted, n = 1'][0]:.0f}).\n"
        )
    md.append("![cap fixed-stir](figures/cap_fixedstir.png)\n")
    md.append("![cap Exp 17](figures/cap_exp17.png)\n")
    chi_c = cap_ref["EMM (C), n = 1"][0]
    chi_e = cap_ref["EMM refitted, n fitted (G)"][0]
    chi_one = cap_ref["MEML refitted, n = 1"][0]
    chi_fit, x_ref, our_ref = cap_ref["MEML refitted, n fitted (F)"]
    md.append(
        "- **EMM's cap form does not rescue MEML.** With n = 1 the refit "
        f"reaches chi2 {chi_one:.0f} against {chi_c:.0f} for EMM, with "
        f"x_max pushed to {c_one.values['x_max']:.1f} kg/m3, above any "
        "measured biomass. That cap slows growth from the start of the "
        "batch. EMM absorbs this in a free y_xs; MEML cannot, because its "
        "growth per glucose comes from the surrogate.\n"
        f"- **A sharper cap does.** With n fitted ({c_fit.values['cap_n']:.1f}"
        f") MEML reaches chi2 {chi_fit:.0f} with one x_max of "
        f"{c_fit.values['x_max']:.2f} kg/m3 for both reactors, and it "
        "follows the rise and fall of the Exp 17 oxygen uptake.\n"
        f"- **EMM gains from the same change** ({chi_c:.0f} to "
        f"{chi_e:.0f}, n {e_cap.values['cap_n']:.1f}, x_max "
        f"{e_cap.values['x_max']:.2f}). On Exp 17 alone it takes EMM from "
        f"chi2 {r_17.chi2:.0f} to {e_cap_17.chi2:.0f} for "
        f"{e_cap_17.n_residuals} residuals; the fixed-stir runs alone "
        f"hardly tell the two forms apart ({r_fs.chi2:.0f} and "
        f"{e_cap_fs.chi2:.0f}).\n"
        "- **Most of the gap in x_max between the experiments came from "
        "the n = 1 form.** Fitted separately, EMM's x_max is "
        f"{r_fs.values['x_max']:.1f} (fixed-stir) and "
        f"{r_17.values['x_max']:.1f} (Exp 17) with n = 1, and "
        f"{e_cap_fs.values['x_max']:.1f} and "
        f"{e_cap_17.values['x_max']:.1f} with n fitted.\n"
        + (
            "- **The exact sharpness matters little above about 4** (table "
            "above); n = 30, close to a hard stop at x_max, costs "
            f"{prof_cap[-1][1] - chi_fit:.0f} in chi2.\n"
            if prof_cap
            else ""
        )
        + "- **Oxygen demand in the CFD.** For the reference update capped "
        f"MEML gives X {x_ref:.2f} and a mean oxygen sink of {our_ref:.1f} "
        f"mol/m3/h, against {Cc['emm'][2]:.1f} for EMM set C and "
        f"{E['meml'][2]:.1f} for uncapped MEML set E.\n"
        "- **In the solver.** putidaFoam reads an optional `cap_n` entry in "
        "`microbeUpdateProperties`, a plain number, and applies the cap to "
        "growth in both models (`biomass_cap` in `microbeModel.C`). "
        "Without the entry nothing changes: EMM keeps n = 1 and MEML has "
        "no cap. The tests compile `microbeModel.C` on its own and compare "
        "it with this Python model. On 2026-11-25 the change was also "
        "built with OpenFOAM 9 on Kestrel and one bio update of the 1000 L "
        "reacting case was run for EMM and MEML, each with and without "
        "the entry. The liquid averages after the update agree with this "
        "Python model to 2e-5 when it is given the model file of that "
        "case.\n"
        "- The data do not say what stops growth (question 4 below), so "
        "x_max and n describe these batches; a run with a different "
        "medium or feed may plateau elsewhere.\n"
    )
    md.append(
        "Parameter blocks for the two capped sets (they need a putidaFoam "
        "build that reads `cap_n`):\n"
    )
    for name, r, model in (
        ("F MEML with the cap", c_fit, "meml_cap"),
        ("G EMM with the sharper cap", e_cap, "emm"),
    ):
        p = r.params
        block = p.to_foam(kLa=p.y_os / r.values["theta"], model=model)
        (out / "candidates" / f"{name.split()[0]}.txt").write_text(
            block + "\n"
        )
        md.append(f"**{name}**\n\n```\n{block}\n```\n")

    md.append("## 7. Candidate parameter sets\n")
    md.append(
        table(
            ["set"] + list(fk.KINETIC) + ["implied fixed-stir kLa"],
            [
                [c["name"]]
                + vals_row(c["v"], fk.KINETIC)
                + [f"{c['kla_fs']:.0f}"]
                for c in cand_rows
            ],
        )
        + "\n"
    )
    md.append(
        table(
            [
                "set",
                "chi2 all data, own model",
                "chi2 MEML, growth phase",
                "electron recovery",
                f"qs at {do_pct} % DO",
                "surrogate qs range",
                "CFD ref. update EMM: X, OUR",
                "MEML: X, OUR",
            ],
            [
                [
                    c["name"].split()[0],
                    f"{c['chi_meml' if 'MEML' in c['name'] else 'chi_emm']:.0f}",
                    f"{c['chi_meml_growth']:.0f}",
                    f"{c['recovery'] * 100:.0f} %",
                    f(c["qs_set"]),
                    f"{c['hull'][0]:.2f}-{c['hull'][1]:.2f}",
                    f"{c['emm'][0]:.2f}, {c['emm'][2]:.1f}",
                    f"{c['meml'][0]:.2f}, {c['meml'][2]:.1f}",
                ]
                for c in cand_rows
            ],
        )
        + "\n"
    )
    md.append("![qs](figures/qs_vs_o2.png)\n")
    md.append(
        "- **CFD ref. update** is one 4 h putidaFoam bio update from the "
        "state of the 1000 L runs (X 7.33, O2 0.2528, glucose 228.8): "
        "biomass after it and the mean oxygen sink in mol/m3/h.\n"
        "- **Set C** is the one I would start from for EMM: fitted to both "
        f"experiments, K_o at the value the CFD already uses as its O2-sink "
        f"floor ({K_O_CFD} mol/m3 = 0.1 mg/L). The data prefer an even "
        "smaller K_o (set B, where it sits on the search bound); fixing it "
        f"raises chi2 by {(Cc['chi_emm'] / B['chi_emm'] - 1) * 100:.0f} % "
        "and keeps a finite O2 shut-off for the CFD. The lab data only "
        f"cover {o2_set:.3f}-0.18 mol/m3 dissolved O2; below that every "
        "set is an extrapolation.\n"
        "- **Set E** is the MEML counterpart. Use it knowing that MEML "
        "overshoots biomass once growth should stop. Set F in section 6 "
        "removes the overshoot; it needs a putidaFoam build that reads "
        "`cap_n`.\n"
        "- x_max in B-D is a compromise between the two experiments; pick "
        "it for the reactor being simulated.\n"
        "- All refits roughly double the oxygen demand per glucose "
        "relative to production and remove the strong O2 limitation of "
        "uptake, so CFD oxygen sinks are higher at the same state "
        f"({A['emm'][2]:.0f} -> {Cc['emm'][2]:.0f} mol/m3/h for the "
        "reference update with EMM).\n"
    )
    for name, block in cand_blocks:
        md.append(f"**{name}**\n\n```\n{block}\n```\n")

    md.append("## 8. Questions for the experimentalists\n")
    md.append(
        "1. Stir rate and air flow of the fixed-stir runs BR1-BR3. With "
        f"y_os measured, they imply kLa of {kla_lo:.0f}-{kla_hi:.0f} 1/h; "
        "Exp 17 reaches that in its upper stir-rate bins (section 3).\n"
        "2. BR1 and BR2 have identical DO % at all six samples in "
        "`0325ODDCW_conversion.xlsx`. Is BR2's DO a copy?\n"
        "3. Which concentration is 100 % DO: 0.188 mol/m3 (6.02 mg/L) is "
        "used here for all runs; the fixedstir2/3 files used 0.214.\n"
        "4. What stops growth at 7-9 kg/m3 while glucose remains (nitrogen, "
        "another nutrient, product)? It sets x_max and it is what MEML "
        "lacks.\n"
        f"5. About {(1 - carbon) * 100:.0f} % of the carbon in Exp 17 is "
        "unaccounted for. By-product time series (gluconate, "
        "2-ketogluconate) would close it.\n"
    )

    md.append("## 9. Note for the CFD session\n")
    md.append(
        "```\n"
        "Kinetics fit rework (branch kinetics-fit-rework, "
        "scripts/ma_reactions/kinetics_fit):\n"
        "- The production set (y_os 1.195, K_o 0.1134, q_max 8.26) is a "
        "fit to fixedstir1 only at an assumed kLa of 336 1/h; the kLa 428 "
        "in the case files is a mislabel of an unused entry.\n"
        f"- Off-gas (Exp 17, 4 reactors): y_os = {mean['y_os']:.2f} mol "
        f"O2/mol glucose over the run, ~{y_grow:.1f} in growth, "
        f"~{y_end:.1f} at the end.\n"
        f"- Refit to both experiments: q_max {r_ko.values['q_max']:.2f}, "
        f"K_s {r_ko.values['K_s']:.2f}, y_os {r_ko.values['y_os']:.2f}, "
        f"y_xs {r_ko.values['y_xs']:.4f}, y_ms {r_ko.values['y_ms']:.4f}, "
        f"x_max {r_ko.values['x_max']:.2f}, K_o fixed at {K_O_CFD} "
        "(data only bound it from above).\n"
        f"- No O2 limitation of uptake is visible between {do_pct} % and "
        "95 % DO, so the 450B pressure-step result (flat uptake) is "
        "expected.\n"
        f"- MEML refit (growth phase): q_max {m_grow.values['q_max']:.2f}, "
        f"y_os {m_grow.values['y_os']:.2f}. MEML has no biomass cap and "
        "overshoots after growth should stop.\n"
        "- New optional entry cap_n in putidaFoam on this branch (plain "
        "number; growth is multiplied by 1 - (X/x_max)^cap_n; without the "
        "entry EMM and MEML behave as before). With it MEML fits both "
        f"experiments over the whole run: cap_n {c_fit.values['cap_n']:.1f}"
        f", x_max {c_fit.values['x_max']:.2f}, q_max "
        f"{c_fit.values['q_max']:.2f}, K_s {c_fit.values['K_s']:.1f}, y_os "
        f"{c_fit.values['y_os']:.2f} (chi2 {chi_fit:.0f} against "
        f"{chi_c:.0f} for EMM set C). cap_n 1, EMM's old form, does not "
        f"work for MEML. Reference update with it: X {x_ref:.2f}, mean OUR "
        f"{our_ref:.1f} mol/m3/h. Built and smoke-tested on Kestrel on "
        "2026-11-25 in /scratch/nwimer/agent/accelerate/capn-check (see "
        "the README there); build/two-phase-reset was not touched.\n"
        f"- Reference 4 h update (X 7.33, O2 0.2528, G 228.8): set C EMM "
        f"mean OUR {Cc['emm'][2]:.1f} mol/m3/h (production 50.3); set E "
        f"MEML {E['meml'][2]:.1f} (production 56.1). The stage-2 notes "
        "put the tank's maximum transfer at 64 mol/m3/h (CFD kLa 254 "
        "1/h), so with the refits demand exceeds supply at that state.\n"
        f"- Small-tank kLa measured during Exp 17: "
        + ", ".join(
            f"{lo}-{lo + RPM_BIN} rpm {med:.0f}"
            for lo, _, med, _, _ in kla_bins
        )
        + " 1/h.\n"
        "```\n"
    )

    md.append("## Appendix\n")
    md.append(
        "Noise levels (standard deviation across replicates, RMS over the "
        "samples):\n"
    )
    md.append(
        table(
            ["experiment", "X", "G", "M", "DO fraction", "OUR"],
            [
                [g] + [f(sigma[g][q]) for q in ("X", "G", "M", "DO", "OUR")]
                for g in sigma
            ],
        )
        + "\n"
    )
    md.append(
        "Multi-start agreement (starts within 1 % of the best chi2): "
        + "; ".join(
            f"{label} {sum(c < 1.01 * r.start_chi2[0] + 0.1 for c in r.start_chi2)}"
            f"/{len(r.start_chi2)}"
            for label, r in (
                ("fixed-stir", r_fs),
                ("Exp 17", r_17),
                ("B", r_all),
                ("C", r_ko),
                ("D", r_kk),
                ("E", m_grow),
                ("MEML cap n = 1", c_one),
                ("MEML cap n fitted", c_fit),
                ("EMM n fitted", e_cap),
            )
        )
        + ".\n"
    )
    md.append(
        "Regenerate with `python make_report.py`; tests with "
        "`python -m pytest scripts/ma_reactions/kinetics_fit/tests`.\n"
    )
    (out / "REPORT.md").write_text("\n".join(md))
    cap_names = names + ("cap_n",)
    summary = {
        "candidates": {
            c["name"]: {q: float(c["v"][q]) for q in names} for c in cand_rows
        },
        "y_os_measured": mean["y_os"],
        "theta_fixedstir": th_fs,
        # These need a putidaFoam build that reads cap_n (section 6).
        "capped": {
            "F MEML with the cap": {
                q: float(c_fit.values[q]) for q in cap_names
            },
            "G EMM with the sharper cap": {
                q: float(e_cap.values[q]) for q in cap_names
            },
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    say("wrote", out / "REPORT.md")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true")
    run(ap.parse_args().quick)
