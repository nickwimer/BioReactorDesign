"""Before/after of the standard EMM and MEML runs on Kestrel.

Usage: python -I before_after.py <download root> <output dir>
Reads the liquidO2 and tankVolumes monitors of the old and new runs (as
downloaded under <download root>/scratch/nwimer/agent/accelerate/runs) and
writes do_time.png and a table of the tank-mean dissolved oxygen.
"""

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

root = Path(sys.argv[1]) / "scratch/nwimer/agent/accelerate/runs"
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
CSAT = 8.091  # mg/L, saturated tank mean of the spinup (spin0 at 60 s)
C = {
    "blue": "#2a78d6",
    "orange": "#eb6834",
    "aqua": "#1baf7a",
    "ink": "#0b0b0b",
    "ink2": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "surface": "#fcfcfb",
}


def monitor(case, name):
    """All rows of a monitor, restarts joined in time order."""
    rows = []
    base = root / case / "postProcessing" / name
    if not base.exists():
        return np.empty((0, 2))
    for d in sorted(base.iterdir(), key=lambda p: float(p.name)):
        for f in d.glob("*.dat"):
            for line in f.read_text().splitlines():
                if line and not line.startswith("#"):
                    rows.append([float(x) for x in line.split()])
    a = np.array(rows)
    return (
        a[np.argsort(a[:, 0], kind="stable")] if len(a) else np.empty((0, 2))
    )


def do_series(*cases):
    """(t, DO % of saturation) for cases chained in time."""
    ts, ds = [], []
    for case in cases:
        o2 = monitor(case, "liquidO2")
        vol = monitor(case, "tankVolumes")
        if not len(o2):
            continue
        v = np.interp(o2[:, 0], vol[:, 0], vol[:, 1])
        keep = o2[:, 0] > (ts[-1][-1] if ts else -1.0)
        ts.append(o2[keep, 0])
        ds.append(100.0 * 1e6 * o2[keep, 1] / v[keep] / CSAT)
    if not ts:
        return np.array([]), np.array([])
    return np.concatenate(ts), np.concatenate(ds)


panels = [
    (
        "EMM",
        ("old parameters", C["blue"], ("rx_emm", "rx_emm_k01")),
        ("new parameters (set G)", C["orange"], ("rx_emm_setG",)),
    ),
    (
        "MEML",
        ("old parameters", C["blue"], ("rx_meml", "rx_meml_k01")),
        ("new parameters (set F)", C["aqua"], ("rx_meml_setF",)),
    ),
]
# rx_meml (no floor) only supplies 60-60.2 s; rx_meml_k01 has the floor from
# 60.2 s. For EMM the floor was added at 100 s (rx_emm_k01).
series = {}
for title, *lines in panels:
    for label, col, cases in lines:
        if cases == ("rx_meml", "rx_meml_k01"):
            t0, d0 = do_series("rx_meml")
            m = t0 <= 60.2
            t1, d1 = do_series("rx_meml_k01")
            t, d = np.concatenate([t0[m], t1]), np.concatenate([d0[m], d1])
        else:
            t, d = do_series(*cases)
        series[(title, label)] = (t, d, col)

plt.rcParams.update(
    {
        "figure.facecolor": C["surface"],
        "axes.facecolor": C["surface"],
        "savefig.facecolor": C["surface"],
        "axes.edgecolor": C["axis"],
        "axes.labelcolor": C["ink2"],
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "xtick.color": C["muted"],
        "ytick.color": C["muted"],
        "axes.grid": True,
        "grid.color": C["grid"],
        "grid.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 2.0,
        "font.size": 10,
        "legend.frameon": False,
        "figure.dpi": 110,
    }
)
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
for ax, (title, *lines) in zip(axes, panels):
    for label, col, _ in lines:
        t, d, _ = series[(title, label)]
        if not len(t):
            continue
        ax.plot(t, d, color=col, label=label)
        ax.annotate(
            f"{d[-1]:.1f} %",
            (t[-1], d[-1]),
            textcoords="offset points",
            xytext=(6, 0),
            va="center",
            color=C["ink2"],
            fontsize=9,
        )
    ax.set_title(title, loc="left")
    ax.set_xlabel("time (s); biology update at 60 s")
    ax.set_xlim(58, 150)
    ax.set_ylim(0, 105)
    ax.legend(loc="upper right")
axes[0].set_ylabel("tank-mean dissolved O2 (% of saturation)")
axes[0].axvline(100, color=C["axis"], lw=1, ls="--")
axes[0].text(
    101,
    60,
    "old run: O2 floor\non the sink from here",
    color=C["muted"],
    fontsize=8.5,
    va="center",
)
fig.suptitle(
    "Dissolved oxygen after one 4 h biology update, 1000 L tank at 166 rpm, 410 SLPM",
    x=0.01,
    ha="left",
    fontsize=11,
    color=C["ink"],
)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(out / "do_time.png")
plt.close(fig)

lines_out = ["tank-mean dissolved O2, % of saturation (8.091 mg/L)", ""]
times = [60.2, 65, 70, 75, 80, 90, 100, 120, 140]
head = "case".ljust(34) + "".join(f"{t:>8g}" for t in times) + "    last (t)"
lines_out.append(head)
for (title, label), (t, d, _) in series.items():
    if not len(t):
        lines_out.append(f"{title} {label}".ljust(34) + "  no data yet")
        continue
    cells = [
        f"{np.interp(x, t, d):8.1f}" if t[0] <= x <= t[-1] else "       -"
        for x in times
    ]
    lines_out.append(
        f"{title} {label}".ljust(34)
        + "".join(cells)
        + f"    {d[-1]:.1f} ({t[-1]:.1f} s)"
    )
(out / "do_table.txt").write_text("\n".join(lines_out) + "\n")
print("\n".join(lines_out))
print("\nwrote", out / "do_time.png")
