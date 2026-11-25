"""Fit the putidaFoam kinetic parameters to 0D batch data.

The same code fits EMM and MEML to the same data with the same weighting.
Dissolved-O2 data only determine theta = y_os / kLa, so theta is fitted and
any kLa or y_os assumption is arithmetic afterwards (kLa = y_os / theta).

Example::

    python fit_kinetics.py --datasets fixedstir1 fixedstir2 fixedstir3 \\
        --free q_max K_o x_max y_xs y_ms theta --fix K_s=0.4 y_os=2.4
"""

import argparse
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from kinetics0d import (
    MEML_MODELS,
    Params,
    Surrogate,
    load_datasets,
    simulate,
    simulate_o2_ode,
)
from scipy.optimize import least_squares, minimize

HERE = Path(__file__).resolve().parent
KINETIC = ("q_max", "K_o", "K_s", "y_xs", "y_ms", "y_os", "x_max")
NAMES = KINETIC + ("theta",)
# Optional in putidaFoam; fixed at these values unless a fit frees them.
OPTIONAL = {"cap_n": 1.0}
BOUNDS = {
    "cap_n": (1.0, 30.0),
    "q_max": (0.5, 40.0),
    "K_o": (1e-5, 1.0),
    "K_s": (1e-3, 30.0),
    "y_xs": (0.005, 0.2),
    "y_ms": (0.05, 1.0),
    "y_os": (0.5, 4.0),
    "x_max": (7.0, 12.0),
    "theta": (1e-4, 0.1),
}
# Starting point of every search: the production set, theta for kLa 336.
START = {
    "q_max": 8.26,
    "K_o": 0.1134,
    "K_s": 0.4008,
    "y_xs": 0.04498,
    "y_ms": 0.3189,
    "y_os": 1.195,
    "x_max": 7.6,
    "theta": 1.195 / 336.0,
}
START.update(OPTIONAL)
# Floors for the noise levels, in data units.
SIGMA_FLOOR = {"X": 0.1, "G": 5.0, "M": 2.0, "DO": 0.03, "OUR": 2.0}


def params_of(values):
    """Params from a dict of fitted values (theta and extras ignored)."""
    kw = {k: values[k] for k in KINETIC}
    kw.update({k: values.get(k, d) for k, d in OPTIONAL.items()})
    return Params(**kw)


def local_config(path=None):
    """Read local.yaml (paths to files that are not in the repository)."""
    path = Path(path) if path else HERE / "local.yaml"
    if not path.exists():
        return {}
    from ruamel.yaml import YAML

    return YAML(typ="safe").load(path.read_text()) or {}


def noise_levels(datasets):
    """Noise level per group and quantity from the replicate scatter.

    For each quantity: standard deviation across the replicates of a group
    at each sample time, root-mean-squared over the times after t = 0.
    Dissolved O2 flagged as suspect is left out.
    """
    out = {}
    for group in sorted({d.group for d in datasets}):
        members = [d for d in datasets if d.group == group]
        level = dict(SIGMA_FLOOR)
        if len(members) > 1:
            series = {
                "X": [d.X for d in members],
                "G": [d.G for d in members],
                "M": [d.M for d in members],
            }
            do = [
                d.do
                for d in members
                if d.do is not None and not d.notes.get("do_suspect")
            ]
            if len(do) > 1:
                series["DO"] = do
            our = [d.our for d in members if d.our is not None]
            if len(our) > 1 and len({len(o) for o in our}) == 1:
                series["OUR"] = our
            for key, rows in series.items():
                sd = np.std(np.array(rows), axis=0, ddof=1)
                first = 0 if key in ("DO", "OUR") else 1
                rms = float(np.sqrt(np.mean(sd[first:] ** 2)))
                level[key] = max(rms, SIGMA_FLOOR[key])
        out[group] = level
    return out


def meml_bounds(sur):
    """Search bounds of q_max and y_os for MEML, from the surrogate.

    They keep a fit inside the range of the surrogate's training inputs:
    q_max below the largest glucose uptake, y_os within the central 98 %
    of the O2/glucose ratios.
    """
    q_hi = float(np.floor(-10.0 * sur.points[:, 0].min()) / 10.0)
    lo, hi = np.percentile(sur.points[:, 1] / sur.points[:, 0], [1, 99])
    y_lo = float(np.ceil(10.0 * lo) / 10.0)
    y_hi = float(np.floor(10.0 * hi) / 10.0)
    return (1.0, q_hi), (y_lo, y_hi)


@dataclass
class Problem:
    datasets: list
    free: tuple
    fixed: dict
    model: str = "emm"
    sur: Surrogate = None
    sigma: dict = None
    use_do: bool = True
    use_our: bool = True
    hull_weight: float = 50.0
    bounds: dict = field(default_factory=lambda: dict(BOUNDS))

    def __post_init__(self):
        self.free = tuple(self.free)
        missing = set(NAMES) - set(self.free) - set(self.fixed)
        if missing:
            raise ValueError(f"neither free nor fixed: {sorted(missing)}")
        defaults = {k: v for k, v in OPTIONAL.items() if k not in self.free}
        self.fixed = {**defaults, **self.fixed}
        if self.sigma is None:
            self.sigma = noise_levels(self.datasets)
        if self.model in MEML_MODELS:
            q_b, y_b = meml_bounds(self.sur)
            self.bounds["q_max"], self.bounds["y_os"] = q_b, y_b

    def values(self, x_log):
        v = dict(self.fixed)
        v.update(zip(self.free, np.exp(x_log)))
        return v

    def split(self, values):
        return params_of(values), values["theta"]

    def blocks(self, values):
        """Scaled residual blocks, keyed by (dataset, quantity)."""
        p, theta = self.split(values)
        out = {}
        for ds in self.datasets:
            sg = self.sigma[ds.group]
            sim = simulate(p, ds, self.model, self.sur, theta)
            for key in ("X", "G", "M"):
                obs = getattr(ds, key)[1:]
                ok = np.isfinite(obs)
                out[(ds.name, key)] = (sim[key][1:] - obs)[ok] / sg[key]
            if (
                self.use_do
                and ds.do is not None
                and not ds.notes.get("do_suspect")
            ):
                out[(ds.name, "DO")] = (sim["O2"] / ds.o2sat - ds.do) / sg[
                    "DO"
                ]
            if self.use_our and ds.our is not None:
                out[(ds.name, "OUR")] = (sim["OUR_at_our_t"] - ds.our) / sg[
                    "OUR"
                ]
            if self.model in MEML_MODELS:
                out[(ds.name, "hull")] = np.array(
                    [self.hull_weight * max(0.0, sim["hull_distance"])]
                )
        return out

    def residuals(self, x_log):
        try:
            blocks = self.blocks(self.values(x_log))
        except (RuntimeError, FloatingPointError, ZeroDivisionError):
            return np.full(self.n_residuals, 1e3)
        return np.concatenate(list(blocks.values()))

    @property
    def n_residuals(self):
        if not hasattr(self, "_n"):
            start = {**START, **self.fixed}
            self._n = sum(len(b) for b in self.blocks(start).values())
        return self._n

    def log_bounds(self):
        lo = np.log([self.bounds[n][0] for n in self.free])
        hi = np.log([self.bounds[n][1] for n in self.free])
        return lo, hi


@dataclass
class FitResult:
    values: dict
    free: tuple
    chi2: float
    n_residuals: int
    rmse: dict  # (dataset, quantity) -> RMSE in units of the noise level
    start_chi2: list
    singular_values: list
    correlation: list
    at_bound: list
    t_outside: float = 0.0

    @property
    def params(self):
        return params_of(self.values)

    def summary(self):
        return {
            "values": {k: float(v) for k, v in self.values.items()},
            "free": list(self.free),
            "chi2": self.chi2,
            "n_residuals": self.n_residuals,
            "at_bound": self.at_bound,
            "start_chi2": [float(c) for c in self.start_chi2],
            "singular_values": [float(s) for s in self.singular_values],
            "t_outside": self.t_outside,
        }


def _solve(problem, x0):
    lo, hi = problem.log_bounds()
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
    return least_squares(
        problem.residuals,
        x0,
        jac="3-point",
        bounds=(lo, hi),
        method="trf",
        x_scale="jac",
        diff_step=1e-3,
        max_nfev=120,
    )


def fit(problem, start=None, n_start=20, seed=0):
    """Multi-start least squares in log-parameters."""
    lo, hi = problem.log_bounds()
    base = {**START, **(start or {})}
    x_base = np.log([base[n] for n in problem.free])
    rng = np.random.default_rng(seed)
    starts = [x_base]
    for _ in range(max(0, n_start - 1)):
        mix = rng.uniform(0.1, 0.9, size=len(lo))
        starts.append(lo + mix * (hi - lo))
    sols = [_solve(problem, x0) for x0 in starts]
    best = min(sols, key=lambda s: s.cost)
    values = problem.values(best.x)
    blocks = problem.blocks(values)
    rmse = {k: float(np.sqrt(np.mean(b**2))) for k, b in blocks.items()}
    sv = np.linalg.svd(best.jac, compute_uv=False)
    corr = []
    if sv[-1] > 1e-8 * sv[0]:
        cov = np.linalg.inv(best.jac.T @ best.jac)
        d = np.sqrt(np.diag(cov))
        corr = (cov / np.outer(d, d)).tolist()
    at_bound = [
        n
        for n, x, a, b in zip(problem.free, best.x, lo, hi)
        if x - a < 0.05 or b - x < 0.05  # within 5 % of a bound
    ]
    t_out = 0.0
    if problem.model in MEML_MODELS:
        p, theta = problem.split(values)
        t_out = sum(
            simulate(p, ds, problem.model, problem.sur, theta)["t_outside"]
            for ds in problem.datasets
        )
    return FitResult(
        values=values,
        free=problem.free,
        chi2=float(2.0 * best.cost),
        n_residuals=problem.n_residuals,
        rmse=rmse,
        start_chi2=sorted(2.0 * s.cost for s in sols),
        singular_values=sv.tolist(),
        correlation=corr,
        at_bound=at_bound,
        t_outside=t_out,
    )


def profile(problem, result, name, grid):
    """Profile of chi2 over one parameter, the others refitted.

    Each grid point starts from its neighbour's solution, walking outwards
    from the optimum, which is added to the grid. Returns a list of
    (value, chi2, values) sorted by value.
    """
    centre = result.values[name] if name in result.values else grid[0]
    grid = np.asarray(sorted({*map(float, grid), float(centre)}))
    others = tuple(n for n in problem.free if n != name)
    order = np.argsort(np.abs(np.log(grid / centre)))
    done = {}
    for i in order:
        near = [j for j in done if abs(j - i) == 1]
        warm = done[near[0]][2] if near else result.values
        sub = Problem(
            problem.datasets,
            others,
            {**problem.fixed, name: float(grid[i])},
            problem.model,
            problem.sur,
            problem.sigma,
            problem.use_do,
            problem.use_our,
            problem.hull_weight,
        )
        sol = _solve(sub, np.log([warm[n] for n in others]))
        alt = _solve(sub, np.log([result.values[n] for n in others]))
        if alt.cost < sol.cost:
            sol = alt
        done[i] = (float(grid[i]), float(2.0 * sol.cost), sub.values(sol.x))
    return [done[i] for i in range(len(grid))]


# ------------------------------------------------------------ legacy metric
def legacy_loss(p, ds, kLa, o2sat=None):
    """Loss of the original fit scripts for one fixed-stir dataset.

    Sum over biomass, glucose, muconate and dissolved O2 of
    MSE / mean(obs^2), with the stiff O2 balance integrated by Radau at
    default tolerances.
    """
    o2sat = ds.o2sat if o2sat is None else o2sat
    obs = np.column_stack([ds.X, ds.G, ds.M, ds.do * o2sat])
    sim = simulate_o2_ode(p, kLa, o2sat, ds.t, obs[0])
    if sim.shape != obs.shape:
        return 1e6
    return float(
        np.sum(np.mean((obs - sim) ** 2, axis=0) / np.mean(obs**2, axis=0))
    )


LEGACY_ORDER = ("q_max", "K_o", "y_xs", "y_ms", "y_os", "K_s")


def legacy_fit(ds, start, kLa, x_max, o2sat=None):
    """The original procedure: scipy ``minimize`` on the scalar loss."""

    def loss(v):
        p = Params(**dict(zip(LEGACY_ORDER, v)), x_max=x_max)
        return legacy_loss(p, ds, kLa, o2sat)

    res = minimize(loss, start, bounds=[(0, None)] * 6)
    return Params(**dict(zip(LEGACY_ORDER, res.x)), x_max=x_max), res.fun


# ---------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Fit putidaFoam kinetics to 0D batch data."
    )
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument(
        "--model",
        choices=("emm",) + MEML_MODELS,
        default="emm",
        help="meml_cap is MEML with the biomass cap (cap_n in putidaFoam)",
    )
    every = NAMES + tuple(OPTIONAL)
    ap.add_argument("--free", nargs="+", required=True, choices=every)
    ap.add_argument(
        "--fix", nargs="*", default=[], metavar="NAME=VALUE", help="fixed"
    )
    ap.add_argument("--starts", type=int, default=20)
    ap.add_argument("--profile", choices=every)
    ap.add_argument("--no-do", action="store_true")
    ap.add_argument("--no-our", action="store_true")
    ap.add_argument("--surrogate-dir")
    ap.add_argument("--exp17-dir")
    args = ap.parse_args(argv)

    cfg = local_config()
    exp17_dir = args.exp17_dir or cfg.get("exp17_dir")
    datasets = load_datasets(
        HERE / "datasets.yaml", args.datasets, local_dir=exp17_dir
    )
    sur = None
    if args.model in MEML_MODELS:
        sur = Surrogate(args.surrogate_dir or cfg["surrogate_dir"])
    fixed = {k: float(v) for k, v in (item.split("=") for item in args.fix)}
    for name in NAMES:
        if name not in args.free and name not in fixed:
            fixed[name] = START[name]
    problem = Problem(
        datasets,
        args.free,
        fixed,
        args.model,
        sur,
        use_do=not args.no_do,
        use_our=not args.no_our,
    )
    result = fit(problem, n_start=args.starts)
    print(json.dumps(result.summary(), indent=2))
    for key, value in result.rmse.items():
        print(f"  {key[0]:12s} {key[1]:4s} rmse / sigma = {value:.2f}")
    theta = result.values["theta"]
    y_os = result.values["y_os"]
    print(f"\ntheta = {theta:.5f} h -> kLa = {y_os / theta:.0f} 1/h at y_os")
    block = result.params.to_foam(model=args.model)
    print(block)
    if "cap_n" in block:
        print("// cap_n needs a putidaFoam build that reads it")
    if args.profile:
        lo, hi = problem.bounds[args.profile]
        grid = np.geomspace(lo, hi, 13)
        for value, chi2, _ in profile(problem, result, args.profile, grid):
            print(f"  {args.profile} = {value:10.4g}   chi2 = {chi2:.3f}")


if __name__ == "__main__":
    main()
