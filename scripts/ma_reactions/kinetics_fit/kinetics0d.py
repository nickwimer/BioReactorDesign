"""0D P. putida kinetics as implemented in putidaFoam.

EMM is the Michaelis-Menten model (``get_rhs`` in
OFsolvers/birdmultiphaseEulerFoam/putidaFoam/microbeModel.C). MEML keeps the
same glucose and oxygen uptake and takes growth, H+ and muconate rates from
the ME-model surrogate (``get_rhs_ml``).

Units: time h, biomass X kg/m3, glucose G, muconate M and oxygen O2 mol/m3,
specific rates mol/kg/h (= mmol/gDCW/h).
"""

from dataclasses import asdict, dataclass, field, replace
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.spatial import ConvexHull

MW_GLUCOSE = 180.156  # g/mol
MW_MUCONIC_ACID = 142.11  # g/mol
MW_BIOMASS_CMOL = 24.6  # g per C-mol, CH1.8 O0.5 N0.2
GAMMA_GLUCOSE = 24.0  # electrons per mol
GAMMA_MUCONATE = 22.0  # C6H6O4
GAMMA_BIOMASS_CMOL = 4.2  # per C-mol, ammonia as N source


@dataclass(frozen=True)
class Params:
    """Kinetic parameters, named as in constant/microbeUpdateProperties."""

    q_max: float  # mol glucose / kg biomass / h
    K_o: float  # mol O2 / m3
    K_s: float  # mol glucose / m3
    y_xs: float  # kg biomass / mol glucose (EMM only)
    y_ms: float  # mol muconate / mol glucose (EMM only)
    y_os: float  # mol O2 / mol glucose
    x_max: float  # kg/m3, logistic cap (EMM only)

    def with_(self, **changes):
        return replace(self, **changes)

    def as_dict(self):
        return asdict(self)

    def to_foam(self, o2_max=0.25906, kLa=0.0, use_ml=False):
        """Entries for constant/microbeUpdateProperties.

        ``o2_max`` and ``kLa`` are read by putidaFoam but not used.
        """
        lines = [
            f"use_ml          {'true' if use_ml else 'false'};",
            f"y_xs            [ 0 0 0 0 0 0 0 ] {self.y_xs:.5g};",
            f"y_ms            [ 0 0 0 0 0 0 0 ] {self.y_ms:.5g};",
            f"y_os            [ 0 0 0 0 0 0 0 ] {self.y_os:.5g};",
            f"x_max           [ -1 0 0 0 1 0 0 ] {self.x_max:.5g};",
            f"q_max           [ -1 0 -1 0 1 0 0 ] {self.q_max:.5g};",
            f"o2_max          [ 0 -3 0 0 1 0 0 ] {o2_max:.5g};",
            f"K_o             [ 0 -3 0 0 1 0 0 ] {self.K_o:.5g};",
            f"K_s             [ 0 -3 0 0 1 0 0 ] {self.K_s:.5g};",
            f"kLa             [ 0 0 -1 0 0 0 0 ] {kLa:.5g};",
        ]
        return "\n".join(lines)


# Parameter set used in the putidaFoam cases (fit of 2026-03-30).
PRODUCTION = Params(
    q_max=8.26,
    K_o=0.1134,
    K_s=0.4008,
    y_xs=0.04498,
    y_ms=0.3189,
    y_os=1.195,
    x_max=7.6,
)


def uptake(G, O2, p, clamp=True):
    """Specific glucose uptake qs, mol/kg/h."""
    if clamp:
        G = np.maximum(G, 0.0)
        O2 = np.maximum(O2, 0.0)
    return p.q_max * G / (G + p.K_s) * O2 / (O2 + p.K_o)


class Surrogate:
    """ME-model surrogate as putidaFoam evaluates it.

    Reads the files a putidaFoam case directory holds: the TorchScript
    model, the four scaler text files and the Delaunay triangulation of the
    training inputs. Inputs are (glucose exchange, O2 exchange), negative for
    uptake. Outside the triangulation the solver returns zeros.

    The network is a small ReLU MLP. Its weights are read once through torch
    and evaluated with numpy in float64, which is fast and smooth enough for
    finite-difference fitting. ``torch_eval`` runs the original module.
    """

    def __init__(
        self, directory, model_file="scripted_model_biomass_h_muconate.pt"
    ):
        import torch

        d = Path(directory)
        self.directory = d
        self.x_mean = np.loadtxt(d / "scaler_x_mean.txt")
        self.x_scale = np.loadtxt(d / "scaler_x_scale.txt")
        self.y_mean = np.loadtxt(d / "scaler_y_mean.txt")
        self.y_scale = np.loadtxt(d / "scaler_y_scale.txt")
        self.points = np.loadtxt(d / "tri_points.txt")
        self.simplices = np.loadtxt(d / "tri_simplices.txt").astype(int)
        self._torch = torch
        self._module = torch.jit.load(str(d / model_file), map_location="cpu")
        self._module.eval()
        tensors = [
            t.detach().double().numpy()
            for _, t in self._module.named_parameters()
        ]
        self.layers = list(zip(tensors[0::2], tensors[1::2]))
        for w, b in self.layers:
            if w.ndim != 2 or b.shape != (w.shape[0],):
                raise ValueError("surrogate is not a plain MLP")
        hull = ConvexHull(self.points)
        self.hull_normals = hull.equations[:, :2]
        self.hull_offsets = hull.equations[:, 2]
        self.hull_vertices = self.points[hull.vertices]
        check = self.points[:: max(1, len(self.points) // 200)]
        err = np.max(np.abs(self._raw(check) - self.torch_eval(check)))
        if err > 1e-4:
            raise ValueError(f"numpy and torch surrogate differ by {err}")

    def _raw(self, inputs):
        """Network output in physical units, no hull test."""
        h = (np.atleast_2d(inputs) - self.x_mean) / self.x_scale
        for w, b in self.layers[:-1]:
            h = np.maximum(h @ w.T + b, 0.0)
        w, b = self.layers[-1]
        return (h @ w.T + b) * self.y_scale + self.y_mean

    def torch_eval(self, inputs, dtype="float32"):
        """Original TorchScript module (float32 is what the solver uses)."""
        torch = self._torch
        x = (np.atleast_2d(inputs) - self.x_mean) / self.x_scale
        module = self._module.double() if dtype == "float64" else self._module
        with torch.no_grad():
            y = module(torch.tensor(x, dtype=getattr(torch, dtype)))
        if dtype == "float64":
            self._module.float()
        return y.numpy().astype(float) * self.y_scale + self.y_mean

    def signed_distance(self, inputs):
        """Distance to the hull boundary; negative or zero inside."""
        pts = np.atleast_2d(inputs)
        return np.max(pts @ self.hull_normals.T + self.hull_offsets, axis=1)

    def inside(self, inputs, tol=1e-12):
        return self.signed_distance(inputs) <= tol

    def inside_simplices(self, inputs):
        """The solver's own test: inside any triangle. Slow; for checks."""
        a = self.points[self.simplices[:, 0]]
        b = self.points[self.simplices[:, 1]]
        c = self.points[self.simplices[:, 2]]
        det = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (
            b[:, 1] - a[:, 1]
        ) * (c[:, 0] - a[:, 0])
        out = []
        with np.errstate(divide="ignore", invalid="ignore"):
            for p in np.atleast_2d(inputs):
                al = (
                    (b[:, 0] - p[0]) * (c[:, 1] - p[1])
                    - (b[:, 1] - p[1]) * (c[:, 0] - p[0])
                ) / det
                be = (
                    (c[:, 0] - p[0]) * (a[:, 1] - p[1])
                    - (c[:, 1] - p[1]) * (a[:, 0] - p[0])
                ) / det
                out.append(
                    bool(np.any((al >= 0) & (be >= 0) & (1 - al - be >= 0)))
                )
        return np.array(out)

    def __call__(self, qs, our_spec):
        """(mu 1/h, r_H, r_muc mol/kg/h) for specific uptakes in mol/kg/h.

        ``qs`` and ``our_spec`` are positive uptakes; the model is fed
        (-qs, -our_spec) as in ``get_rhs_ml``. Zeros outside the hull.
        """
        qs = np.asarray(qs, dtype=float)
        shape = qs.shape
        inputs = np.column_stack(
            [-qs.ravel(), -np.asarray(our_spec, dtype=float).ravel()]
        )
        out = self._raw(inputs)
        out[~self.inside(inputs)] = 0.0
        return out.reshape(shape + (3,))

    def point(self, qs, our_spec):
        """Scalar version of ``__call__`` for use inside an ODE."""
        x = np.array([-qs, -our_spec])
        if np.max(self.hull_normals @ x + self.hull_offsets) > 1e-12:
            return 0.0, 0.0, 0.0
        h = (x - self.x_mean) / self.x_scale
        for w, b in self.layers[:-1]:
            h = np.maximum(w @ h + b, 0.0)
        w, b = self.layers[-1]
        out = (w @ h + b) * self.y_scale + self.y_mean
        return float(out[0]), float(out[1]), float(out[2])

    def qs_range(self, y_os):
        """Interval of qs inside the hull along O2 uptake = y_os * qs."""
        qs = np.linspace(0.0, 12.0, 4801)
        ok = self.inside(np.column_stack([-qs, -y_os * qs]))
        if not ok.any():
            return (np.nan, np.nan)
        return (float(qs[ok][0]), float(qs[ok][-1]))


def rates(X, G, O2, p, model="emm", sur=None, clamp=True):
    """Return (dX, dG, dM, dH, OUR) in volumetric units per hour.

    OUR is the oxygen uptake rate, mol/m3/h; putidaFoam solves O2 in the
    CFD, so it is returned separately instead of as an O2 derivative.
    """
    qs = uptake(G, O2, p, clamp)
    our = p.y_os * qs * X
    if model == "emm":
        dX = p.y_xs * qs * X * (1.0 - X / p.x_max)
        dM = p.y_ms * qs * X
        dH = np.zeros_like(dX)
    elif model == "meml":
        out = sur(qs, p.y_os * qs)
        dX = out[..., 0] * X
        dH = out[..., 1] * X
        dM = out[..., 2] * X
    else:
        raise ValueError(f"unknown model {model!r}")
    return dX, -qs * X, dM, dH, our


def advance_solver(
    state, p, t_adv=4.0, dt=0.01, t_now=0.0, model="emm", sur=None
):
    """Copy of ``microbemodel::advance`` plus the OUR average of bioReact.H.

    ``state`` is (X, O2, G, M, H). O2 stays frozen. The loop counter
    accumulates ``dt`` exactly as the C++ does, so 4 h at dt 0.01 takes 401
    steps. Returns (new state, mean of OUR before and after, steps).
    """
    X, O2, G, M, H = (float(v) for v in state)

    def f(y):
        dX, dG, dM, dH, _ = rates(y[0], y[1], O2, p, model, sur, clamp=False)
        return np.array([dX, dG, dM, dH], dtype=float).reshape(4)

    def our(y):
        return float(p.y_os * uptake(y[1], O2, p, clamp=False) * y[0])

    y = np.array([X, G, M, H])
    our_before = our(y)
    current = t_now
    final = t_now + t_adv
    steps = 0
    while current < final:
        current += dt
        y_mid = y + 0.5 * dt * f(y)
        y = y + dt * f(y_mid)
        steps += 1
    new_state = (y[0], O2, y[1], y[2], y[3])
    return new_state, 0.5 * (our_before + our(y)), steps


def o2_qss(X, G, p, theta, o2sat):
    """Quasi-steady dissolved O2 for kLa (o2sat - O2) = y_os qs X.

    ``theta`` = y_os / kLa in hours: dissolved-O2 data only determine this
    ratio, so it is the parameter that is fitted.
    """
    G = np.maximum(G, 0.0)
    a = theta * p.q_max * G / (G + p.K_s) * X
    b = p.K_o + a - o2sat
    disc = np.sqrt(b * b + 4.0 * o2sat * p.K_o)
    with np.errstate(divide="ignore", invalid="ignore"):
        pos = 2.0 * o2sat * p.K_o / (b + disc)
    return np.where(b >= 0, pos, 0.5 * (disc - b))


@dataclass
class Dataset:
    """One batch run.

    Fixed-stir runs carry dissolved O2 at the sample times as a fraction of
    saturation (``do``). Exp 17 runs carry a dissolved-O2 series
    (``do_t``, ``do_frac``) that is prescribed, and the measured oxygen
    uptake rate (``our_t``, ``our``, mol/m3/h).
    """

    name: str
    group: str
    t: np.ndarray
    X: np.ndarray
    G: np.ndarray
    M: np.ndarray
    o2sat: float
    do: np.ndarray = None
    do_t: np.ndarray = None
    do_frac: np.ndarray = None
    our_t: np.ndarray = None
    our: np.ndarray = None
    notes: dict = field(default_factory=dict)

    @property
    def o2_mode(self):
        return "prescribed" if self.do_t is not None else "qss"


def _scalar_rhs(p, ds, model, sur, theta):
    """Right-hand side on plain floats, for speed inside the integrator."""
    q_max, K_o, K_s = float(p.q_max), float(p.K_o), float(p.K_s)
    y_xs, y_ms, y_os = float(p.y_xs), float(p.y_ms), float(p.y_os)
    x_max, o2sat = float(p.x_max), float(ds.o2sat)
    prescribed = ds.o2_mode == "prescribed"
    do_t = ds.do_t
    do_c = None if not prescribed else ds.do_frac * o2sat
    emm = model == "emm"

    def f(y, t):
        X = y[0] if y[0] > 1e-12 else 1e-12
        G = y[1] if y[1] > 0.0 else 0.0
        fs = G / (G + K_s)
        if prescribed:
            O2 = float(np.interp(t, do_t, do_c))
        else:
            b = K_o + theta * q_max * fs * X - o2sat
            disc = sqrt(b * b + 4.0 * o2sat * K_o)
            O2 = 2.0 * o2sat * K_o / (b + disc) if b >= 0 else 0.5 * (disc - b)
        qs = q_max * fs * O2 / (O2 + K_o)
        if emm:
            return (
                y_xs * qs * X * (1.0 - X / x_max),
                -qs * X,
                y_ms * qs * X,
            )
        mu, _, r_muc = sur.point(qs, y_os * qs)
        return (mu * X, -qs * X, r_muc * X)

    return f


def simulate(
    p, ds, model="emm", sur=None, theta=None, t_dense=None, rtol=1e-8
):
    """Integrate one dataset from its first sample.

    Dissolved O2 is either prescribed (Exp 17) or quasi-steady with
    ``theta`` = y_os / kLa (fixed-stir). Returns a dict with X, G, M, O2,
    qs and OUR at ``ds.t`` (or ``t_dense``), OUR at ``ds.our_t`` when
    present, and for MEML the largest signed distance to the surrogate hull
    along the trajectory and the time spent outside it.
    """
    t_out = ds.t if t_dense is None else np.asarray(t_dense, dtype=float)
    t_eval = t_out
    if ds.our_t is not None and t_dense is None:
        t_eval = np.union1d(t_out, ds.our_t)
    t0 = ds.t[0]
    t_eval = t_eval[t_eval >= t0]
    if model == "meml" and sur is None:
        raise ValueError("MEML needs a Surrogate")
    if ds.o2_mode == "qss" and theta is None:
        raise ValueError(f"{ds.name} needs theta = y_os / kLa")

    y, info = odeint(
        _scalar_rhs(p, ds, model, sur, theta),
        [ds.X[0], ds.G[0], ds.M[0]],
        t_eval,
        rtol=rtol,
        atol=1e-10,
        hmax=0.25 if ds.o2_mode == "prescribed" else 0.0,
        mxstep=50000,
        full_output=True,
    )
    if info["message"] != "Integration successful." or not np.all(
        np.isfinite(y)
    ):
        raise RuntimeError(f"integration failed for {ds.name}")
    X = np.maximum(y[:, 0], 1e-12)
    G = np.maximum(y[:, 1], 0.0)
    M = y[:, 2]
    if ds.o2_mode == "prescribed":
        O2 = ds.o2sat * np.interp(t_eval, ds.do_t, ds.do_frac)
    else:
        O2 = o2_qss(X, G, p, theta, ds.o2sat)
    qs = uptake(G, O2, p)
    our = p.y_os * qs * X
    pick = np.searchsorted(t_eval, t_out)
    out = {
        "t": t_out,
        "X": X[pick],
        "G": G[pick],
        "M": M[pick],
        "O2": O2[pick],
        "qs": qs[pick],
        "OUR": our[pick],
    }
    if ds.our_t is not None and t_dense is None:
        out["OUR_at_our_t"] = our[np.searchsorted(t_eval, ds.our_t)]
    if model == "meml":
        dist = sur.signed_distance(np.column_stack([-qs, -p.y_os * qs]))
        live = G > 1e-6
        out["hull_distance"] = float(np.max(dist[live])) if live.any() else 0.0
        out["t_outside"] = float(
            np.sum(np.diff(t_eval)[(dist[1:] > 0) & live[1:]])
        )
    return out


def simulate_o2_ode(p, kLa, o2sat, t, y0, rtol=1e-3, atol=1e-6):
    """EMM with the dissolved-O2 balance as a fourth, stiff equation.

    This is the model of the original fit scripts, integrated with Radau at
    scipy's default tolerances as they did. ``y0`` is (X, G, M, O2).
    """

    def f(_, y):
        X, G, M, O2 = y
        dX, dG, dM, _, our = rates(X, G, O2, p, "emm", clamp=False)
        return [dX, dG, dM, kLa * (o2sat - O2) - our]

    sol = solve_ivp(
        f, (t[0], t[-1]), y0, method="Radau", t_eval=t, rtol=rtol, atol=atol
    )
    return sol.y.T


def electron_balance(y_x, y_ms, y_os=None):
    """Degree-of-reduction and carbon balance per mol of glucose.

    ``y_x`` is the realised biomass yield in kg/mol (for EMM that is y_xs
    times the logistic factor, or simply dX / -dG from data). Returns the
    O2 yield that closes the electron balance when everything else leaves as
    CO2, the matching CO2 yield and RQ and, if ``y_os`` is given, the
    electrons left unaccounted for.
    """
    cmol_x = y_x * 1000.0 / MW_BIOMASS_CMOL
    e_products = GAMMA_MUCONATE * y_ms + GAMMA_BIOMASS_CMOL * cmol_x
    y_os_closed = (GAMMA_GLUCOSE - e_products) / 4.0
    co2 = 6.0 - 6.0 * y_ms - cmol_x
    out = {
        "y_os_closed": y_os_closed,
        "co2_closed": co2,
        "rq_closed": co2 / y_os_closed,
    }
    if y_os is not None:
        missing = GAMMA_GLUCOSE - e_products - 4.0 * y_os
        out["electrons_unaccounted"] = missing
        out["electron_recovery"] = 1.0 - missing / GAMMA_GLUCOSE
    return out


def load_datasets(yaml_path, names=None, local_dir=None):
    """Read datasets.yaml. Exp 17 runs need the extracted CSVs."""
    from ruamel.yaml import YAML

    yaml_path = Path(yaml_path)
    cfg = YAML(typ="safe").load(yaml_path.read_text())
    base = yaml_path.parent
    o2sat = float(cfg["o2sat_lab"])
    out = []
    for name, spec in cfg["datasets"].items():
        if names is not None and name not in names:
            continue
        if spec["group"] == "fixedstir":
            d = np.loadtxt(base / spec["file"])
            out.append(
                Dataset(
                    name=name,
                    group="fixedstir",
                    t=d[:, 0],
                    X=d[:, 1],
                    G=d[:, 2],
                    M=d[:, 3],
                    do=d[:, 4] / float(spec["do_basis"]),
                    o2sat=o2sat,
                    notes=dict(spec.get("notes", {})),
                )
            )
        elif spec["group"] == "exp17":
            if local_dir is None:
                raise ValueError("Exp 17 datasets need local_dir")
            from exp17 import load_reactor

            out.append(
                load_reactor(Path(local_dir), int(spec["reactor"]), o2sat)
            )
        else:
            raise ValueError(f"unknown group {spec['group']!r}")
    return out
