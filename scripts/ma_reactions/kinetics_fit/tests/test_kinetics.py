import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from fit_kinetics import LEGACY_ORDER, Problem, fit, legacy_loss, noise_levels
from kinetics0d import (
    PRODUCTION,
    Dataset,
    Params,
    advance_solver,
    electron_balance,
    load_datasets,
    o2_qss,
    rates,
    simulate,
    simulate_o2_ode,
    uptake,
)

HERE = Path(__file__).resolve().parent.parent
OLD = HERE.parent / "0dWillComp"
# First putidaFoam bio update of the reacting 1000 L runs (solver log).
REF_STATE = (7.33, 0.2528, 228.8, 0.0, 0.0)


@pytest.fixture(scope="module")
def fixedstir():
    names = ["fixedstir1", "fixedstir2", "fixedstir3"]
    return load_datasets(HERE / "datasets.yaml", names)


# ------------------------------------------------------------------- model
def test_emm_rates_by_hand():
    p = Params(
        q_max=8.0, K_o=0.1, K_s=1.0, y_xs=0.05, y_ms=0.3, y_os=2.0, x_max=10.0
    )
    dX, dG, dM, dH, our = rates(2.0, 9.0, 0.1, p)
    qs = 8.0 * 0.9 * 0.5
    assert uptake(9.0, 0.1, p) == pytest.approx(qs)
    assert dX == pytest.approx(0.05 * qs * 2.0 * 0.8)
    assert dG == pytest.approx(-qs * 2.0)
    assert dM == pytest.approx(0.3 * qs * 2.0)
    assert dH == 0.0
    assert our == pytest.approx(2.0 * qs * 2.0)


def test_advance_solver_emm_reference():
    state, our, steps = advance_solver(REF_STATE, PRODUCTION)
    assert steps == 401  # the C++ loop overshoots by one step
    assert state[0] == pytest.approx(7.5009, abs=2e-4)
    assert our == pytest.approx(50.31, abs=0.01)
    assert state[1] == REF_STATE[1]  # O2 is frozen


def test_o2_qss_closes_the_balance():
    rng = np.random.default_rng(1)
    p = PRODUCTION
    for _ in range(50):
        X, G = rng.uniform(0.3, 8.0), rng.uniform(0.1, 600.0)
        kLa, o2sat = rng.uniform(50.0, 1000.0), 0.188125
        c = float(o2_qss(X, G, p, p.y_os / kLa, o2sat))
        supply = kLa * (o2sat - c)
        demand = p.y_os * uptake(G, c, p) * X
        assert 0.0 < c < o2sat
        assert supply == pytest.approx(demand, rel=1e-9)


def test_qss_matches_stiff_o2_ode(fixedstir):
    ds = fixedstir[0]
    p, kLa = PRODUCTION, 336.0
    y0 = [ds.X[0], ds.G[0], ds.M[0], ds.do[0] * ds.o2sat]
    full = simulate_o2_ode(p, kLa, ds.o2sat, ds.t, y0, rtol=1e-10, atol=1e-12)
    qss = simulate(p, ds, theta=p.y_os / kLa)
    assert qss["X"] == pytest.approx(full[:, 0], rel=1e-3)
    assert qss["M"][1:] == pytest.approx(full[1:, 2], rel=1e-3)
    assert qss["G"] == pytest.approx(full[:, 1], rel=1e-3, abs=0.05)
    assert qss["O2"][1:] == pytest.approx(full[1:, 3], rel=1e-3)


def test_only_y_os_over_kla_matters(fixedstir):
    """Doubling kLa and y_os together leaves the trajectories unchanged."""
    ds = fixedstir[0]
    p = PRODUCTION
    y0 = [ds.X[0], ds.G[0], ds.M[0], ds.do[0] * ds.o2sat]
    kw = dict(rtol=1e-10, atol=1e-12)
    a = simulate_o2_ode(p, 336.0, ds.o2sat, ds.t, y0, **kw)
    b = simulate_o2_ode(
        p.with_(y_os=2 * p.y_os), 672.0, ds.o2sat, ds.t, y0, **kw
    )
    # glucose is nearly exhausted at the last sample: allow 0.02 mol/m3
    assert a[1:] == pytest.approx(b[1:], rel=2e-3, abs=0.02)


def test_electron_balance_arithmetic():
    full = electron_balance(0.0, 0.0)
    assert full["y_os_closed"] == pytest.approx(6.0)
    assert full["co2_closed"] == pytest.approx(6.0)
    assert full["rq_closed"] == pytest.approx(1.0)
    b = electron_balance(0.0246, 0.5, y_os=2.0)
    # 1 C-mol biomass (4.2 e) and 0.5 mol muconate (11 e) per mol glucose
    assert b["y_os_closed"] == pytest.approx((24 - 11 - 4.2) / 4)
    assert b["co2_closed"] == pytest.approx(6 - 3 - 1)
    assert b["electrons_unaccounted"] == pytest.approx(24 - 15.2 - 8.0)


def test_to_foam_block():
    text = PRODUCTION.to_foam(kLa=336, use_ml=True)
    assert "use_ml          true;" in text
    assert "q_max           [ -1 0 -1 0 1 0 0 ] 8.26;" in text
    assert "K_o             [ 0 -3 0 0 1 0 0 ] 0.1134;" in text
    assert "x_max           [ -1 0 0 0 1 0 0 ] 7.6;" in text
    # cap_n is written only where it differs from the solver's default.
    assert "cap_n" not in text
    assert "cap_n" not in PRODUCTION.to_foam()
    sharp = PRODUCTION.with_(cap_n=5.0)
    assert "cap_n           5;" in sharp.to_foam()
    assert "cap_n           5;" in sharp.to_foam(model="meml_cap")
    assert "use_ml          true;" in sharp.to_foam(model="meml_cap")
    assert "cap_n" not in sharp.to_foam(model="meml")
    assert "cap_n           1;" in PRODUCTION.to_foam(model="meml_cap")


# -------------------------------------------------------------------- data
def test_fixedstir_do_is_a_fraction(fixedstir):
    for ds in fixedstir:
        assert 0.94 < ds.do[0] < 0.97
        assert ds.o2sat == pytest.approx(0.188125)
    # Documented data problem: BR2 carries BR1's DO readings.
    assert fixedstir[1].do == pytest.approx(fixedstir[0].do, abs=1e-4)
    assert fixedstir[1].notes["do_suspect"]


def test_stored_sets_reproduce_their_solutions(fixedstir):
    """Both stored parameter sets regenerate their soln.dat and loss."""
    ds = fixedstir[0]
    cases = [
        (OLD / "gDCW_conv" / "fixedstir1", 540.0, 7.65, 0.0068067),
        (OLD / "production_kla336" / "fixedstir1", 336.0, 7.6, 0.007114),
    ]
    for folder, kLa, x_max, loss in cases:
        v = np.loadtxt(folder / "fitted_param.dat")
        p = Params(**dict(zip(LEGACY_ORDER, v)), x_max=x_max)
        assert legacy_loss(p, ds, kLa) == pytest.approx(loss, rel=5e-3)
        stored = np.loadtxt(folder / "soln.dat")
        y0 = [ds.X[0], ds.G[0], ds.M[0], ds.do[0] * ds.o2sat]
        sim = simulate_o2_ode(p, kLa, ds.o2sat, ds.t, y0)
        assert sim[:, 0] == pytest.approx(stored[:, 1], rel=5e-3)
        assert sim[:, 2] == pytest.approx(stored[:, 3], rel=5e-3, abs=0.05)
        assert sim[:, 1] == pytest.approx(stored[:, 2], rel=5e-3, abs=1.0)


def test_production_set_is_not_a_kla_428_fit(fixedstir):
    v = np.loadtxt(
        OLD / "production_kla336" / "fixedstir1" / "fitted_param.dat"
    )
    p = Params(**dict(zip(LEGACY_ORDER, v)), x_max=7.6)
    assert legacy_loss(p, fixedstir[0], 428.0) > 10 * legacy_loss(
        p, fixedstir[0], 336.0
    )


# --------------------------------------------------------------------- fit
def test_fit_recovers_synthetic_parameters():
    truth = Params(
        q_max=4.5,
        K_o=0.002,
        K_s=0.4,
        y_xs=0.05,
        y_ms=0.33,
        y_os=2.4,
        x_max=7.8,
    )
    theta = 2.4 / 650.0
    t = np.array([0.0, 6.0, 12.0, 18.0, 23.0, 27.0, 30.0, 32.0])
    data = []
    for i, x0 in enumerate((0.45, 0.5, 0.55)):
        seed = Dataset(
            name=f"syn{i}",
            group="syn",
            t=t,
            X=np.full(len(t), x0),
            G=np.full(len(t), 590.0),
            M=np.full(len(t), 0.2),
            o2sat=0.188125,
        )
        sim = simulate(truth, seed, theta=theta)
        data.append(
            Dataset(
                name=f"syn{i}",
                group="syn",
                t=t,
                X=sim["X"],
                G=sim["G"],
                M=sim["M"],
                do=sim["O2"] / 0.188125,
                o2sat=0.188125,
            )
        )
    sigma = {"syn": {"X": 0.1, "G": 5.0, "M": 2.0, "DO": 0.03, "OUR": 2.0}}
    free = ("q_max", "y_xs", "y_ms", "x_max", "theta")
    fixed = {"K_o": 0.002, "K_s": 0.4, "y_os": 2.4}
    result = fit(Problem(data, free, fixed, sigma=sigma), n_start=3)
    assert result.chi2 < 1e-6
    for name in ("q_max", "y_xs", "y_ms", "x_max"):
        assert result.values[name] == pytest.approx(
            getattr(truth, name), rel=1e-3
        )
    assert result.values["theta"] == pytest.approx(theta, rel=1e-3)


def test_new_search_beats_the_production_fit(fixedstir):
    """A real search improves on the production set in its own metric."""
    ds = fixedstir[:1]
    sigma = noise_levels(fixedstir)
    free = ("q_max", "K_o", "y_xs", "y_ms", "theta")
    fixed = {"K_s": 0.4008, "x_max": 7.6, "y_os": 1.195}
    result = fit(Problem(ds, free, fixed, sigma=sigma), n_start=4)
    kLa = 1.195 / result.values["theta"]
    assert legacy_loss(result.params, ds[0], kLa) < 0.0055
    assert result.values["q_max"] < 6.0


# --------------------------------------------------------------- surrogate
def test_surrogate_numpy_matches_torch(sur):
    rng = np.random.default_rng(0)
    pts = sur.points[rng.choice(len(sur.points), 400, replace=False)]
    assert sur._raw(pts) == pytest.approx(
        sur.torch_eval(pts, "float32"), abs=1e-5
    )
    assert sur._raw(pts) == pytest.approx(
        sur.torch_eval(pts, "float64"), abs=1e-9
    )


def test_hull_test_matches_the_solver(sur):
    rng = np.random.default_rng(0)
    box = np.column_stack(
        [rng.uniform(-10.5, 0.5, 1500), rng.uniform(-21.0, 0.5, 1500)]
    )
    v = sur.hull_vertices
    mid = 0.5 * (v + np.roll(v, -1, axis=0))
    centre = sur.points.mean(axis=0)
    near = np.vstack(
        [centre + (mid - centre) * s for s in (1 - 1e-6, 1 + 1e-6)]
    )
    pts = np.vstack([box, near])
    assert np.array_equal(sur.inside(pts), sur.inside_simplices(pts))
    assert not sur.inside([[4.0, 8.0]])[0]  # positive inputs are outside
    assert np.all(sur(4.0, 30.0) == 0.0)  # zeros outside


def test_surrogate_point_matches_array(sur):
    assert sur.point(4.0, 8.0) == pytest.approx(tuple(sur(4.0, 8.0)))
    assert sur.point(4.0, 30.0) == (0.0, 0.0, 0.0)  # outside the hull


def test_surrogate_expected_values(sur, expected):
    """Values noted in local.yaml for the surrogate files in use."""
    if "surrogate_growth" not in expected:
        pytest.skip("no expected surrogate values in local.yaml")
    for qs, our_spec, mu in expected["surrogate_growth"]:
        assert sur(qs, our_spec)[0] == pytest.approx(mu, abs=1e-3)
    lo, hi = sur.qs_range(PRODUCTION.y_os)
    assert lo == pytest.approx(expected["surrogate_qs_range"][0], abs=0.02)
    assert hi == pytest.approx(expected["surrogate_qs_range"][1], abs=0.05)


def test_meml_bounds_stay_inside_the_training_inputs(sur):
    from fit_kinetics import meml_bounds

    (q_lo, q_hi), (y_lo, y_hi) = meml_bounds(sur)
    ratio = sur.points[:, 1] / sur.points[:, 0]
    assert 0.0 < q_lo < q_hi <= -sur.points[:, 0].min()
    assert ratio.min() <= y_lo < y_hi <= ratio.max()


def test_advance_solver_meml_reference(sur):
    state, our, steps = advance_solver(
        REF_STATE, PRODUCTION, model="meml", sur=sur
    )
    assert steps == 401
    assert state[0] == pytest.approx(9.245, abs=2e-3)
    assert state[2] == pytest.approx(40.86, abs=0.02)
    assert our == pytest.approx(56.12, abs=0.02)
    # putidaFoam log of the same update: X 9.26, G 40.7, mean OUR 56.16
    assert state[0] == pytest.approx(9.26, rel=0.01)
    assert state[2] == pytest.approx(40.7, rel=0.01)
    assert our == pytest.approx(56.16, rel=0.01)


def test_reproduces_solver_wellmixed_output(sur, surrogate_dir):
    """wellMixed.H of putidaFoam, against the wellmixed.dat it wrote."""
    path = surrogate_dir / "wellmixed.dat"
    if not path.exists():
        pytest.skip("no wellmixed.dat next to the surrogate files")
    rows = [
        [float(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line and line[0].isdigit()
    ]
    ref = np.array([r[:6] for r in rows])  # t, X, O2, G, M, H
    state = tuple(ref[0, 1:6])
    state = (0.350147052723467, 0.214, 545.9946375, 0.157763, 0.0)
    now = 0.0
    for row in ref[1:]:
        now += 5.0
        state, _, _ = advance_solver(
            state, PRODUCTION, 5.0, 0.1, now, "meml", sur
        )
        assert now == row[0]
        assert np.array(state) == pytest.approx(row[1:6], rel=2e-5)


def test_meml_cap_scales_only_growth(sur):
    p = PRODUCTION.with_(q_max=5.0, K_o=0.003, y_os=2.4, x_max=8.0)
    free = rates(6.0, 200.0, 0.1, p, "meml", sur)
    capped = rates(6.0, 200.0, 0.1, p, "meml_cap", sur)
    assert free[0] > 0.0
    assert capped[0] == pytest.approx(free[0] * (1.0 - 6.0 / 8.0))
    assert [float(v) for v in capped[1:]] == [float(v) for v in free[1:]]
    at_cap = rates(8.0, 200.0, 0.1, p, "meml_cap", sur)
    assert at_cap[0] == pytest.approx(0.0, abs=1e-12)
    assert at_cap[1] < 0.0 and at_cap[2] > 0.0  # uptake and product go on
    # A sharper cap leaves growth far below x_max almost untouched.
    sharp = rates(6.0, 200.0, 0.1, p.with_(cap_n=4.0), "meml_cap", sur)
    assert sharp[0] == pytest.approx(free[0] * (1.0 - 0.75**4))
    # The same sharpness applies to EMM's cap; 1 is the solver's EMM.
    emm = rates(6.0, 200.0, 0.1, p, "emm")[0]
    emm_sharp = rates(6.0, 200.0, 0.1, p.with_(cap_n=4.0), "emm")[0]
    assert emm_sharp == pytest.approx(emm * (1.0 - 0.75**4) / 0.25)
    assert rates(6.0, 200.0, 0.1, p.with_(cap_n=1.0), "emm")[0] == emm


def test_meml_cap_holds_biomass_below_x_max(sur, fixedstir):
    ds = fixedstir[0]
    p = PRODUCTION.with_(q_max=5.0, K_o=0.003, K_s=4.0, y_os=2.4, x_max=8.0)
    theta = 2.4 / 700.0
    free = simulate(p, ds, "meml", sur, theta)
    capped = simulate(p, ds, "meml_cap", sur, theta)
    assert free["X"][-1] > p.x_max
    assert capped["X"][-1] < p.x_max
    assert np.all(np.diff(capped["X"]) >= 0.0)
    # With the cap far away the two models are the same.
    far = simulate(p.with_(x_max=1e9), ds, "meml_cap", sur, theta)
    assert far["X"] == pytest.approx(free["X"], rel=1e-6)
    assert far["M"] == pytest.approx(free["M"], rel=1e-6, abs=1e-6)


def test_unknown_model_is_refused(fixedstir):
    with pytest.raises(ValueError):
        simulate(PRODUCTION, fixedstir[0], "memlcap", theta=0.0036)


# ------------------------------------------------------------- solver code
SOLVER = HERE.parents[2] / "OFsolvers/birdmultiphaseEulerFoam/putidaFoam"
SURROGATE_FILES = (
    "scaler_x_mean.txt",
    "scaler_x_scale.txt",
    "scaler_y_mean.txt",
    "scaler_y_scale.txt",
    "tri_points.txt",
    "tri_simplices.txt",
    "scripted_model_biomass_h_muconate.pt",
)
CAPPED = PRODUCTION.with_(
    q_max=4.7, K_o=0.003125, K_s=13.5, y_os=2.5, x_max=7.8, cap_n=5.0
)


@pytest.fixture(scope="module")
def solver_driver(surrogate_dir, tmp_path_factory):
    """putidaFoam's microbeModel.C built with a small driver, no OpenFOAM."""
    torch = pytest.importorskip("torch")
    cxx = shutil.which(os.environ.get("CXX", "c++"))
    tdir = Path(torch.__file__).parent
    inc = tdir / "include"
    if cxx is None or not (inc / "torch" / "script.h").exists():
        pytest.skip("no C++ compiler or no libtorch headers")
    work = tmp_path_factory.mktemp("solver")
    # microbeModel.H includes OpenFOAM's scalar.H and uses nothing from it.
    (work / "scalar.H").write_text("\n")
    for name in SURROGATE_FILES:  # the solver loads them from its cwd
        shutil.copy(surrogate_dir / name, work / name)
    cmd = [
        cxx,
        "-std=c++17",
        "-O1",
        f"-I{work}",
        f"-I{SOLVER}",
        f"-I{inc}",
        f"-I{inc / 'torch' / 'csrc' / 'api' / 'include'}",
        str(HERE / "tests" / "solver_driver.cpp"),
        str(SOLVER / "microbeModel.C"),
        f"-L{tdir / 'lib'}",
        "-ltorch",
        "-ltorch_cpu",
        "-lc10",
        f"-Wl,-rpath,{tdir / 'lib'}",
        "-o",
        str(work / "driver"),
    ]
    built = subprocess.run(cmd, capture_output=True, text=True)
    assert built.returncode == 0, built.stderr[-4000:]
    return work / "driver"


@pytest.mark.parametrize(
    "use_ml, solver_cap_n, model, p",
    [
        # Without a cap_n entry putidaFoam passes 1 (EMM) or 0 (MEML).
        (0, 1.0, "emm", PRODUCTION),
        (1, 0.0, "meml", PRODUCTION),
        (1, 5.0, "meml_cap", CAPPED),
        (1, 1.0, "meml_cap", CAPPED.with_(cap_n=1.0)),
        (0, 5.0, "emm", CAPPED),
        # cap_n 0 switches the cap off; slow uptake so glucose lasts.
        (0, 0.0, "emm", CAPPED.with_(q_max=2.0, x_max=1e30, cap_n=1.0)),
    ],
)
def test_compiled_solver_matches_python(
    solver_driver, sur, use_ml, solver_cap_n, model, p
):
    """Seven 5 h updates as in wellMixed.H, solver code against Python."""
    state = (0.350147052723467, 0.214, 545.9946375, 0.157763, 0.0)
    values = [use_ml, solver_cap_n, p.y_xs, p.y_ms, p.y_os, p.x_max]
    values += [p.q_max, p.K_o, p.K_s, *state[:4], 5.0, 0.1]
    out = subprocess.run(
        [str(solver_driver)] + [repr(float(v)) for v in values] + ["7"],
        cwd=solver_driver.parent,
        capture_output=True,
        text=True,
        check=True,
    )
    solver = np.array([line.split() for line in out.stdout.splitlines()])
    now, rows = 0.0, []
    for _ in range(7):
        state, our, _ = advance_solver(state, p, 5.0, 0.1, now, model, sur)
        now += 5.0
        rows.append([now, *state, our])
    # EMM agrees to round-off; MEML runs the network in float32 in the
    # solver and in float64 here.
    assert solver.astype(float) == pytest.approx(
        np.array(rows), rel=1e-6, abs=1e-9
    )


# ------------------------------------------------------------------ exp 17
def test_exp17_oxygen_yield(exp17_dir, expected):
    """Whole-run yields against the ranges noted in local.yaml."""
    import exp17

    if "exp17_y_os" not in expected:
        pytest.skip("no expected Exp 17 yields in local.yaml")
    y_lo, y_hi = expected["exp17_y_os"]
    rq_lo, rq_hi = expected["exp17_rq"]
    for reactor in exp17.REACTORS:
        n = len(exp17.samples(exp17_dir, reactor))
        row = exp17.yields(exp17_dir, reactor, intervals=[(0, n - 1)])[0]
        assert y_lo < row["y_os"] < y_hi
        assert rq_lo < row["rq"] < rq_hi
