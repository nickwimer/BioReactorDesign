"""Tests for apply_fixes.py and check_run.py.

    python -m pytest scripts/putida_stir/two_phase_reset

The case edits are exercised on the two-phase stirred-tank cases that are in
this repository. They check the text of the dictionaries, not that OpenFOAM
accepts them: that needs a run.
"""

import filecmp
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))

import apply_fixes as af  # noqa: E402
import check_run as cr  # noqa: E402

PUTIDA = REPO / "scripts" / "putida_stir" / "nonreact"
TUTORIAL = REPO / "tutorial_cases" / "stirred_tank"


def run(*args):
    return af.main([str(a) for a in args])


def snapshot(folder):
    return {
        p.relative_to(folder): p.read_bytes()
        for p in folder.rglob("*")
        if p.is_file()
    }


def changed(src, dst):
    cmp = filecmp.dircmp(src, dst)
    out = set()

    def walk(c, prefix=Path()):
        out.update(prefix / n for n in c.diff_files + c.right_only)
        for name, sub in c.subdirs.items():
            walk(sub, prefix / name)

    walk(cmp)
    return {str(p) for p in out}


def balanced(text):
    m = af.mask(text)
    return all(m.count(a) == m.count(b) for a, b in ("{}", "()", "[]"))


# ------------------------------------------------------------ dictionary parsing

SAMPLE = """\
FoamFile { version 2.0; format ascii; }
#include "$FOAM_CASE/constant/globalTransport"
// phases (wrong one);
phases (gas liquid);   /* referencePhase liquid; */
D_O2 #codeStream { code #{ os << 1; } #}; };
note "a } string; with { braces";
gas
{
    type            multiComponentPhaseModel; // model
    residualAlpha   1e-6;
}
drag
(
    (gas in liquid) { type Grace; residualRe 1e-3; }
);
"alpha.*" { nAlphaCorr 1; }
blending { default { type linear; } massTransfer { $default; } }
"""


def test_entries_skip_comments_strings_and_code():
    keys = [e[0] for e in af.entries(SAMPLE, af.mask(SAMPLE))]
    assert keys == [
        "FoamFile",
        "phases",
        "D_O2",
        "note",
        "gas",
        "drag",
        '"alpha.*"',
        "blending",
    ]


def test_values_are_read_through_nested_dictionaries():
    assert af.value(SAMPLE, ["phases"]) == "(gas liquid)"
    assert af.value(SAMPLE, ["referencePhase"]) is None
    assert af.value(SAMPLE, ["gas", "type"]) == "multiComponentPhaseModel"
    assert af.value(SAMPLE, ["alpha.*", "nAlphaCorr"]) == "1"
    assert af.value(SAMPLE, ["blending", "default", "type"]) == "linear"
    assert af.value(SAMPLE, ["gas", "missing"]) is None


def test_set_value_replaces_adds_and_is_repeatable():
    once = af.set_value(SAMPLE, ["referencePhase"], "gas", after="phases")
    assert "phases (gas liquid);\nreferencePhase gas;" in once
    again = af.set_value(once, ["referencePhase"], "gas", after="phases")
    assert again == once
    swapped = af.set_value(once, ["referencePhase"], "liquid", after="phases")
    assert af.value(swapped, ["referencePhase"]) == "liquid"
    assert (
        swapped.count("referencePhase") == SAMPLE.count("referencePhase") + 1
    )

    nested = af.set_value(SAMPLE, ["alpha.*", "meanFluxReference"], "true")
    assert af.value(nested, ["alpha.*", "meanFluxReference"]) == "true"
    assert af.value(nested, ["alpha.*", "nAlphaCorr"]) == "1"
    typed = af.set_value(
        SAMPLE, ["gas", "type"], "multiComponentIsothermalPhaseModel"
    )
    assert (
        "type            multiComponentIsothermalPhaseModel; // model" in typed
    )
    assert balanced(once) and balanced(nested) and balanced(typed)


def test_set_value_needs_a_place_for_a_new_keyword():
    with pytest.raises(af.CaseError):
        af.set_value(SAMPLE, ["referencePhase"], "gas", after="nothing")
    with pytest.raises(af.CaseError):
        af.set_value(SAMPLE, ["absent", "key"], "1")


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("10.0", 10.0),
        ("17.38", 17.38),
        ("constant 17.38", 17.38),
        ("table ((0 0) (0.1 17.38))", 17.38),
        ("{ type table; values ((0 0) (1e-1 1.738e1)); }", 17.38),
    ],
)
def test_final_omega(raw, expected):
    assert af.final_omega(raw) == pytest.approx(expected)


def test_final_omega_refuses_what_it_cannot_read():
    with pytest.raises(af.CaseError):
        af.final_omega("$impellerSpeed")
    with pytest.raises(af.CaseError):
        af.final_omega("{ type linearRamp; duration 2; }")


# ------------------------------------------------------------ whole cases


@pytest.mark.parametrize(
    "src, start", [(PUTIDA, "0.org"), (TUTORIAL, "orig0")]
)
def test_recommended_set_on_repository_cases(tmp_path, src, start, capsys):
    before = snapshot(src)
    dst = tmp_path / "case"
    assert run("apply", src, dst, "--recommended") == 0
    assert snapshot(src) == before  # the source case is never touched

    pp = (dst / "constant" / "phaseProperties").read_text()
    assert af.value(pp, ["referencePhase"]) == "gas"
    assert (
        af.value(pp, ["gas", "type"]) == "multiComponentIsothermalPhaseModel"
    )
    assert (
        af.value(pp, ["liquid", "type"])
        == "multiComponentIsothermalPhaseModel"
    )

    fs = (dst / "system" / "fvSolution").read_text()
    assert af.value(fs, ["solvers", "alpha.*", "meanFluxReference"]) == "true"
    assert af.value(fs, ["solvers", "alpha.*", "nAlphaSubCycles"]) == af.value(
        (src / "system" / "fvSolution").read_text(),
        ["solvers", "alpha.*", "nAlphaSubCycles"],
    )

    u = (dst / start / "U.liquid").read_text()
    a = (dst / start / "alpha.liquid").read_text()
    assert af.value(u, ["boundaryField", "outlet"]) == "{ type slip; }"
    assert af.value(a, ["boundaryField", "outlet"]) == "{ type zeroGradient; }"
    assert af.value(u, ["boundaryField", "walls", "type"]) == "noSlip"
    assert af.value(a, ["boundaryField", "walls", "type"]) == "zeroGradient"

    rate = af.final_omega(
        af.value(
            (src / "constant" / "MRFProperties").read_text(), ["MRF1", "omega"]
        )
    )
    mrf = (dst / "constant" / "MRFProperties").read_text()
    assert af.value(mrf, ["MRF1", "omega"]) == f"table ((0 0) (2 {rate:.10g}))"

    cd = (dst / "system" / "controlDict").read_text()
    assert af.value(cd, ["maxCo"]) == "0.3"
    assert af.value(cd, ["maxDeltaT"]) == "0.0005"
    assert '#include "tankMonitors"' in cd[af.locate(cd, ["functions"])[2] :]
    assert (dst / "system" / "tankMonitors").read_bytes() == (
        HERE / "tankMonitors"
    ).read_bytes()

    assert changed(src, dst) == {
        f"{start}/U.liquid",
        f"{start}/alpha.liquid",
        "constant/MRFProperties",
        "constant/phaseProperties",
        "system/controlDict",
        "system/fvSolution",
        "system/tankMonitors",
        "two_phase_reset_changes.txt",
    }
    for path in dst.rglob("*"):
        if path.is_file() and path.parent.name in (
            "constant",
            "system",
            start,
        ):
            assert balanced(path.read_text(errors="replace")), path

    # Applying the same set to the edited copy changes no dictionary again;
    # only the change record differs, because it names its source.
    again = tmp_path / "again"
    assert run("apply", dst, again, "--recommended") == 0
    assert changed(dst, again) == {"two_phase_reset_changes.txt"}


@pytest.mark.parametrize(
    "flags, files",
    [
        (["--reference-phase", "gas"], {"constant/phaseProperties"}),
        (["--reference-phase", "liquid"], {"constant/phaseProperties"}),
        (["--mean-flux-reference"], {"system/fvSolution"}),
        (["--liquid-tight-top"], {"0.org/U.liquid", "0.org/alpha.liquid"}),
        (["--mrf-ramp", "1.5"], {"constant/MRFProperties"}),
        (["--max-co", "0.25"], {"system/controlDict"}),
        (["--isothermal"], {"constant/phaseProperties"}),
        (["--monitors"], {"system/controlDict", "system/tankMonitors"}),
    ],
)
def test_one_flag_changes_only_its_own_files(tmp_path, flags, files):
    dst = tmp_path / "case"
    assert run("apply", PUTIDA, dst, *flags) == 0
    assert changed(PUTIDA, dst) == files | {"two_phase_reset_changes.txt"}


def test_run_outputs_are_not_copied(tmp_path):
    src = tmp_path / "src"
    af.copy_case(PUTIDA, src)
    for name in (
        "processor0",
        "processors52",
        "postProcessing",
        "0.1",
        "31.199982372783058",
        "0",
    ):
        (src / name).mkdir()
        (src / name / "alpha.liquid").write_text("x")
    (src / "log.putidaFoam").write_text("x")
    dst = tmp_path / "dst"
    af.copy_case(src, dst)
    names = {p.name for p in dst.iterdir()}
    assert {"0", "0.org", "constant", "system"} <= names
    assert not names & {
        "processor0",
        "processors52",
        "postProcessing",
        "0.1",
        "31.199982372783058",
        "log.putidaFoam",
    }


def test_failures_leave_no_copy_behind(tmp_path, capsys):
    dst = tmp_path / "case"
    assert (
        run("apply", PUTIDA, dst, "--liquid-tight-top", "--top-patch", "lid")
        == 2
    )
    assert not dst.exists()
    assert "No patch 'lid'" in capsys.readouterr().err
    assert run("apply", PUTIDA, dst, "--reference-phase", "water") == 2
    assert not dst.exists()
    assert run("apply", PUTIDA, dst) == 2  # nothing requested
    assert not dst.exists()
    dst.mkdir()
    assert run("apply", PUTIDA, dst, "--monitors") == 2  # never overwrite


def test_binary_start_fields_are_refused(tmp_path):
    src = tmp_path / "src"
    af.copy_case(PUTIDA, src)
    path = src / "0.org" / "U.liquid"
    text = path.read_text()
    assert "format      ascii" in text and "uniform (0 0 0)" in text
    text = text.replace("format      ascii", "format      binary")
    path.write_text(
        text.replace(
            "uniform (0 0 0)", "nonuniform List<vector> 2((0 0 0)(0 0 0))", 1
        )
    )
    assert run("apply", src, tmp_path / "dst", "--liquid-tight-top") == 2
    assert not (tmp_path / "dst").exists()


def test_report_names_the_settings_that_matter(capsys):
    assert run("report", PUTIDA) == 0
    out = capsys.readouterr().out
    lines = dict(
        line.split(None, 1)
        for line in out.splitlines()
        if "  " in line and line.split()[0] in ("referencePhase", "maxCo")
    )
    assert lines["referencePhase"].strip() == "(not set)"
    assert lines["maxCo"].strip() == "0.5"
    assert (
        "0.org/U.liquid outlet" in out and "pressureInletOutletVelocity" in out
    )


# ------------------------------------------------------------ log checker


def log_text(fractions, with_sum=None, liquid="liquid"):
    """A putidaFoam-like log: two time lines, several fraction lines per step."""
    lines = ["Create time", "liquid fraction is not on this line"]
    for n, frac in enumerate(fractions, start=1):
        t = n * 1e-3
        lines += [
            f"Time = {(n - 1) * 1e-3:g}",
            "bio update time: 14400",
            "Courant Number mean: 0.01 max: 0.2",
        ]
        lines += [f"deltaT = {1e-3 / n:g}", f"Time = {t:.6f}", ""]
        for sub in range(3):
            last = sub == 2
            lines.append(f"gas fraction, min, max = {1 - frac:g} 0 1")
            lines.append(
                f"{liquid} fraction, min, max = {frac if last else 0.5:.6g} {-1e-9 if last else 0:g} 1"
            )
            if with_sum:
                lines.append(
                    f"Phase-sum volume fraction, min, max = 1 {with_sum[0]:g} {with_sum[1]:g}"
                )
        lines.append("ExecutionTime = 1 s  ClockTime = 1 s")
    return [line + "\n" for line in lines]


def test_log_gives_one_record_per_step_with_the_last_values():
    records = cr.read_log(log_text([0.609, 0.609, 0.61]))
    assert [r["time"] for r in records] == pytest.approx([1e-3, 2e-3, 3e-3])
    assert [r["fraction"] for r in records] == [0.609, 0.609, 0.61]
    assert [r["deltaT"] for r in records] == pytest.approx(
        [1e-3, 5e-4, 1e-3 / 3]
    )
    assert all("sum_min" not in r for r in records)
    assert records[0]["alpha_min"] == -1e-9


def test_summary_reports_gain_in_litres_and_the_phase_sum():
    records = cr.read_log(
        log_text([0.609, 0.6181, 0.6272], with_sum=(0.97, 1.0002))
    )
    s = cr.summarise(records, volume=1.642)
    assert s["rescaled"] and s["sum_min"] == 0.97 and s["sum_max"] == 1.0002
    assert s["relative_change"] == pytest.approx(0.0182 / 0.609)
    assert s["change_L"] == pytest.approx(0.0182 * 1642)
    assert s["rate_L_per_s"] == pytest.approx(0.0182 * 1642 / 2e-3)


def test_checker_exit_status_and_csv(tmp_path, capsys):
    good, bad, empty = (
        tmp_path / n for n in ("good.log", "bad.log", "empty.log")
    )
    good.write_text("".join(log_text([0.609, 0.609001, 0.609])))
    bad.write_text(
        "".join(log_text([0.609, 0.615, 0.62], with_sum=(0.98, 1.0)))
    )
    empty.write_text("Time = 1\nEnd\n")
    out = tmp_path / "history.csv"

    assert cr.main([str(good), "--csv", str(out)]) == 0
    text = capsys.readouterr().out
    assert "PASS" in text and "a reference phase is in use" in text
    assert (
        out.read_text().splitlines()[0]
        == "time,fraction,alpha_min,alpha_max,sum_min,sum_max,deltaT"
    )
    assert len(out.read_text().splitlines()) == 4

    assert cr.main([str(bad), "--volume", "1.642"]) == 1
    text = capsys.readouterr().out
    assert "FAIL" in text and "no reference phase" in text and "L/s" in text
    assert cr.main([str(bad), "--tol", "0.05"]) == 0
    assert cr.main([str(empty)]) == 2


def test_checker_follows_a_differently_named_liquid():
    records = cr.read_log(log_text([0.4, 0.4], liquid="broth"), liquid="broth")
    assert len(records) == 2
    assert cr.read_log(log_text([0.4, 0.4], liquid="broth")) == []
