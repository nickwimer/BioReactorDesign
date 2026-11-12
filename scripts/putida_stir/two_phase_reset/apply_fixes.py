#!/usr/bin/env python3
"""Copy a putidaFoam stirred-tank case and apply the two-phase reset changes.

    apply_fixes.py report CASE
    apply_fixes.py apply SRC DST [--recommended | individual flags]

The source case is never modified: `apply` copies its inputs to a new
directory and edits the copy. Every change has its own flag, so a variant
with exactly one change can be made for a comparison run. An entry that is
expected but missing is an error, never a silent skip.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

START_DIRS = ("0", "0.org", "0.orig", "orig0")
MONITORS = "tankMonitors"
SKIP_DIR = re.compile(
    r"^(processors?\d.*|postProcessing|dynamicCode|VTK|\.git)$"
)
SKIP_FILE = re.compile(r"^(log.*|.*\.log|slurm-.*|core(\.\d+)?|.*\.foam)$")
ISOTHERMAL = {
    "multiComponentPhaseModel": "multiComponentIsothermalPhaseModel",
    "purePhaseModel": "pureIsothermalPhaseModel",
}


class CaseError(Exception):
    pass


# ---------------------------------------------------------------- parsing


def mask(text):
    """Blank comments, string contents and #{ #} code, keeping all offsets."""
    out = list(text)
    n = len(text)

    def blank(a, b):
        for k in range(a, b):
            if out[k] != "\n":
                out[k] = " "

    i = 0
    while i < n:
        two = text[i : i + 2]
        if two == "//":
            j = text.find("\n", i)
            j = n if j < 0 else j
            blank(i, j)
        elif two == "/*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            blank(i, j)
        elif two == "#{":
            j = text.find("#}", i + 2)
            j = n if j < 0 else j + 2
            blank(i, j)
        elif text[i] == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            j = min(j + 1, n)
            blank(i + 1, j - 1)  # the quotes stay, so a quoted key is a token
        else:
            j = i + 1
        i = j
    return "".join(out)


def _close(m, i):
    """Index of the bracket that closes the one at m[i]."""
    depth = 0
    for j in range(i, len(m)):
        if m[j] in "({[":
            depth += 1
        elif m[j] in ")}]":
            depth -= 1
            if depth == 0:
                return j
    raise CaseError("Unbalanced brackets")


def _stream_end(m, i, end):
    """Index of the ';' that ends a keyword's value starting at m[i]."""
    depth = 0
    for j in range(i, end):
        if m[j] in "({[":
            depth += 1
        elif m[j] in ")}]":
            depth -= 1
        elif m[j] == ";" and depth == 0:
            return j
    raise CaseError("Entry without a terminating ';'")


def entries(text, m, start=0, end=None):
    """Yield (key, key_start, value_start, entry_end, is_dict) at one level."""
    end = len(text) if end is None else end
    i = start
    while i < end:
        c = m[i]
        if c.isspace() or c in ";}":
            i += 1
        elif c == "#":  # #include and friends end with the line
            j = m.find("\n", i, end)
            i = end if j < 0 else j
        else:
            if c == '"':
                j = m.find('"', i + 1, end) + 1
                if j <= i:
                    raise CaseError("Unterminated quoted keyword")
            else:
                j = i
                depth = 0
                while j < end and (
                    depth or not (m[j].isspace() or m[j] in "{;")
                ):
                    depth += (m[j] == "(") - (m[j] == ")")
                    j += 1
            k = j
            while k < end and m[k].isspace():
                k += 1
            if k < end and m[k] == "{":
                e = _close(m, k) + 1
                yield text[i:j], i, k, e, True
            else:
                e = _stream_end(m, k, end) + 1
                yield text[i:j], i, k, e, False
            i = e


def find(text, key, start=0, end=None, m=None):
    m = mask(text) if m is None else m
    for ent in entries(text, m, start, end):
        if ent[0].strip('"') == key:
            return ent
    return None


def locate(text, path):
    """Follow a list of keys through nested dictionaries."""
    m = mask(text)
    start, end, ent = 0, len(text), None
    for depth, key in enumerate(path):
        ent = find(text, key, start, end, m)
        if ent is None:
            return None
        if depth < len(path) - 1:
            if not ent[4]:
                return None
            start, end = ent[2] + 1, ent[3] - 1
    return ent


def value(text, path):
    """Value of a keyword as one line of text, or None when absent."""
    ent = locate(text, path)
    if ent is None:
        return None
    raw = text[ent[2] : ent[3]]
    if ent[4]:  # a dictionary: show it without its comments
        return " ".join(mask(raw).split())
    return " ".join(raw[:-1].split())


def _indent(text, pos):
    line = text[text.rfind("\n", 0, pos) + 1 : pos]
    return line[: len(line) - len(line.lstrip())]


# ---------------------------------------------------------------- editing


def set_value(text, path, new, after=None):
    """Set `key new;` at `path`, adding the keyword when it is not there.

    A new top-level keyword goes after the entry named `after`; a new nested
    keyword goes at the end of its parent dictionary.
    """
    ent = locate(text, path)
    key = path[-1]
    if ent is not None:
        lead = text[ent[1] : ent[2]]
        lead = lead if lead[-1:].isspace() else lead + " "
        return text[: ent[1]] + f"{lead}{new};" + text[ent[3] :]
    if len(path) == 1:
        anchor = locate(text, [after]) if after else None
        if anchor is None:
            raise CaseError(f"Cannot place '{key}': no '{after}' entry")
        pad = _indent(text, anchor[1])
        return text[: anchor[3]] + f"\n{pad}{key} {new};" + text[anchor[3] :]
    parent = locate(text, path[:-1])
    if parent is None or not parent[4]:
        raise CaseError(f"No dictionary '{'/'.join(path[:-1])}'")
    pad = _indent(text, parent[1]) + "    "
    close = parent[3] - 1
    head = text[:close].rstrip()
    return (
        f"{head}\n{pad}{key} {new};\n{_indent(text, parent[1])}" + text[close:]
    )


def set_patch(text, patch, body):
    """Replace one patch's boundary condition in a field file."""
    field = locate(text, ["boundaryField"])
    if field is None:
        raise CaseError("No boundaryField")
    ent = locate(text, ["boundaryField", patch])
    if ent is None or not ent[4]:
        names = [
            e[0] for e in entries(text, mask(text), field[2] + 1, field[3] - 1)
        ]
        raise CaseError(f"No patch '{patch}' (found: {', '.join(names)})")
    pad = _indent(text, ent[1])
    lines = "".join(f"{pad}    {line}\n" for line in body)
    return (
        text[: ent[1]] + f"{patch}\n{pad}{{\n{lines}{pad}}}" + text[ent[3] :]
    )


def final_omega(raw):
    """Final rotation rate of a constant or tabulated omega entry."""
    number = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    numbers = re.findall(number, raw)
    words = re.findall(r"[A-Za-z_$]\w*", re.sub(number, " ", raw))
    known = {
        "table",
        "constant",
        "type",
        "values",
        "value",
        "outOfBounds",
        "clamp",
    }
    if not numbers or any(w not in known for w in words):
        raise CaseError(f"Cannot read omega '{raw}'; pass --omega")
    return float(numbers[-1])


def ramp_mrf(text, seconds, omega=None):
    """Make every MRF zone reach its rotation rate linearly over `seconds`."""
    zones = [
        e for e in entries(text, mask(text)) if e[4] and e[0] != "FoamFile"
    ]
    done = []
    for name in [z[0] for z in zones]:
        raw = value(text, [name, "omega"])
        if raw is None:
            continue
        rate = final_omega(raw) if omega is None else omega
        text = set_value(
            text, [name, "omega"], f"table ((0 0) ({seconds:g} {rate:.10g}))"
        )
        done.append((name, rate))
    if not done:
        raise CaseError("No MRF zone with an omega entry")
    return text, done


def add_monitors(text):
    """Include the monitor function objects in controlDict."""
    include = f'#include "{MONITORS}"'
    ent = locate(text, ["functions"])
    if ent is None:
        return text.rstrip() + f"\n\nfunctions\n{{\n    {include}\n}}\n"
    if include in text[ent[2] : ent[3]]:
        return text
    close = ent[3] - 1
    return text[:close].rstrip() + f"\n    {include}\n" + text[close:]


# ---------------------------------------------------------------- the case


def start_dirs(case):
    return [case / d for d in START_DIRS if (case / d).is_dir()]


def phase_names(phase_properties):
    raw = value(phase_properties, ["phases"])
    if raw is None:
        raise CaseError("No 'phases' entry in phaseProperties")
    return raw.strip("()").split()


def alpha_solver_key(fv_solution):
    solvers = locate(fv_solution, ["solvers"])
    if solvers is None:
        raise CaseError("No 'solvers' in fvSolution")
    m = mask(fv_solution)
    for ent in entries(fv_solution, m, solvers[2] + 1, solvers[3] - 1):
        if ent[4] and ent[0].strip('"').startswith("alpha"):
            return ent[0].strip('"')
    raise CaseError("No alpha solver entry in fvSolution")


def read(path):
    """Text of a dictionary or field file; binary field data is refused."""
    data = path.read_bytes()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        text = None
    # A header may say binary while the file only holds uniform values.
    if text is None or (
        value(text, ["FoamFile", "format"]) == "binary"
        and "nonuniform" in text
    ):
        raise CaseError(f"{path} holds binary data; use ASCII start fields")
    return text


def copy_case(src, dst):
    if dst.exists():
        raise CaseError(f"{dst} already exists")
    if not (src / "system" / "controlDict").is_file():
        raise CaseError(f"{src} is not an OpenFOAM case")

    def ignore(folder, names):
        if Path(folder) != src:
            return []
        skip = []
        for name in names:
            is_dir = (src / name).is_dir()
            is_time = is_dir and re.fullmatch(
                r"\d+(\.\d+)?([eE][-+]?\d+)?", name
            )
            if is_dir and (
                SKIP_DIR.match(name) or (is_time and name not in START_DIRS)
            ):
                skip.append(name)
            elif not is_dir and SKIP_FILE.match(name):
                skip.append(name)
        return skip

    shutil.copytree(src, dst, ignore=ignore, symlinks=True)


def apply(case, opts):
    """Edit the case in place (it is already a copy) and return a change log."""
    log = []
    liquid, top = opts.liquid, opts.top_patch
    pp = case / "constant" / "phaseProperties"

    if opts.reference_phase:
        text = read(pp)
        if opts.reference_phase not in phase_names(text):
            raise CaseError(
                f"'{opts.reference_phase}' is not one of the phases"
            )
        pp.write_text(
            set_value(
                text, ["referencePhase"], opts.reference_phase, after="phases"
            )
        )
        log.append(f"phaseProperties: referencePhase {opts.reference_phase}")

    if opts.mean_flux_reference:
        fs = case / "system" / "fvSolution"
        text = read(fs)
        key = alpha_solver_key(text)
        fs.write_text(
            set_value(text, ["solvers", key, "meanFluxReference"], "true")
        )
        log.append(f'fvSolution: solvers/"{key}"/meanFluxReference true')

    if opts.liquid_tight_top:
        edits = {
            f"U.{liquid}": ["type            slip;"],
            f"alpha.{liquid}": ["type            zeroGradient;"],
        }
        dirs = start_dirs(case)
        if not dirs:
            raise CaseError("No start-time field directory")
        for folder in dirs:
            for name, body in edits.items():
                path = folder / name
                if not path.is_file():
                    raise CaseError(f"Missing {path}")
                path.write_text(set_patch(read(path), top, body))
                log.append(
                    f"{folder.name}/{name}: {top} -> {body[0].split()[1].rstrip(';')}"
                )

    if opts.mrf_ramp:
        mrf = case / "constant" / "MRFProperties"
        text, zones = ramp_mrf(read(mrf), opts.mrf_ramp, opts.omega)
        mrf.write_text(text)
        for name, rate in zones:
            log.append(
                f"MRFProperties: {name} ramps 0 -> {rate:g} rad/s over {opts.mrf_ramp:g} s"
            )

    cd = case / "system" / "controlDict"
    for key, new in (("maxCo", opts.max_co), ("maxDeltaT", opts.max_delta_t)):
        if new is not None:
            cd.write_text(
                set_value(read(cd), [key], f"{new:g}", after="deltaT")
            )
            log.append(f"controlDict: {key} {new:g}")

    if opts.isothermal:
        text = read(pp)
        for phase in phase_names(text):
            old = value(text, [phase, "type"])
            if old in ISOTHERMAL.values():
                continue
            if old not in ISOTHERMAL:
                raise CaseError(f"No isothermal form of phase model '{old}'")
            text = set_value(text, [phase, "type"], ISOTHERMAL[old])
            log.append(f"phaseProperties: {phase} type {ISOTHERMAL[old]}")
        pp.write_text(text)

    if opts.monitors:
        shutil.copy(
            Path(__file__).with_name(MONITORS), case / "system" / MONITORS
        )
        cd.write_text(add_monitors(read(cd)))
        log.append(f"controlDict: functions include system/{MONITORS}")

    return log


def report(case, top="outlet"):
    """One line per setting that bears on conservation and the free surface."""
    pp = read(case / "constant" / "phaseProperties")
    fs = read(case / "system" / "fvSolution")
    sch = read(case / "system" / "fvSchemes")
    cd = read(case / "system" / "controlDict")
    rows = [
        ("phases", value(pp, ["phases"])),
        ("referencePhase", value(pp, ["referencePhase"])),
    ]
    for phase in phase_names(pp):
        rows.append((f"{phase} model", value(pp, [phase, "type"])))
    rows.append(("interfaceCompression", value(pp, ["interfaceCompression"])))
    key = alpha_solver_key(fs)
    for name in ("nAlphaSubCycles", "nAlphaCorr", "meanFluxReference"):
        rows.append((name, value(fs, ["solvers", key, name])))
    for name in (
        "nOuterCorrectors",
        "nCorrectors",
        "nNonOrthogonalCorrectors",
        "faceMomentum",
    ):
        rows.append((name, value(fs, ["PIMPLE", name])))
    rows.append(("snGrad default", value(sch, ["snGradSchemes", "default"])))
    rows.append(
        ("laplacian default", value(sch, ["laplacianSchemes", "default"]))
    )
    for name in ("application", "deltaT", "maxCo", "maxDeltaT"):
        rows.append((name, value(cd, [name])))
    mrf = case / "constant" / "MRFProperties"
    if mrf.is_file():
        text = read(mrf)
        for ent in entries(text, mask(text)):
            if ent[4] and ent[0] != "FoamFile":
                rows.append(
                    (f"MRF {ent[0]} omega", value(text, [ent[0], "omega"]))
                )
    rows.append(
        ("hRef", "set" if (case / "constant" / "hRef").is_file() else None)
    )
    for folder in start_dirs(case)[:1]:
        for path in (
            sorted(folder.glob("alpha.*"))
            + sorted(folder.glob("U.*"))
            + [folder / "p_rgh"]
        ):
            if path.is_file():
                bc = value(read(path), ["boundaryField", top, "type"])
                rows.append((f"{folder.name}/{path.name} {top}", bc))
    return rows


# ---------------------------------------------------------------- command line


def parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("report", help="print the relevant settings of a case")
    r.add_argument("case", type=Path)
    r.add_argument("--top-patch", default="outlet")
    a = sub.add_parser(
        "apply", help="copy a case and apply changes to the copy"
    )
    a.add_argument("src", type=Path)
    a.add_argument("dst", type=Path)
    a.add_argument(
        "--recommended",
        action="store_true",
        help="the full set, with default values",
    )
    a.add_argument(
        "--reference-phase",
        metavar="PHASE",
        help="phase obtained as the remainder",
    )
    a.add_argument(
        "--mean-flux-reference",
        action="store_true",
        help="bounded MULES flux reference",
    )
    a.add_argument(
        "--liquid-tight-top",
        action="store_true",
        help="no liquid flux through the top patch",
    )
    a.add_argument(
        "--mrf-ramp",
        type=float,
        metavar="SECONDS",
        help="linear impeller start-up",
    )
    a.add_argument(
        "--omega",
        type=float,
        help="final rad/s, when it cannot be read from the case",
    )
    a.add_argument("--max-co", type=float)
    a.add_argument("--max-delta-t", type=float)
    a.add_argument(
        "--isothermal", action="store_true", help="isothermal phase models"
    )
    a.add_argument(
        "--monitors",
        action="store_true",
        help="inventory, oxygen and outlet monitors",
    )
    a.add_argument("--liquid", default="liquid", help="liquid phase name")
    a.add_argument("--top-patch", default="outlet")
    return p


def main(argv=None):
    opts = parser().parse_args(argv)
    try:
        if opts.command == "report":
            for name, val in report(opts.case, opts.top_patch):
                print(f"{name:28s} {'(not set)' if val is None else val}")
            return 0
        if opts.recommended:
            opts.reference_phase = opts.reference_phase or "gas"
            opts.mean_flux_reference = opts.liquid_tight_top = True
            opts.isothermal = opts.monitors = True
            opts.mrf_ramp = opts.mrf_ramp or 2.0
            opts.max_co = opts.max_co or 0.3
            opts.max_delta_t = opts.max_delta_t or 5e-4
        wanted = (
            opts.reference_phase,
            opts.mean_flux_reference,
            opts.liquid_tight_top,
            opts.mrf_ramp,
            opts.max_co,
            opts.max_delta_t,
            opts.isothermal,
            opts.monitors,
        )
        if not any(w for w in wanted):
            raise CaseError("No change requested")
        copy_case(opts.src, opts.dst)
        try:
            log = apply(opts.dst, opts)
        except CaseError:
            shutil.rmtree(opts.dst)  # never leave a half-edited copy
            raise
        (opts.dst / "two_phase_reset_changes.txt").write_text(
            f"Copied from {opts.src.resolve()}\n"
            + "".join(f"{line}\n" for line in log)
        )
        print("\n".join(log))
        return 0
    except CaseError as err:
        print(f"apply_fixes: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
