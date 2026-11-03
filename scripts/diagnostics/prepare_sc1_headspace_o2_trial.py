#!/usr/bin/env python3
"""Create a scratch-only SC1 trial case with bounded headspace-liquid O2 transfer.

Default source is the existing scratch bubble-audit continuation, NEVER the
production case. Copies only the most recent collated time and mesh and
adds the Henry / diffusive / heat-transfer pair to a scratch phaseProperties.
Does not compile OpenFOAM, change production files, or start CFD.

Run after compiling the separate libsc1HeadspaceFilmModels.so in scratch.
"""

import argparse
import re
import shutil
from pathlib import Path

SOURCE = Path("/scratch/nwimer/sc1_bubble_audit_20261029/run")
TRIAL = Path("/scratch/nwimer/sc1_headspace_o2_trial_20261103/run")
LIBRARY = Path(
    "/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/"
    "libsc1HeadspaceFilmModels.so"
)


def require_scratch(path: Path) -> Path:
    real = path.resolve()
    scratch = Path("/scratch/nwimer").resolve()
    if real == scratch or not real.is_relative_to(scratch):
        raise ValueError(f"Must remain inside {scratch}: {real}")
    return real


def find_keyword_list(text: str, keyword: str) -> tuple[int, int] | None:
    """Locate one OpenFOAM top-level keyword (...) list, ignoring comments."""
    rex = re.compile(
        r"(?m)^[ \t]*" + re.escape(keyword)
        + r"[ \t]*(?:\n[ \t]*)?\("
    )
    found = list(rex.finditer(text))
    if not found:
        return None
    if len(found) != 1:
        raise ValueError(f"Expected one list for {keyword}, found {len(found)}")

    start = found[0].end() - 1
    depth = 0
    state = "plain"
    i = start
    while i < len(text):
        c = text[i]
        n = text[i + 1] if i + 1 < len(text) else ""
        if state == "line":
            if c == "\n":
                state = "plain"
        elif state == "block":
            if c == "*" and n == "/":
                state = "plain"
                i += 1
        elif state == "quote":
            if c == "\\":
                i += 1
            elif c == '"':
                state = "plain"
        else:
            if c == "/" and n == "/":
                state = "line"
                i += 1
            elif c == "/" and n == "*":
                state = "block"
                i += 1
            elif c == '"':
                state = "quote"
            elif c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    return (start, i)
        i += 1
    raise ValueError(f"Unterminated list in {keyword}")


def append_list_entry(text: str, keyword: str, entry: str) -> str:
    position = find_keyword_list(text, keyword)
    if position is None:
        return text.rstrip() + f"\n\n{keyword}\n(\n{entry}\n);\n"
    _, closing = position
    return text[:closing] + "\n" + entry + "\n" + text[closing:]


def section_body(text: str, keyword: str) -> str:
    bounds = find_keyword_list(text, keyword)
    if bounds is None:
        raise ValueError(f"Missing source dictionary list {keyword}")
    return text[bounds[0] + 1 : bounds[1]]


def pair_exists(section: str, phase1: str, direction: str, phase2: str) -> bool:
    return bool(
        re.search(
            r"\(\s*" + re.escape(phase1)
            + r"\s+" + re.escape(direction)
            + r"\s+" + re.escape(phase2) + r"\s*\)",
            section,
        )
    )


def model_dictionary(model: str, report: bool) -> str:
    return """    {
        type resolvedHeadspaceFilm;
        headspacePhase headspaceGas;
        liquidPhase liquid;
        bubblePhase bubbleGas;

        // kL = D_O2 / filmThickness. Pilot value; NOT calibrated.
        filmThickness 1e-5;

        // Geometrically resolved macro-interface only by default.
        // Explore bounded foam enhancement after baseline smoke test.
        foamAreaMultiplier 0;
        foamDiameter 0.003;
        foamYMin 0.65;
        foamYMax 1.05;
        foamRMax 0.30;
        foamAxisX 0;
        foamAxisZ 0;

        minInterfaceAlpha 0.02;
        fullInterfaceAlpha 0.10;
        maxBubbleAlpha 0.10;
        maxAreaDensity 2000;
        reportArea %s;
    }""" % ("true" if report else "false")


def install_phase_models(original: str) -> str:
    phase_system = re.search(
        r"(?m)^\s*type\s+([^;]+);", original
    )
    if not phase_system or "interfaceCompositionPhaseChange" not in phase_system[1]:
        raise ValueError("Case is not the expected interface-composition phase system.")

    if "bubbleGas" not in original or "headspaceGas" not in original:
        raise ValueError("Expected bubbleGas and headspaceGas phases missing.")

    if "( bubbleGas in headspaceGas )" not in original and not re.search(
        r"\(\s*bubbleGas\s+in\s+headspaceGas\s*\)", original
    ):
        raise ValueError("Expected existing bubble-to-headspace deposition not found.")

    # Preserve the exact solubility and species constants of the active
    # bubble-liquid Henry model, instead of inventing He_O2 numbers.
    old_liquid = section_body(original, "interfaceComposition.liquid")
    m = re.search(
        r"\(\s*liquid\s+and\s+bubbleGas\s*\)\s*\{",
        old_liquid,
    )
    if not m:
        raise ValueError("Missing existing (liquid and bubbleGas) Henry model.")
    start = m.start()
    depth = 0
    close = -1
    for i in range(m.end() - 1, len(old_liquid)):
        if old_liquid[i] == "{":
            depth += 1
        elif old_liquid[i] == "}":
            depth -= 1
            if depth == 0:
                close = i + 1
                break
    if close < 0:
        raise ValueError("Unclosed bubble-liquid Henry model.")

    henry_entry = old_liquid[start:close]
    if not re.search(r"\bspecies\s*\(\s*O2\s*\)", henry_entry):
        raise ValueError("Existing bubble-liquid Henry model is not O2-only.")
    if not re.search(r"\btype\s+Henry\s*;", henry_entry):
        raise ValueError("Expected unmodified Henry O2 model.")

    if pair_exists(old_liquid, "liquid", "and", "headspaceGas"):
        raise ValueError("A headspace-liquid composition model already exists; abort.")

    updated = append_list_entry(
        original, "interfaceComposition.liquid",
        henry_entry.replace("bubbleGas", "headspaceGas"),
    )

    mass_key = "diffusiveMassTransfer.liquid"
    for dispersed, continuous, report in (
        ("headspaceGas", "liquid", True),
        ("liquid", "headspaceGas", False),
    ):
        body = section_body(updated, mass_key)
        if pair_exists(body, dispersed, "in", continuous):
            raise ValueError(
                f"Model already exists: ({dispersed} in {continuous})"
            )
        entry = (
            f"    ({dispersed} in {continuous})\n"
            + model_dictionary("resolvedHeadspaceFilm", report)
        )
        updated = append_list_entry(updated, mass_key, entry)

    # TwoResistanceHeatTransferPhaseSystem requires heat-transfer models
    # on BOTH phase sides of every interface-composition pair.
    for phase, model1, model2 in (
        ("headspaceGas", "spherical", "RanzMarshall"),
        ("liquid", "RanzMarshall", "spherical"),
    ):
        key = f"heatTransfer.{phase}"
        for direction, model in (
            (("headspaceGas", "liquid"), model1),
            (("liquid", "headspaceGas"), model2),
        ):
            body = section_body(updated, key) if find_keyword_list(updated, key) else ""
            if pair_exists(body, direction[0], "in", direction[1]):
                continue
            entry = (
                f"    ({direction[0]} in {direction[1]})\n"
                "    {\n"
                f"        type {model};\n"
                "        residualAlpha 1e-4;\n"
                "    }"
            )
            updated = append_list_entry(updated, key, entry)

    return updated


def set_control_entry(text: str, key: str, value: str) -> str:
    pat = re.compile(r"(?m)^[ \t]*" + re.escape(key) + r"\s+[^;]*;")
    matches = list(pat.finditer(text))
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {key} in controlDict")
    return pat.sub(f"{key}    {value};", text, count=1)


def add_library_to_control(text: str, library: Path) -> str:
    lib = str(library.resolve())
    if lib in text:
        raise ValueError("controlDict already references trial library.")
    pair = find_keyword_list(text, "libs")
    if pair is None:
        return text.rstrip() + '\n\nlibs ("' + lib + '");\n'
    _, closing = pair
    return text[:closing] + '\n    "' + lib + '"\n' + text[closing:]


def latest_collated_time(root: Path) -> tuple[str, float]:
    times = []
    for item in (root / "processors52").iterdir():
        if item.is_dir():
            try:
                value = float(item.name)
            except ValueError:
                continue
            times.append((item.name, value))
    if not times:
        raise ValueError(f"No numeric collated times in {root}")
    return max(times, key=lambda pair: pair[1])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out", type=Path, default=TRIAL)
    parser.add_argument("--library", type=Path, default=LIBRARY)
    parser.add_argument(
        "--duration", type=float, default=0.2,
        help="Very short diagnostic-only CFD advance (s), not equilibrium.",
    )
    args = parser.parse_args()

    try:
        src = require_scratch(args.source)
        dst = require_scratch(args.out)
        lib = require_scratch(args.library)

        if src == dst or src in dst.parents or dst in src.parents:
            raise ValueError("Source and trial paths may not overlap.")
        if dst.exists():
            raise FileExistsError(f"Trial output exists; refusing overwrite: {dst}")
        if args.duration <= 0 or args.duration > 0.5:
            raise ValueError("Diagnostic duration must be >0 and <=0.5 s.")

        start_name, start_value = latest_collated_time(src)
        for field in ("O2.headspaceGas", "O2.liquid", "alpha.headspaceGas"):
            if not (src / "processors52" / start_name / field).is_file():
                raise FileNotFoundError(
                    f"Missing critical field at selected collated time: {field}"
                )

        original = (src / "constant/phaseProperties").read_text()
        updated_phase = install_phase_models(original)

        original_control = (src / "system/controlDict").read_text()
        control = set_control_entry(
            original_control, "startFrom", "latestTime"
        )
        control = set_control_entry(
            control, "endTime", f"{start_value + args.duration:.10f}"
        )
        control = set_control_entry(
            control, "writeControl", "runTime"
        )
        control = set_control_entry(
            control, "writeInterval", "0.1"
        )
        control = add_library_to_control(control, lib)

        # Create only after all model validation has succeeded.
        dst.mkdir(parents=True, exist_ok=False)
        shutil.copytree(src / "constant", dst / "constant")
        shutil.copytree(src / "system", dst / "system")
        (dst / "processors52").mkdir()
        shutil.copytree(
            src / "processors52/constant",
            dst / "processors52/constant",
        )
        shutil.copytree(
            src / "processors52" / start_name,
            dst / "processors52" / start_name,
        )

        (dst / "constant/phaseProperties").write_text(updated_phase)
        (dst / "system/controlDict").write_text(control)

        print(f"Trial case prepared: {dst}")
        print(f"  Starting collated checkpoint: {start_name}")
        print(f"  End time: {start_value + args.duration:.10f}")
        print(f"  Experimental dynamic library: {lib}")
        print("  Existing bubble-liquid Henry/transfer retained unchanged.")
        print("  Headspace-liquid Henry + two film orientations added.")
        print("  Foam area multiplier: zero (geometric-area baseline).")
        print("  Production case and source scratch continuation untouched.")
        print("  NOT YET compiled or runtime verified on Kestrel.")
    except (OSError, ValueError, FileNotFoundError) as err:
        print(f"Trial preparation stopped: {err}")
        print("No production case was modified.")


if __name__ == "__main__":
    main()
