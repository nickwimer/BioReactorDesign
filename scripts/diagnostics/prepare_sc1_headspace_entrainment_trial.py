#!/usr/bin/env python3
"""Prepare a DIFFERENT scratch-only SC1 headspace re-entrainment trial.

Uses the existing headspace-O2 helper to retain Henry/film transfer, then
substitutes ONLY the scratch phaseTransfer deposition with the optional
deposition+entrainment model and adds a matching population-balance
nucleation source. Keeps the ORIGINAL production case untouched.

Required OpenFOAM shared libraries must already exist in scratch.
No production edits, no background work, no automatic CFD launch.
"""

import argparse
import importlib.util
import math
import re
import shutil
from pathlib import Path

ROOT = Path("/scratch/nwimer/sc1_headspace_entrainment_20261104")
HEADSPACE_O2_LIBRARY = Path(
    "/scratch/nwimer/sc1_headspace_o2_trial_20261103/lib/"
    "libsc1HeadspaceFilmModels.so"
)
ENTRAINMENT_LIBRARY = ROOT / "lib/libsc1HeadspaceEntrainmentModels.so"
SOURCE = Path("/scratch/nwimer/sc1_bubble_audit_20261029/run")


def load_headspace_helper():
    path = Path(__file__).with_name("prepare_sc1_headspace_o2_trial.py")
    spec = importlib.util.spec_from_file_location("sc1_headspace_base", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import existing scratch helper {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def enclosing_brace_end(text: str, opening: int) -> int:
    """Find a matching '}' for a known '{' in a phase-transfer dictionary."""
    if text[opening] != "{":
        raise ValueError("Expected a sub-dictionary opening brace")
    depth = 0
    state = "plain"
    i = opening
    while i < len(text):
        c = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if state == "line":
            if c == "\n":
                state = "plain"
        elif state == "block":
            if c == "*" and nxt == "/":
                state = "plain"
                i += 1
        elif state == "quote":
            if c == "\\":
                i += 1
            elif c == '"':
                state = "plain"
        else:
            if c == "/" and nxt == "/":
                state = "line"
                i += 1
            elif c == "/" and nxt == "*":
                state = "block"
                i += 1
            elif c == '"':
                state = "quote"
            elif c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    return i
        i += 1
    raise ValueError("Unterminated phase transfer model dictionary")


def install_entrainment_models(
    original: str,
    *,
    efficiency: float,
    maximum_rate_kg_s: float,
    nucleation_diameter_m: float,
    foam_multiplier: float = 0,
    velocity_criterion: str = "legacyDownwardLiquid",
) -> str:
    """Patch source string; do not touch any production/scratch file."""
    for name, value, maximum in (
        ("efficiency", efficiency, 1),
        ("maximum_rate_kg_s", maximum_rate_kg_s, 0.02),
        ("nucleation_diameter_m", nucleation_diameter_m, 0.011),
    ):
        if not math.isfinite(value) or not 0 <= value <= maximum:
            raise ValueError(f"{name} must be in [0,{maximum}]")

    if maximum_rate_kg_s <= 0 or nucleation_diameter_m <= 0:
        raise ValueError("Maximum rate and bubble diameter must be positive")
    if velocity_criterion not in ("legacyDownwardLiquid", "normalGasIntoLiquid"):
        raise ValueError("Unknown entrainment velocity criterion")

    base = load_headspace_helper()
    updated = base.install_phase_models(original, foam_multiplier)

    start, end = base.find_keyword_list(updated, "phaseTransfer") or (None, None)
    if start is None:
        raise ValueError("phaseTransfer list missing")
    transfer_body = updated[start + 1 : end]

    pattern = re.compile(
        r"\(\s*bubbleGas\s+in\s+headspaceGas\s*\)\s*\{"
    )
    matches = list(pattern.finditer(transfer_body))
    if len(matches) != 1:
        raise ValueError("Expected exactly one original bubbleGas deposition")
    opening = matches[0].end() - 1
    closing = enclosing_brace_end(transfer_body, opening)
    current = transfer_body[opening + 1:closing]
    type_match = re.compile(r"\btype\s+deposition\s*;")
    if len(type_match.findall(current)) != 1:
        raise ValueError("Original pair is not the expected deposition model")
    if "entrainmentEfficiency" in current:
        raise ValueError("Entrainment already configured")
    if not re.search(r"\bdroplet\s+bubbleGas\s*;", current):
        raise ValueError("Existing bubble phase configuration missing")
    if not re.search(r"\bsurface\s+headspaceGas\s*;", current):
        raise ValueError("Existing headspace phase configuration missing")

    new_body = type_match.sub("type headspaceEntrainingDeposition;", current)
    new_body += f"""
        // EXPERIMENTAL, separate from the validated original production run.
        entrainmentEfficiency {efficiency:.8g};
        liquid liquid;
        entrainmentYMin 0.68;
        entrainmentYMax 0.96;
        entrainmentRMin 0.05;
        entrainmentRMax 0.30;
        entrainmentAxisX 0;
        entrainmentAxisZ 0;
        minInterfaceAlpha 0.02;
        fullInterfaceAlpha 0.10;
        maxBubbleAlpha 0.10;
        maxSlipSpeed 5;
        donorTimeScale 0.05;
        maxGlobalEntrainedKgS {maximum_rate_kg_s:.8g};
        entrainmentVelocityCriterion {velocity_criterion};
        requireDownwardLiquid {'true' if velocity_criterion == 'legacyDownwardLiquid' else 'false'};
        reportRate true;
    """
    transfer_body = (
        transfer_body[:opening + 1]
        + new_body
        + transfer_body[closing:]
    )
    updated = updated[:start + 1] + transfer_body + updated[end:]

    if "nucleationModels" not in updated:
        raise ValueError("Population balance nucleationModels entry missing")

    sizes = [
        float(s) for s in re.findall(
            r"\bdSph\s+([+\-.\deE]+)\s*;", updated
        )
    ]
    if not sizes:
        raise ValueError("No size group diameters were found")
    if not min(sizes) <= nucleation_diameter_m <= max(sizes):
        raise ValueError("Entrainment bubble diameter outside size group range")

    # BiRD's PB sources() loops each nucleationModel over all size groups.
    # A companion model allocates POSITIVE NET bubbleGas creation, after
    # deposition and phaseTransfer blending, to adjacent size bins.
    nucleation_entry = f"""
    headspaceEntrainmentNucleation
    {{
        velocityGroup bubbleGas;
        headspacePhase headspaceGas;
        nucleationDiameter {nucleation_diameter_m:.8g};
    }}
"""
    if "headspaceEntrainmentNucleation" in updated:
        raise ValueError("Companion PB nucleation is already configured")
    updated = base.append_list_entry(
        updated, "nucleationModels", nucleation_entry
    )

    return updated


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--out", type=Path, default=ROOT / "run_entrain_001")
    parser.add_argument("--duration", type=float, default=0.5)
    parser.add_argument("--foam-multiplier", type=float, default=0)
    parser.add_argument(
        "--entrainment-efficiency", type=float, default=0,
        help="0 reproduces original deposition (safe control)."
    )
    parser.add_argument(
        "--velocity-criterion",
        choices=("legacyDownwardLiquid", "normalGasIntoLiquid"),
        default="legacyDownwardLiquid",
        help="Explicit optional kinematic eligibility; legacy is default.",
    )
    parser.add_argument(
        "--max-entrained-kg-s", type=float, default=0.002,
        help="Gross UNBLENDED gas entrainment cap, kg/s."
    )
    parser.add_argument(
        "--nucleation-diameter", type=float, default=0.002,
        help="PB nucleus diameter in metres, within existing size groups."
    )
    parser.add_argument(
        "--entrainment-library", type=Path, default=ENTRAINMENT_LIBRARY
    )
    parser.add_argument(
        "--headspace-library", type=Path, default=HEADSPACE_O2_LIBRARY
    )
    args = parser.parse_args()

    try:
        base = load_headspace_helper()
        src = base.require_scratch(args.source)
        out = base.require_scratch(args.out)
        film_lib = base.require_scratch(args.headspace_library)
        entrain_lib = base.require_scratch(args.entrainment_library)
        if not film_lib.is_file() or not entrain_lib.is_file():
            raise FileNotFoundError(
                "Compile BOTH scratch libraries before creating the trial."
            )
        if out.exists():
            raise FileExistsError(f"Trial exists; refusing overwrite: {out}")
        if src == out or src in out.parents or out in src.parents:
            raise ValueError("Scratch input/output paths must not overlap")
        if not 0 < args.duration <= 0.5:
            raise ValueError("Diagnostic duration must be >0 and <=0.5 s")

        tname, tstart = base.latest_collated_time(src)
        for field in (
            "O2.headspaceGas", "O2.bubbleGas",
            "O2.liquid", "alpha.headspaceGas", "alpha.bubbleGas",
            "f.bubbleGas",
        ):
            if not (src / "processors52" / tname / field).is_file():
                raise FileNotFoundError(
                    f"Missing collated starting field: {field}"
                )

        new_phase = install_entrainment_models(
            (src / "constant/phaseProperties").read_text(),
            efficiency=args.entrainment_efficiency,
            maximum_rate_kg_s=args.max_entrained_kg_s,
            nucleation_diameter_m=args.nucleation_diameter,
            foam_multiplier=args.foam_multiplier,
            velocity_criterion=args.velocity_criterion,
        )

        control = (src / "system/controlDict").read_text()
        control = base.set_control_entry(
            control, "startFrom", "latestTime"
        )
        control = base.set_control_entry(
            control, "endTime", f"{tstart + args.duration:.10f}"
        )
        control = base.set_control_entry(
            control, "writeControl", "runTime", insert_if_missing=True
        )
        control = base.set_control_entry(
            control, "writeInterval", "0.1"
        )
        control = base.add_library_to_control(control, film_lib)
        control = base.add_library_to_control(control, entrain_lib)

        # Copy only after everything above has passed. Never overwrite
        # an existing trial, and never edit a production case.
        out.mkdir(parents=True, exist_ok=False)
        shutil.copytree(src / "constant", out / "constant")
        shutil.copytree(src / "system", out / "system")
        (out / "processors52").mkdir()
        shutil.copytree(
            src / "processors52/constant",
            out / "processors52/constant"
        )
        shutil.copytree(
            src / "processors52" / tname,
            out / "processors52" / tname
        )
        (out / "constant/phaseProperties").write_text(new_phase)
        (out / "system/controlDict").write_text(control)

        print(f"Prepared isolated SC1 entrainment trial: {out}")
        print(f"  Starting checkpoint: {tname}; end={tstart + args.duration:.6f}")
        print(f"  Bubble entrainment efficiency: {args.entrainment_efficiency}")
        print(f"  Entrainment velocity criterion: {args.velocity_criterion}")
        print(f"  Raw entrainment cap: {args.max_entrained_kg_s} kg/s")
        print(f"  Nucleated bubble diameter: {args.nucleation_diameter} m")
        print(f"  Headspace O2 foam multiplier: {args.foam_multiplier}")
        print("  Existing bubble deposition correlation retained.")
        print("  PB headspaceEntrainmentNucleation added to bubble size groups.")
        print("  No production or source scratch files modified; no CFD run.")
        return 0

    except (OSError, ValueError, RuntimeError) as error:
        print(f"SC1 entrainment trial preparation failed: {error}")
        print("No production case modified. No CFD launched.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
