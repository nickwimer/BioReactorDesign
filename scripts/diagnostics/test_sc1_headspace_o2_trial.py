#!/usr/bin/env python3
"""Pure-Python checks of the isolated SC1 headspace-liquid patcher and closure.

These tests do not replace a RHEL9 C++ compile and runtime solver test.
They require no OpenFOAM installation and write no files.
"""

import importlib.util
from pathlib import Path


def load_installer():
    path = Path(__file__).with_name("prepare_sc1_headspace_o2_trial.py")
    spec = importlib.util.spec_from_file_location("sc1_prepare_headspace", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load installer")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def synthetic_phase_properties() -> str:
    return """
type interfaceCompositionPhaseChangePopulationBalanceMultiphaseSystem;
phases (bubbleGas headspaceGas liquid);
phaseTransfer
(
    (bubbleGas in headspaceGas)
    {
        type deposition;
        droplet bubbleGas;
        surface headspaceGas;
    }
);
interfaceComposition.bubbleGas ();
interfaceComposition.headspaceGas ();
interfaceComposition.liquid
(
    (liquid and bubbleGas)
    {
        type Henry;
        species (O2);
        k ($He_O2);
        Le $LeLiqMix;
    }
);
diffusiveMassTransfer.bubbleGas ();
diffusiveMassTransfer.headspaceGas ();
diffusiveMassTransfer.liquid
(
    (bubbleGas in liquid)
    {
        type Higbie;
        Le $LeLiqMix;
    }
    (liquid in bubbleGas)
    {
        type spherical;
        Le 1;
    }
);
heatTransfer.bubbleGas
(
    (bubbleGas in liquid) { type spherical; residualAlpha 1e-4; }
    (liquid in bubbleGas) { type RanzMarshall; residualAlpha 1e-4; }
);
heatTransfer.headspaceGas
(
);
heatTransfer.liquid
(
    (bubbleGas in liquid) { type RanzMarshall; residualAlpha 1e-4; }
    (liquid in bubbleGas) { type spherical; residualAlpha 1e-4; }
);
"""


def test_installer() -> None:
    patcher = load_installer()
    source = synthetic_phase_properties()
    changed = patcher.install_phase_models(source)
    assert changed.count("type resolvedHeadspaceFilm;") == 2
    assert changed.count("type Henry;") == 2
    assert changed.count("species (O2);") == 2
    assert changed.count("foamAreaMultiplier 0;") == 2
    assert changed.count("reportArea true;") == 1
    assert changed.count("reportArea false;") == 1
    assert "(liquid and headspaceGas)" in changed
    assert "(headspaceGas in liquid)" in changed
    assert "(liquid in headspaceGas)" in changed
    assert "k ($He_O2)" in changed

    for key in (
        "interfaceComposition.liquid",
        "diffusiveMassTransfer.liquid",
        "heatTransfer.headspaceGas",
        "heatTransfer.liquid",
    ):
        bounds = patcher.find_keyword_list(changed, key)
        assert bounds is not None
        assert changed[bounds[0]] == "("
        assert changed[bounds[1]] == ")"

    # Parameter selection is made per isolated scenario, not by hand-
    # editing the production phaseProperties file.
    with_foam = patcher.install_phase_models(source, foam_multiplier=0.25)
    assert with_foam.count("foamAreaMultiplier 0.25;") == 2
    try:
        patcher.install_phase_models(source, foam_multiplier=-1)
    except ValueError:
        pass
    else:
        raise AssertionError("Negative foam-area multiplier was accepted")

    # Reapplication must refuse a duplicate pair, not silently double it.
    try:
        patcher.install_phase_models(changed)
    except ValueError:
        pass
    else:
        raise AssertionError("Installer did not refuse existing pair")

    before = """startFrom startTime;
endTime 31.2;
writeControl timeStep;
writeInterval 10;
libs ("libbirdphaseSystem.so");
"""
    modified = patcher.add_library_to_control(
        before, Path("/scratch/nwimer/libsc1HeadspaceFilm.so")
    )
    assert modified.count("libsc1HeadspaceFilm.so") == 1
    assert modified.count("libbirdphaseSystem.so") == 1
    assert patcher.set_control_entry(modified, "endTime", "31.4").count(
        "endTime    31.4;"
    ) == 1

    # The real Kestrel scratch controlDict omitted writeControl.
    # The installer must add an explicit runTime setting to the *copy*.
    missing_write_control = """startFrom latestTime;
endTime 31.2;
writeInterval 0.1;
// writeControl timeStep;  commented out, not an active control
"""
    appended = patcher.set_control_entry(
        missing_write_control,
        "writeControl", "runTime", insert_if_missing=True
    )
    assert appended.count("writeControl    runTime;") == 1
    assert missing_write_control.count("writeControl    runTime;") == 0
    assert appended.count("writeInterval 0.1;") == 1

    # If explicitly present, replace rather than insert a second entry.
    explicit = patcher.set_control_entry(
        before, "writeControl", "runTime", insert_if_missing=True
    )
    assert explicit.count("writeControl    runTime;") == 1
    assert "writeControl timeStep;" not in explicit

    # Never silently accept a duplicated setting.
    duplicate = before + "writeControl adjustableRunTime;\n"
    try:
        patcher.set_control_entry(
            duplicate, "writeControl", "runTime", insert_if_missing=True
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Duplicate writeControl was not rejected")

    # Actual SC1 failure: 15 nested function-object writeControl entries,
    # but no top-level writeControl. Do not mistake these for 15 duplicates.
    nested_blocks = "\n".join(
        f"""
    sample{i}
    {{
        type volFieldValue;
        writeControl timeStep;
        writeInterval 2;
        // the nested settings must not change
    }}
"""
        for i in range(15)
    )
    complex_control = (
        "application putidaFoam;\n"
        "startFrom latestTime;\n"
        "endTime 31.2;\n"
        "writeInterval 0.1;\n"
        "functions\n{\n"
        + nested_blocks
        + "\n}\n"
    )
    assert len(patcher.top_level_control_matches(
        complex_control, "writeControl"
    )) == 0
    assert len(patcher.top_level_control_matches(
        complex_control, "writeInterval"
    )) == 1
    inserted = patcher.set_control_entry(
        complex_control,
        "writeControl", "runTime", insert_if_missing=True
    )
    assert len(patcher.top_level_control_matches(
        inserted, "writeControl"
    )) == 1
    assert inserted.count("writeControl timeStep;") == 15
    assert inserted.count("writeControl    runTime;") == 1
    assert inserted.count("writeInterval 2;") == 15
    updated_interval = patcher.set_control_entry(
        inserted, "writeInterval", "0.1"
    )
    assert updated_interval.count("writeInterval 2;") == 15
    assert len(patcher.top_level_control_matches(
        updated_interval, "writeInterval"
    )) == 1

    # Changing a root entry must not overwrite nested entries or a
    # commented-out root entry.
    with_root = "writeControl adjustableRunTime;\n" + complex_control
    updated_root = patcher.set_control_entry(
        with_root, "writeControl", "runTime", insert_if_missing=True
    )
    assert len(patcher.top_level_control_matches(
        updated_root, "writeControl"
    )) == 1
    assert updated_root.count("writeControl timeStep;") == 15
    assert updated_root.count("writeControl    runTime;") == 1

    commented = "/*\nwriteControl badSetting;\n*/\n" + complex_control
    assert len(patcher.top_level_control_matches(
        commented, "writeControl"
    )) == 0
    assert patcher.set_control_entry(
        commented, "writeControl", "runTime", insert_if_missing=True
    ).count("writeControl    runTime;") == 1

    # Other critical settings must remain required.
    try:
        patcher.set_control_entry(
            missing_write_control, "startTime", "31.2"
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Missing required setting was accepted")

    print("PASS: phase-pair insertion, preservation and duplicate checks")


def closure_area_density(
    h, l, b, grad_h, x, y, z,
    foam_area_multiplier, foam_diameter=0.003,
    min_alpha=0.02, full_alpha=0.10,
    max_bubble_alpha=0.10, max_area_density=2000
):
    if b >= max_bubble_alpha:
        return 0
    gate = (
        max(0, min(1, (h-min_alpha)/(full_alpha-min_alpha)))
      * max(0, min(1, (l-min_alpha)/(full_alpha-min_alpha)))
    )
    geom = gate * abs(grad_h)/(1 - min_alpha - full_alpha)
    foam = 0
    if (
        foam_area_multiplier > 0
        and 0.65 <= y <= 1.05
        and (x*x + z*z) ** 0.5 <= 0.30
    ):
        foam = gate * foam_area_multiplier * 6*h*l / foam_diameter
    return min(max(geom, foam), max_area_density)


def test_area_closure() -> None:
    # No transfer in pure liquid, pure headspace or bubble-rich regions.
    assert closure_area_density(0, 1, 0, 100, 0, 0.8, 0, 1) == 0
    assert closure_area_density(1, 0, 0, 100, 0, 0.8, 0, 1) == 0
    assert closure_area_density(.5, .3, .2, 100, 0, 0.8, 0, 1) == 0

    # Resolved-only baseline contains no area from homogeneous foam.
    assert closure_area_density(.5, .5, 0, 0, 0, 0.8, 0, 0) == 0
    # Foam geometry is strictly bounded to the upper impeller region.
    assert closure_area_density(.5, .5, 0, 0, 0, 0.8, 0, 1) > 0
    assert closure_area_density(.5, .5, 0, 0, 0, 0.3, 0, 1) == 0
    assert closure_area_density(.5, .5, 0, 0, 0.5, 0.8, 0, 1) == 0

    # MAX rather than a sum means no resolved/foam area double count.
    assert abs(
        closure_area_density(.5, .5, 0, 500, 0, 0.8, 0, 1) - 500/0.88
    ) < 1e-9
    assert closure_area_density(.5, .5, 0, 100, 0, 0.8, 0, 1) == 500
    assert closure_area_density(.5, .5, 0, 999999, 0, 0.8, 0, 1) == 2000

    # For a transition alphaH=z/L, integral |grad alpha|dz ~ unity.
    # The smooth gates reduce this slightly, in a reproducible way.
    samples = 1000
    dz = 0.1/samples
    integrated_area = sum(
        closure_area_density(
            (i+.5)/samples,
            1-(i+.5)/samples,
            0, 10, 0, 0.8, 0, 0
        )*dz
        for i in range(samples)
    )
    assert abs(integrated_area - 1) < 0.002
    print(f"PASS: area gates, bounded foam enhancement; planar area={integrated_area:.5f}")


if __name__ == "__main__":
    test_installer()
    test_area_closure()
    print("PASS: Python preflight complete; C++ build/runtime validation still required")
