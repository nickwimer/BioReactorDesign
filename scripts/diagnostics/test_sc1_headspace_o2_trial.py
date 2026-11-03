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
    geom = gate * abs(grad_h)
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
    assert closure_area_density(.5, .5, 0, 500, 0, 0.8, 0, 1) == 500
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
    assert 0.75 < integrated_area < 1
    print(f"PASS: area gates, bounded foam enhancement; planar area={integrated_area:.5f}")


if __name__ == "__main__":
    test_installer()
    test_area_closure()
    print("PASS: Python preflight complete; C++ build/runtime validation still required")
