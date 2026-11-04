#!/usr/bin/env python3
"""No-write regression tests for SCRATCH-only entrainment+PB integration.

These tests do not certify the C++ model or replace the required Kestrel
compile/scratch runtime conservation and phase-fraction checks.
"""

import importlib.util
import math
from pathlib import Path


def load(relative: str, alias: str):
    path = Path(__file__).with_name(relative)
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not find test file {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def case_fixture() -> str:
    original_test = load(
        "test_sc1_headspace_o2_trial.py", "sc1_existing_film_test"
    )
    original = original_test.synthetic_phase_properties().replace(
        "surface headspaceGas;", "surface headspaceGas; efficiency 1;"
    )
    return original + """
bubbleGas
{
    diameterModel velocityGroup;
    velocityGroupCoeffs
    {
        populationBalance bubbles;
        sizeGroups
        (
            f1 { dSph 0.001; value 0; }
            f2 { dSph 0.0015; value 0; }
            f3 { dSph 0.002; value 1; }
            f4 { dSph 0.0025; value 0; }
            f5 { dSph 0.011; value 0; }
        );
    }
}
populationBalanceCoeffs
{
    bubbles
    {
        nucleationModels
        (
        );
    }
}
"""


def test_installer() -> None:
    trial = load(
        "prepare_sc1_headspace_entrainment_trial.py", "sc1_reentrain"
    )
    base = trial.load_headspace_helper()
    original = case_fixture()

    control = trial.install_entrainment_models(
        original, efficiency=0, maximum_rate_kg_s=0.002,
        nucleation_diameter_m=0.002,
    )
    active = trial.install_entrainment_models(
        original, efficiency=0.01, maximum_rate_kg_s=0.002,
        nucleation_diameter_m=0.002,
    )

    for output in (control, active):
        assert output.count("type headspaceEntrainingDeposition;") == 1
        assert output.count("headspaceEntrainmentNucleation") == 1
        assert output.count("type deposition;") == 0
        assert output.count("type resolvedHeadspaceFilm;") == 2
        assert output.count("type Henry;") == 2
        assert output.count("k ($He_O2);") == 2
        assert output.count("species (O2);") == 2
        assert output.count("droplet bubbleGas;") == 1
        assert output.count("surface headspaceGas;") == 1
        assert output.count("efficiency 1;") == 1
        assert "nucleationDiameter 0.002;" in output
        assert "velocityGroup bubbleGas;" in output
        assert "headspacePhase headspaceGas;" in output
        assert "maxGlobalEntrainedKgS 0.002;" in output
        assert output.count("foamAreaMultiplier 0;") == 2

        for key in (
            "phaseTransfer", "nucleationModels",
            "interfaceComposition.liquid",
            "diffusiveMassTransfer.liquid",
        ):
            bounds = base.find_keyword_list(output, key)
            assert bounds is not None

    assert "entrainmentEfficiency 0;" in control
    assert "entrainmentEfficiency 0.01;" in active

    # Reject existing transfer, unphysical parameters and size targets.
    bad_inputs = [
        dict(efficiency=-0.01, maximum_rate_kg_s=0.002,
             nucleation_diameter_m=0.002),
        dict(efficiency=2, maximum_rate_kg_s=0.002,
             nucleation_diameter_m=0.002),
        dict(efficiency=0.01, maximum_rate_kg_s=0,
             nucleation_diameter_m=0.002),
        dict(efficiency=0.01, maximum_rate_kg_s=1,
             nucleation_diameter_m=0.002),
        dict(efficiency=0.01, maximum_rate_kg_s=0.002,
             nucleation_diameter_m=0.015),
        dict(efficiency=math.nan, maximum_rate_kg_s=0.002,
             nucleation_diameter_m=0.002),
    ]
    for kwargs in bad_inputs:
        try:
            trial.install_entrainment_models(original, **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Invalid trial input accepted: {kwargs}")

    try:
        trial.install_entrainment_models(
            active, efficiency=0.01,
            maximum_rate_kg_s=0.002, nucleation_diameter_m=0.002,
        )
    except ValueError:
        pass
    else:
        raise AssertionError("Repeated model installation not rejected")

    print("PASS: separate phase-transfer and PB size-injection configuration")


def gate(alpha: float, min_alpha=0.02, full_alpha=0.10) -> float:
    return min(max((alpha-min_alpha)/(full_alpha-min_alpha), 0), 1)


def candidate(
    h, l, b, gradient, relative_speed, liquid_y, y, r,
    entrainment_efficiency=0.01, rho_h=1.2,
    max_slip=5, donor_time=0.05, down_only=True
) -> float:
    """Pure Python mirror of the C++ per-cell bounded source."""
    if (y < 0.68 or y > 0.96 or r < 0.05 or r > 0.30
            or b >= 0.1 or (down_only and liquid_y >= 0)):
        return 0
    h = min(max(h, 0), 1)
    l = min(max(l, 0), 1)
    interfacial_area = gate(h)*gate(l)*abs(gradient)
    raw = rho_h*entrainment_efficiency*interfacial_area*min(
        abs(relative_speed), max_slip
    )
    donor = rho_h*max(h-0.02, 0)/donor_time
    return min(raw, donor)


def test_bounds() -> None:
    kw = dict(h=0.5, l=0.5, b=0, gradient=10,
              relative_speed=2, liquid_y=-1, y=0.8, r=0.15)
    assert candidate(**kw) > 0
    assert candidate(**dict(kw, y=1.2)) == 0
    assert candidate(**dict(kw, r=0.40)) == 0
    assert candidate(**dict(kw, b=0.2)) == 0
    assert candidate(**dict(kw, h=0)) == 0
    assert candidate(**dict(kw, l=0)) == 0
    assert candidate(**dict(kw, liquid_y=0.1)) == 0
    assert candidate(**dict(kw, entrainment_efficiency=0)) == 0
    assert candidate(**dict(kw, relative_speed=500)) <= candidate(
        **dict(kw, relative_speed=5)
    ) + 1e-12
    assert candidate(**dict(kw, entrainment_efficiency=1)) <= (
        1.2*(0.5-0.02)/0.05
    ) + 1e-12

    # In C++ the global cap rescales all positive candidate cells by
    # min(1, global_cap/sum(R_raw V)). The net bubble gain after local
    # deposition and BiRD blending is then assigned to PB size bins.
    rates = [0.004, 0.002, 0.008]
    cap = 0.002
    scale = min(1.0, cap/sum(rates))
    assert abs(sum(r*scale for r in rates)-cap) < 1e-12

    # Positive bubble phase dmdtf is entrainment and cannot also be
    # counted as deposition; only max(0, net bubble gain) enters PB.
    for net in (-0.001, 0, 0.002):
        pb = max(0, net)
        assert pb >= 0
        assert pb == 0 if net <= 0 else pb == net

    print("PASS: bounded interface/downflow entrainment and global cap")


def test_kinematic_diagnostics() -> None:
    """Check the sign of the diagnostic normal-penetration proxy.

    The unit normal points from liquid toward headspace: grad(alphaH).
    Relative gas penetration into liquid corresponds to
    (U_liquid - U_headspace) dot n > 0. This is NOT a validated
    physical breakup/entrainment law.
    """
    def dot(a, b):
        return sum(x*y for x, y in zip(a, b))

    n_liquid_to_gas = (0, 1, 0)
    liquid_static = (0, 0, 0)
    gas_downward = (0, -1, 0)
    gas_upward = (0, 1, 0)

    gas_to_liquid = dot(
        tuple(a-b for a, b in zip(liquid_static, gas_downward)),
        n_liquid_to_gas,
    )
    gas_away_from_liquid = dot(
        tuple(a-b for a, b in zip(liquid_static, gas_upward)),
        n_liquid_to_gas,
    )

    assert gas_to_liquid > 0
    assert gas_away_from_liquid < 0

    # A downward-moving LIQUID is neither a necessary nor a
    # sufficient condition for normal gas penetration: interfacial
    # geometry and relative gas velocity matter.
    moving_up_liquid = (0, 0.5, 0)
    faster_upward_gas = (0, 0.75, 0)
    assert moving_up_liquid[1] > 0
    assert dot(
        tuple(a-b for a, b in zip(moving_up_liquid, gas_downward)),
        n_liquid_to_gas,
    ) > 0
    assert dot(
        tuple(a-b for a, b in zip(moving_up_liquid, faster_upward_gas)),
        n_liquid_to_gas,
    ) < 0

    print("PASS: candidate gas-into-liquid normal direction and velocity gates")


if __name__ == "__main__":
    test_installer()
    test_bounds()
    test_kinematic_diagnostics()
    print("PASS: isolated Python preflight; C++ and 52-rank runtime still require separate validation")
