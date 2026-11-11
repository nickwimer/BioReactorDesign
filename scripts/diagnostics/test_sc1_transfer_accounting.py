#!/usr/bin/env python3
"""Analytic SC1 transfer-accounting counterexamples, NOT CFD validation.

These tests evaluate small mathematical examples motivated by the reviewed
BiRD/source equations. They do not import, compile, execute or validate an
OpenFOAM closure, its discretization, its MPI behavior or the Kestrel cases.

E and D below are hypothetical GROSS POST-BLENDING entrainment/deposition
mass rates. They are not the 2 g/s raw entrainment cap or measured SC1 rates.
"""

import math
import unittest


def gross_species_gain(entrainment, deposition, headspace_y, bubble_y):
    """O2 rate into bubbles if both gross transfer channels are retained."""
    return entrainment * headspace_y - deposition * bubble_y


def net_species_gain(entrainment, deposition, headspace_y, bubble_y):
    """O2 rate if only the net mass transfer selects the donor phase."""
    net = entrainment - deposition
    return max(net, 0) * headspace_y + min(net, 0) * bubble_y


def pivot_weights(pivots, nucleus_volume):
    """In-range triangular eta algebra; independent mathematical example.

    Motivated by populationBalanceModel.C eta(), not a call into that code.
    Out-of-range behavior is deliberately excluded from this example.
    """
    if not pivots[0] <= nucleus_volume <= pivots[-1]:
        raise ValueError("This analytic example requires an in-range nucleus")
    weights = [0.0] * len(pivots)
    for i, pivot in enumerate(pivots):
        if nucleus_volume == pivot:
            weights[i] = 1.0
            return weights
    for i, (lower, upper) in enumerate(zip(pivots, pivots[1:])):
        if lower < nucleus_volume < upper:
            weights[i] = (upper - nucleus_volume) / (upper - lower)
            weights[i + 1] = (nucleus_volume - lower) / (upper - lower)
            return weights
    raise ValueError("Pivots must be strictly increasing")


class AnalyticTransferAccountingTests(unittest.TestCase):
    def assert_counterexchange_identity(self, entrainment, deposition, yh, yb):
        actual = gross_species_gain(entrainment, deposition, yh, yb)
        net_only = net_species_gain(entrainment, deposition, yh, yb)
        omitted = min(entrainment, deposition) * (yh - yb)
        self.assertAlmostEqual(actual - net_only, omitted, places=15)

    def test_net_entrainment_omits_counterexchange(self):
        self.assert_counterexchange_identity(0.002, 0.001, 0.23, 0.20)

    def test_net_deposition_omits_counterexchange(self):
        self.assert_counterexchange_identity(0.001, 0.002, 0.23, 0.20)

    def test_equal_gross_rates_can_exchange_oxygen_at_zero_net_mass(self):
        e = d = 0.001  # Hypothetical 1 g/s each direction.
        self.assertEqual(e - d, 0)
        self.assertEqual(net_species_gain(e, d, 0.23, 0.20), 0)
        self.assertAlmostEqual(gross_species_gain(e, d, 0.23, 0.20), 30e-6)
        self.assert_counterexchange_identity(e, d, 0.23, 0.20)

    def test_equal_compositions_remove_counterexchange_difference(self):
        for e, d in ((0.002, 0.001), (0.001, 0.002), (0.001, 0.001)):
            with self.subTest(entrainment=e, deposition=d):
                self.assertAlmostEqual(
                    gross_species_gain(e, d, 0.23, 0.23),
                    net_species_gain(e, d, 0.23, 0.23),
                )

    def test_one_way_transfer_is_unchanged(self):
        for e, d in ((0.002, 0), (0, 0.002), (0, 0)):
            with self.subTest(entrainment=e, deposition=d):
                self.assertAlmostEqual(
                    gross_species_gain(e, d, 0.23, 0.20),
                    net_species_gain(e, d, 0.23, 0.20),
                )

    def test_counterexchange_changes_sign_with_composition_contrast(self):
        self.assert_counterexchange_identity(0.001, 0.002, 0.18, 0.23)
        self.assertLess(
            gross_species_gain(0.001, 0.002, 0.18, 0.23)
            - net_species_gain(0.001, 0.002, 0.18, 0.23), 0
        )

    def test_both_models_conserve_phase_summed_oxygen_algebraically(self):
        # Omitting counterexchange is a model-form error; it need not produce
        # an overall oxygen conservation residual. Compute both donor/receiver
        # contributions explicitly rather than simply negating their sum.
        e, d, yh, yb = 0.002, 0.001, 0.23, 0.20
        b_gross = e * yh - d * yb
        h_gross = d * yb - e * yh
        net = e - d
        b_net = max(net, 0) * yh + min(net, 0) * yb
        h_net = -max(net, 0) * yh - min(net, 0) * yb
        self.assertAlmostEqual(b_gross + h_gross, 0, places=15)
        self.assertAlmostEqual(b_net + h_net, 0, places=15)
        self.assertNotAlmostEqual(b_gross, b_net, places=8)

    def test_deposition_only_normalization_supplies_bin_removal(self):
        # Homogeneous, fixed-density, no transport, breakup, coalescence,
        # drift, residual-alpha stabilization or explicit bin-removal source.
        # Phase alpha falls while d(alpha*f_i)/dt=0 for each raw PB class.
        alpha_old, alpha_new = 0.10, 0.09
        f_old = (0.25, 0.75)
        raw = [alpha_old * f / alpha_new for f in f_old]
        self.assertGreater(sum(raw), 1)
        raw_volume = alpha_new * sum(raw)
        self.assertAlmostEqual(raw_volume, alpha_old)
        normalized = [f / sum(raw) for f in raw]
        corrected_volume = alpha_new * sum(normalized)
        self.assertAlmostEqual(corrected_volume, alpha_new)
        self.assertAlmostEqual(corrected_volume - raw_volume, -0.01)
        for old, corrected in zip(f_old, normalized):
            self.assertAlmostEqual(old, corrected)

    def test_exact_pivot_birth_matches_volume_and_number(self):
        pivots = [math.pi * d**3 / 6 for d in (0.0015, 0.002, 0.0025)]
        nucleus = pivots[1]  # The documented 2 mm nucleus is a class pivot.
        weights = pivot_weights(pivots, nucleus)
        volume_rate = 1e-4  # Hypothetical S_positive/rho_B, m3/m3/s.
        number_rates = [eta * volume_rate / x for eta, x in zip(weights, pivots)]
        self.assertAlmostEqual(sum(x * j for x, j in zip(pivots, number_rates)), volume_rate)
        self.assertAlmostEqual(sum(number_rates) / (volume_rate / nucleus), 1)

    def test_off_pivot_volume_allocation_does_not_preserve_nucleus_number(self):
        pivots = [math.pi * d**3 / 6 for d in (0.0015, 0.002)]
        nucleus = math.pi * 0.00175**3 / 6
        weights = pivot_weights(pivots, nucleus)
        volume_rate = 1e-4
        number_rates = [eta * volume_rate / x for eta, x in zip(weights, pivots)]
        self.assertAlmostEqual(sum(weights), 1)
        self.assertAlmostEqual(sum(x * j for x, j in zip(pivots, number_rates)), volume_rate)
        ratio = sum(number_rates) / (volume_rate / nucleus)
        self.assertAlmostEqual(ratio, 1.1940737123842593)
        # A number-based adjacent-pivot allocation can preserve both moments.
        alternative = [eta * volume_rate / nucleus for eta in weights]
        self.assertAlmostEqual(sum(alternative) / (volume_rate / nucleus), 1)
        self.assertAlmostEqual(sum(x * j for x, j in zip(pivots, alternative)), volume_rate)


if __name__ == "__main__":
    print("Analytic counterexamples only; no OpenFOAM/CFD implementation validation.", flush=True)
    unittest.main(verbosity=2)
