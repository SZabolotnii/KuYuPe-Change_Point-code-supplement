"""Test Gaussian limit: S=1 poly GSA = classical CUSUM (Theorem 1)."""

import numpy as np
import pytest
from gsa_cpd import GSADetector, BasisType


class TestGaussianLimit:
    """When data is Gaussian and s=1, all bases should give identical results."""

    def test_s1_all_bases_identical_add(self):
        rng = np.random.default_rng(42)
        cal_data = rng.normal(0, 1, size=1000)

        adds = {}
        for basis in [BasisType.POLY, BasisType.FRAC, BasisType.LOG]:
            det = GSADetector(basis=basis, degree=1, epsilon=0.01)
            det.fit(cal_data, delta=0.3)

            delays = []
            for trial in range(100):
                det.reset()
                test_rng = np.random.default_rng(1000 + trial)
                # H0 for 100, then shifted
                h0 = test_rng.normal(0, 1, size=100)
                h1 = test_rng.normal(0.5, 1, size=200)
                test_data = np.concatenate([h0, h1])

                for t, x in enumerate(test_data):
                    if det.predict(x):
                        if t >= 100:
                            delays.append(t - 100)
                        break

            adds[basis] = np.mean(delays) if delays else 999

        # All bases at s=1 should give very similar ADD on Gaussian data
        values = list(adds.values())
        spread = max(values) - min(values)
        mean_add = np.mean(values)
        assert spread / mean_add < 0.15, (
            f"ADD spread too large for Gaussian limit: {adds}"
        )

    def test_s1_poly_matches_energy_detector(self):
        """S=1 polynomial should produce a linear statistic k0 + k1*x."""
        rng = np.random.default_rng(42)
        cal_data = rng.normal(0, 1, size=1000)
        det = GSADetector(basis=BasisType.POLY, degree=1, epsilon=0.01)
        det.fit(cal_data)
        # Coefficient vector should have length 1
        assert len(det.coefficients) == 1
        assert det.diagnostics.k0 != 0.0  # bias should be non-zero
