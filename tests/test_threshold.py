"""Tests for threshold computation."""

import numpy as np
import pytest
from gsa_cpd.core.threshold import ThresholdType, compute_threshold


class TestThresholdFormulas:
    def test_chebyshev_formula(self):
        h = compute_threshold(0.0, 1.0, 0.01, ThresholdType.CHEBYSHEV)
        assert h == pytest.approx(10.0, rel=0.01)  # 0 + 1/sqrt(0.01)

    def test_vp_formula(self):
        h = compute_threshold(0.0, 1.0, 0.01, ThresholdType.VP)
        expected = (2.0 / 3.0) / np.sqrt(0.01)
        assert h == pytest.approx(expected, rel=0.01)

    def test_cantelli_formula(self):
        h = compute_threshold(0.0, 1.0, 0.01, ThresholdType.CANTELLI)
        expected = np.sqrt(1.0 / 0.01 - 1.0)
        assert h == pytest.approx(expected, rel=0.01)

    def test_chebyshev_most_conservative(self):
        E, V, eps = 1.0, 2.0, 0.01
        h_pe = compute_threshold(E, V, eps, ThresholdType.CHEBYSHEV)
        h_vp = compute_threshold(E, V, eps, ThresholdType.VP)
        assert h_pe > h_vp

    def test_threshold_scale(self):
        h1 = compute_threshold(0, 1, 0.01, ThresholdType.CHEBYSHEV, threshold_scale=1.0)
        h2 = compute_threshold(0, 1, 0.01, ThresholdType.CHEBYSHEV, threshold_scale=2.0)
        assert h2 == pytest.approx(2 * h1)

    def test_simulation_raises(self):
        with pytest.raises(ValueError):
            compute_threshold(0, 1, 0.01, ThresholdType.SIMULATION)
