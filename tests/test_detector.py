"""Tests for GSADetector: fit/predict, all basis types."""

import numpy as np
import pytest
from gsa_cpd import GSADetector, BasisType, ThresholdType


@pytest.fixture
def gaussian_data():
    rng = np.random.default_rng(42)
    return rng.normal(0, 1, size=1000)


class TestGSADetectorBasic:
    def test_fit_returns_self(self, gaussian_data):
        det = GSADetector(degree=2, epsilon=0.01)
        result = det.fit(gaussian_data)
        assert result is det

    def test_predict_before_fit_raises(self):
        det = GSADetector()
        with pytest.raises(RuntimeError):
            det.predict(0.0)

    def test_diagnostics_populated(self, gaussian_data):
        det = GSADetector(degree=2).fit(gaussian_data)
        d = det.diagnostics
        assert d.condition_number > 0
        assert d.threshold > 0
        assert d.coeffs is not None
        assert len(d.coeffs) == 2
        assert d.solver_method in ("direct", "ridge", "svd")

    def test_no_alarm_on_short_h0(self, gaussian_data):
        """On a short H0 sequence, PE threshold should prevent alarm."""
        det = GSADetector(
            degree=2, epsilon=0.001, threshold_scale=2.0
        ).fit(gaussian_data[:500])
        rng = np.random.default_rng(123)
        test_data = rng.normal(0, 1, size=30)
        detected = any(det.predict(x) for x in test_data)
        assert not detected, "Alarm should not trigger on 30 H0 samples"

    def test_alarm_on_shift(self, gaussian_data):
        det = GSADetector(degree=2, epsilon=0.01).fit(gaussian_data[:500])
        # Strong mean shift should trigger alarm
        rng = np.random.default_rng(99)
        shifted = rng.normal(3.0, 1, size=200)
        detected = any(det.predict(x) for x in shifted)
        assert detected


class TestAllBasisTypes:
    @pytest.mark.parametrize("basis", list(BasisType))
    def test_fit_predict_all_bases(self, gaussian_data, basis):
        det = GSADetector(basis=basis, degree=2, epsilon=0.01)
        det.fit(gaussian_data[:500])
        assert det.diagnostics.threshold > 0
        # Should not crash on predict
        det.predict(0.5)
        det.predict(-1.2)
        det.predict(5.0)


class TestThresholdTypes:
    @pytest.mark.parametrize("ttype", [ThresholdType.CHEBYSHEV,
                                        ThresholdType.VP,
                                        ThresholdType.CANTELLI])
    def test_analytic_thresholds(self, gaussian_data, ttype):
        det = GSADetector(degree=2, epsilon=0.01, threshold_type=ttype)
        det.fit(gaussian_data[:500])
        assert det.diagnostics.threshold > 0

    def test_chebyshev_most_conservative(self, gaussian_data):
        thresholds = {}
        for ttype in [ThresholdType.CHEBYSHEV, ThresholdType.VP, ThresholdType.CANTELLI]:
            det = GSADetector(degree=2, epsilon=0.01, threshold_type=ttype)
            det.fit(gaussian_data[:500])
            thresholds[ttype] = det.diagnostics.threshold
        # Chebyshev should give the highest (most conservative) threshold
        assert thresholds[ThresholdType.CHEBYSHEV] >= thresholds[ThresholdType.VP]


class TestReset:
    def test_reset_clears_state(self, gaussian_data):
        det = GSADetector(degree=2).fit(gaussian_data[:500])
        det.predict(10.0)  # likely triggers internal state change
        det.reset()
        assert det.alarm_time is None
