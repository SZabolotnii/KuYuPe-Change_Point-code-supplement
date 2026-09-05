"""The clip at +-phi_max makes an unstandardised detector inert.

`evaluate_basis_*` clips every basis value to [-phi_max, phi_max].  On a series
living outside that window the clipped dictionary is a constant, so it carries
no information and the CUSUM statistic can never leave zero: the detector cannot
alarm at any threshold.  That is silence from inability, not from confidence,
and the two are indistinguishable from the outside unless the detector says so.

These tests pin the defect, the warning and the fix.
"""

import numpy as np
import pytest

from gsa_cpd import GSADetector, BasisType


def _off_scale_series(offset=850.0, scale=40.0, n=600, seed=5):
    rng = np.random.default_rng(seed)
    return offset + scale * rng.standard_t(df=6, size=n)


def test_unstandardised_off_scale_series_is_flagged_and_warns():
    data = _off_scale_series()
    detector = GSADetector(basis=BasisType.POLY, degree=1, standardize=None)

    with pytest.warns(RuntimeWarning, match="clipped"):
        detector.fit(data, delta=0.5)

    assert detector.diagnostics.basis_degenerate
    llr = {detector._compute_llr(float(x)) for x in data}
    assert len(llr) == 1, "expected a constant LLR from a constant dictionary"


def test_standardisation_restores_the_dictionary():
    data = _off_scale_series()
    detector = GSADetector(basis=BasisType.POLY, degree=1, standardize="robust")
    detector.fit(data, delta=0.5)

    assert not detector.diagnostics.basis_degenerate
    assert detector.diagnostics.standardize_scale == pytest.approx(
        1.4826 * np.median(np.abs(data - np.median(data))), rel=1e-12
    )
    llr = np.array([detector._compute_llr(float(x)) for x in data])
    assert np.ptp(llr) > 0.0


def test_standardised_detector_can_alarm_where_the_raw_one_cannot():
    """The witness is the peak of the statistic, not the alarm itself."""
    calib = _off_scale_series()
    rng = np.random.default_rng(11)
    post = 850.0 + 40.0 * (2.5 + rng.standard_t(df=6, size=400))

    def peak(standardize):
        det = GSADetector(basis=BasisType.POLY, degree=1, standardize=standardize)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                det.fit(calib, delta=0.5)
        g = m = 0.0
        for x in post:
            g = max(0.0, g + det._compute_llr(float(x)))
            m = max(m, g)
        return m

    assert peak(None) == 0.0, (
        "the unstandardised detector should be inert on this series; if it is "
        "not, the defect this test documents has changed shape"
    )
    assert peak("robust") > 0.0


def test_standardize_none_is_the_identity():
    """The default is unchanged, so the Monte-Carlo study still reproduces."""
    detector = GSADetector(standardize=None)
    assert detector.standardize is None
    detector.fit(np.random.default_rng(3).normal(size=500))
    assert detector._loc == 0.0 and detector._scale == 1.0
    assert not detector.diagnostics.basis_degenerate


def test_rejects_unknown_standardiser():
    with pytest.raises(ValueError, match="standardize"):
        GSADetector(standardize="whiten")
