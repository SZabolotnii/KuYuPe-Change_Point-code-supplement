"""Tests for the competitor baselines (SCUSUM, KQTEWMA).

These tests are web-independent and fast (target total runtime < ~60s). They
check the closed-form Hyvarinen score, the CUSUM/EWMA recursion invariants, the
package API (fit/predict/reset/alarm_time), and that each detector calibrated to
a target ARL_0 on i.i.d. H0 data achieves an empirical ARL_0 in the right
ballpark (loose tolerance, small N for speed).
"""

import numpy as np
import scipy.stats as stats
import pytest

from gsa_cpd.baselines import SCUSUM, KQTEWMA, gaussian_hyvarinen_score


# --------------------------------------------------------------------------- #
# SCUSUM: Hyvarinen score                                                      #
# --------------------------------------------------------------------------- #
class TestGaussianHyvarinenScore:
    def test_closed_form_matches_finite_difference(self):
        """S_H(x, N(mu,sigma^2)) = 0.5*(x-mu)^2/sigma^4 - 1/sigma^2 numerically."""
        mu, sigma = 0.7, 1.3
        xs = np.array([-1.0, 0.0, 0.5, 2.0, 3.5])

        def logp(x):
            return stats.norm.logpdf(x, mu, sigma)

        h = 1e-4
        d1 = (logp(xs + h) - logp(xs - h)) / (2 * h)
        d2 = (logp(xs + h) - 2 * logp(xs) + logp(xs - h)) / (h * h)
        fd = 0.5 * d1**2 + d2
        closed = gaussian_hyvarinen_score(xs, mu, sigma)
        assert np.allclose(fd, closed, atol=1e-5)

    def test_explicit_formula(self):
        """Match the explicit algebraic form directly."""
        mu, sigma = -0.3, 0.8
        xs = np.array([-2.0, 0.0, 1.0])
        expected = 0.5 * (xs - mu) ** 2 / sigma**4 - 1.0 / sigma**2
        assert np.allclose(gaussian_hyvarinen_score(xs, mu, sigma), expected)

    def test_score_difference_drift_signs(self):
        """Score difference has negative drift under H0, positive under H1.

        With pre-change N(0,1) and post-change N(0.5,1), the SCUSUM instantaneous
        statistic S_H(.,p_inf) - S_H(.,p_1) must have E_inf[.] < 0 and E_1[.] > 0
        for the CUSUM to behave correctly. Moreover, for equal variances the
        Hyvarinen/Fisher relation gives E_inf[.] = -0.5*((mu1-mu0)/sigma^2)^2.
        """
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.5, "sigma": 1.0})
        x0 = stats.norm.rvs(0.0, 1.0, size=200000, random_state=1)
        x1 = stats.norm.rvs(0.5, 1.0, size=200000, random_state=2)
        e0 = float(np.mean(det.score_difference(x0)))
        e1 = float(np.mean(det.score_difference(x1)))
        assert e0 < 0.0
        assert e1 > 0.0
        # Fisher-divergence identity (equal sigma): E_inf = -0.5*(Delta/sigma^2)^2.
        assert e0 == pytest.approx(-0.5 * (0.5) ** 2, abs=0.01)


# --------------------------------------------------------------------------- #
# SCUSUM: recursion invariants and lambda selection                            #
# --------------------------------------------------------------------------- #
class TestSCUSUMRecursion:
    def test_fit_returns_self(self):
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.6, "sigma": 1.0})
        h0 = stats.norm.rvs(0.0, 1.0, size=5000, random_state=0)
        assert det.fit(h0) is det

    def test_lambda_satisfies_mgf_constraint(self):
        """Selected lambda keeps the empirical H0 MGF E_inf[exp(z_lam)] <= 1."""
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.6, "sigma": 1.0})
        h0 = stats.norm.rvs(0.0, 1.0, size=20000, random_state=3)
        det.fit(h0)
        assert det.lam > 0.0
        mgf = float(np.mean(np.exp(det.lam * det.score_difference(h0))))
        assert mgf <= 1.0 + 1e-6

    def test_statistic_nonnegative_and_resets(self):
        """Z(n) = max(0, .) is non-negative; reset returns it to zero."""
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.6, "sigma": 1.0}, lam=0.2)
        det.fit(np.zeros(10))  # lam already set; just resets state
        # threshold defaults to inf -> never alarm; just track the statistic.
        rng = np.random.default_rng(0)
        for x in rng.normal(size=500):
            det.predict(float(x))
            assert det._g >= 0.0
        det.reset()
        assert det._g == 0.0
        assert det.alarm_time is None

    def test_alarm_time_set_on_shift(self):
        """A clear post-change mean shift triggers an alarm and records the time."""
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 1.0, "sigma": 1.0}, lam=0.3)
        det.fit(np.zeros(10))
        det._threshold = 5.0
        rng = np.random.default_rng(7)
        detected = any(det.predict(float(x)) for x in rng.normal(2.0, 1.0, size=300))
        assert detected
        assert det.alarm_time is not None

    def test_calibrated_arl0_ballpark(self):
        """SCUSUM calibrated to ARL_0=300 on N(0,1) achieves it within ~30%."""
        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.6, "sigma": 1.0})

        def h0_sampler(n):
            return stats.norm.rvs(0.0, 1.0, size=n)

        det.fit(h0_sampler(20000))
        target = 300.0
        det.calibrate(h0_sampler, target_arl0=target, stream_len=2500,
                      n_runs=600, tol=0.05)

        # Independent verification of the achieved ARL_0.
        rng = np.random.default_rng(123)
        rls = []
        for _ in range(600):
            det.reset()
            stream = rng.normal(0.0, 1.0, size=2500)
            rl = 2500
            for t, x in enumerate(stream):
                if det.predict(float(x)):
                    rl = t + 1
                    break
            rls.append(rl)
        emp = float(np.mean(rls))
        assert 0.7 * target <= emp <= 1.3 * target


# --------------------------------------------------------------------------- #
# KQT-EWMA                                                                     #
# --------------------------------------------------------------------------- #
class TestKQTEWMA:
    def test_fit_returns_self(self):
        det = KQTEWMA(n_bins=8, lam=0.05)
        assert det.fit(stats.norm.rvs(size=3000, random_state=0)) is det

    def test_statistic_nonnegative_and_in_range(self):
        """T_t = sum (Z_j - pi_hat_j)^2 / pi_hat_j is finite and >= 0."""
        det = KQTEWMA(n_bins=8, lam=0.05)
        det.fit(stats.norm.rvs(size=5000, random_state=0))
        rng = np.random.default_rng(1)
        for x in rng.normal(size=500):
            det.predict(float(x))
            T = det.statistic()
            assert T >= 0.0 and np.isfinite(T)
        # Z is a convex combination of indicators and the initial pi_hat, so it
        # stays within [0, 1] componentwise.
        assert np.all(det.Z >= -1e-12) and np.all(det.Z <= 1.0 + 1e-12)

    def test_bins_equiprobable_under_h0(self):
        """Empirical bin membership under H0 is approximately 1/K per bin."""
        K = 8
        det = KQTEWMA(n_bins=K, lam=0.05)
        det.fit(stats.norm.rvs(size=20000, random_state=2))
        test = stats.norm.rvs(size=40000, random_state=3)
        counts = np.zeros(K)
        for x in test:
            counts += det._bin_indicators(float(x))
        freqs = counts / counts.sum()
        assert np.allclose(freqs, 1.0 / K, atol=0.03)

    def test_reset_restores_initial_state(self):
        det = KQTEWMA(n_bins=6, lam=0.1)
        det.fit(stats.norm.rvs(size=3000, random_state=4))
        for x in stats.norm.rvs(size=100, random_state=5):
            det.predict(float(x))
        det.reset()
        assert np.allclose(det.Z, det.pi_hat)
        assert det.statistic() == pytest.approx(0.0, abs=1e-12)
        assert det.alarm_time is None

    def test_calibrated_arl0_ballpark(self):
        """KQT-EWMA calibrated to ARL_0=300 on N(0,1) achieves it within ~30%."""
        det = KQTEWMA(n_bins=8, lam=0.05)

        def h0_sampler(n):
            return stats.norm.rvs(0.0, 1.0, size=n)

        det.fit(h0_sampler(10000))
        target = 300.0
        det.calibrate(h0_sampler, target_arl0=target, stream_len=2500,
                      n_runs=600, tol=0.05)

        rng = np.random.default_rng(321)
        rls = []
        for _ in range(600):
            stream = rng.normal(0.0, 1.0, size=2500)
            rls.append(det._run_length(stream, det._threshold, 2500))
        emp = float(np.mean(rls))
        assert 0.7 * target <= emp <= 1.3 * target
