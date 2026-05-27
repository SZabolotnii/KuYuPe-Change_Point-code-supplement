"""Score-based CUSUM (SCUSUM) baseline detector.

Faithful implementation of the Score-based CUSUM detector of

    Suya Wu, Enmao Diao, Taposh Banerjee, Jie Ding, Vahid Tarokh,
    "Score-based Quickest Change Detection for Unnormalized Models,"
    Proc. 26th Int. Conf. on Artificial Intelligence and Statistics
    (AISTATS), PMLR vol. 206, pp. 10546-10565, 2023.
    https://proceedings.mlr.press/v206/wu23b.html

The detector is a CUSUM recursion driven by the *Hyvarinen score difference*
(score matching, Hyvarinen 2005), which avoids the (often intractable)
normalizing constants of the pre- and post-change densities.

Notation follows the paper (Def. 3, Eqs. 7-8 and Thm. 3):

  Hyvarinen score of a density p (1-D form used here):

      S_H(x, p) = 1/2 * (d/dx log p(x))^2 + d^2/dx^2 log p(x)

  Instantaneous statistic with pre-change density p_inf and post-change p_1:

      z_lambda(x) = lambda * ( S_H(x, p_inf) - S_H(x, p_1) ),   lambda > 0

  chosen so that the pre-change drift is negative; the paper requires
  E_inf[exp(z_lambda(X))] <= 1 (Eq. 10) for the false-alarm bound to hold.

  CUSUM recursion and stopping rule:

      Z(0) = 0,   Z(n) = max(0, Z(n-1) + z_lambda(X_n)),
      T = inf{ n >= 1 : Z(n) >= b }.

  False-alarm guarantee (Thm. 3): E_inf[T] >= exp(b). The paper's design rule
  for a target ARL_0 = gamma is therefore b = log(gamma); this implementation
  also supports direct Monte Carlo calibration of b to a target ARL_0, which is
  the convention used for the operating-characteristic comparison in this
  package (all methods calibrated to the same ARL_0 before comparing ADD).

This module provides the univariate Gaussian-model instantiation of the score
difference (the closed-form Hyvarinen score for a normal density), which is the
standard worked example in the paper.

API note: to match the other ``gsa_cpd`` baselines (``SignCUSUM``, ``EWMA``),
the online interface is ``fit(calibration_data) -> self``, ``predict(x) -> bool``,
``reset()`` and an ``alarm_time`` property. The pre-/post-change Gaussian
working model is set in the constructor; ``fit`` selects the score-difference
scale ``lambda`` from the calibration (H0) data via the MGF constraint.
"""

import numpy as np
from typing import Callable, Optional


def gaussian_hyvarinen_score(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """Hyvarinen score S_H(x, p) for a Gaussian density p = N(mu, sigma^2).

    For p(x) = N(mu, sigma^2):
        d/dx log p(x)   = -(x - mu) / sigma^2
        d^2/dx^2 log p  = -1 / sigma^2
    so
        S_H(x, p) = 1/2 * (x - mu)^2 / sigma^4 - 1 / sigma^2.

    Args:
        x: Observation(s).
        mu: Mean of the Gaussian density.
        sigma: Standard deviation (> 0).

    Returns:
        Hyvarinen score evaluated at x.
    """
    x = np.asarray(x, dtype=float)
    grad = -(x - mu) / sigma**2
    lap = -1.0 / sigma**2
    return 0.5 * grad**2 + lap


class SCUSUM:
    """Score-based CUSUM detector (Wu et al., AISTATS 2023).

    Maintains the CUSUM statistic Z(n) driven by the Hyvarinen score difference
    between the pre- and post-change Gaussian working models and raises an alarm
    when Z(n) crosses the threshold b. The score difference uses a Gaussian
    working model (a valid score-matching surrogate even when the data are
    non-Gaussian: the detector stays distribution-robust because lambda and the
    threshold are calibrated on the actual H0 data).

    Example::

        det = SCUSUM(pre_params={"mu": 0.0, "sigma": 1.0},
                     post_params={"mu": 0.6, "sigma": 1.0})
        det.fit(calibration_data)              # selects lambda from H0 data
        det.calibrate(h0_sampler, target_arl0=500, stream_len=2500)
        for x in stream:
            if det.predict(x):
                print(f"Change at t={det.alarm_time}")
                break

    Args:
        pre_params: Pre-change Gaussian params {"mu", "sigma"}.
        post_params: Post-change Gaussian params {"mu", "sigma"}.
        lam: Score-difference scaling lambda > 0. If None, selected by ``fit``
            via the MGF constraint on calibration data.
    """

    def __init__(
        self,
        pre_params: dict,
        post_params: dict,
        lam: Optional[float] = None,
    ):
        self._mu0 = float(pre_params["mu"])
        self._sig0 = float(pre_params["sigma"])
        self._mu1 = float(post_params["mu"])
        self._sig1 = float(post_params["sigma"])

        self.lam = lam
        self._threshold: float = np.inf

        # Runtime state.
        self._g: float = 0.0
        self._t: int = 0
        self._alarm_time: Optional[int] = None

    @property
    def alarm_time(self) -> Optional[int]:
        """Time index of the alarm (None if no alarm)."""
        return self._alarm_time

    def score_difference(self, x: np.ndarray) -> np.ndarray:
        """Score difference S_H(x, p_inf) - S_H(x, p_1) (without lambda)."""
        x = np.asarray(x, dtype=float)
        return (
            gaussian_hyvarinen_score(x, self._mu0, self._sig0)
            - gaussian_hyvarinen_score(x, self._mu1, self._sig1)
        )

    def instantaneous_statistic(self, x: np.ndarray) -> np.ndarray:
        """z_lambda(x) = lambda * (S_H(x, p_inf) - S_H(x, p_1))."""
        if self.lam is None:
            raise ValueError("lambda is not set; call fit() or pass lam.")
        return self.lam * self.score_difference(x)

    def fit(
        self,
        calibration_data: np.ndarray,
        lam_grid: Optional[np.ndarray] = None,
    ) -> "SCUSUM":
        """Select lambda > 0 satisfying the MGF constraint E_inf[exp(z_lambda)] <= 1.

        Per Eq. (10) of the paper, lambda must keep the pre-change drift negative
        and the exponentiated statistic a non-positive supermartingale increment.
        We estimate E_inf[exp(lambda * d(X))] on H0 calibration data and pick the
        largest lambda on a grid for which the empirical MGF stays <= 1 (larger
        lambda gives a stronger post-change drift and thus shorter delay). This
        mirrors the data-driven choice discussed in the paper when the constant
        is not known in closed form. If ``lam`` was given to the constructor it
        is kept and only the runtime state is reset.

        Args:
            calibration_data: Pre-change (H0) samples used to estimate the MGF.
            lam_grid: Candidate lambda values (default: geometric grid).

        Returns:
            self (for chaining).
        """
        if self.lam is None:
            d = self.score_difference(calibration_data)
            if lam_grid is None:
                lam_grid = np.geomspace(1e-3, 5.0, 60)

            best = None
            for lam in lam_grid:
                # Empirical MGF E_inf[exp(lambda * d(X))]; numerically guarded.
                mgf = float(np.mean(np.exp(lam * d)))
                if np.isfinite(mgf) and mgf <= 1.0:
                    best = lam
            # Fall back to the smallest grid value if none satisfy the constraint
            # (e.g. degenerate H0); a small positive lambda is always safe-ish.
            self.lam = float(best) if best is not None else float(lam_grid[0])

        self.reset()
        return self

    def predict(self, x: float) -> bool:
        """Advance the CUSUM recursion by one sample and test for an alarm.

        Z(n) = max(0, Z(n-1) + z_lambda(x));  alarm iff Z(n) >= threshold.

        Args:
            x: New observation.

        Returns:
            True if the CUSUM statistic crosses the threshold (alarm).
        """
        self._t += 1
        z = float(self.instantaneous_statistic(np.array([x]))[0])
        self._g = max(0.0, self._g + z)
        if self._g >= self._threshold:
            if self._alarm_time is None:
                self._alarm_time = self._t
            return True
        return False

    def reset(self) -> None:
        """Reset the CUSUM statistic to zero (Z(0) = 0) and runtime state."""
        self._g = 0.0
        self._t = 0
        self._alarm_time = None

    def calibrate(
        self,
        h0_sampler: Callable[[int], np.ndarray],
        target_arl0: float,
        stream_len: int,
        n_runs: int = 2000,
        tol: float = 0.05,
        max_iter: int = 30,
    ) -> float:
        """Monte Carlo calibrate the threshold b to a target ARL_0 by bisection.

        Runs the detector on i.i.d. H0 streams and measures the empirical mean
        run length to (false) alarm; bisects b until the empirical ARL_0 matches
        the target within ``tol`` (relative). Runs that never alarm within
        ``stream_len`` are right-censored at ``stream_len``, so ``stream_len``
        should be a few times the target ARL_0.

        Args:
            h0_sampler: Function n -> array of n i.i.d. H0 samples.
            target_arl0: Desired in-control average run length.
            stream_len: Length of each H0 stream (right-censoring horizon).
            n_runs: Number of H0 streams.
            tol: Relative tolerance on the achieved ARL_0.
            max_iter: Maximum bisection iterations.

        Returns:
            The calibrated threshold b (also stored on the instance).
        """
        if self.lam is None:
            raise ValueError("Set lambda before calibrating (call fit()).")

        streams = [h0_sampler(stream_len) for _ in range(n_runs)]

        def empirical_arl0(b: float) -> float:
            run_lengths = []
            for stream in streams:
                # Vectorize the score difference, scan the recursion.
                z = self.lam * self.score_difference(stream)
                g, rl = 0.0, stream_len
                for t, zt in enumerate(z):
                    g = max(0.0, g + zt)
                    if g >= b:
                        rl = t + 1
                        break
                run_lengths.append(rl)
            return float(np.mean(run_lengths))

        # Bracket: ARL_0 increases monotonically with b.
        b_low, b_high = 0.0, max(np.log(target_arl0), 1.0)
        for _ in range(40):
            if empirical_arl0(b_high) >= target_arl0:
                break
            b_high *= 1.6

        best_b, best_gap = b_high, float("inf")
        for _ in range(max_iter):
            b_mid = 0.5 * (b_low + b_high)
            arl = empirical_arl0(b_mid)
            gap = abs(arl - target_arl0) / target_arl0
            if gap < best_gap:
                best_gap, best_b = gap, b_mid
            if arl < target_arl0:
                b_low = b_mid
            else:
                b_high = b_mid
            if gap < tol:
                break

        self._threshold = best_b
        self.reset()
        return best_b
