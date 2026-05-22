"""Shiryaev-Roberts Procedure (SRP) stopping rule."""

import numpy as np
from typing import Optional, Callable


class SRPRule:
    """Shiryaev-Roberts stopping rule.

    R_t = (1 + R_{t-1}) * exp(lambda_t)
    Alarm when R_t > H.

    Quasi-minimax optimal. Requires separate threshold calibration
    via ARL matching (see calibrate_arl).
    """

    def __init__(self, threshold: float):
        self.threshold = threshold
        self._R: float = 0.0

    def update(self, lambda_t: float) -> bool:
        """Process one LLR increment. Returns True if alarm."""
        self._R = (1.0 + self._R) * np.exp(lambda_t)
        return self._R > self.threshold

    def reset(self) -> None:
        """Reset statistic to zero."""
        self._R = 0.0

    @property
    def statistic(self) -> float:
        """Current value of the Shiryaev-Roberts statistic."""
        return self._R

    @staticmethod
    def calibrate_arl(
        compute_llr: Callable[[float], float],
        sample_h0: Callable[[int], np.ndarray],
        target_arl: float,
        max_iter: int = 12,
        tol: float = 0.05,
        n_runs: int = 200,
        max_run_length: int = 5000,
        rng: Optional[np.random.Generator] = None,
    ) -> float:
        """Calibrate SRP threshold via binary search on ARL_0.

        Args:
            compute_llr: Function that computes Lambda(x) for one observation.
            sample_h0: Function that generates n samples under H0.
            target_arl: Target Average Run Length under H0.
            max_iter: Maximum binary search iterations.
            tol: Relative tolerance for ARL match.
            n_runs: Number of MC runs to estimate ARL.
            max_run_length: Max observations per run.
            rng: Random generator.

        Returns:
            Calibrated threshold H.
        """
        if rng is None:
            rng = np.random.default_rng()

        H_low = 1.0
        H_high = np.exp(10.0)

        def estimate_arl(H: float) -> float:
            run_lengths = []
            for _ in range(n_runs):
                data = sample_h0(max_run_length)
                R = 0.0
                for t, x in enumerate(data):
                    R = (1.0 + R) * np.exp(compute_llr(x))
                    if R > H:
                        run_lengths.append(t + 1)
                        break
                else:
                    run_lengths.append(max_run_length)
            return float(np.mean(run_lengths))

        for _ in range(max_iter):
            H_mid = (H_low + H_high) / 2
            arl = estimate_arl(H_mid)
            if abs(arl - target_arl) / target_arl < tol:
                return H_mid
            if arl < target_arl:
                H_low = H_mid
            else:
                H_high = H_mid

        return (H_low + H_high) / 2
