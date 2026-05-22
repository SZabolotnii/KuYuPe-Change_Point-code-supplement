"""Exponentially Weighted Moving Average (EWMA) baseline detector."""

import numpy as np
from typing import Optional


class EWMA:
    """EWMA control chart detector.

    Tracks exponentially weighted moving average and signals when
    it deviates from the calibration mean by more than L*sigma.
    """

    def __init__(self, lam: float = 0.1, L: float = 3.5):
        self.lam = lam
        self.L = L
        self._mu: float = 0.0
        self._sigma: float = 1.0
        self._z: float = 0.0
        self._alarm_time: Optional[int] = None
        self._t: int = 0

    def fit(self, calibration_data: np.ndarray) -> "EWMA":
        self._mu = float(np.mean(calibration_data))
        self._sigma = float(np.std(calibration_data))
        if self._sigma < 1e-10:
            self._sigma = 1.0
        self.reset()
        return self

    def predict(self, x: float) -> bool:
        self._t += 1
        self._z = self.lam * x + (1 - self.lam) * self._z
        # EWMA control limit (asymptotic)
        sigma_z = self._sigma * np.sqrt(self.lam / (2 - self.lam))
        if abs(self._z - self._mu) > self.L * sigma_z:
            if self._alarm_time is None:
                self._alarm_time = self._t
            return True
        return False

    def reset(self) -> None:
        self._z = self._mu
        self._t = 0
        self._alarm_time = None

    @property
    def alarm_time(self) -> Optional[int]:
        return self._alarm_time
