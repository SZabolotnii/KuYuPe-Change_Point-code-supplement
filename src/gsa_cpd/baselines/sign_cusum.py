"""Nonparametric Sign-CUSUM baseline detector."""

import numpy as np
from typing import Optional


class SignCUSUM:
    """Sign-based CUSUM detector.

    Uses the sign of (x - median) as the test statistic.
    Fully nonparametric: no distributional assumptions.
    """

    def __init__(self, epsilon: float = 0.01):
        self.epsilon = epsilon
        self._median: float = 0.0
        self._threshold: float = 0.0
        self._g: float = 0.0
        self._alarm_time: Optional[int] = None
        self._t: int = 0

    def fit(self, calibration_data: np.ndarray) -> "SignCUSUM":
        self._median = float(np.median(calibration_data))
        self._threshold = np.sqrt(1.0 / self.epsilon)
        self.reset()
        return self

    def predict(self, x: float) -> bool:
        self._t += 1
        sign_val = 1.0 if x > self._median else -1.0
        self._g = max(0.0, self._g + sign_val)
        if self._g > self._threshold:
            if self._alarm_time is None:
                self._alarm_time = self._t
            return True
        return False

    def reset(self) -> None:
        self._g = 0.0
        self._t = 0
        self._alarm_time = None

    @property
    def alarm_time(self) -> Optional[int]:
        return self._alarm_time
