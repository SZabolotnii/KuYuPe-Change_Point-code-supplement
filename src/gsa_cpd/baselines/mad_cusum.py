"""MAD-normalized robust CUSUM baseline detector."""

import numpy as np
from typing import Optional


class MADCUSUM:
    """MAD-CUSUM detector.

    Normalizes observations using Median Absolute Deviation (MAD),
    then applies CUSUM on the squared z-scores.
    Robust to outliers in the calibration data.
    """

    def __init__(self, epsilon: float = 0.01):
        self.epsilon = epsilon
        self._median: float = 0.0
        self._mad: float = 1.0
        self._threshold: float = 0.0
        self._g: float = 0.0
        self._alarm_time: Optional[int] = None
        self._t: int = 0

    def fit(self, calibration_data: np.ndarray) -> "MADCUSUM":
        self._median = float(np.median(calibration_data))
        self._mad = float(np.median(np.abs(calibration_data - self._median)))
        if self._mad < 1e-10:
            self._mad = float(np.std(calibration_data)) or 1.0
        self._threshold = 1.0 / self.epsilon
        self.reset()
        return self

    def predict(self, x: float) -> bool:
        self._t += 1
        z = (x - self._median) / (1.4826 * self._mad)  # scale to ~N(0,1)
        llr = z ** 2 - 1.0  # centered chi-squared
        self._g = max(0.0, self._g + llr)
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
