"""Oracle CUSUM baseline (requires known true LLR)."""

import numpy as np
from typing import Optional, Callable


class OracleCUSUM:
    """Oracle CUSUM detector using true analytical LLR.

    Used as an upper-bound baseline: shows the best possible CUSUM
    performance when both f0 and f1 are fully known.
    """

    def __init__(self, true_llr_func: Callable[[float], float], epsilon: float = 0.01):
        self.true_llr_func = true_llr_func
        self.epsilon = epsilon
        self._threshold: float = 0.0
        self._E_L_H0: float = 0.0
        self._g: float = 0.0
        self._alarm_time: Optional[int] = None
        self._t: int = 0

    def fit(self, calibration_data: np.ndarray) -> "OracleCUSUM":
        llr_values = np.array([self.true_llr_func(x) for x in calibration_data])
        self._E_L_H0 = float(np.mean(llr_values))
        Var_L_H0 = float(np.var(llr_values))
        sigma = np.sqrt(max(Var_L_H0, 1e-12))
        self._threshold = self._E_L_H0 + sigma / np.sqrt(self.epsilon)
        self.reset()
        return self

    def predict(self, x: float) -> bool:
        self._t += 1
        llr = self.true_llr_func(x)
        self._g = max(0.0, self._g + llr - self._E_L_H0)
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
