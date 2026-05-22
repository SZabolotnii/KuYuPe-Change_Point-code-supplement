"""CUSUM stopping rule (Lorden minimax)."""


class CUSUMRule:
    """Cumulative Sum (CUSUM) stopping rule.

    g_t = max(0, g_{t-1} + lambda_t)
    Alarm when g_t > threshold.

    Minimax-optimal: minimizes worst-case ADD over all change times.
    """

    def __init__(self, threshold: float):
        self.threshold = threshold
        self._g: float = 0.0

    def update(self, lambda_t: float) -> bool:
        """Process one LLR increment. Returns True if alarm."""
        self._g = max(0.0, self._g + lambda_t)
        return self._g > self.threshold

    def reset(self) -> None:
        """Reset statistic to zero."""
        self._g = 0.0

    @property
    def statistic(self) -> float:
        """Current value of the CUSUM statistic."""
        return self._g
