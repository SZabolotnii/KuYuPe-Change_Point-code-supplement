"""GRSh-like Bayesian stopping rule (Girshick-Rubin-Shiryaev)."""


class GRShRule:
    """Additive Bayesian stopping rule (GRSh-like).

    S_t = S_{t-1} + lambda_t  (no reset)
    Alarm when S_t > threshold.

    Bayesian-optimal for geometric prior on change time.
    More sensitive to drift, but less robust to past false spikes.
    """

    def __init__(self, threshold: float):
        self.threshold = threshold
        self._S: float = 0.0

    def update(self, lambda_t: float) -> bool:
        """Process one LLR increment. Returns True if alarm."""
        self._S += lambda_t
        return self._S > self.threshold

    def reset(self) -> None:
        """Reset statistic to zero."""
        self._S = 0.0

    @property
    def statistic(self) -> float:
        """Current value of the accumulated sum."""
        return self._S
