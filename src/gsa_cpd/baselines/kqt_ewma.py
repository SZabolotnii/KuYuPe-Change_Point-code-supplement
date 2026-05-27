"""KQT-EWMA baseline (univariate Quantile-Tree EWMA core).

Faithful implementation of the QT-EWMA detector that underlies the KQT-EWMA
method of

    Michelangelo Olmo Nogara Notarianni, Filippo Leveni, Diego Stucchi,
    Luca Frittoli, Giacomo Boracchi, "Change Detection in Multivariate Data
    Streams: Online Analysis with Kernel-QuantTree," AALTD @ ECML-PKDD 2024
    (Springer LNCS 15433; arXiv:2410.13778).

QT-EWMA itself is the distribution-free, online change detector of
Frittoli, Carrera and Boracchi (QuantTree-EWMA). It builds a histogram with K
bins of (target) probabilities pi_j from a stationary training set, then tracks
exponentially weighted moving averages of the per-bin membership indicators and
monitors their squared deviation from the target probabilities. Thresholds are
calibrated by Monte Carlo to operate at a pre-determined ARL_0, independently of
the data distribution.

Recursion (Eqs. 3-8 of the KQT-EWMA paper):

  Per-bin indicator for sample x_t falling in bin S_j:
      y_{j,t} = 1(x_t in S_j)

  Bias-corrected target bin probabilities (finite training set of size N):
      pi_hat_j = N*pi_j / (N + 1),         j < K
      pi_hat_K = (N*pi_K + 1) / (N + 1)
  (these are E[y_{j,t}] under H0 for a uniformly-at-random tie-break of the
  last training point; they also serve as the EWMA initial values).

  EWMA recursion:
      Z_{j,0} = pi_hat_j
      Z_{j,t} = (1 - lambda) * Z_{j,t-1} + lambda * y_{j,t}

  Monitoring statistic (Pearson-style normalized squared deviation):
      T_t = sum_j (Z_{j,t} - pi_hat_j)^2 / pi_hat_j

  Stopping rule with time-varying thresholds {h_t}:
      T = inf{ t >= 1 : T_t > h_t },
  where the h_t are chosen so that P(T_t > h_t | T_k <= h_k, k < t) = alpha for
  all t, giving ARL_0 = 1 / alpha (Eq. 8).

SIMPLIFICATION (documented): the full KQT method uses a Kernel-QuantTree
histogram for multivariate streams. Here the data are univariate, so the
QuantTree partition reduces to its 1-D special case: K equiprobable bins whose
edges are the empirical (j/K)-quantiles of the training set (target pi_j = 1/K).
This is the exact univariate QuantTree construction; only the multivariate
kernel partition is omitted, which is irrelevant for scalar series. The EWMA
recursion, statistic and ARL_0 calibration below follow the paper verbatim.

Threshold calibration: the paper precomputes the {h_t} sequence by Monte Carlo
on synthetic H0 data so that the per-step conditional false-alarm probability is
a constant alpha. We follow the same Monte Carlo route. For a compact, robust
implementation we calibrate a *single* constant threshold h (the stationary
regime of {h_t}, reached quickly because the EWMA statistic is stationary under
H0) by bisection to the target ARL_0. The per-step-constant-alpha variant is
provided via ``calibrate_sequential_thresholds`` for completeness.

API note: to match the other ``gsa_cpd`` baselines (``SignCUSUM``, ``EWMA``),
the online interface is ``fit(calibration_data) -> self``, ``predict(x) -> bool``,
``reset()`` and an ``alarm_time`` property. ``fit`` builds the QuantTree
partition from the calibration (H0) data; ``calibrate`` tunes the threshold.
"""

import numpy as np
from typing import Callable, Optional, List


class KQTEWMA:
    """Univariate QT-EWMA detector (Frittoli/Carrera/Boracchi; KQT-EWMA core).

    Example::

        det = KQTEWMA(n_bins=8, lam=0.05)
        det.fit(calibration_data)                 # builds the QuantTree bins
        det.calibrate(h0_sampler, target_arl0=500, stream_len=2500)
        for x in stream:
            if det.predict(x):
                print(f"Change at t={det.alarm_time}")
                break

    Args:
        n_bins: Number of equiprobable QuantTree bins K (target pi_j = 1/K).
        lam: EWMA forgetting factor lambda in (0, 1].
    """

    def __init__(self, n_bins: int = 8, lam: float = 0.03):
        if not (0.0 < lam <= 1.0):
            raise ValueError("lam must be in (0, 1].")
        self.n_bins = int(n_bins)
        self.lam = float(lam)
        self._threshold: float = np.inf

        # Set by fit(): bin edges and bias-corrected target probabilities.
        self.edges: Optional[np.ndarray] = None
        self.pi_hat: Optional[np.ndarray] = None

        # EWMA state Z_{j,t} and runtime state.
        self.Z: Optional[np.ndarray] = None
        self._t: int = 0
        self._alarm_time: Optional[int] = None

    @property
    def alarm_time(self) -> Optional[int]:
        """Time index of the alarm (None if no alarm)."""
        return self._alarm_time

    def fit(self, calibration_data: np.ndarray) -> "KQTEWMA":
        """Build the equiprobable 1-D QuantTree partition from a training set.

        Edges are the empirical (j/K) quantiles; target probabilities pi_j = 1/K
        are bias-corrected to pi_hat per Eq. (4) with N = len(calibration_data).

        Args:
            calibration_data: Stationary (H0) training samples.

        Returns:
            self (for chaining).
        """
        data = np.asarray(calibration_data, dtype=float)
        N = len(data)
        K = self.n_bins

        # Interior edges at the (1/K, 2/K, ..., (K-1)/K) empirical quantiles.
        qs = np.linspace(0.0, 1.0, K + 1)[1:-1]
        interior = np.quantile(data, qs)
        self.edges = np.concatenate(([-np.inf], interior, [np.inf]))

        # Target pi_j = 1/K, bias-corrected (Eq. 4): the last bin carries the
        # +1 correction for the tie-break of the held-out training point.
        pi = np.full(K, 1.0 / K)
        pi_hat = N * pi / (N + 1)
        pi_hat[-1] = (N * pi[-1] + 1.0) / (N + 1)
        self.pi_hat = pi_hat

        self.reset()
        return self

    def _bin_indicators(self, x: float) -> np.ndarray:
        """One-hot indicator vector y_{.,t} = 1(x in S_j) over the K bins."""
        # np.searchsorted on the interior edges maps x to a bin index in [0, K).
        idx = int(np.searchsorted(self.edges[1:-1], x, side="right"))
        idx = min(idx, self.n_bins - 1)
        y = np.zeros(self.n_bins)
        y[idx] = 1.0
        return y

    def statistic(self) -> float:
        """Current monitoring statistic T_t = sum_j (Z_j - pi_hat_j)^2 / pi_hat_j."""
        return float(np.sum((self.Z - self.pi_hat) ** 2 / self.pi_hat))

    def predict(self, x: float) -> bool:
        """Advance the EWMA recursion by one sample and test for an alarm.

        Z_{j,t} = (1 - lambda) Z_{j,t-1} + lambda * 1(x in S_j);
        alarm iff T_t > threshold.

        Args:
            x: New observation.

        Returns:
            True if the monitoring statistic crosses the threshold (alarm).
        """
        self._t += 1
        y = self._bin_indicators(float(x))
        self.Z = (1.0 - self.lam) * self.Z + self.lam * y
        if self.statistic() > self._threshold:
            if self._alarm_time is None:
                self._alarm_time = self._t
            return True
        return False

    def reset(self) -> None:
        """Reset the EWMA state to its H0 init Z_{j,0} = pi_hat_j and runtime state."""
        if self.pi_hat is None:
            raise ValueError("Call fit() before reset().")
        self.Z = self.pi_hat.copy()
        self._t = 0
        self._alarm_time = None

    def _run_length(self, stream: np.ndarray, b: float, horizon: int) -> int:
        """Run length to first crossing of constant threshold b on a stream."""
        Z = self.pi_hat.copy()
        edges_int = self.edges[1:-1]
        lam = self.lam
        pi_hat = self.pi_hat
        for t, x in enumerate(stream):
            idx = int(np.searchsorted(edges_int, x, side="right"))
            if idx >= self.n_bins:
                idx = self.n_bins - 1
            Z *= (1.0 - lam)
            Z[idx] += lam
            T = np.sum((Z - pi_hat) ** 2 / pi_hat)
            if T > b:
                return t + 1
        return horizon

    def calibrate(
        self,
        h0_sampler: Callable[[int], np.ndarray],
        target_arl0: float,
        stream_len: int,
        n_runs: int = 2000,
        tol: float = 0.05,
        max_iter: int = 30,
    ) -> float:
        """Monte Carlo calibrate a constant threshold h to a target ARL_0.

        Bisects the constant threshold until the empirical mean run length to a
        false alarm on i.i.d. H0 streams matches ``target_arl0`` within ``tol``.
        Runs that never alarm are right-censored at ``stream_len`` (choose
        ``stream_len`` a few times the target ARL_0). Calibrating a constant
        threshold targets the stationary regime of the paper's {h_t} sequence,
        which the stationary H0 EWMA statistic reaches quickly.

        Args:
            h0_sampler: Function n -> array of n i.i.d. H0 samples.
            target_arl0: Desired in-control average run length.
            stream_len: Length of each H0 stream (right-censoring horizon).
            n_runs: Number of H0 streams.
            tol: Relative tolerance on the achieved ARL_0.
            max_iter: Maximum bisection iterations.

        Returns:
            The calibrated constant threshold (also stored on the instance).
        """
        if self.pi_hat is None:
            raise ValueError("Call fit() before calibrate().")

        streams = [h0_sampler(stream_len) for _ in range(n_runs)]

        def empirical_arl0(b: float) -> float:
            rls = [self._run_length(s, b, stream_len) for s in streams]
            return float(np.mean(rls))

        # T_t is bounded; bracket b in [0, b_max]. ARL_0 increases with b.
        b_low, b_high = 0.0, max(1e-3, 4.0 / self.n_bins)
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

    def calibrate_sequential_thresholds(
        self,
        h0_sampler: Callable[[int], np.ndarray],
        target_arl0: float,
        horizon: int,
        n_runs: int = 5000,
    ) -> List[float]:
        """Calibrate the time-varying thresholds {h_t} (Eq. 8) by Monte Carlo.

        For each step t, h_t is the (1 - alpha) quantile of T_t over the H0 runs
        that have survived (not yet alarmed) up to t, with alpha = 1/ARL_0. This
        enforces P(T_t > h_t | T_k <= h_k, k < t) = alpha, the paper's defining
        property of QT-EWMA thresholds.

        Args:
            h0_sampler: Function n -> array of n i.i.d. H0 samples.
            target_arl0: Desired ARL_0; sets alpha = 1 / target_arl0.
            horizon: Number of steps t for which to compute thresholds.
            n_runs: Number of H0 runs (more is needed for tail quantiles).

        Returns:
            List of thresholds [h_1, ..., h_horizon].
        """
        if self.pi_hat is None:
            raise ValueError("Call fit() before calibrating thresholds.")

        alpha = 1.0 / target_arl0
        edges_int = self.edges[1:-1]
        lam = self.lam
        pi_hat = self.pi_hat

        # Precompute the T_t trajectories for all runs.
        Z = np.tile(pi_hat, (n_runs, 1))
        traj = np.empty((n_runs, horizon))
        streams = np.stack([h0_sampler(horizon) for _ in range(n_runs)])
        for t in range(horizon):
            xs = streams[:, t]
            idx = np.searchsorted(edges_int, xs, side="right")
            idx = np.minimum(idx, self.n_bins - 1)
            Z *= (1.0 - lam)
            Z[np.arange(n_runs), idx] += lam
            traj[:, t] = np.sum((Z - pi_hat) ** 2 / pi_hat, axis=1)

        thresholds: List[float] = []
        alive = np.ones(n_runs, dtype=bool)
        for t in range(horizon):
            Tt = traj[alive, t]
            if Tt.size == 0:
                thresholds.append(thresholds[-1] if thresholds else np.inf)
                continue
            h_t = float(np.quantile(Tt, 1.0 - alpha))
            thresholds.append(h_t)
            # Survivors are those still below their threshold at step t.
            alive_idx = np.where(alive)[0]
            crossed = traj[alive_idx, t] > h_t
            alive[alive_idx[crossed]] = False
        return thresholds
