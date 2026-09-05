"""Unified GSA Detector for sequential change-point detection."""

import warnings

import numpy as np
from typing import Optional

from gsa_cpd.core.basis import BasisType, evaluate_basis_vector, evaluate_basis_matrix
from gsa_cpd.core.solver import solve_fk_y
from gsa_cpd.core.threshold import (
    ThresholdType, compute_threshold, calibrate_threshold_mc,
)
from gsa_cpd.core.diagnostics import GSADiagnostics


class GSADetector:
    """Generalized Stochastic Approximation detector for sequential change-point detection.

    Approximates the log-likelihood ratio (LLR) using a finite set of basis
    functions and moment-based optimization (Kunchenko's KU1 criterion).
    Supports polynomial, logarithmic, fractional, and Hermite bases.

    Example::

        detector = GSADetector(basis=BasisType.FRAC, degree=2, epsilon=0.01)
        detector.fit(calibration_data, delta=0.3)
        for x in stream:
            if detector.predict(x):
                print(f"Change at t={detector.alarm_time}")
                break

    Args:
        basis: Basis function type.
        degree: Order of approximation s (number of basis functions).
        epsilon: Target false alarm rate.
        threshold_type: Method for threshold computation.
        threshold_scale: Safety multiplier for threshold (default 1.0).
        ridge_lambda: Tikhonov regularization for FK=Y solver.
        phi_max: Clip basis function values to [-phi_max, phi_max].
        winsor_pct: Winsorization percentage for each tail (0 to disable).
        standardize: Location-scale reduction applied to the observation before
            the dictionary is evaluated, estimated on the calibration sample
            alone. None (the default, and what the Monte-Carlo study of the
            paper uses) feeds the raw observation to the dictionary; "robust"
            uses the median and 1.4826 * MAD; "zscore" the mean and the standard
            deviation.

            **Use it on any series that is not already on a unit scale.** The
            dictionary is a set of powers and every value is clipped to
            [-phi_max, phi_max], so a series living outside that window has
            every basis value clipped to the same constant: the dictionary then
            carries no information, and the detector is inert at every threshold
            rather than merely conservative. `fit` warns when that happens and
            records it in `diagnostics.basis_degenerate`.
    """

    def __init__(
        self,
        basis: BasisType = BasisType.POLY,
        degree: int = 2,
        epsilon: float = 0.01,
        threshold_type: ThresholdType = ThresholdType.CHEBYSHEV,
        threshold_scale: float = 1.0,
        ridge_lambda: float = 1e-6,
        phi_max: float = 10.0,
        winsor_pct: float = 0.05,
        standardize: Optional[str] = None,
    ):
        self.basis = basis if isinstance(basis, BasisType) else BasisType(basis)
        self.degree = degree
        self.epsilon = epsilon
        self.threshold_type = (threshold_type if isinstance(threshold_type, ThresholdType)
                               else ThresholdType(threshold_type))
        self.threshold_scale = threshold_scale
        self.ridge_lambda = ridge_lambda
        self.phi_max = phi_max
        self.winsor_pct = winsor_pct
        if standardize not in (None, "robust", "zscore"):
            raise ValueError(
                f"standardize must be None, 'robust' or 'zscore', got {standardize!r}"
            )
        self.standardize = standardize

        # Location and scale of the standardising map, fitted in fit().
        # (0.0, 1.0) is the identity.
        self._loc: float = 0.0
        self._scale: float = 1.0
        self._basis_degenerate: bool = False

        # Fitted parameters
        self._coeffs: Optional[np.ndarray] = None
        self._k0: float = 0.0
        self._threshold: float = 0.0

        # Runtime state
        self._g_stat: float = 0.0
        self._t: int = 0
        self._alarm_time: Optional[int] = None

        # Diagnostics
        self._diagnostics = GSADiagnostics()
        self._fitted = False

    @property
    def diagnostics(self) -> GSADiagnostics:
        """Diagnostic information from calibration."""
        return self._diagnostics

    @property
    def coefficients(self) -> Optional[np.ndarray]:
        """Coefficient vector K from FK=Y."""
        return self._coeffs

    @property
    def alarm_time(self) -> Optional[int]:
        """Time index of the alarm (None if no alarm)."""
        return self._alarm_time

    def fit(self, calibration_data: np.ndarray, delta: float = 0.2,
            h1_data: Optional[np.ndarray] = None) -> "GSADetector":
        """Calibrate the detector on H0 data.

        Args:
            calibration_data: Array of observations under normal regime (H0).
            delta: MDE coefficient — expected relative change in moments.
            h1_data: Optional post-change (H1) sample. When given, the reference
                anomaly is taken from its empirical basis moments (a general
                change, e.g. a shape/skew/kurtosis change at matched mean and
                variance), overriding the MDE heuristic.

        Returns:
            self (for chaining).
        """
        self._fit_standardizer(calibration_data)
        data = self._standardize(self._winsorize(calibration_data))
        s = self.degree

        # 1. Compute basis function values
        B = evaluate_basis_matrix(data, s, self.basis, self.phi_max)
        self._check_basis_degenerate(B)

        # 2. Moments under H0
        u = np.mean(B, axis=0)
        Cov0 = np.cov(B, rowvar=False)
        if s == 1:
            Cov0 = np.array([[np.var(B)]])

        # 3. Reference hypothesis H1
        if h1_data is not None:
            # General reference anomaly: empirical H1 moments from actual
            # post-change data (e.g. a pure shape change at matched mean/var).
            Bh1 = evaluate_basis_matrix(
                self._standardize(self._winsorize(h1_data)),
                s, self.basis, self.phi_max)
            m = np.mean(Bh1, axis=0)
            Cov1 = np.cov(Bh1, rowvar=False) if s > 1 else np.array([[np.var(Bh1)]])
        else:
            # Reference hypothesis H1 (MDE strategy)
            m = u.copy()
            std_x = np.std(calibration_data)
            for i in range(s):
                power = i + 1
                if power % 2 == 0:
                    m[i] = u[i] * (1.0 + delta)
                else:
                    m[i] = u[i] + 0.1 * std_x
            Cov1 = Cov0 * (1.0 + delta)

        # 4. Build and solve FK=Y with the paper's normalisation
        # F = 0.5 * (C0 + C1) (eq:normal-system). Under it the surrogate
        # Lambda = K^T (phi - 0.5 * (u + m)) is asymptotically on the
        # log-likelihood-ratio scale: for local alternatives the exponential
        # tilt root of E_0 exp(theta Lambda) = 1 satisfies theta_0 -> 1.
        # Detection is unaffected (the threshold is calibrated to the realised
        # H0 moments of the statistic); only the absolute value of J(s) changes.
        F = 0.5 * (Cov0 + Cov1)
        Y = m - u

        K, cond_F, method = solve_fk_y(F, Y, ridge_lambda=self.ridge_lambda)
        self._coeffs = K

        # 5. Bias term
        self._k0 = -0.5 * float(np.dot(K, m + u))

        # 6. Statistics of Lambda under H0 and H1
        E_L_H0 = self._k0 + float(np.dot(K, u))
        E_L_H1 = self._k0 + float(np.dot(K, m))
        Var_L_H0 = max(float(K @ Cov0 @ K), 1e-12)

        # 7. Threshold
        if self.threshold_type == ThresholdType.SIMULATION:
            self._threshold = calibrate_threshold_mc(
                lambda x: self._compute_llr(x),
                calibration_data, self.epsilon,
            )
        else:
            self._threshold = compute_threshold(
                E_L_H0, Var_L_H0, self.epsilon,
                self.threshold_type, self.threshold_scale,
            )

        # 8. Diagnostics
        sigma_0 = np.sqrt(Var_L_H0)
        eta = (E_L_H1 - E_L_H0) / sigma_0 if sigma_0 > 0 else 0.0

        self._diagnostics = GSADiagnostics(
            condition_number=cond_F,
            J_s=float(np.dot(K, Y)),
            J_s_lean=float(np.dot(K, K)),
            E_L_H0=E_L_H0,
            E_L_H1=E_L_H1,
            Var_L_H0=Var_L_H0,
            threshold=self._threshold,
            coeffs=K.copy(),
            k0=self._k0,
            eta=eta,
            solver_method=method,
            basis_degenerate=self._basis_degenerate,
            standardize_loc=self._loc,
            standardize_scale=self._scale,
        )

        self._fitted = True
        self.reset()
        return self

    def predict(self, x: float) -> bool:
        """Process one observation and return True if change detected.

        Args:
            x: New observation.

        Returns:
            True if the CUSUM statistic exceeds the threshold (alarm).
        """
        if not self._fitted:
            raise RuntimeError("Detector not fitted. Call fit() first.")

        self._t += 1
        llr = self._compute_llr(x)
        self._g_stat = max(0.0, self._g_stat + llr)

        if self._g_stat > self._threshold:
            if self._alarm_time is None:
                self._alarm_time = self._t
            return True
        return False

    def reset(self) -> None:
        """Reset the runtime state for a new monitoring session."""
        self._g_stat = 0.0
        self._t = 0
        self._alarm_time = None

    def _compute_llr(self, x: float) -> float:
        """Compute approximate LLR increment Lambda(x)."""
        z = (x - self._loc) / self._scale
        v = evaluate_basis_vector(z, self.degree, self.basis, self.phi_max)
        return self._k0 + float(np.dot(self._coeffs, v))

    def _standardize(self, data: np.ndarray) -> np.ndarray:
        """Apply the fitted location-scale map. Identity when standardize=None."""
        return (data - self._loc) / self._scale

    def _fit_standardizer(self, calibration_data: np.ndarray) -> None:
        """Fit the location-scale map from the calibration sample alone.

        A non-positive or non-finite scale estimate falls back to the next
        cruder one and finally to the identity: a constant calibration sample
        carries no scale, and inventing one would be worse than leaving the
        observation alone.
        """
        if self.standardize is None:
            self._loc, self._scale = 0.0, 1.0
        elif self.standardize == "zscore":
            self._loc = float(np.mean(calibration_data))
            scale = float(np.std(calibration_data))
            self._scale = scale if np.isfinite(scale) and scale > 0 else 1.0
        else:  # "robust"
            self._loc = float(np.median(calibration_data))
            scale = 1.4826 * float(np.median(np.abs(calibration_data - self._loc)))
            if not (np.isfinite(scale) and scale > 0):
                scale = float(np.std(calibration_data))
            self._scale = scale if np.isfinite(scale) and scale > 0 else 1.0

        self._basis_degenerate = False

    def _check_basis_degenerate(self, B: np.ndarray) -> None:
        """Warn when the clip at +-phi_max has flattened a basis column.

        A constant column carries no information at any threshold, so the
        detector cannot alarm however the threshold is set.  Left silent this is
        indistinguishable from a confident non-detection, which is exactly how
        it goes unnoticed.
        """
        degenerate = bool(np.any(np.ptp(B, axis=0) < 1e-12))
        self._basis_degenerate = degenerate
        if degenerate:
            warnings.warn(
                "The dictionary is constant on the calibration sample: every "
                f"basis value was clipped to +-phi_max={self.phi_max}. The "
                "detector cannot raise an alarm at any threshold. Pass "
                'standardize="robust" (or rescale the series yourself) so the '
                "dictionary sees an observation on a unit scale.",
                RuntimeWarning,
                stacklevel=3,
            )

    def _winsorize(self, data: np.ndarray) -> np.ndarray:
        """Apply Winsorization to calibration data."""
        if self.winsor_pct <= 0:
            return data.copy()
        lower = np.percentile(data, self.winsor_pct * 100)
        upper = np.percentile(data, (1 - self.winsor_pct) * 100)
        return np.clip(data, lower, upper)
