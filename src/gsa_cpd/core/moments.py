"""Moment estimation for GSA detector calibration.

Provides both empirical (from data) and theoretical (analytical) moment
computation for basis functions under H0 and H1.
"""

import numpy as np
from typing import Tuple

from gsa_cpd.core.basis import BasisType, evaluate_basis_matrix


def estimate_empirical_moments(
    data: np.ndarray,
    degree: int,
    basis_type: BasisType,
    phi_max: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Estimate mean vector and covariance matrix of basis functions.

    Args:
        data: Calibration data under H0, shape (n,).
        degree: Number of basis functions.
        basis_type: Type of basis.
        phi_max: Clipping bound for basis functions.

    Returns:
        Tuple of (u, Cov) where:
            u: Mean vector, shape (degree,).
            Cov: Covariance matrix, shape (degree, degree).
    """
    B = evaluate_basis_matrix(data, degree, basis_type, phi_max)
    u = np.mean(B, axis=0)
    if degree == 1:
        Cov = np.array([[np.var(B)]])
    else:
        Cov = np.cov(B, rowvar=False)
    return u, Cov


def compute_mde_hypothesis(
    u: np.ndarray,
    Cov0: np.ndarray,
    delta: float,
    std_x: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute H1 moments using Minimal Detectable Effect (MDE) strategy.

    Args:
        u: Mean vector under H0.
        Cov0: Covariance matrix under H0.
        delta: Expected relative change in moments.
        std_x: Standard deviation of raw data.

    Returns:
        Tuple of (m, Cov1) where:
            m: Expected mean vector under H1.
            Cov1: Expected covariance under H1.
    """
    s = len(u)
    m = u.copy()
    for i in range(s):
        power = i + 1
        if power % 2 == 0:
            m[i] = u[i] * (1.0 + delta)
        else:
            m[i] = u[i] + 0.1 * std_x
    Cov1 = Cov0 * (1.0 + delta)
    return m, Cov1


def winsorize(data: np.ndarray, pct: float = 0.05) -> np.ndarray:
    """Apply Winsorization to data.

    Clips extreme values at the given percentile on each tail.

    Args:
        data: Input array.
        pct: Fraction to clip on each tail (e.g., 0.05 = 5%).

    Returns:
        Winsorized copy of data.
    """
    if pct <= 0:
        return data.copy()
    lower = np.percentile(data, pct * 100)
    upper = np.percentile(data, (1 - pct) * 100)
    return np.clip(data, lower, upper)
