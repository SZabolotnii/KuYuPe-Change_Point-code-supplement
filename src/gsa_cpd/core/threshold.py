"""Threshold computation for GSA detector (PE, VP, Cantelli, simulation)."""

import numpy as np
from enum import Enum
from typing import Optional


class ThresholdType(Enum):
    """Available threshold computation methods."""
    CHEBYSHEV = "chebyshev"   # PE: h = E + sqrt(Var/eps) — universal
    VP = "vp"                 # Vysochansky-Petunin: h = E + (2/3)*sqrt(Var/eps) — unimodal
    CANTELLI = "cantelli"     # One-sided: h = E + sigma*sqrt(1/eps - 1)
    SIMULATION = "simulation"  # MC-calibrated via binary search on ARL


def compute_threshold(
    E_L_H0: float,
    Var_L_H0: float,
    epsilon: float,
    threshold_type: ThresholdType = ThresholdType.CHEBYSHEV,
    threshold_scale: float = 1.0,
) -> float:
    """Compute detection threshold analytically.

    Args:
        E_L_H0: Expected value of Lambda under H0.
        Var_L_H0: Variance of Lambda under H0.
        epsilon: Target FAR level.
        threshold_type: Threshold computation method.
        threshold_scale: Safety multiplier (default 1.0).

    Returns:
        Threshold value h.
    """
    sigma_0 = np.sqrt(max(Var_L_H0, 1e-12))

    if threshold_type == ThresholdType.CHEBYSHEV:
        h = E_L_H0 + sigma_0 / np.sqrt(epsilon)

    elif threshold_type == ThresholdType.VP:
        h = E_L_H0 + (2.0 / 3.0) * sigma_0 / np.sqrt(epsilon)

    elif threshold_type == ThresholdType.CANTELLI:
        h = E_L_H0 + sigma_0 * np.sqrt(1.0 / epsilon - 1.0)

    else:
        raise ValueError(f"Use calibrate_threshold_mc for {threshold_type}")

    return h * threshold_scale


def calibrate_threshold_mc(
    compute_llr_func,
    calibration_data: np.ndarray,
    epsilon: float,
    n_runs: int = 50,
    block_length: int = 500,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """Calibrate threshold via Monte Carlo simulation.

    Estimates the (1-epsilon)-quantile of max CUSUM under H0.

    Args:
        compute_llr_func: Callable(x) -> float, computes Lambda(x).
        calibration_data: H0 calibration data.
        epsilon: Target FAR level.
        n_runs: Number of bootstrap runs.
        block_length: Length of each simulation block.
        rng: Random number generator.

    Returns:
        Calibrated threshold.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(calibration_data)
    max_cusums = []

    for _ in range(n_runs):
        indices = rng.choice(n, size=block_length, replace=True)
        block = calibration_data[indices]

        g = 0.0
        g_max = 0.0
        for x in block:
            llr = compute_llr_func(x)
            g = max(0.0, g + llr)
            g_max = max(g_max, g)
        max_cusums.append(g_max)

    return float(np.quantile(max_cusums, 1.0 - epsilon))
