"""Preprocessing utilities for time series data."""

import numpy as np
from typing import Optional


def winsorize(data: np.ndarray, pct: float = 0.05) -> np.ndarray:
    """Winsorize data at given percentile."""
    if pct <= 0:
        return data.copy()
    lower = np.percentile(data, pct * 100)
    upper = np.percentile(data, (1 - pct) * 100)
    return np.clip(data, lower, upper)


def zscore(data: np.ndarray) -> np.ndarray:
    """Standardize to zero mean and unit variance."""
    mu = np.mean(data)
    sigma = np.std(data)
    if sigma < 1e-10:
        return data - mu
    return (data - mu) / sigma


def log_transform(data: np.ndarray, offset: float = 1.0) -> np.ndarray:
    """Apply log(1 + |x|) * sign(x) transform."""
    return np.sign(data) * np.log1p(np.abs(data) + offset - 1.0)


def difference(data: np.ndarray, order: int = 1) -> np.ndarray:
    """Compute n-th order differences."""
    result = data.copy()
    for _ in range(order):
        result = np.diff(result)
    return result


def ar1_decorrelate(data: np.ndarray) -> np.ndarray:
    """Remove AR(1) autocorrelation via residuals."""
    if len(data) < 3:
        return data.copy()
    phi = np.corrcoef(data[:-1], data[1:])[0, 1]
    residuals = data[1:] - phi * data[:-1]
    return residuals
