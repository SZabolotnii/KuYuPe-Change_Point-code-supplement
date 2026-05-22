"""Numerical solver for the FK=Y linear system with stability control."""

import numpy as np
from typing import Tuple


def solve_fk_y(
    F: np.ndarray,
    Y: np.ndarray,
    ridge_lambda: float = 1e-6,
    cond_threshold_ridge: float = 1e6,
    cond_threshold_svd: float = 1e8,
    svd_cutoff: float = 1e-10,
) -> Tuple[np.ndarray, float, str]:
    """Solve FK=Y with 3-level numerical stability strategy.

    Level 1 (cond < cond_threshold_ridge): Direct LU solve.
    Level 2 (cond < cond_threshold_svd): Ridge regularization F + lambda*I.
    Level 3 (cond >= cond_threshold_svd): Truncated SVD pseudo-inverse.

    Args:
        F: System matrix (s x s), symmetric positive semi-definite.
        Y: Right-hand side vector (s,).
        ridge_lambda: Tikhonov regularization parameter.
        cond_threshold_ridge: Condition number threshold for ridge.
        cond_threshold_svd: Condition number threshold for SVD.
        svd_cutoff: Singular value cutoff for SVD pseudo-inverse.

    Returns:
        Tuple of (K, cond_F, method) where:
            K: Solution vector (s,).
            cond_F: Condition number of F.
            method: "direct", "ridge", or "svd".
    """
    cond_F = float(np.linalg.cond(F))

    # Level 1: Direct solve
    if cond_F < cond_threshold_ridge:
        try:
            K = np.linalg.solve(F, Y)
            return K, cond_F, "direct"
        except np.linalg.LinAlgError:
            pass  # fall through to ridge

    # Level 2: Ridge regularization
    if cond_F < cond_threshold_svd:
        try:
            F_reg = F + ridge_lambda * np.eye(F.shape[0])
            K = np.linalg.solve(F_reg, Y)
            return K, cond_F, "ridge"
        except np.linalg.LinAlgError:
            pass  # fall through to SVD

    # Level 3: Truncated SVD
    U, sigma, Vt = np.linalg.svd(F, full_matrices=False)
    sigma_inv = np.where(sigma > svd_cutoff, 1.0 / sigma, 0.0)
    K = Vt.T @ np.diag(sigma_inv) @ U.T @ Y

    return K, cond_F, "svd"
