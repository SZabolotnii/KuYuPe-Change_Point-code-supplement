"""Diagnostic information from GSA detector calibration."""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class GSADiagnostics:
    """Diagnostic information collected during GSA detector calibration.

    Attributes:
        condition_number: Condition number of system matrix F.
        J_s: Information functional J(s) = K^T Y (Kunchenko definition).
        J_s_lean: J(s) = sum(K_i^2) (Lean/Parseval definition).
        E_L_H0: Expected value of Lambda under H0.
        E_L_H1: Expected value of Lambda under H1.
        Var_L_H0: Variance of Lambda under H0.
        threshold: Detection threshold h.
        coeffs: Coefficient vector K.
        k0: Bias term.
        eta: Efficiency coefficient (E[L|H1] - E[L|H0]) / sqrt(Var[L|H0]).
        solver_method: Method used to solve FK=Y ("direct", "ridge", "svd").
        basis_degenerate: True when the clip at +-phi_max leaves a basis column
            constant on the calibration sample. That column carries no
            information at any threshold, so the detector is inert rather than
            conservative. Almost always means the series was not standardised.
        standardize_loc: Location of the fitted standardising map (0.0 = none).
        standardize_scale: Scale of the fitted standardising map (1.0 = none).
    """
    condition_number: float = 0.0
    J_s: float = 0.0
    J_s_lean: float = 0.0
    E_L_H0: float = 0.0
    E_L_H1: float = 0.0
    Var_L_H0: float = 0.0
    threshold: float = 0.0
    coeffs: Optional[np.ndarray] = None
    k0: float = 0.0
    eta: float = 0.0
    solver_method: str = ""
    basis_degenerate: bool = False
    standardize_loc: float = 0.0
    standardize_scale: float = 1.0
