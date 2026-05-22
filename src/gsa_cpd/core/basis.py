"""Basis function types and evaluation for GSA-LLR approximation."""

import numpy as np
from enum import Enum


class BasisType(Enum):
    """Available basis function families."""
    POLY = "poly"       # Polynomial: {x, x^2, ..., x^s}
    LOG = "log"         # Logarithmic: {x, ln|x|, x*ln|x|, (ln|x|)^2, ...}
    FRAC = "frac"       # Fractional: {sgn(x)|x|^a1, sgn(x)|x|^a2, ...}
    HERMITE = "hermite"  # Probabilist's Hermite polynomials


# Predefined fractional exponents (theoretically motivated)
FRAC_EXPONENTS = [1.0, 0.5, 1 / 3, 2 / 3, 0.25, 0.75]


def evaluate_basis_single(x: float, power: int, basis_type: BasisType) -> float:
    """Evaluate a single basis function phi_{power}(x).

    Args:
        x: Input value.
        power: Basis function index (1-based).
        basis_type: Type of basis.

    Returns:
        Basis function value (unclipped).
    """
    if basis_type == BasisType.POLY:
        return x ** power

    elif basis_type == BasisType.LOG:
        x_abs = abs(x) + 1e-8
        if power == 1:
            return x
        elif power == 2:
            return np.log(x_abs)
        elif power == 3:
            return x * np.log(x_abs)
        elif power == 4:
            return np.log(x_abs) ** 2
        else:
            return x * (np.log(x_abs) ** max(0, power - 3))

    elif basis_type == BasisType.FRAC:
        idx = min(power - 1, len(FRAC_EXPONENTS) - 1)
        alpha = FRAC_EXPONENTS[idx]
        return np.sign(x) * (abs(x) ** alpha)

    elif basis_type == BasisType.HERMITE:
        return _hermite_poly(x, power)

    return x ** power  # fallback


def evaluate_basis_vector(x: float, degree: int, basis_type: BasisType,
                          phi_max: float = 10.0) -> np.ndarray:
    """Evaluate all basis functions phi_1(x), ..., phi_s(x) with clipping.

    Args:
        x: Input value.
        degree: Number of basis functions (s).
        basis_type: Type of basis.
        phi_max: Clip values to [-phi_max, phi_max].

    Returns:
        Array of shape (degree,) with clipped basis function values.
    """
    values = np.array([evaluate_basis_single(x, i + 1, basis_type)
                       for i in range(degree)])
    return np.clip(values, -phi_max, phi_max)


def evaluate_basis_matrix(data: np.ndarray, degree: int, basis_type: BasisType,
                          phi_max: float = 10.0) -> np.ndarray:
    """Evaluate basis functions for an array of data points.

    Args:
        data: Array of shape (n,).
        degree: Number of basis functions.
        basis_type: Type of basis.
        phi_max: Clip values to [-phi_max, phi_max].

    Returns:
        Matrix of shape (n, degree).
    """
    n = len(data)
    B = np.zeros((n, degree))
    for i in range(degree):
        power = i + 1
        B[:, i] = np.array([evaluate_basis_single(x, power, basis_type)
                            for x in data])
    return np.clip(B, -phi_max, phi_max)


def _hermite_poly(x: float, n: int) -> float:
    """Probabilist's Hermite polynomial He_n(x).

    Orthogonal w.r.t. standard Gaussian measure N(0,1).
    He_0=1, He_1=x, He_2=x^2-1, He_3=x^3-3x, He_4=x^4-6x^2+3, ...
    """
    if n == 0:
        return 1.0
    elif n == 1:
        return x
    elif n == 2:
        return x ** 2 - 1
    elif n == 3:
        return x ** 3 - 3 * x
    elif n == 4:
        return x ** 4 - 6 * x ** 2 + 3
    elif n == 5:
        return x ** 5 - 10 * x ** 3 + 15 * x
    else:
        he_prev = x ** 4 - 6 * x ** 2 + 3
        he_curr = x ** 5 - 10 * x ** 3 + 15 * x
        for k in range(5, n):
            he_next = x * he_curr - k * he_prev
            he_prev = he_curr
            he_curr = he_next
        return he_curr
