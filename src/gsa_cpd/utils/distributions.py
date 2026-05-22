"""Distribution factory for Monte Carlo experiments."""

import numpy as np
import scipy.stats as stats
from typing import Dict, Any, Callable, Optional


SUPPORTED_DISTRIBUTIONS = [
    "normal", "pearson3", "student_t", "laplace",
    "pareto", "lognormal", "exponential",
]


def create_distribution(
    name: str,
    params: Optional[Dict[str, Any]] = None,
) -> stats.rv_continuous:
    """Create a scipy distribution object.

    Args:
        name: Distribution name (normal, pearson3, student_t, etc.).
        params: Distribution parameters.

    Returns:
        Frozen scipy distribution.
    """
    params = params or {}

    if name == "normal":
        return stats.norm(
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "pearson3":
        return stats.pearson3(
            skew=params.get("skew", 0.0),
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "student_t":
        return stats.t(
            df=params.get("df", 5),
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "laplace":
        return stats.laplace(
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "pareto":
        return stats.pareto(
            b=params.get("b", 2.5),
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "lognormal":
        return stats.lognorm(
            s=params.get("s", 1.0),
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    elif name == "exponential":
        return stats.expon(
            loc=params.get("loc", 0.0),
            scale=params.get("scale", 1.0),
        )
    else:
        raise ValueError(f"Unknown distribution: {name}. "
                         f"Supported: {SUPPORTED_DISTRIBUTIONS}")


def compute_true_llr(
    x: float,
    dist_h0: stats.rv_continuous,
    dist_h1: stats.rv_continuous,
    eps: float = 1e-10,
) -> float:
    """Compute true log-likelihood ratio ln(f1(x)/f0(x)).

    Args:
        x: Observation value.
        dist_h0: Distribution under H0.
        dist_h1: Distribution under H1.
        eps: Regularization to avoid log(0).

    Returns:
        LLR value.
    """
    f0 = max(dist_h0.pdf(x), eps)
    f1 = max(dist_h1.pdf(x), eps)
    return float(np.log(f1 / f0))
