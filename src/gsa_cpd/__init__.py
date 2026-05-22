"""
GSA Change-Point Detection.

Generalized Stochastic Approximation of Log-Likelihood Ratio
for Robust Sequential Change-Point Detection.
"""

from gsa_cpd.core.detector import GSADetector
from gsa_cpd.core.basis import BasisType
from gsa_cpd.core.threshold import ThresholdType
from gsa_cpd.core.diagnostics import GSADiagnostics

__version__ = "1.0.0"

__all__ = [
    "GSADetector",
    "BasisType",
    "ThresholdType",
    "GSADiagnostics",
]
