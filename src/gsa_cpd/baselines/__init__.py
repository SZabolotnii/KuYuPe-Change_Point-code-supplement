"""Baseline detectors for comparison with GSA."""

from gsa_cpd.baselines.sign_cusum import SignCUSUM
from gsa_cpd.baselines.mad_cusum import MADCUSUM
from gsa_cpd.baselines.ewma import EWMA
from gsa_cpd.baselines.oracle_cusum import OracleCUSUM
from gsa_cpd.baselines.scusum import SCUSUM, gaussian_hyvarinen_score
from gsa_cpd.baselines.kqt_ewma import KQTEWMA

__all__ = [
    "SignCUSUM",
    "MADCUSUM",
    "EWMA",
    "OracleCUSUM",
    "SCUSUM",
    "gaussian_hyvarinen_score",
    "KQTEWMA",
]
