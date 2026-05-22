"""Baseline detectors for comparison with GSA."""

from gsa_cpd.baselines.sign_cusum import SignCUSUM
from gsa_cpd.baselines.mad_cusum import MADCUSUM
from gsa_cpd.baselines.ewma import EWMA
from gsa_cpd.baselines.oracle_cusum import OracleCUSUM

__all__ = ["SignCUSUM", "MADCUSUM", "EWMA", "OracleCUSUM"]
