"""Modular stopping rules for sequential change-point detection."""

from gsa_cpd.stopping_rules.cusum import CUSUMRule
from gsa_cpd.stopping_rules.grsh import GRShRule
from gsa_cpd.stopping_rules.srp import SRPRule

__all__ = ["CUSUMRule", "GRShRule", "SRPRule"]
