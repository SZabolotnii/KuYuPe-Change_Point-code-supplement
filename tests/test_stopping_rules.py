"""Tests for CUSUM, GRSh, and SRP stopping rules."""

import numpy as np
import pytest
from gsa_cpd.stopping_rules import CUSUMRule, GRShRule, SRPRule


class TestCUSUM:
    def test_no_alarm_on_negative_input(self):
        rule = CUSUMRule(threshold=10.0)
        for _ in range(100):
            assert not rule.update(-1.0)
        assert rule.statistic == 0.0  # reset to 0

    def test_alarm_on_positive_input(self):
        rule = CUSUMRule(threshold=5.0)
        for i in range(10):
            if rule.update(1.0):
                assert rule.statistic > 5.0
                return
        pytest.fail("Should have triggered alarm")

    def test_reset(self):
        rule = CUSUMRule(threshold=100.0)
        rule.update(50.0)
        assert rule.statistic > 0
        rule.reset()
        assert rule.statistic == 0.0


class TestGRSh:
    def test_accumulates_without_reset(self):
        rule = GRShRule(threshold=100.0)
        for _ in range(5):
            rule.update(1.0)
        assert rule.statistic == pytest.approx(5.0)
        # Feed negative — should NOT reset
        rule.update(-2.0)
        assert rule.statistic == pytest.approx(3.0)

    def test_alarm(self):
        rule = GRShRule(threshold=5.0)
        alarmed = False
        for _ in range(10):
            if rule.update(1.0):
                alarmed = True
                break
        assert alarmed


class TestSRP:
    def test_multiplicative_update(self):
        rule = SRPRule(threshold=1e10)
        rule.update(0.0)  # R = (1+0)*exp(0) = 1
        assert rule.statistic == pytest.approx(1.0)
        rule.update(0.0)  # R = (1+1)*exp(0) = 2
        assert rule.statistic == pytest.approx(2.0)

    def test_alarm_on_large_input(self):
        rule = SRPRule(threshold=100.0)
        assert rule.update(10.0)  # (1+0)*exp(10) >> 100
