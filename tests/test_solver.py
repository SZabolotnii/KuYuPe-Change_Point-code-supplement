"""Tests for FK=Y solver with stability control."""

import numpy as np
import pytest
from gsa_cpd.core.solver import solve_fk_y


class TestDirectSolve:
    def test_well_conditioned_system(self):
        F = np.array([[2.0, 0.5], [0.5, 3.0]])
        Y = np.array([1.0, 2.0])
        K, cond, method = solve_fk_y(F, Y)
        assert method == "direct"
        np.testing.assert_allclose(F @ K, Y, atol=1e-10)

    def test_identity_system(self):
        F = np.eye(3)
        Y = np.array([1.0, 2.0, 3.0])
        K, cond, method = solve_fk_y(F, Y)
        np.testing.assert_allclose(K, Y, atol=1e-10)


class TestRidgeFallback:
    def test_moderately_ill_conditioned(self):
        # Create a system with cond ~ 1e7
        F = np.diag([1.0, 1e-7])
        Y = np.array([1.0, 1.0])
        K, cond, method = solve_fk_y(F, Y, cond_threshold_ridge=1e6)
        assert method in ("ridge", "svd")
        assert K is not None


class TestSVDFallback:
    def test_singular_system(self):
        F = np.array([[1.0, 1.0], [1.0, 1.0]])  # rank 1
        Y = np.array([1.0, 1.0])
        K, cond, method = solve_fk_y(F, Y)
        assert method == "svd"
        assert not np.any(np.isnan(K))

    def test_near_singular_preserves_direction(self):
        F = np.diag([1.0, 1e-15])
        Y = np.array([2.0, 0.0])
        K, _, method = solve_fk_y(F, Y)
        assert method == "svd"
        assert abs(K[0] - 2.0) < 0.1  # first component should be ~2


class TestConditionNumber:
    def test_condition_number_returned(self):
        F = np.array([[1.0, 0.0], [0.0, 0.01]])
        Y = np.ones(2)
        _, cond, _ = solve_fk_y(F, Y)
        assert cond == pytest.approx(100.0, rel=0.01)
