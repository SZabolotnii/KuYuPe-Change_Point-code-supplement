"""Metrics for evaluating change-point detection performance."""

import numpy as np
from typing import List, Optional


def compute_add(
    detection_times: List[int],
    change_point: int,
    max_delay: int = 1000,
) -> float:
    """Compute Average Detection Delay (ADD).

    Args:
        detection_times: List of detection times across MC runs.
        change_point: True change-point time tau.
        max_delay: Penalty for missed detections.

    Returns:
        Average detection delay.
    """
    delays = []
    for t in detection_times:
        if t >= change_point:
            delays.append(t - change_point)
        else:
            delays.append(max_delay)  # penalize false alarm before change
    return float(np.mean(delays)) if delays else max_delay


def compute_far(
    detection_times: List[Optional[int]],
    change_point: int,
) -> float:
    """Compute False Alarm Rate (FAR).

    Args:
        detection_times: List of detection times (None = no detection).
        change_point: True change-point time tau.

    Returns:
        Fraction of runs with detection before change_point.
    """
    if not detection_times:
        return 0.0
    false_alarms = sum(1 for t in detection_times
                       if t is not None and t < change_point)
    return false_alarms / len(detection_times)


def compute_detection_rate(
    detection_times: List[Optional[int]],
    change_point: int,
    max_time: int = 1000,
) -> float:
    """Compute detection rate (fraction of changes detected).

    Args:
        detection_times: List of detection times (None = no detection).
        change_point: True change-point time tau.
        max_time: Maximum observation time.

    Returns:
        Fraction of runs where change was detected after tau.
    """
    if not detection_times:
        return 0.0
    detected = sum(1 for t in detection_times
                   if t is not None and change_point <= t <= max_time)
    return detected / len(detection_times)


def efficiency_coefficient(
    E_L_H1: float,
    E_L_H0: float,
    Var_L_H0: float,
) -> float:
    """Compute efficiency coefficient eta = (E[L|H1] - E[L|H0]) / sqrt(Var[L|H0]).

    Args:
        E_L_H1: Expected Lambda under H1.
        E_L_H0: Expected Lambda under H0.
        Var_L_H0: Variance of Lambda under H0.

    Returns:
        Efficiency coefficient eta.
    """
    sigma = np.sqrt(max(Var_L_H0, 1e-12))
    return (E_L_H1 - E_L_H0) / sigma
