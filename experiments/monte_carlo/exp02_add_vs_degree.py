"""Experiment 02: ADD vs Degree for Pearson III.

Tests how Average Detection Delay changes with approximation degree
s=1..4 for Pearson III distribution at different skewness levels
(gamma_3 = 0, 2, 10) and mean shift delta=0.3.
"""

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.utils.distributions import create_distribution
from gsa_cpd.utils.metrics import compute_add, compute_far


def run_single_trial(detector, h0_stream, h1_stream, change_point):
    """Run a single detection trial and return detection time."""
    detector.reset()
    for t, x in enumerate(h0_stream):
        if detector.predict(x):
            return t
    for t, x in enumerate(h1_stream):
        if detector.predict(x):
            return change_point + t
    return None


def run_experiment(degrees=None, skewness_values=None, n_cal=1000,
                   n_test=500, n_trials=500, delta=0.3, seed=42,
                   results_dir=None):
    """Run ADD vs degree experiment.

    Args:
        degrees: List of approximation degrees to test.
        skewness_values: List of gamma_3 values for Pearson III.
        n_cal: Calibration sample size.
        n_test: Test stream length per regime.
        n_trials: Number of Monte Carlo trials.
        delta: Mean shift magnitude.
        seed: Random seed.
        results_dir: Directory for saving results.
    """
    if degrees is None:
        degrees = [1, 2, 3, 4]
    if skewness_values is None:
        skewness_values = [0.0, 2.0, 10.0]

    rng = np.random.default_rng(seed)
    change_point = n_test

    results = {}

    for skew in skewness_values:
        print(f"\n--- Pearson III (gamma_3={skew}) ---")
        dist_h0 = create_distribution("pearson3", {"skew": skew})
        dist_h1 = create_distribution("pearson3",
                                      {"skew": skew, "loc": delta})

        results[f"skew={skew}"] = {}

        for s in degrees:
            print(f"  degree s={s}...", end="", flush=True)
            det_times = []

            # Use higher ridge for s>=3 to handle ill-conditioning
            ridge = 1e-4 if s >= 3 else 1e-6

            for trial in range(n_trials):
                cal_data = dist_h0.rvs(n_cal, random_state=rng)
                detector = GSADetector(
                    basis=BasisType.POLY,
                    degree=s,
                    epsilon=0.01,
                    threshold_type=ThresholdType.CHEBYSHEV,
                    ridge_lambda=ridge,
                )
                detector.fit(cal_data, delta=delta)

                h0_stream = dist_h0.rvs(n_test, random_state=rng)
                h1_stream = dist_h1.rvs(n_test, random_state=rng)

                det_time = run_single_trial(detector, h0_stream, h1_stream,
                                            change_point)
                det_times.append(det_time)

            valid = [t for t in det_times
                     if t is not None and t >= change_point]
            add = compute_add(valid, change_point) if valid else float("nan")
            far = compute_far(det_times, change_point)
            det_rate = len(valid) / n_trials

            results[f"skew={skew}"][f"s={s}"] = {
                "ADD": add, "FAR": far, "detection_rate": det_rate,
            }
            add_str = f"{add:.2f}" if not np.isnan(add) else "N/A"
            print(f" ADD={add_str}, FAR={far:.4f}, det_rate={det_rate:.2f}")

    # --- Print summary table ---
    print("\n" + "=" * 75)
    print("Experiment 02: ADD vs Degree (Pearson III, poly basis)")
    print(f"  N_cal={n_cal}, N_test={n_test}, trials={n_trials}, delta={delta}")
    print("=" * 75)

    header = f"{'gamma_3':<12}"
    for s in degrees:
        header += f"{'s=' + str(s):>14}"
    print(header)
    print("-" * 75)

    for skew in skewness_values:
        row = f"{skew:<12.1f}"
        for s in degrees:
            entry = results[f"skew={skew}"][f"s={s}"]
            add = entry["ADD"]
            add_str = f"{add:.2f}" if not np.isnan(add) else "N/A"
            row += f"{add_str:>14}"
        print(row)
    print("-" * 75)
    print("(Values are ADD; lower is better)")

    # --- Save results ---
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp02_add_vs_degree",
            "params": {"degrees": degrees, "skewness_values": skewness_values,
                       "n_cal": n_cal, "n_test": n_test,
                       "n_trials": n_trials, "delta": delta},
            "results": results,
        }
        path = os.path.join(results_dir, "exp02_add_vs_degree.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nResults saved to {path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Exp02: ADD vs degree for Pearson III")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with fewer trials")
    parser.add_argument("--n_cal", type=int, default=1000)
    parser.add_argument("--n_test", type=int, default=500)
    parser.add_argument("--n_trials", type=int, default=500)
    parser.add_argument("--delta", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_cal = 500
        args.n_test = 200
        args.n_trials = 50

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp02_{ts}")

    run_experiment(
        n_cal=args.n_cal,
        n_test=args.n_test,
        n_trials=args.n_trials,
        delta=args.delta,
        seed=args.seed,
        results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
