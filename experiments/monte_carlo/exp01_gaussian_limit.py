"""Experiment 01: Gaussian Limit Verification (Theorem 1).

Verifies that at S=1 polynomial basis, the GSA-LLR approximation
reproduces classical CUSUM on Gaussian data. Compares ADD across
poly/frac/log bases at degree s=1 on Normal(0,1) -> Normal(delta,1).
"""

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.utils.distributions import create_distribution, compute_true_llr
from gsa_cpd.utils.metrics import compute_add, compute_far
from gsa_cpd.baselines import OracleCUSUM


def run_single_trial(
    detector,
    cal_data: np.ndarray,
    h0_stream: np.ndarray,
    h1_stream: np.ndarray,
    change_point: int,
):
    """Run a single detection trial and return detection time."""
    detector.reset()
    for t, x in enumerate(h0_stream):
        if detector.predict(x):
            return t  # false alarm before change
    for t, x in enumerate(h1_stream):
        if detector.predict(x):
            return change_point + t
    return None


def run_experiment(n_cal=1000, n_test=500, n_trials=100, delta=0.3,
                   seed=42, results_dir=None):
    """Run Gaussian limit experiment.

    Args:
        n_cal: Calibration sample size (H0 data for fitting).
        n_test: Test stream length per regime (H0 + H1 segments).
        n_trials: Number of Monte Carlo trials.
        delta: Mean shift magnitude (H0: mu=0, H1: mu=delta).
        seed: Random seed.
        results_dir: Directory for saving results.
    """
    rng = np.random.default_rng(seed)

    dist_h0 = create_distribution("normal", {"loc": 0.0, "scale": 1.0})
    dist_h1 = create_distribution("normal", {"loc": delta, "scale": 1.0})

    change_point = n_test  # change happens after n_test H0 observations

    # Configurations to compare
    configs = {
        "GSA poly s=1": {"basis": BasisType.POLY, "degree": 1},
        "GSA frac s=1": {"basis": BasisType.FRAC, "degree": 1},
        "GSA log  s=1": {"basis": BasisType.LOG, "degree": 1},
    }

    results = {}

    # --- Oracle CUSUM baseline ---
    print("Running Oracle CUSUM baseline...")
    oracle_llr = lambda x: compute_true_llr(x, dist_h0, dist_h1)
    oracle_det_times = []
    for trial in range(n_trials):
        cal_data = dist_h0.rvs(n_cal, random_state=rng)
        oracle = OracleCUSUM(oracle_llr, epsilon=0.01)
        oracle.fit(cal_data)

        h0_stream = dist_h0.rvs(n_test, random_state=rng)
        h1_stream = dist_h1.rvs(n_test, random_state=rng)

        det_time = run_single_trial(oracle, cal_data, h0_stream, h1_stream,
                                    change_point)
        oracle_det_times.append(det_time)

    valid_oracle = [t for t in oracle_det_times if t is not None and t >= change_point]
    oracle_add = compute_add(valid_oracle, change_point) if valid_oracle else float("nan")
    oracle_far = compute_far(oracle_det_times, change_point)
    results["Oracle CUSUM"] = {"ADD": oracle_add, "FAR": oracle_far,
                               "detected": len(valid_oracle)}

    # --- GSA detectors ---
    for name, cfg in configs.items():
        print(f"Running {name}...")
        det_times = []

        for trial in range(n_trials):
            cal_data = dist_h0.rvs(n_cal, random_state=rng)
            detector = GSADetector(
                basis=cfg["basis"],
                degree=cfg["degree"],
                epsilon=0.01,
                threshold_type=ThresholdType.CHEBYSHEV,
            )
            detector.fit(cal_data, delta=delta)

            h0_stream = dist_h0.rvs(n_test, random_state=rng)
            h1_stream = dist_h1.rvs(n_test, random_state=rng)

            det_time = run_single_trial(detector, cal_data, h0_stream,
                                        h1_stream, change_point)
            det_times.append(det_time)

        valid = [t for t in det_times if t is not None and t >= change_point]
        add = compute_add(valid, change_point) if valid else float("nan")
        far = compute_far(det_times, change_point)
        results[name] = {"ADD": add, "FAR": far, "detected": len(valid)}

    # --- Print results table ---
    print("\n" + "=" * 65)
    print("Experiment 01: Gaussian Limit (Theorem 1)")
    print(f"  N_cal={n_cal}, N_test={n_test}, trials={n_trials}, delta={delta}")
    print("=" * 65)
    print(f"{'Method':<20} {'ADD':>8} {'FAR':>8} {'Detected':>10}")
    print("-" * 65)
    for name, r in results.items():
        add_str = f"{r['ADD']:.2f}" if not np.isnan(r["ADD"]) else "N/A"
        print(f"{name:<20} {add_str:>8} {r['FAR']:>8.4f} {r['detected']:>10}")
    print("-" * 65)

    # Check Theorem 1: poly s=1 should be close to Oracle on Gaussian
    poly_add = results["GSA poly s=1"]["ADD"]
    if not np.isnan(poly_add) and not np.isnan(oracle_add):
        ratio = poly_add / oracle_add if oracle_add > 0 else float("inf")
        status = "PASS" if ratio < 2.0 else "FAIL"
        print(f"\nTheorem 1 check: ADD(poly_s1)/ADD(oracle) = {ratio:.3f}  [{status}]")
    else:
        print("\nTheorem 1 check: insufficient detections to compare.")

    # --- Save results ---
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp01_gaussian_limit",
            "params": {"n_cal": n_cal, "n_test": n_test,
                       "n_trials": n_trials, "delta": delta},
            "results": results,
        }
        path = os.path.join(results_dir, "exp01_gaussian_limit.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nResults saved to {path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Exp01: Gaussian limit verification (Theorem 1)")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with fewer trials")
    parser.add_argument("--n_cal", type=int, default=1000)
    parser.add_argument("--n_test", type=int, default=500)
    parser.add_argument("--n_trials", type=int, default=100)
    parser.add_argument("--delta", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_cal = 500
        args.n_test = 200
        args.n_trials = 20

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp01_{ts}")

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
