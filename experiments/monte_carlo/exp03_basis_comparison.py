"""Experiment 03: Basis Function Comparison.

Compares poly/frac/log bases at degree s=2 on heavy-tailed and
skewed distributions: Student-t(df=5) and Pearson III(skew=10).
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


def run_experiment(degree=2, n_cal=1000, n_test=500, n_trials=500,
                   delta=0.3, seed=42, results_dir=None):
    """Run basis comparison experiment.

    Args:
        degree: Approximation degree s (fixed at 2 for this experiment).
        n_cal: Calibration sample size.
        n_test: Test stream length per regime.
        n_trials: Number of Monte Carlo trials.
        delta: Mean shift magnitude.
        seed: Random seed.
        results_dir: Directory for saving results.
    """
    rng = np.random.default_rng(seed)
    change_point = n_test

    distributions = {
        "Student-t(df=5)": {
            "h0": ("student_t", {"df": 5}),
            "h1": ("student_t", {"df": 5, "loc": delta}),
        },
        "Pearson III(skew=10)": {
            "h0": ("pearson3", {"skew": 10.0}),
            "h1": ("pearson3", {"skew": 10.0, "loc": delta}),
        },
    }

    bases = [BasisType.POLY, BasisType.FRAC, BasisType.LOG]

    results = {}

    for dist_name, dist_cfg in distributions.items():
        print(f"\n--- {dist_name} ---")
        dist_h0 = create_distribution(*dist_cfg["h0"])
        dist_h1 = create_distribution(*dist_cfg["h1"])

        results[dist_name] = {}

        for basis in bases:
            basis_name = basis.value
            print(f"  basis={basis_name}...", end="", flush=True)
            det_times = []

            for trial in range(n_trials):
                cal_data = dist_h0.rvs(n_cal, random_state=rng)
                detector = GSADetector(
                    basis=basis,
                    degree=degree,
                    epsilon=0.01,
                    threshold_type=ThresholdType.CHEBYSHEV,
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

            # Collect diagnostics from last trial
            diag = detector.diagnostics
            cond_num = diag.condition_number if diag else None
            eta = diag.eta if diag else None

            results[dist_name][basis_name] = {
                "ADD": add, "FAR": far, "detection_rate": det_rate,
                "cond_F": cond_num, "eta": eta,
            }
            add_str = f"{add:.2f}" if not np.isnan(add) else "N/A"
            print(f" ADD={add_str}, FAR={far:.4f}, eta={eta:.3f}")

    # --- Print summary table ---
    print("\n" + "=" * 75)
    print(f"Experiment 03: Basis Comparison (degree s={degree})")
    print(f"  N_cal={n_cal}, N_test={n_test}, trials={n_trials}, delta={delta}")
    print("=" * 75)

    for dist_name in distributions:
        print(f"\n  {dist_name}:")
        print(f"  {'Basis':<10} {'ADD':>8} {'FAR':>8} {'DetRate':>8} "
              f"{'cond(F)':>12} {'eta':>8}")
        print("  " + "-" * 60)
        for basis in bases:
            r = results[dist_name][basis.value]
            add_str = f"{r['ADD']:.2f}" if not np.isnan(r["ADD"]) else "N/A"
            cond_str = (f"{r['cond_F']:.1e}" if r["cond_F"] is not None
                        else "N/A")
            eta_str = f"{r['eta']:.3f}" if r["eta"] is not None else "N/A"
            print(f"  {basis.value:<10} {add_str:>8} {r['FAR']:>8.4f} "
                  f"{r['detection_rate']:>8.2f} {cond_str:>12} {eta_str:>8}")
        print("  " + "-" * 60)

    # --- Save results ---
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp03_basis_comparison",
            "params": {"degree": degree, "n_cal": n_cal, "n_test": n_test,
                       "n_trials": n_trials, "delta": delta},
            "results": results,
        }
        path = os.path.join(results_dir, "exp03_basis_comparison.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nResults saved to {path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Exp03: Basis comparison at s=2")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with fewer trials")
    parser.add_argument("--degree", type=int, default=2)
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
        args.results_dir = os.path.join("results", f"exp03_{ts}")

    run_experiment(
        degree=args.degree,
        n_cal=args.n_cal,
        n_test=args.n_test,
        n_trials=args.n_trials,
        delta=args.delta,
        seed=args.seed,
        results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
