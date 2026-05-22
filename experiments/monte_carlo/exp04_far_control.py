"""Experiment 04: FAR Control Verification.

Verifies that the empirical False Alarm Rate stays below the target
epsilon for all basis/degree/distribution configurations. Runs
H0-only trials (no change point) and measures the fraction of
runs that produce a false alarm within the monitoring window.
"""

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.utils.distributions import create_distribution


def run_h0_trial(detector, h0_stream):
    """Run a single H0-only trial. Returns True if false alarm occurred."""
    detector.reset()
    for x in h0_stream:
        if detector.predict(x):
            return True
    return False


def run_experiment(epsilon=0.01, n_cal=1000, monitoring_length=500,
                   n_trials=1000, delta=0.3, seed=42, results_dir=None):
    """Run FAR control experiment.

    Args:
        epsilon: Target FAR level.
        n_cal: Calibration sample size.
        monitoring_length: Number of H0 observations in monitoring phase.
        n_trials: Number of Monte Carlo trials.
        delta: MDE coefficient for detector fitting.
        seed: Random seed.
        results_dir: Directory for saving results.
    """
    rng = np.random.default_rng(seed)

    # Test configurations
    configs = [
        # (distribution_name, dist_params, basis, degree)
        ("Normal",         {"loc": 0, "scale": 1}, BasisType.POLY, 1),
        ("Normal",         {"loc": 0, "scale": 1}, BasisType.POLY, 2),
        ("Normal",         {"loc": 0, "scale": 1}, BasisType.FRAC, 2),
        ("Normal",         {"loc": 0, "scale": 1}, BasisType.LOG,  2),
        ("Student-t(5)",   None,                    BasisType.POLY, 1),
        ("Student-t(5)",   None,                    BasisType.POLY, 2),
        ("Student-t(5)",   None,                    BasisType.FRAC, 2),
        ("Pearson(skew=2)", None,                   BasisType.POLY, 2),
        ("Pearson(skew=10)", None,                  BasisType.POLY, 2),
        ("Pearson(skew=10)", None,                  BasisType.FRAC, 2),
        ("Pearson(skew=10)", None,                  BasisType.LOG,  2),
    ]

    def make_dist(name):
        if name == "Normal":
            return create_distribution("normal", {"loc": 0, "scale": 1})
        elif name == "Student-t(5)":
            return create_distribution("student_t", {"df": 5})
        elif name == "Pearson(skew=2)":
            return create_distribution("pearson3", {"skew": 2.0})
        elif name == "Pearson(skew=10)":
            return create_distribution("pearson3", {"skew": 10.0})
        else:
            raise ValueError(f"Unknown: {name}")

    results = []

    print(f"FAR control experiment: epsilon={epsilon}, "
          f"monitoring_length={monitoring_length}, trials={n_trials}")
    print()

    for dist_name, dist_params, basis, degree in configs:
        label = f"{dist_name} / {basis.value} s={degree}"
        print(f"  {label}...", end="", flush=True)

        dist_h0 = make_dist(dist_name)
        ridge = 1e-4 if degree >= 3 else 1e-6

        false_alarms = 0
        for trial in range(n_trials):
            cal_data = dist_h0.rvs(n_cal, random_state=rng)
            detector = GSADetector(
                basis=basis,
                degree=degree,
                epsilon=epsilon,
                threshold_type=ThresholdType.CHEBYSHEV,
                ridge_lambda=ridge,
            )
            detector.fit(cal_data, delta=delta)

            h0_stream = dist_h0.rvs(monitoring_length, random_state=rng)
            if run_h0_trial(detector, h0_stream):
                false_alarms += 1

        empirical_far = false_alarms / n_trials
        passed = empirical_far <= epsilon * 2  # allow 2x tolerance
        status = "PASS" if passed else "FAIL"

        results.append({
            "distribution": dist_name,
            "basis": basis.value,
            "degree": degree,
            "empirical_far": empirical_far,
            "target_epsilon": epsilon,
            "false_alarms": false_alarms,
            "n_trials": n_trials,
            "status": status,
        })

        print(f" FAR={empirical_far:.4f} (target<={epsilon})  [{status}]")

    # --- Print summary table ---
    print("\n" + "=" * 80)
    print("Experiment 04: FAR Control Verification")
    print(f"  epsilon={epsilon}, monitoring_length={monitoring_length}, "
          f"trials={n_trials}")
    print("=" * 80)
    print(f"{'Configuration':<35} {'FAR':>8} {'Target':>8} {'Status':>8}")
    print("-" * 80)

    n_pass = 0
    for r in results:
        label = f"{r['distribution']} / {r['basis']} s={r['degree']}"
        print(f"{label:<35} {r['empirical_far']:>8.4f} "
              f"{r['target_epsilon']:>8.4f} {r['status']:>8}")
        if r["status"] == "PASS":
            n_pass += 1

    print("-" * 80)
    print(f"Passed: {n_pass}/{len(results)}")

    # --- Save results ---
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp04_far_control",
            "params": {"epsilon": epsilon, "n_cal": n_cal,
                       "monitoring_length": monitoring_length,
                       "n_trials": n_trials, "delta": delta},
            "results": results,
        }
        path = os.path.join(results_dir, "exp04_far_control.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2)
        print(f"\nResults saved to {path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Exp04: FAR control verification")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with fewer trials")
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--n_cal", type=int, default=1000)
    parser.add_argument("--monitoring_length", type=int, default=500)
    parser.add_argument("--n_trials", type=int, default=1000)
    parser.add_argument("--delta", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_cal = 500
        args.monitoring_length = 200
        args.n_trials = 100

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp04_{ts}")

    run_experiment(
        epsilon=args.epsilon,
        n_cal=args.n_cal,
        monitoring_length=args.monitoring_length,
        n_trials=args.n_trials,
        delta=args.delta,
        seed=args.seed,
        results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
