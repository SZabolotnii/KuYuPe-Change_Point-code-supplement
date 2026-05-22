"""Experiment 05: Stopping Rule Comparison.

Compares CUSUM, GRSh, and SRP stopping rules using the same
GSA-LLR approximation on Pearson III(skew=10) data. Thresholds
are calibrated to match the same target ARL under H0 for a fair
comparison.
"""

import argparse
import json
import os
import sys
from datetime import datetime

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.stopping_rules import CUSUMRule, GRShRule, SRPRule
from gsa_cpd.utils.distributions import create_distribution
from gsa_cpd.utils.metrics import compute_add, compute_far


def calibrate_threshold_by_arl(compute_llr_fn, dist_h0, target_arl,
                               rule_class, n_runs=200,
                               max_run_length=5000, rng=None):
    """Calibrate threshold via binary search on ARL for a given rule.

    Args:
        compute_llr_fn: Function x -> Lambda(x).
        dist_h0: H0 distribution for sampling.
        target_arl: Target average run length under H0.
        rule_class: One of CUSUMRule, GRShRule, SRPRule.
        n_runs: MC runs for ARL estimation.
        max_run_length: Max run length per trial.
        rng: Random generator.

    Returns:
        Calibrated threshold.
    """
    if rng is None:
        rng = np.random.default_rng()

    h_low = 0.1
    h_high = 1000.0

    # For SRP, thresholds are typically much larger
    if rule_class is SRPRule:
        h_high = np.exp(12.0)

    def estimate_arl(h):
        run_lengths = []
        for _ in range(n_runs):
            data = dist_h0.rvs(max_run_length, random_state=rng)
            rule = rule_class(threshold=h)
            for t, x in enumerate(data):
                llr = compute_llr_fn(x)
                if rule.update(llr):
                    run_lengths.append(t + 1)
                    break
            else:
                run_lengths.append(max_run_length)
        return float(np.mean(run_lengths))

    for _ in range(15):
        h_mid = (h_low + h_high) / 2.0
        if rule_class is SRPRule:
            h_mid = np.exp((np.log(h_low) + np.log(h_high)) / 2.0)

        arl = estimate_arl(h_mid)
        if abs(arl - target_arl) / target_arl < 0.1:
            return h_mid
        if arl < target_arl:
            h_low = h_mid
        else:
            h_high = h_mid

    return (h_low + h_high) / 2.0


def run_experiment(skew=10.0, degree=2, n_cal=1000, n_test=500,
                   n_trials=500, delta=0.3, target_arl=500,
                   seed=42, results_dir=None):
    """Run stopping rule comparison experiment.

    Args:
        skew: Skewness parameter for Pearson III.
        degree: Approximation degree for GSA.
        n_cal: Calibration sample size.
        n_test: Test stream length per regime.
        n_trials: Number of MC trials for ADD measurement.
        delta: Mean shift magnitude.
        target_arl: Target ARL for threshold calibration.
        seed: Random seed.
        results_dir: Directory for saving results.
    """
    rng = np.random.default_rng(seed)
    change_point = n_test

    dist_h0 = create_distribution("pearson3", {"skew": skew})
    dist_h1 = create_distribution("pearson3", {"skew": skew, "loc": delta})

    # Step 1: Fit GSA detector to get LLR function
    print(f"Fitting GSA detector (basis=poly, s={degree}, "
          f"Pearson III skew={skew})...")
    cal_data = dist_h0.rvs(n_cal * 2, random_state=rng)
    gsa = GSADetector(
        basis=BasisType.POLY,
        degree=degree,
        epsilon=0.01,
        threshold_type=ThresholdType.CHEBYSHEV,
    )
    gsa.fit(cal_data, delta=delta)

    # Extract LLR function from fitted detector
    def compute_llr(x):
        return gsa._compute_llr(x)

    # Step 2: Calibrate thresholds for each stopping rule
    rules_info = {
        "CUSUM": CUSUMRule,
        "GRSh": GRShRule,
        "SRP": SRPRule,
    }

    thresholds = {}
    print(f"\nCalibrating thresholds (target ARL={target_arl})...")
    for rule_name, rule_class in rules_info.items():
        print(f"  {rule_name}...", end="", flush=True)
        h = calibrate_threshold_by_arl(
            compute_llr, dist_h0, target_arl, rule_class,
            n_runs=100, max_run_length=target_arl * 3, rng=rng,
        )
        thresholds[rule_name] = h
        print(f" h={h:.4f}")

    # Step 3: Run detection trials
    results = {}

    for rule_name, rule_class in rules_info.items():
        print(f"\nRunning {rule_name} trials...", end="", flush=True)
        h = thresholds[rule_name]
        det_times = []

        for trial in range(n_trials):
            rule = rule_class(threshold=h)

            # H0 phase
            h0_stream = dist_h0.rvs(n_test, random_state=rng)
            alarm = False
            alarm_time = None
            for t, x in enumerate(h0_stream):
                llr = compute_llr(x)
                if rule.update(llr):
                    alarm_time = t
                    alarm = True
                    break

            if not alarm:
                # H1 phase
                h1_stream = dist_h1.rvs(n_test, random_state=rng)
                for t, x in enumerate(h1_stream):
                    llr = compute_llr(x)
                    if rule.update(llr):
                        alarm_time = change_point + t
                        break

            det_times.append(alarm_time)

        valid = [t for t in det_times if t is not None and t >= change_point]
        add = compute_add(valid, change_point) if valid else float("nan")
        far = compute_far(det_times, change_point)
        det_rate = len(valid) / n_trials

        results[rule_name] = {
            "ADD": add, "FAR": far, "detection_rate": det_rate,
            "threshold": h, "n_detected": len(valid),
        }

        add_str = f"{add:.2f}" if not np.isnan(add) else "N/A"
        print(f" ADD={add_str}, FAR={far:.4f}")

    # --- Print summary table ---
    print("\n" + "=" * 70)
    print("Experiment 05: Stopping Rule Comparison")
    print(f"  Pearson III(skew={skew}), GSA poly s={degree}, "
          f"delta={delta}, ARL={target_arl}")
    print(f"  N_test={n_test}, trials={n_trials}")
    print("=" * 70)
    print(f"{'Rule':<10} {'ADD':>8} {'FAR':>8} {'DetRate':>8} "
          f"{'Threshold':>12}")
    print("-" * 70)

    for rule_name in rules_info:
        r = results[rule_name]
        add_str = f"{r['ADD']:.2f}" if not np.isnan(r["ADD"]) else "N/A"
        print(f"{rule_name:<10} {add_str:>8} {r['FAR']:>8.4f} "
              f"{r['detection_rate']:>8.2f} {r['threshold']:>12.4f}")
    print("-" * 70)

    # --- Save results ---
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp05_stopping_rules",
            "params": {"skew": skew, "degree": degree, "n_cal": n_cal,
                       "n_test": n_test, "n_trials": n_trials,
                       "delta": delta, "target_arl": target_arl},
            "results": results,
        }
        path = os.path.join(results_dir, "exp05_stopping_rules.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nResults saved to {path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Exp05: Stopping rule comparison (CUSUM/GRSh/SRP)")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with fewer trials")
    parser.add_argument("--skew", type=float, default=10.0)
    parser.add_argument("--degree", type=int, default=2)
    parser.add_argument("--n_cal", type=int, default=1000)
    parser.add_argument("--n_test", type=int, default=500)
    parser.add_argument("--n_trials", type=int, default=500)
    parser.add_argument("--delta", type=float, default=0.3)
    parser.add_argument("--target_arl", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_cal = 500
        args.n_test = 200
        args.n_trials = 50
        args.target_arl = 200

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp05_{ts}")

    run_experiment(
        skew=args.skew,
        degree=args.degree,
        n_cal=args.n_cal,
        n_test=args.n_test,
        n_trials=args.n_trials,
        delta=args.delta,
        target_arl=args.target_arl,
        seed=args.seed,
        results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
