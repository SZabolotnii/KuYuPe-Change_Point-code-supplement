"""Run all Monte Carlo experiments (exp01-exp08).

Usage:
    python run_all.py           # full run
    python run_all.py --quick   # quick run with fewer trials
"""

import argparse
import os
import sys
import time
from datetime import datetime

from experiments.monte_carlo.exp01_gaussian_limit import (
    run_experiment as run_exp01,
)
from experiments.monte_carlo.exp02_add_vs_degree import (
    run_experiment as run_exp02,
)
from experiments.monte_carlo.exp03_basis_comparison import (
    run_experiment as run_exp03,
)
from experiments.monte_carlo.exp04_far_control import (
    run_experiment as run_exp04,
)
from experiments.monte_carlo.exp05_stopping_rules import (
    run_experiment as run_exp05,
)
from experiments.monte_carlo.exp06_oc_curves import (
    run_experiment as run_exp06,
)
from experiments.monte_carlo.exp07_drift_verification import (
    run_experiment as run_exp07,
)
from experiments.monte_carlo.exp08_f_misspecification import (
    run_experiment as run_exp08,
)


def main():
    parser = argparse.ArgumentParser(
        description="Run all Monte Carlo experiments")
    parser.add_argument("--quick", action="store_true",
                        help="Quick run with reduced trial counts")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.results_dir is None:
        args.results_dir = os.path.join("results", f"mc_all_{ts}")
    os.makedirs(args.results_dir, exist_ok=True)

    if args.quick:
        n_cal, n_test = 500, 200
        n_trials_small, n_trials_large = 20, 50
        n_trials_far = 100
        target_arl = 200
        # exp06 OC curves (reduced smoke grid).
        oc_dists = ["normal", "student_t5"]
        oc_arl0s = [300.0, 600.0]
        oc_kwargs = dict(n_fit=3000, n_calib=150, n_verify=200, n_add=600)
        # exp07 drift verification.
        drift_kwargs = dict(n_samples=10_000, oracle_n=40_000)
        # exp08 F-misspecification.
        fmis_kwargs = dict(n_trials=1000, n_calibration=1500,
                           oracle_n=30000, calib_runs=250)
    else:
        n_cal, n_test = 1000, 500
        n_trials_small, n_trials_large = 100, 500
        n_trials_far = 1000
        target_arl = 500
        oc_dists = ["normal", "student_t5", "laplace"]
        oc_arl0s = [500.0, 1000.0]
        oc_kwargs = dict(n_fit=5000, n_calib=400, n_verify=400, n_add=2000)
        drift_kwargs = dict(n_samples=100_000, oracle_n=200_000)
        fmis_kwargs = dict(n_trials=2000, n_calibration=2000,
                           oracle_n=40000, calib_runs=350)

    experiments = [
        (
            "Exp01: Gaussian Limit",
            lambda: run_exp01(
                n_cal=n_cal, n_test=n_test,
                n_trials=n_trials_small, seed=args.seed,
                results_dir=args.results_dir,
            ),
        ),
        (
            "Exp02: ADD vs Degree",
            lambda: run_exp02(
                n_cal=n_cal, n_test=n_test,
                n_trials=n_trials_large, seed=args.seed,
                results_dir=args.results_dir,
            ),
        ),
        (
            "Exp03: Basis Comparison",
            lambda: run_exp03(
                n_cal=n_cal, n_test=n_test,
                n_trials=n_trials_large, seed=args.seed,
                results_dir=args.results_dir,
            ),
        ),
        (
            "Exp04: FAR Control",
            lambda: run_exp04(
                n_cal=n_cal, monitoring_length=n_test,
                n_trials=n_trials_far, seed=args.seed,
                results_dir=args.results_dir,
            ),
        ),
        (
            "Exp05: Stopping Rules",
            lambda: run_exp05(
                n_cal=n_cal, n_test=n_test,
                n_trials=n_trials_large, target_arl=target_arl,
                seed=args.seed, results_dir=args.results_dir,
            ),
        ),
        (
            "Exp06: OC Curves",
            lambda: run_exp06(
                dists=oc_dists, arl0s=oc_arl0s, seed=args.seed,
                results_dir=args.results_dir, **oc_kwargs,
            ),
        ),
        (
            "Exp07: Drift Verification",
            lambda: run_exp07(
                seed=args.seed, results_dir=args.results_dir, **drift_kwargs,
            ),
        ),
        (
            "Exp08: F-Misspecification",
            lambda: run_exp08(
                seed=args.seed, results_dir=args.results_dir, **fmis_kwargs,
            ),
        ),
    ]

    total_start = time.time()
    print("=" * 70)
    mode = "QUICK" if args.quick else "FULL"
    print(f"Running all Monte Carlo experiments ({mode} mode)")
    print(f"Results directory: {args.results_dir}")
    print("=" * 70)

    for name, run_fn in experiments:
        print(f"\n{'#' * 70}")
        print(f"# {name}")
        print(f"{'#' * 70}")
        start = time.time()
        try:
            run_fn()
            elapsed = time.time() - start
            print(f"\n[{name}] completed in {elapsed:.1f}s")
        except Exception as e:
            elapsed = time.time() - start
            print(f"\n[{name}] FAILED after {elapsed:.1f}s: {e}")

    total_elapsed = time.time() - total_start
    print(f"\n{'=' * 70}")
    print(f"All experiments completed in {total_elapsed:.1f}s")
    print(f"Results saved to: {args.results_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
