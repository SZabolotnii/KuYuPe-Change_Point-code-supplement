"""Experiment 08: Covariance-matrix (F) misspecification ablation.

The GSA detector solves the linear system F K = Y for its LLR coefficients,
where F = Cov0 + Cov1 is the (pooled) covariance of the basis functions and
Y = E[phi | H1] - E[phi | H0] is the moment-difference direction. This study
quantifies how the *source* used to estimate F affects detection delay (ADD):

    (a) training mixture : F from a blend of pre-change (H0) and post-change
                           (H1) basis observations (realistic labelled training).
    (b) H0 only          : F from H0 data only; Cov1 approximated from the H0
                           covariance scaled by the MDE delta (the standard
                           online setting, no post-change data).
    (c) oracle           : F from large-sample Cov0 and Cov1 (basis covariances
                           under each hypothesis from large i.i.d. samples).

Across the three sources we hold the moment-difference Y and the basis fixed,
so only F differs. Following the repo invariant ("equalize FAR before comparing
ADD"), the detection threshold is recalibrated per F-source on independent H0
streams to hit a common target FAR; ADD is then read at matched FAR.

Output (under --results_dir, default results/exp08_<ts>/):
    exp08_f_misspecification.json   raw per-source rows + degradation summary

Usage:
    python exp08_f_misspecification.py --quick
    python exp08_f_misspecification.py
"""

import argparse
import json
import os
from datetime import datetime
from typing import Dict, Any, List

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.core.basis import evaluate_basis_matrix
from gsa_cpd.core.solver import solve_fk_y
from gsa_cpd.core.threshold import compute_threshold
from gsa_cpd.core.moments import compute_mde_hypothesis
from gsa_cpd.utils.distributions import create_distribution
from gsa_cpd.utils.metrics import compute_add, compute_far


# Distributions span tail/asymmetry regimes. mean_shift drives the change; the
# poly basis is used throughout so F is well conditioned.
SCENARIOS = [
    {
        "id": "normal",
        "label": "Normal",
        "dist_type": "normal",
        "H0": {"loc": 0.0, "scale": 1.0},
        "H1": {"loc": 0.6, "scale": 1.0},
        "mean_shift": 0.6,
    },
    {
        "id": "student_t5",
        "label": "Student-t (df=5)",
        "dist_type": "student_t",
        "H0": {"df": 5, "loc": 0.0, "scale": 1.0},
        "H1": {"df": 5, "loc": 0.6, "scale": 1.0},
        "mean_shift": 0.6,
    },
    {
        "id": "laplace",
        "label": "Laplace",
        "dist_type": "laplace",
        "H0": {"loc": 0.0, "scale": 1.0},
        "H1": {"loc": 0.6, "scale": 1.0},
        "mean_shift": 0.6,
    },
]

F_SOURCES = ["mixture", "h0_only", "oracle"]
F_SOURCE_LABELS = {
    "mixture": "Training mixture",
    "h0_only": "H0 only",
    "oracle": "Oracle (true F)",
}

DEGREE = 2
BASIS = BasisType.POLY
PHI_MAX = 10.0
DELTA = 0.3  # MDE shift used by the H0-only source to approximate Cov1


def _cov(phi: np.ndarray, s: int) -> np.ndarray:
    if s == 1:
        return np.array([[np.var(phi)]])
    return np.cov(phi, rowvar=False)


def compute_F(f_source: str, s: int, phi0: np.ndarray, phi1: np.ndarray,
              Cov0_emp: np.ndarray, u0: np.ndarray, std_x: float,
              scenario: Dict[str, Any], oracle_n: int, seed: int,
              ridge_lambda: float) -> np.ndarray:
    """Covariance matrix F = Cov0 + Cov1 from the chosen source.

    Only F varies across sources; the moment-difference Y and the bias are held
    fixed by the caller, so this isolates the effect of F misspecification.

        mixture : F = 2 * Cov(pooled H0 and H1 basis observations).
        h0_only : F = Cov0(H0) + Cov1, Cov1 = MDE-scaled H0 covariance.
        oracle  : F = Cov0(true) + Cov1(true) from large i.i.d. samples.
    """
    if f_source == "mixture":
        pooled = np.vstack([phi0, phi1])
        cov_mix = _cov(pooled, s)
        F = cov_mix + cov_mix
    elif f_source == "oracle":
        dist_h0 = create_distribution(scenario["dist_type"], scenario["H0"])
        dist_h1 = create_distribution(scenario["dist_type"], scenario["H1"])
        big_H0 = dist_h0.rvs(oracle_n, random_state=np.random.default_rng(seed + 778))
        big_H1 = dist_h1.rvs(oracle_n, random_state=np.random.default_rng(seed + 777))
        from gsa_cpd.core.moments import winsorize
        phi0_big = evaluate_basis_matrix(winsorize(big_H0), s, BASIS, PHI_MAX)
        phi1_big = evaluate_basis_matrix(winsorize(big_H1), s, BASIS, PHI_MAX)
        F = _cov(phi0_big, s) + _cov(phi1_big, s)
    else:  # h0_only
        _, Cov1 = compute_mde_hypothesis(u0, Cov0_emp, DELTA, std_x)
        F = Cov0_emp + Cov1
    return F + np.eye(s) * ridge_lambda


def fit_with_F_source(scenario: Dict[str, Any], f_source: str,
                      calib_H0: np.ndarray, calib_H1: np.ndarray,
                      epsilon: float, oracle_n: int, seed: int) -> GSADetector:
    """Build a GSA detector whose covariance F comes from a chosen source.

    The basis is poly. The moment-difference Y = E[phi|H1] - E[phi|H0] (and the
    bias k0) are estimated ONCE from the labelled H0/H1 calibration data and held
    fixed across all F-sources, so the only thing that changes between rows is
    the covariance matrix F. We then solve F K = Y and set the bias/threshold
    exactly as GSADetector.fit does.
    """
    det = GSADetector(
        basis=BASIS, degree=DEGREE, epsilon=epsilon,
        threshold_type=ThresholdType.CHEBYSHEV, threshold_scale=2.0,
        phi_max=PHI_MAX,
    )
    s = det.degree
    ridge_lambda = det.ridge_lambda

    # Basis matrices under H0 / H1 (winsorised, matching the detector's fit()).
    h0 = det._winsorize(calib_H0)
    h1 = det._winsorize(calib_H1)
    phi0 = evaluate_basis_matrix(h0, s, BASIS, PHI_MAX)
    phi1 = evaluate_basis_matrix(h1, s, BASIS, PHI_MAX)
    u0 = np.mean(phi0, axis=0)
    m1 = np.mean(phi1, axis=0)
    Cov0_emp = _cov(phi0, s)
    std_x = float(np.std(calib_H0))

    # Fixed moment-difference direction (same for every F-source).
    Y = m1 - u0

    # F depends on the source.
    F = compute_F(f_source, s, phi0, phi1, Cov0_emp, u0, std_x,
                  scenario, oracle_n, seed, ridge_lambda)

    K, cond_F, method = solve_fk_y(F, Y, ridge_lambda=ridge_lambda)
    det._coeffs = K

    # Bias term and LLR statistics use the empirical H0 moments (the data the
    # detector actually sees online), exactly as in GSADetector.fit.
    det._k0 = -0.5 * float(np.dot(K, m1 + u0))

    E_L_H0 = det._k0 + float(np.dot(K, u0))
    Var_L_H0 = max(float(np.dot(K, Cov0_emp @ K)), 1e-12)

    det._threshold = compute_threshold(
        E_L_H0, Var_L_H0, epsilon, det.threshold_type, det.threshold_scale,
    )
    det._fitted = True
    det._cond_F = cond_F
    det.reset()
    return det


def _alarm_time(detector: GSADetector, stream: np.ndarray):
    """First alarm index on a stream (None if no alarm)."""
    detector.reset()
    for t, x in enumerate(stream):
        if detector.predict(float(x)):
            return t
    return None


def calibrate_threshold(detector: GSADetector, scenario: Dict[str, Any],
                        epsilon: float, n_runs: int, stream_len: int,
                        seed: int) -> float:
    """Binary-search the threshold so empirical FAR ~ epsilon on H0 streams.

    Equalises FAR across F-sources before ADD is compared. The CUSUM
    coefficients (det._coeffs, det._k0) are fixed; only the scalar threshold is
    tuned.
    """
    dist_h0 = create_distribution(scenario["dist_type"], scenario["H0"])
    streams = [dist_h0.rvs(stream_len, random_state=np.random.default_rng(seed + 9000 + r))
               for r in range(n_runs)]

    def empirical_far(h: float) -> float:
        alarms = 0
        for stream in streams:
            detector._threshold = h
            if _alarm_time(detector, stream) is not None:
                alarms += 1
        return alarms / n_runs

    h_init = max(detector._threshold, 1e-3)
    h_low, h_high = h_init * 0.02, h_init * 5.0

    # Expand upper bound until FAR <= target (conservative side).
    for _ in range(10):
        if empirical_far(h_high) <= epsilon:
            break
        h_high *= 2.0
    # Lower bound should sit above target so the bracket is valid.
    for _ in range(10):
        if empirical_far(h_low) >= epsilon:
            break
        h_low *= 0.5

    # Bisection: empirical_far is a deterministic monotone step function of h
    # (the H0 streams are fixed), so it converges to the threshold band whose
    # FAR is closest to the target. Run a fixed, generous number of steps.
    best_h, best_gap = h_init, float("inf")
    for _ in range(30):
        h_mid = 0.5 * (h_low + h_high)
        far_mid = empirical_far(h_mid)
        gap = abs(far_mid - epsilon)
        if gap < best_gap:
            best_gap, best_h = gap, h_mid
        if far_mid > epsilon:
            h_low = h_mid
        else:
            h_high = h_mid
    return best_h


def run_scenario(scenario: Dict[str, Any], f_source: str,
                 params: Dict[str, Any]) -> Dict[str, Any]:
    """Calibrate one (scenario, F-source) detector and measure FAR/ADD."""
    seed = params["seed"]
    dist_h0 = create_distribution(scenario["dist_type"], scenario["H0"])
    dist_h1 = create_distribution(scenario["dist_type"], scenario["H1"])

    # Calibration data: H0 for all; mixture additionally needs H1.
    calib_H0 = dist_h0.rvs(params["n_calibration"],
                           random_state=np.random.default_rng(seed + 1))
    calib_H1 = dist_h1.rvs(params["n_calibration"],
                           random_state=np.random.default_rng(seed + 2))

    detector = fit_with_F_source(
        scenario, f_source, calib_H0, calib_H1,
        epsilon=params["epsilon"], oracle_n=params["oracle_n"], seed=seed,
    )

    # FAR equalization.
    h_cal = calibrate_threshold(
        detector, scenario, params["epsilon"],
        n_runs=params["calib_runs"], stream_len=params["change_point"], seed=seed,
    )
    detector._threshold = h_cal

    # Evaluate ADD/FAR over independent trials (change at tau).
    tau, T_max = params["change_point"], params["t_max"]
    detection_times: List = []
    for trial in range(params["n_trials"]):
        data_H0 = dist_h0.rvs(tau, random_state=np.random.default_rng(seed + 100000 + trial))
        data_H1 = dist_h1.rvs(T_max - tau, random_state=np.random.default_rng(seed + 200000 + trial))
        test = np.concatenate([data_H0, data_H1])
        detection_times.append(_alarm_time(detector, test))

    far = compute_far(detection_times, tau)
    valid = [t for t in detection_times if t is not None and t >= tau]
    add = compute_add(valid, tau) if valid else float("nan")
    det_rate = len(valid) / params["n_trials"]

    return {
        "scenario_id": scenario["id"],
        "label": scenario["label"],
        "f_source": f_source,
        "add": add,
        "far": far,
        "detection_rate": det_rate,
        "threshold": h_cal,
        "condition_number": detector._cond_F,
        "coeffs": detector._coeffs.tolist(),
    }


def run_experiment(n_trials=2000, n_calibration=2000, oracle_n=40000,
                   calib_runs=350, change_point=200, t_max=1200,
                   epsilon=0.02, seed=42, results_dir=None):
    """Run the F-misspecification ablation and write results.

    Args:
        n_trials: MC trials for ADD/FAR measurement.
        n_calibration: Calibration sample size for H0/H1.
        oracle_n: Sample size for the oracle covariance.
        calib_runs: H0 streams for FAR threshold calibration.
        change_point: True change-point tau (and FAR window length).
        t_max: Total stream length per trial.
        epsilon: Target FAR.
        seed: Random seed.
        results_dir: Output directory.

    Returns:
        List of per-(scenario, F-source) result dicts.
    """
    params = {
        "n_trials": n_trials, "n_calibration": n_calibration,
        "oracle_n": oracle_n, "calib_runs": calib_runs,
        "change_point": change_point, "t_max": t_max,
        "epsilon": epsilon, "seed": seed,
    }

    print("=" * 72)
    print("Experiment 08: F-MISSPECIFICATION ABLATION (FAR equalized before ADD)")
    print("=" * 72)
    print(f"trials={n_trials} calib={n_calibration} oracle_n={oracle_n} "
          f"calib_runs={calib_runs} eps={epsilon} tau={change_point}")
    print()

    rows: List[Dict[str, Any]] = []
    for scenario in SCENARIOS:
        print(f"[{scenario['label']}]")
        for f_source in F_SOURCES:
            row = run_scenario(scenario, f_source, params)
            rows.append(row)
            add_str = f"{row['add']:.2f}" if np.isfinite(row["add"]) else "inf"
            print(f"  {F_SOURCE_LABELS[f_source]:<18} ADD={add_str:>8}  "
                  f"FAR={row['far']:.4f}  thr={row['threshold']:.3f}  "
                  f"cond(F)={row['condition_number']:.2e}")
        # Degradation summary vs oracle.
        oracle_add = next(r["add"] for r in rows
                          if r["scenario_id"] == scenario["id"] and r["f_source"] == "oracle")
        for f_source in ("mixture", "h0_only"):
            r = next(r for r in rows
                     if r["scenario_id"] == scenario["id"] and r["f_source"] == f_source)
            if np.isfinite(r["add"]) and np.isfinite(oracle_add) and oracle_add > 0:
                deg = 100.0 * (r["add"] - oracle_add) / oracle_add
                print(f"    -> {F_SOURCE_LABELS[f_source]} degradation vs oracle: {deg:+.1f}%")
        print()

    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp08_f_misspecification",
            "params": params,
            "results": rows,
        }
        path = os.path.join(results_dir, "exp08_f_misspecification.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"Results saved to {path}")

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exp08: F-misspecification ablation")
    parser.add_argument("--quick", action="store_true", help="Fast settings")
    parser.add_argument("--n_trials", type=int, default=2000)
    parser.add_argument("--n_calibration", type=int, default=2000)
    parser.add_argument("--oracle_n", type=int, default=40000)
    parser.add_argument("--calib_runs", type=int, default=350)
    parser.add_argument("--change_point", type=int, default=200)
    parser.add_argument("--t_max", type=int, default=1200)
    parser.add_argument("--epsilon", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_trials = 1000
        args.n_calibration = 1500
        args.oracle_n = 30000
        args.calib_runs = 250

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp08_{ts}")

    run_experiment(
        n_trials=args.n_trials, n_calibration=args.n_calibration,
        oracle_n=args.oracle_n, calib_runs=args.calib_runs,
        change_point=args.change_point, t_max=args.t_max,
        epsilon=args.epsilon, seed=args.seed, results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
