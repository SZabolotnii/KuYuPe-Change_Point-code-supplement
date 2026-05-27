"""Experiment 07: Drift Lemma verification (Lemma 2.bis) via Monte Carlo.

Verifies the drift conditions for the GSA approximate-LLR statistic

    Lambda^(s)(x) = K^T (phi(x) - mu_comb),

with mu_comb = (m_0 + m_1) / 2 (average basis mean under H0/H1) and
K = F^{-1} Y the solution of the normal system F K = Y. The lemma states

    E_0[Lambda^(s)] = -1/2 J(s) < 0,
    E_1[Lambda^(s)] = +1/2 J(s) > 0,
    J(s) = Y^T F^{-1} Y > 0   (positive-definite F).

Design note on F, Y, K
----------------------
F (basis covariance averaged over H0/H1), Y (basis-mean difference) and K are
built from a large reference (``oracle``) sample drawn from the SciPy density and
evaluated with the package basis functions (gsa_cpd.core.basis). This is the
empirical analogue of the analytical basis expectations; for large oracle_n it
converges to the exact moments for every distribution used here. The drift
identity holds by construction whenever (mu_comb, F, Y, K) are mutually
consistent; the Monte Carlo step then validates that the basis moments used to
build them match an *independent* sampled draw of the distribution.

Usage:
    python exp07_drift_verification.py --quick
    python exp07_drift_verification.py
"""

import argparse
import json
import os
from datetime import datetime
from typing import Any, Dict, List

import numpy as np

from gsa_cpd.core.basis import BasisType, evaluate_basis_matrix
from gsa_cpd.utils.distributions import create_distribution


# Configurations: (label, dist_type, H0, H1, degrees).
CONFIGS: List[Dict[str, Any]] = [
    {
        "label": "Normal mean shift (loc: 0.0 -> 0.5)",
        "dist_type": "normal",
        "H0": {"loc": 0.0, "scale": 1.0},
        "H1": {"loc": 0.5, "scale": 1.0},
        "degrees": [1, 2, 3],
    },
    {
        "label": "Student-t(df=5) mean shift (loc: 0.0 -> 0.6)",
        "dist_type": "student_t",
        "H0": {"df": 5, "loc": 0.0, "scale": 1.0},
        "H1": {"df": 5, "loc": 0.6, "scale": 1.0},
        # df=5 -> only moments of order < 5 exist; F at degree s needs E[X^{2s}],
        # so the polynomial basis is restricted to s <= 2 here.
        "degrees": [1, 2],
    },
    {
        "label": "Pearson III loc shift (loc: 0.0 -> 0.4, skew=1.0)",
        "dist_type": "pearson3",
        "H0": {"skew": 1.0, "loc": 0.0, "scale": 1.0},
        "H1": {"skew": 1.0, "loc": 0.4, "scale": 1.0},
        "degrees": [1, 2, 3],
    },
]

BASIS = BasisType.POLY
# Large clip so the polynomial basis is effectively unclipped at the reference
# moments (the drift identity is about the raw polynomial moments).
PHI_MAX = 1e12


def _basis_moments(dist, n: int, degree: int, rng) -> Dict[str, np.ndarray]:
    """Mean vector E[phi] and covariance Cov(phi) from a large sample."""
    data = dist.rvs(n, random_state=rng)
    B = evaluate_basis_matrix(data, degree, BASIS, PHI_MAX)
    E = np.mean(B, axis=0)
    if degree == 1:
        F = np.array([[np.var(B)]])
    else:
        F = np.cov(B, rowvar=False)
    return {"E": E, "F": F}


def build_drift_system(dist_h0, dist_h1, degree: int, oracle_n: int,
                       seed: int) -> Dict[str, Any]:
    """Build F, Y, K, J(s) and mu_comb from large reference samples.

    F = (F_H0 + F_H1) / 2,  Y = E_1[phi] - E_0[phi],  K = F^{-1} Y,
    J(s) = Y^T F^{-1} Y,  mu_comb = (E_0[phi] + E_1[phi]) / 2.

    Returns a dict with keys F, Y, K, J, mu_comb, cond_F.
    """
    res_H0 = _basis_moments(dist_h0, oracle_n, degree,
                            np.random.default_rng(seed + 11))
    res_H1 = _basis_moments(dist_h1, oracle_n, degree,
                            np.random.default_rng(seed + 22))

    E0, E1 = res_H0["E"], res_H1["E"]
    F = 0.5 * (res_H0["F"] + res_H1["F"])
    Y = E1 - E0

    K = np.linalg.solve(F, Y)
    J = float(Y @ K)  # Y^T F^{-1} Y == K^T Y
    mu_comb = 0.5 * (E0 + E1)
    cond_F = float(np.linalg.cond(F))

    return {"F": F, "Y": Y, "K": K, "J": J, "mu_comb": mu_comb, "cond_F": cond_F}


def evaluate_lambda(samples: np.ndarray, K: np.ndarray,
                    mu_comb: np.ndarray) -> np.ndarray:
    """Evaluate Lambda^(s)(x) = K^T (phi(x) - mu_comb) for each sample."""
    degree = len(K)
    phi = evaluate_basis_matrix(samples, degree, BASIS, PHI_MAX)
    return phi @ K - float(K @ mu_comb)


def run_one(dist_type: str, H0: Dict[str, Any], H1: Dict[str, Any],
            degree: int, n_samples: int, oracle_n: int,
            seed: int) -> Dict[str, Any]:
    """Build the drift system and Monte Carlo estimate E_0/E_1 of Lambda."""
    dist_h0 = create_distribution(dist_type, H0)
    dist_h1 = create_distribution(dist_type, H1)

    sys_dict = build_drift_system(dist_h0, dist_h1, degree, oracle_n, seed)
    K, J, mu_comb = sys_dict["K"], sys_dict["J"], sys_dict["mu_comb"]

    x0 = dist_h0.rvs(n_samples, random_state=np.random.default_rng(seed + 101))
    x1 = dist_h1.rvs(n_samples, random_state=np.random.default_rng(seed + 202))

    L0 = evaluate_lambda(x0, K, mu_comb)
    L1 = evaluate_lambda(x1, K, mu_comb)

    e0_emp = float(np.mean(L0))
    e1_emp = float(np.mean(L1))
    se0 = float(np.std(L0, ddof=1) / np.sqrt(n_samples))
    se1 = float(np.std(L1, ddof=1) / np.sqrt(n_samples))

    e0_th = -0.5 * J
    e1_th = +0.5 * J

    # Relative errors against +/- J/2 (J > 0 expected).
    denom = abs(0.5 * J) + 1e-6
    rel0 = abs(e0_emp - e0_th) / denom
    rel1 = abs(e1_emp - e1_th) / denom

    return {
        "dist_type": dist_type,
        "degree": degree,
        "J": J,
        "cond_F": sys_dict["cond_F"],
        "E0_emp": e0_emp,
        "E0_th": e0_th,
        "E0_se": se0,
        "E0_rel_err": rel0,
        "E1_emp": e1_emp,
        "E1_th": e1_th,
        "E1_se": se1,
        "E1_rel_err": rel1,
    }


def _within_tolerance(row: Dict[str, Any], rel_tol: float, se_mult: float) -> bool:
    """A row passes if the drift sign is correct, J > 0, and each empirical
    mean is within rel_tol of +/- J/2 OR within se_mult standard errors."""
    if not (row["J"] > 0):
        return False
    ok0 = (row["E0_rel_err"] < rel_tol) or (
        abs(row["E0_emp"] - row["E0_th"]) < se_mult * row["E0_se"]
    )
    ok1 = (row["E1_rel_err"] < rel_tol) or (
        abs(row["E1_emp"] - row["E1_th"]) < se_mult * row["E1_se"]
    )
    sign_ok = (row["E0_emp"] < 0) and (row["E1_emp"] > 0)
    return ok0 and ok1 and sign_ok


def run_experiment(n_samples=100_000, oracle_n=200_000, rel_tol=0.05,
                   se_mult=3.0, seed=42, results_dir=None):
    """Run the Drift Lemma Monte Carlo verification.

    Args:
        n_samples: Monte Carlo sample size per hypothesis (verification draw).
        oracle_n: Reference sample size used to build F, Y, K.
        rel_tol: Relative tolerance on E[Lambda] vs +/- J/2.
        se_mult: Standard-error multiplier for the alternative pass criterion.
        seed: Random seed.
        results_dir: Output directory.

    Returns:
        List of per-configuration result dicts.
    """
    print("=" * 88)
    print("Experiment 07: DRIFT LEMMA VERIFICATION (Lemma 2.bis) -- Monte Carlo")
    print("=" * 88)
    print("Lambda^(s)(x) = K^T (phi(x) - mu_comb),  K = F^{-1} Y,  J(s) = Y^T F^{-1} Y")
    print("Expected:  E_0[Lambda] = -J/2 < 0,   E_1[Lambda] = +J/2 > 0,   J > 0")
    print(f"Basis: poly | N per hypothesis: {n_samples} | oracle_n: {oracle_n} "
          f"| seed: {seed}")
    print(f"Tolerance: rel < {rel_tol} OR within {se_mult} SE")
    print()

    header = (
        f"{'distribution':<40}{'s':>3}{'J(s)':>10}{'E0_emp':>11}{'E0=-J/2':>11}"
        f"{'E1_emp':>11}{'E1=+J/2':>11}{'rel0':>8}{'rel1':>8}{'pass':>6}"
    )
    print(header)
    print("-" * len(header))

    rows: List[Dict[str, Any]] = []
    all_pass = True
    for cfg in CONFIGS:
        for s in cfg["degrees"]:
            row = run_one(cfg["dist_type"], cfg["H0"], cfg["H1"], s,
                          n_samples, oracle_n, seed)
            ok = _within_tolerance(row, rel_tol, se_mult)
            all_pass = all_pass and ok
            rows.append({**row, "label": cfg["label"], "passed": ok})
            print(
                f"{cfg['label']:<40}{s:>3}{row['J']:>10.4f}"
                f"{row['E0_emp']:>11.5f}{row['E0_th']:>11.5f}"
                f"{row['E1_emp']:>11.5f}{row['E1_th']:>11.5f}"
                f"{row['E0_rel_err']:>8.3f}{row['E1_rel_err']:>8.3f}"
                f"{'OK' if ok else 'FAIL':>6}"
            )

    n_pass = sum(1 for r in rows if r["passed"])
    print()
    print(f"Passed {n_pass}/{len(rows)} configurations within tolerance.")
    print(f"All J(s) > 0: {all(r['J'] > 0 for r in rows)}")

    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
        out = {
            "experiment": "exp07_drift_verification",
            "params": {"n_samples": n_samples, "oracle_n": oracle_n,
                       "rel_tol": rel_tol, "se_mult": se_mult, "seed": seed},
            "n_pass": n_pass,
            "n_total": len(rows),
            "results": rows,
        }
        path = os.path.join(results_dir, "exp07_drift_verification.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nResults saved to {path}")

    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Exp07: Monte Carlo verification of the Drift Lemma.")
    parser.add_argument("--quick", action="store_true",
                        help="Quick mode: smaller sample sizes.")
    parser.add_argument("--n_samples", type=int, default=100_000,
                        help="MC sample size per hypothesis (verification).")
    parser.add_argument("--oracle_n", type=int, default=200_000,
                        help="Reference sample size used to build F, Y, K.")
    parser.add_argument("--rel_tol", type=float, default=0.05,
                        help="Relative tolerance on E[Lambda] vs +/- J/2.")
    parser.add_argument("--se_mult", type=float, default=3.0,
                        help="Standard-error multiplier for the pass criterion.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--results_dir", type=str, default=None)
    args = parser.parse_args()

    if args.quick:
        args.n_samples = 10_000
        args.oracle_n = 40_000

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp07_{ts}")

    rows = run_experiment(
        n_samples=args.n_samples, oracle_n=args.oracle_n,
        rel_tol=args.rel_tol, se_mult=args.se_mult, seed=args.seed,
        results_dir=args.results_dir,
    )
    all_pass = all(r["passed"] for r in rows)
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
