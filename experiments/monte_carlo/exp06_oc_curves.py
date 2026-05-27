"""Experiment 06: Operating-characteristic (OC) curves (paper section 4.6).

Compares GSA-CUSUM against the classical Gaussian CUSUM and two recent
distribution-free baselines -- the Hyvarinen-score SCUSUM (Wu et al., AISTATS
2023) and the univariate QT-EWMA core of KQT-EWMA (Nogara Notarianni et al.,
2024) -- via their OC curves: average detection delay (ADD) versus the
in-control average run length ARL_0.

Methodology (repo invariant: equalize ARL_0 / FAR before comparing ADD):
  For every (method, distribution, target ARL_0):
    1. Calibrate the detector's threshold by bisection on i.i.d. H0 streams so
       the empirical ARL_0 matches the target within tolerance.
    2. With that threshold fixed, measure ADD on H1 streams. The change occurs
       at the start of the H1 stream (t=0), so ADD is the pure post-change
       detection delay at matched ARL_0 -- the quantity an OC curve isolates.

Compute plan (reduced grid -- the full 4x3x4 grid at 1e5 runs is infeasible
locally; see --full):
  DEFAULT : 3 distributions x 2 target ARL_0 x {GSA-CUSUM, classical CUSUM,
            SCUSUM, KQT-EWMA}, modest N so the whole run finishes quickly.
  --quick : even smaller (smoke test).
  --full  : the paper-grade grid (4 distributions x 3 ARL_0 x 4 methods, large
            N). Intended for a cluster.

FULL-GRID REPRODUCTION COMMAND (paper grade; run on a cluster):
  python exp06_oc_curves.py --full \\
      --arl0 200,500,1000 \\
      --dists normal,student_t5,laplace,pearson3 \\
      --methods gsa,cusum,scusum,kqt \\
      --n_calib 100000 --calib_stream_mult 5 --n_add 100000

Outputs (under --results_dir, default results/exp06_<ts>/):
  results.csv          raw per-cell rows
  fig_oc_curves.pdf    the OC figure (requires matplotlib; skipped if absent)
"""

import argparse
import csv
import os
import time
from datetime import datetime
from typing import Callable, Dict, Any, List

import numpy as np

from gsa_cpd import GSADetector, BasisType, ThresholdType
from gsa_cpd.baselines import SCUSUM, KQTEWMA
from gsa_cpd.utils.distributions import create_distribution


# --------------------------------------------------------------------------- #
# Scenarios                                                                    #
# --------------------------------------------------------------------------- #
# Each scenario: a pre-change (H0) and post-change (H1) distribution, with a
# pure mean shift of magnitude `mean_shift` (held common so detectors face the
# same change). poly basis is used for the GSA arm so F is well conditioned.
# Each shift is ~1 standard deviation of the base distribution so the change is
# comparably detectable across distributions (a fair OC comparison); the
# non-Gaussian *shape* then drives the method ranking, not the shift magnitude.
# SD(Student-t(5)) = sqrt(5/3) ~ 1.29; SD(Laplace, scale=1) = sqrt(2) ~ 1.41.
SCENARIOS = {
    "normal": {
        "label": "Normal",
        "dist_type": "normal",
        "H0": {"loc": 0.0, "scale": 1.0},
        "H1": {"loc": 1.0, "scale": 1.0},
        "mean_shift": 1.0,
    },
    "student_t5": {
        "label": "Student-t(5)",
        "dist_type": "student_t",
        "H0": {"df": 5, "loc": 0.0, "scale": 1.0},
        "H1": {"df": 5, "loc": 2.0, "scale": 1.0},
        "mean_shift": 2.0,
    },
    "laplace": {
        "label": "Laplace",
        "dist_type": "laplace",
        "H0": {"loc": 0.0, "scale": 1.0},
        "H1": {"loc": 1.4, "scale": 1.0},
        "mean_shift": 1.4,
    },
    "pearson3": {
        "label": "Pearson III (skew 1.0)",
        "dist_type": "pearson3",
        "H0": {"skew": 1.0, "loc": 0.0, "scale": 1.0},
        "H1": {"skew": 1.0, "loc": 1.0, "scale": 1.0},
        "mean_shift": 1.0,
    },
    # Pure SHAPE-change scenarios (constant mean & variance): a mean-shift
    # detector is near-blind; GSA (tuned to the actual H1 moments via h1_data)
    # and the distribution-free QT-EWMA still detect. `mean_shift` is only the
    # nominal design delta for the mean-based baselines.
    "skew_shape": {
        "label": "Skewness change (skew 0->1.5)",
        "dist_type": "pearson3",
        "H0": {"skew": 0.0, "loc": 0.0, "scale": 1.0},
        "H1": {"skew": 1.5, "loc": 0.0, "scale": 1.0},
        "mean_shift": 1.0,
        "kind": "shape",
        "gsa": {"degree": 3, "basis": "poly", "phi_max": 50.0},
    },
    "kurt_shape": {
        "label": "Kurtosis change (kurt 0->6)",
        "dist_type": "student_t",
        "H0": {"df": 200, "loc": 0.0, "scale": 0.99499},
        "H1": {"df": 5, "loc": 0.0, "scale": 0.774597},
        "mean_shift": 1.0,
        "kind": "shape",
        "gsa": {"degree": 2, "basis": "log", "phi_max": 50.0},
    },
}

METHOD_LABELS = {
    "gsa": "GSA-CUSUM",
    "cusum": "Classical CUSUM",
    "scusum": "SCUSUM",
    "kqt": "QT-EWMA",
}


# --------------------------------------------------------------------------- #
# Classical Gaussian mean-shift CUSUM                                          #
# --------------------------------------------------------------------------- #
class ClassicalCUSUM:
    """Textbook one-sided Gaussian CUSUM for a known-magnitude mean shift.

    The mean mu0 and std sigma0 are estimated from H0 calibration data; the
    shift delta is the assumed post-change mean change. The log-likelihood-ratio
    increment for a Gaussian mean shift is

        s_t = (delta / sigma0^2) * (x_t - mu0) - delta^2 / (2 sigma0^2),

    and the CUSUM recursion is g_t = max(0, g_{t-1} + s_t), alarm iff g_t >= h.
    This is the classical Page CUSUM (not oracle: it plugs in sample mu0/sigma0
    and an assumed delta rather than the true densities).
    """

    def __init__(self, delta: float = 0.6):
        self.delta = delta
        self._mu0 = 0.0
        self._sigma0 = 1.0
        self._threshold = np.inf
        self._g = 0.0

    def fit(self, calibration_data: np.ndarray) -> "ClassicalCUSUM":
        self._mu0 = float(np.mean(calibration_data))
        self._sigma0 = float(np.std(calibration_data)) or 1.0
        self.reset()
        return self

    def increment(self, x: np.ndarray) -> np.ndarray:
        return (self.delta / self._sigma0**2) * (x - self._mu0) \
            - self.delta**2 / (2 * self._sigma0**2)

    def reset(self) -> None:
        self._g = 0.0

    def predict(self, x: float) -> bool:
        s = float(self.increment(np.array([x]))[0])
        self._g = max(0.0, self._g + s)
        return self._g >= self._threshold


# --------------------------------------------------------------------------- #
# Uniform calibration / evaluation                                             #
# --------------------------------------------------------------------------- #
def _gsa_set_threshold(detector, h: float) -> None:
    """Set the scalar threshold on any supported detector type."""
    detector._threshold = h


def calibrate_to_arl0(
    detector,
    h0_sampler: Callable[[int], np.ndarray],
    target_arl0: float,
    stream_len: int,
    n_runs: int,
    tol: float = 0.05,
    max_iter: int = 30,
) -> float:
    """Bisect the detector threshold so empirical ARL_0 ~= target.

    Works for any detector exposing reset()/predict(x)->bool and a `_threshold`
    attribute. Run length is censored at `stream_len`. ARL_0 increases with the
    threshold, giving a clean monotone bisection. The SCUSUM and KQTEWMA classes
    have their own (vectorized) calibrate(); this is the shared fallback used
    for GSA-CUSUM and classical CUSUM.
    """
    streams = [h0_sampler(stream_len) for _ in range(n_runs)]

    def empirical_arl0(h: float) -> float:
        rls = []
        for stream in streams:
            detector.reset()
            _gsa_set_threshold(detector, h)
            rl = stream_len
            for t, x in enumerate(stream):
                if detector.predict(float(x)):
                    rl = t + 1
                    break
            rls.append(rl)
        return float(np.mean(rls))

    h_low, h_high = 0.0, 1.0
    for _ in range(60):
        if empirical_arl0(h_high) >= target_arl0:
            break
        h_high *= 1.6

    best_h, best_gap = h_high, float("inf")
    for _ in range(max_iter):
        h_mid = 0.5 * (h_low + h_high)
        arl = empirical_arl0(h_mid)
        gap = abs(arl - target_arl0) / target_arl0
        if gap < best_gap:
            best_gap, best_h = gap, h_mid
        if arl < target_arl0:
            h_low = h_mid
        else:
            h_high = h_mid
        if gap < tol:
            break
    _gsa_set_threshold(detector, best_h)
    return best_h


def measure_arl0(detector, h0_sampler, stream_len, n_runs) -> float:
    """Empirical ARL_0 on fresh H0 streams (verification at fixed threshold)."""
    rls = []
    for _ in range(n_runs):
        detector.reset()
        rl = stream_len
        for t, x in enumerate(h0_sampler(stream_len)):
            if detector.predict(float(x)):
                rl = t + 1
                break
        rls.append(rl)
    return float(np.mean(rls))


def measure_add(detector, h1_sampler, horizon, n_runs) -> Dict[str, float]:
    """Average detection delay on post-change (H1) streams (change at t=0).

    Returns the mean detection delay over runs that detect within `horizon`, the
    detection rate, and the censored-mean delay (non-detections counted as
    `horizon`). Streams start at the change point, so the delay is the index of
    the first alarm.
    """
    delays = []
    detected = 0
    censored = []
    for _ in range(n_runs):
        detector.reset()
        d = horizon
        alarmed = False
        for t, x in enumerate(h1_sampler(horizon)):
            if detector.predict(float(x)):
                d = t
                alarmed = True
                break
        censored.append(d)
        if alarmed:
            delays.append(d)
            detected += 1
    add = float(np.mean(delays)) if delays else float("inf")
    return {
        "add": add,
        "detection_rate": detected / n_runs,
        "add_censored": float(np.mean(censored)),
    }


# --------------------------------------------------------------------------- #
# Per-cell driver                                                              #
# --------------------------------------------------------------------------- #
def build_detector(method: str, scenario: Dict[str, Any], calib_H0: np.ndarray,
                   epsilon: float, calib_H1=None):
    """Construct and fit a detector for the given method on calibration data."""
    if method == "gsa":
        gcfg = scenario.get("gsa", {})
        det = GSADetector(
            basis=gcfg.get("basis", BasisType.POLY), degree=gcfg.get("degree", 2),
            epsilon=epsilon, threshold_type=ThresholdType.CHEBYSHEV,
            threshold_scale=1.0, phi_max=gcfg.get("phi_max", 10.0),
        )
        if scenario.get("kind") == "shape" and calib_H1 is not None:
            # General reference anomaly: tune to the actual post-change moments.
            det.fit(calib_H0, h1_data=calib_H1)
        else:
            det.fit(calib_H0, delta=0.3)
        return det
    if method == "cusum":
        det = ClassicalCUSUM(delta=scenario["mean_shift"])
        det.fit(calib_H0)
        return det
    if method == "scusum":
        # Pre-/post-change std from calibration data; the score difference uses
        # a Gaussian working model (a valid score-matching surrogate even when
        # the data are non-Gaussian -- the detector remains distribution-robust
        # because lambda and the threshold are calibrated on the actual H0 data).
        mu0 = float(np.mean(calib_H0))
        sig0 = float(np.std(calib_H0)) or 1.0
        det = SCUSUM(
            pre_params={"mu": mu0, "sigma": sig0},
            post_params={"mu": mu0 + scenario["mean_shift"], "sigma": sig0},
        )
        det.fit(calib_H0)
        return det
    if method == "kqt":
        det = KQTEWMA(n_bins=8, lam=0.05)
        det.fit(calib_H0)
        return det
    raise ValueError(f"Unknown method: {method}")


def run_cell(method: str, scenario: Dict[str, Any], target_arl0: float,
             params: Dict[str, Any], seed: int) -> Dict[str, Any]:
    """Calibrate one (method, scenario, ARL_0) detector and measure ADD."""
    dist_type = scenario["dist_type"]
    dist_h0 = create_distribution(dist_type, scenario["H0"])
    dist_h1 = create_distribution(dist_type, scenario["H1"])

    def make_sampler(dist, base_seed):
        counter = {"i": 0}

        def sampler(n):
            rs = np.random.default_rng(base_seed + counter["i"])
            counter["i"] += 1
            return dist.rvs(n, random_state=rs)

        return sampler

    h0_sampler = make_sampler(dist_h0, seed + 10_000)
    h1_sampler = make_sampler(dist_h1, seed + 20_000)

    calib_H0 = dist_h0.rvs(params["n_fit"], random_state=np.random.default_rng(seed))
    calib_H1 = dist_h1.rvs(params["n_fit"], random_state=np.random.default_rng(seed + 1))

    detector = build_detector(method, scenario, calib_H0, params["epsilon"], calib_H1)

    calib_stream_len = int(params["calib_stream_mult"] * target_arl0)
    add_horizon = int(params["add_horizon_mult"] * target_arl0)

    # Threshold calibration to the target ARL_0.
    if method in ("scusum", "kqt"):
        detector.calibrate(h0_sampler, target_arl0, calib_stream_len,
                            n_runs=params["n_calib"], tol=params["tol"])
        thr = detector._threshold
    else:
        thr = calibrate_to_arl0(detector, h0_sampler, target_arl0,
                                calib_stream_len, n_runs=params["n_calib"],
                                tol=params["tol"])

    # Verify achieved ARL_0 on independent streams.
    if method in ("scusum", "kqt"):
        achieved_arl0 = _measure_arl0_fast(detector, method, h0_sampler,
                                            calib_stream_len, params["n_verify"])
    else:
        achieved_arl0 = measure_arl0(detector, h0_sampler, calib_stream_len,
                                     params["n_verify"])

    # Measure ADD at the calibrated threshold.
    add_stats = measure_add(detector, h1_sampler, add_horizon, params["n_add"])

    return {
        "method": method,
        "method_label": METHOD_LABELS[method],
        "scenario": next(k for k, v in SCENARIOS.items() if v is scenario),
        "scenario_label": scenario["label"],
        "target_arl0": target_arl0,
        "achieved_arl0": achieved_arl0,
        "threshold": thr,
        "add": add_stats["add"],
        "add_censored": add_stats["add_censored"],
        "detection_rate": add_stats["detection_rate"],
    }


def _measure_arl0_fast(detector, method, h0_sampler, stream_len, n_runs) -> float:
    """Vectorized ARL_0 verification for SCUSUM / KQTEWMA."""
    rls = []
    if method == "scusum":
        for _ in range(n_runs):
            stream = h0_sampler(stream_len)
            z = detector.lam * detector.score_difference(stream)
            g, rl = 0.0, stream_len
            for t, zt in enumerate(z):
                g = max(0.0, g + zt)
                if g >= detector._threshold:
                    rl = t + 1
                    break
            rls.append(rl)
    else:  # kqt
        for _ in range(n_runs):
            stream = h0_sampler(stream_len)
            rls.append(detector._run_length(stream, detector._threshold, stream_len))
    return float(np.mean(rls))


# --------------------------------------------------------------------------- #
# Figure                                                                       #
# --------------------------------------------------------------------------- #
def make_figure(rows: List[Dict[str, Any]], dists: List[str], methods: List[str],
                out_path: str) -> None:
    """OC plot: ADD vs target ARL_0, one line per method, one panel per dist."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping figure.")
        return

    colors = {"gsa": "#1f77b4", "cusum": "#ff7f0e",
              "scusum": "#2ca02c", "kqt": "#d62728"}
    markers = {"gsa": "o", "cusum": "s", "scusum": "^", "kqt": "D"}

    n = len(dists)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 3.8), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, dist in zip(axes, dists):
        for method in methods:
            cells = sorted(
                [r for r in rows if r["scenario"] == dist and r["method"] == method],
                key=lambda r: r["target_arl0"],
            )
            if not cells:
                continue
            xs = [c["achieved_arl0"] for c in cells]
            ys = [c["add"] for c in cells]
            ax.plot(xs, ys, marker=markers[method], color=colors[method],
                    label=METHOD_LABELS[method], linewidth=1.8, markersize=6)
        ax.set_xscale("log")
        ax.set_xlabel(r"In-control ARL$_0$ (achieved)")
        ax.set_title(SCENARIOS[dist]["label"])
        ax.grid(True, which="both", alpha=0.3)
    axes[0].set_ylabel("Average detection delay (ADD)")
    axes[0].legend(loc="upper left", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Experiment driver                                                            #
# --------------------------------------------------------------------------- #
def run_experiment(dists=None, methods=None, arl0s=None, epsilon=0.02,
                   n_fit=5000, n_calib=400, n_verify=400, n_add=2000,
                   calib_stream_mult=4.0, add_horizon_mult=8.0, tol=0.05,
                   seed=42, results_dir=None):
    """Run the OC-curve comparison and write results.csv + the OC figure.

    Args:
        dists: List of scenario keys (default: normal, student_t5, laplace).
        methods: List of method keys (default: gsa, cusum, scusum, kqt).
        arl0s: List of target ARL_0 values (default: [500, 1000]).
        epsilon: Target FAR for the GSA threshold (analytic seed value).
        n_fit: Calibration sample size.
        n_calib: H0 streams for threshold bisection.
        n_verify: H0 streams for ARL_0 verification.
        n_add: H1 streams for ADD.
        calib_stream_mult: H0 stream length = mult * target ARL_0.
        add_horizon_mult: H1 horizon = mult * target ARL_0.
        tol: Relative ARL_0 tolerance.
        seed: Base random seed.
        results_dir: Output directory.

    Returns:
        List of per-cell result dicts.
    """
    if dists is None:
        dists = ["normal", "student_t5", "laplace"]
    if methods is None:
        methods = ["gsa", "cusum", "scusum", "kqt"]
    if arl0s is None:
        arl0s = [500.0, 1000.0]

    for d in dists:
        if d not in SCENARIOS:
            raise ValueError(f"Unknown distribution '{d}'. Choose from {list(SCENARIOS)}.")
    for m in methods:
        if m not in METHOD_LABELS:
            raise ValueError(f"Unknown method '{m}'. Choose from {list(METHOD_LABELS)}.")

    params = {
        "epsilon": epsilon, "n_fit": n_fit, "n_calib": n_calib,
        "n_verify": n_verify, "n_add": n_add,
        "calib_stream_mult": calib_stream_mult,
        "add_horizon_mult": add_horizon_mult, "tol": tol,
    }

    print("=" * 78)
    print("Experiment 06: OPERATING-CHARACTERISTIC CURVES (ADD vs ARL_0)")
    print("=" * 78)
    print(f"dists={dists}  methods={methods}  ARL_0={arl0s}")
    print(f"n_fit={n_fit} n_calib={n_calib} n_verify={n_verify} "
          f"n_add={n_add} eps={epsilon}")
    print(f"cells = {len(dists) * len(methods) * len(arl0s)}")
    print()

    t0 = time.time()
    rows: List[Dict[str, Any]] = []
    cell_idx = 0
    total_cells = len(dists) * len(methods) * len(arl0s)
    for dist in dists:
        scenario = SCENARIOS[dist]
        print(f"[{scenario['label']}]")
        for arl0 in arl0s:
            for method in methods:
                cell_idx += 1
                cell_seed = seed + 1000 * cell_idx
                row = run_cell(method, scenario, arl0, params, cell_seed)
                rows.append(row)
                add_str = f"{row['add']:7.2f}" if np.isfinite(row["add"]) else "    inf"
                print(f"  ARL0={arl0:>6.0f} {METHOD_LABELS[method]:<16} "
                      f"ADD={add_str}  ARL0_ach={row['achieved_arl0']:8.1f}  "
                      f"det={row['detection_rate']:.3f}  thr={row['threshold']:.4f}  "
                      f"[{cell_idx}/{total_cells}]")
        print()

    wall = time.time() - t0
    print(f"Total wall time: {wall:.1f}s ({wall/60:.2f} min)")

    if results_dir:
        os.makedirs(results_dir, exist_ok=True)

        fig_path = os.path.join(results_dir, "fig_oc_curves.pdf")
        make_figure(rows, dists, methods, fig_path)
        if os.path.exists(fig_path):
            print(f"Wrote figure: {fig_path}")

        csv_path = os.path.join(results_dir, "results.csv")
        fieldnames = ["method", "method_label", "scenario", "scenario_label",
                      "target_arl0", "achieved_arl0", "threshold", "add",
                      "add_censored", "detection_rate"]
        with open(csv_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for r in rows:
                w.writerow({k: r[k] for k in fieldnames})
        print(f"Wrote raw results: {csv_path}")

    return rows


# --------------------------------------------------------------------------- #
# Main                                                                         #
# --------------------------------------------------------------------------- #
def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--quick", action="store_true", help="Even smaller smoke grid.")
    p.add_argument("--full", action="store_true",
                   help="Paper-grade grid (4 dist x 3 ARL_0 x 4 methods, large N).")
    p.add_argument("--dists", type=str, default=None,
                   help="Comma list from: " + ",".join(SCENARIOS))
    p.add_argument("--methods", type=str, default=None,
                   help="Comma list from: gsa,cusum,scusum,kqt")
    p.add_argument("--arl0", type=str, default=None,
                   help="Comma list of target ARL_0 values.")
    p.add_argument("--epsilon", type=float, default=0.02)
    p.add_argument("--n_fit", type=int, default=5000, help="Calibration sample size.")
    p.add_argument("--n_calib", type=int, default=400,
                   help="H0 streams for threshold bisection.")
    p.add_argument("--n_verify", type=int, default=400,
                   help="H0 streams for ARL_0 verification.")
    p.add_argument("--n_add", type=int, default=2000, help="H1 streams for ADD.")
    p.add_argument("--calib_stream_mult", type=float, default=4.0,
                   help="H0 stream length = mult * target ARL_0.")
    p.add_argument("--add_horizon_mult", type=float, default=8.0,
                   help="H1 horizon = mult * target ARL_0.")
    p.add_argument("--tol", type=float, default=0.05, help="Relative ARL_0 tolerance.")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--results_dir", type=str, default=None)
    args = p.parse_args()

    # Grid defaults (reduced) unless overridden.
    if args.full:
        dists = ["normal", "student_t5", "laplace", "pearson3"]
        methods = ["gsa", "cusum", "scusum", "kqt"]
        arl0s = [200.0, 500.0, 1000.0]
        args.n_calib = max(args.n_calib, 100000)
        args.n_verify = max(args.n_verify, 50000)
        args.n_add = max(args.n_add, 100000)
        args.n_fit = max(args.n_fit, 50000)
        args.calib_stream_mult = max(args.calib_stream_mult, 5.0)
    else:
        dists = ["normal", "student_t5", "laplace"]
        methods = ["gsa", "cusum", "scusum", "kqt"]
        arl0s = [500.0, 1000.0]

    if args.quick:
        dists = ["normal", "student_t5"]
        methods = ["gsa", "cusum", "scusum", "kqt"]
        arl0s = [300.0, 600.0]
        args.n_calib = 150
        args.n_verify = 200
        args.n_add = 600
        args.n_fit = 3000

    # CLI overrides.
    if args.dists:
        dists = [d.strip() for d in args.dists.split(",")]
    if args.methods:
        methods = [m.strip() for m in args.methods.split(",")]
    if args.arl0:
        arl0s = [float(a) for a in args.arl0.split(",")]

    if args.results_dir is None:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.results_dir = os.path.join("results", f"exp06_{ts}")

    run_experiment(
        dists=dists, methods=methods, arl0s=arl0s, epsilon=args.epsilon,
        n_fit=args.n_fit, n_calib=args.n_calib, n_verify=args.n_verify,
        n_add=args.n_add, calib_stream_mult=args.calib_stream_mult,
        add_horizon_mult=args.add_horizon_mult, tol=args.tol,
        seed=args.seed, results_dir=args.results_dir,
    )


if __name__ == "__main__":
    main()
