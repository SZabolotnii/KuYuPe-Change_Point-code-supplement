# Changelog

All notable changes to the public release of `gsa_cpd` (the public mirror of the
research repository `KuYuPe-Change_Point`) are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

For research-side history (papers, internal experiments, Lean development),
see `../VERSION.md` and `../reports/lean_theory_verification_report.md`.

## [Unreleased]

Targeted for v1.0.0 alongside the arxiv preprint v1 (Ukrainian).

### Added (Sequential Analysis revision, 2026-05-25)
- SCUSUM baseline (`src/gsa_cpd/baselines/scusum.py`) — Hyvärinen-score
  CUSUM (Wu et al., AISTATS 2023) with MGF-constrained `lambda` selection
  and Monte-Carlo ARL₀ calibration.
- KQT-EWMA baseline (`src/gsa_cpd/baselines/kqt_ewma.py`) — distribution-free
  QT-EWMA core (univariate special case of Nogara Notarianni et al., 2024).
- Operating-characteristic experiment
  `experiments/monte_carlo/exp06_oc_curves.py` (ADD vs ARL₀ at matched ARL₀:
  GSA-CUSUM, classical CUSUM, SCUSUM, KQT-EWMA), reproducing the revised
  paper §4.6.
- Drift-lemma and F-misspecification studies
  `exp07_drift_verification.py`, `exp08_f_misspecification.py`; new
  experiments registered in `run_all.py`.
- Baseline unit tests `tests/test_baselines.py` (13 tests).
- Lean: local-Fisher / χ²-projection interpretive doc-comment on `J(s)` in
  `Lean/GSA/Part2/InfoFunctional.lean` (paper §2.7 / C2 reformulation);
  identifiers unchanged, `lake build GSA` remains green.

### Changed (claim audit + shape-change OC, 2026-05-25)
- `GSADetector.fit(..., h1_data=...)` — fit the reference anomaly to the actual
  post-change moments, enabling detection of pure shape changes (skewness /
  kurtosis at constant mean and variance), not just mean shifts.
- `exp06_oc_curves.py` — added `skew_shape` (γ₃: 0→1.5) and `kurt_shape`
  (γ₄: 0→6) scenarios; GSA detects ~2–5× faster than mean-based CUSUM/SCUSUM
  and faster than the QT-EWMA core, whereas on mean shifts all moment/score
  detectors coincide.
- Docs: replaced "guaranteed FAR control" with the accurate "per-step FAR bound
  (proven at s=1) + Monte Carlo calibration for s≥2" (README, `docs/API.md`);
  `J(s)` described as a Hilbert-space quantity `Yᵀ F⁻¹ Y`, the maximum of the
  Rayleigh quotient over the dictionary span.
  **CORRECTION 2026-08-23:** earlier releases described `J(s)` as converging to the
  Jeffreys divergence. That is FALSE. The complete-basis limit is
  `(1/c)·2Δ/(2−Δ)`, where `Δ` is the triangular (Vincze–Le Cam) discrimination and
  `c` is the constant in `F = c(C₀+C₁)`; it stands in no fixed order with respect to
  the Jeffreys divergence. Theorem 2(a) `J(s) ≤ J_Jeffreys` is false as well; only
  part (b), monotonicity, survives. Exact counterexample: Bernoulli(0.2, 0.8) gives
  `J_∞ = 2.25` against a Jeffreys divergence of `1.2 ln 4 ≈ 1.664`. What the Lean
  development proves is a correct abstract Parseval/projection theorem; the bridge to
  the statistical identification was never formalised — see
  `Lean/GSA/Part2/BridgeGap.lean` for the true bridge and for the explicit
  `NOT FORMALISED` inventory.

### Added (Etap 2 — Q1-Q2 closure, 2026-05-04)
- Bootstrap CI utilities (`src/benchmarks/bootstrap_ci.py`,
  `src/benchmarks/compute_tier1_ci.py`) — percentile method, n_boot=1000.
- Wilcoxon signed-rank tests with Holm-Bonferroni for GSA vs every
  baseline (`src/benchmarks/wilcoxon_test.py`); per-trial paired delays.
- Ablation study `src/experiments/exp8_ablation.py` — sweeps
  `winsor_pct`, `max_abs_phi`, `threshold_scale` from baseline; results
  feed into Part5 §5.7.4 + Fig 8.
- Bayesian Online CPD (Adams-MacKay 2007) with Student-t likelihood
  and kernel-CPD (Harchaoui et al.) baselines — `make_bocpd()` and
  `make_kernel_cpd()` factories in `runner.py`. PELT remained from
  prior phase. Re-ran every Tier 1 dataset with the augmented detector
  set and updated `paper/shared/results_manifest.json` accordingly.
- Computational complexity analysis
  (`src/benchmarks/compute_complexity.py`) — per-call wall-time and
  throughput; results in Part5 §5.5.4 (chesna picture: GSA ~4400×
  slower than Sign-CUSUM, ~30-60× slower than BOCPD/PELT — trade-off
  for PE-FAR-control).
- GitHub Actions workflows in `.github/workflows/test.yml`:
  pytest matrix on Python 3.10/3.11/3.12 plus a Lean job that builds
  `GSA` and asserts zero `sorry` in committed sources.

### Planned
- Pinned dependency versions in `requirements.txt`.
- pyproject.toml `[dev]` extras to include `bayesian-changepoint-detection`
  for full reproduction of Etap 2.5 baselines.

## [1.0.0-rc1] — 2026-04-08

First public release candidate. Mirrors research repository at the
`v1.5.2` calibration milestone.

### Added
- Core detector module `src/gsa_cpd/core/` (detector, basis, solver,
  threshold, diagnostics, moments) extracted and cleaned up from the
  research-side `src/experiments/core/gsa_detector_v2.py`.
- Stopping rules `src/gsa_cpd/stopping_rules/` (CUSUM, GRSh, SRP).
- Baseline detectors `src/gsa_cpd/baselines/` (sign-CUSUM, MAD-CUSUM,
  EWMA, oracle CUSUM).
- Data utilities `src/gsa_cpd/data/` (schema, preprocessing).
- Reproducible experiments `experiments/monte_carlo/exp01`–`exp05`,
  with `run_all.py` orchestrator.
- Real-data benchmark scripts `experiments/real_data/`.
- Test suite `tests/` (`test_detector.py`, `test_solver.py`,
  `test_threshold.py`, `test_stopping_rules.py`, `test_gaussian_limit.py`).
- Lean 4 formalization mirrored to `Lean/GSA/Part2/` (Lean 4.26.0 +
  mathlib 4.26, 0 `sorry`, 1 `axiom`).
- `pyproject.toml` (PEP 621), MIT `LICENSE`, bilingual `README.md` /
  `README_UA.md`.
- `CITATION.cff` with author ORCID.

### Known limitations
- `src/gsa_cpd/data/feature_extractors.py` and `loaders.py` not yet
  migrated from the research tree.
- No CI workflow; tests run locally only.
- Theoretical-moments support limited to polynomial basis.
- Lean Theorems 4c (rate O(s^{-2r})) and 5 (criterion Yu) listed as
  TODO in `appendix_lean.tex`.
