# Reproducibility Guide

Step-by-step instructions for reproducing all results from the paper.

## 1. Environment Setup

### System Requirements

- Python 3.9, 3.10, 3.11, or 3.12
- 4 GB RAM minimum (8 GB recommended for full Monte Carlo runs)
- Lean 4 / elan (only for formal proofs)

### Python Environment

```bash
# Clone the repository
git clone https://github.com/SZabolotnii/KuYuPe-Change_Point-code-supplement.git
cd KuYuPe-Change_Point-code-supplement

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
# .venv\Scripts\activate    # Windows

# Install with all dependencies
pip install -e ".[all]"
```

### Verify Installation

```bash
pytest tests/ -v
```

All 5 test modules should pass:

| Test File | What It Verifies |
|---|---|
| `test_detector.py` | GSADetector fit/predict lifecycle, all basis types |
| `test_solver.py` | FK=Y solver (direct, ridge, SVD fallback) |
| `test_threshold.py` | Chebyshev, VP, Cantelli, simulation thresholds |
| `test_stopping_rules.py` | CUSUM, GRSh, SRP rule correctness |
| `test_gaussian_limit.py` | S=1 polynomial basis = classical CUSUM |

### Lean 4 Setup (Optional)

```bash
# Install elan (Lean version manager)
curl https://elan-init.github.io/elan/elan-init.sh -sSf | sh

# Build proofs
cd Lean
lake build GSA
```

## 2. Monte Carlo Simulations (Section 4)

> Run experiments as modules from the repository root (after `pip install -e .`),
> so both the installed `gsa_cpd` package and the local `experiments` package
> are importable.

### Quick Run (~5 minutes)

```bash
python -m experiments.monte_carlo.run_all --quick
```

This runs a reduced set of Monte Carlo trials (fewer samples, fewer parameter combinations) to verify the pipeline works.

### Full Run (~2 hours)

```bash
python -m experiments.monte_carlo.run_all
```

This reproduces all Monte Carlo experiments from Section 4 with the full sample sizes reported in the paper.

### Expected Results

Results are saved to `results/` in timestamped subdirectories.

**Experiment 1 -- Gaussian Limit (Theorem 1):**
- S=1 polynomial GSA should match classical CUSUM within numerical tolerance
- Expected relative error < 1e-6

**Experiment 2 -- Convergence (Theorem 4):**
- J(s) should increase monotonically with s (for well-conditioned systems)
- For s > 5 with polynomial basis, condition numbers may exceed 1e6

**Experiment 3 -- Information Functional (Theorem 2):**
- J(s) values for s = 1, 2, 3, 5, 7, 10 across distributions
- Upper bound: J(s) <= ||z||^2 (KL divergence)

**Experiment 4 -- FAR and ADD (Theorem 6):**
- FAR should be bounded by epsilon for all Chebyshev thresholds
- ADD should decrease as s increases (diminishing returns after s=3)

**Experiment 5 -- Oracle Comparison:**
- GSA-LLR vs. OracleCUSUM efficiency ratio
- Expected 85--95% efficiency for Gaussian, 60--80% for heavy-tailed

## 3. Real-Data Benchmarks (Section 5)

The real-data benchmarks in Section 5 use large external datasets, several of
which require registration or per-dataset manual download. **The benchmark
scripts for these datasets are not bundled in this supplement** — this
repository is scoped to the fully self-contained, immediately reproducible
parts of the paper (Monte Carlo, formal proofs, and the package/tests).

For each dataset, source URLs, versions, and citations are documented in
[docs/DATASETS.md](DATASETS.md). The expected GSA behavior reported in the
paper is summarized below for reference:

| Dataset | Reported GSA Behavior (paper Section 5) |
|---|---|
| US RealInt | Detects 3 structural breaks (matching Bai-Perron 2003) |
| NASA IMS Bearing | Detects bearing degradation; only method working at kurtosis > 20 |
| NSL-KDD | FAR = 0%, DetRate = 100% on test set |
| SKAB | Competitive with ruptures PELT on valve/pump anomalies |
| TCPD | Matches or exceeds median F1 across 42 time series |
| FTSE 100 | Detects volatility regime changes (2008 crisis, COVID-19) |
| FEDFUNDS | Detects interest rate regime shifts (Volcker era) |
| PhysioNet 2019 | Detects sepsis onset with competitive AUROC |

## 5. Troubleshooting

### Numerical Issues

**Problem:** `cond(F) > 1e6` warning during calibration.
**Solution:** This is expected for polynomial basis with s > 5. The solver automatically falls back to SVD. Use `ridge_lambda=1e-4` for additional stability, or switch to `BasisType.FRAC` which is better conditioned.

**Problem:** `J(s)` exceeds `||z||^2` (Theorem 2a violation).
**Solution:** This can occur when orthonormalization changes the J(s) definition. For verification experiments, compare both `J_s` (Kunchenko) and `J_s_lean` (Parseval) in diagnostics.

### Memory Issues

**Problem:** Out-of-memory on full Monte Carlo run.
**Solution:** Use `--quick` flag or reduce `n_samples` parameter. The full run generates ~100k samples per distribution per experiment.

### Dataset Download Issues

**Problem:** Download script fails for a specific dataset.
**Solution:** See [docs/DATASETS.md](DATASETS.md) for manual download URLs and instructions. Some datasets require registration.

### Lean Build Issues

**Problem:** `lake build` fails with version errors.
**Solution:** Ensure you have the correct Lean toolchain. Run `elan default leanprover/lean4:v4.x.0` (check `lean-toolchain` file for the exact version).

## 6. Platform Notes

- **macOS / Linux:** Fully tested. No special configuration needed.
- **Windows:** Works with standard Python. Use PowerShell or WSL for shell commands.
- **Docker:** A Dockerfile is planned for future releases.

## 7. Random Seeds

All experiments use fixed random seeds for reproducibility. The default seed is set in each experiment script. To reproduce exact figures from the paper, do not modify the seed values. Running with `--quick` uses the same seeds but fewer samples, so statistical estimates will differ slightly.
