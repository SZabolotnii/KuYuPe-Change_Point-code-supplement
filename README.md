# GSA Change-Point Detection

Companion code for the paper:

> **"Generalized Stochastic Approximation of Log-Likelihood Ratio for Robust Sequential Change-Point Detection"**
>
> Serhii Zabolotnii

## Highlights

- **Moment-based LLR approximation** -- no PDF knowledge required, only moments up to order 2s
- **Three basis types:** polynomial, fractional, logarithmic (plus Hermite)
- **Kunchenko PE-criterion** for guaranteed FAR control via Chebyshev / Vysochansky-Petunin / Cantelli inequalities
- **Works on heavy-tailed data** where classical methods fail (kurtosis > 20)
- **O(s) per sample** -- suitable for edge/embedded devices
- **Formal proofs in Lean 4** (Mathlib) for core theorems
- **Modular stopping rules:** CUSUM, GRSh (Bayesian), Shiryaev-Roberts

## Installation

```bash
# From source (recommended for experiments)
git clone https://github.com/SZabolotnii/KuYuPe-Change_Point-code-supplement.git
cd KuYuPe-Change_Point-code-supplement
pip install -e "."

# With experiment dependencies (matplotlib, pandas, ruptures)
pip install -e ".[experiments]"

# Everything (experiments + dev/testing)
pip install -e ".[all]"
```

**Requirements:** Python >= 3.9, NumPy >= 1.21, SciPy >= 1.7.

## Quick Start

```python
import numpy as np
from gsa_cpd import GSADetector, BasisType, ThresholdType

# 1. Calibrate on normal-regime data
calibration = np.random.normal(0, 1, size=2000)

detector = GSADetector(
    basis=BasisType.FRAC,       # fractional basis -- best for heavy tails
    degree=2,                    # approximation order s=2
    epsilon=0.01,                # target FAR
    threshold_type=ThresholdType.CHEBYSHEV,
)
detector.fit(calibration, delta=0.3)

# 2. Monitor a live stream
stream = np.concatenate([
    np.random.normal(0, 1, 500),   # normal regime
    np.random.normal(1, 1, 200),   # post-change regime (mean shift)
])

for t, x in enumerate(stream):
    if detector.predict(x):
        print(f"Change detected at t={detector.alarm_time}")
        break

# 3. Inspect diagnostics
diag = detector.diagnostics
print(f"Threshold: {diag.threshold:.3f}")
print(f"J(s): {diag.J_s:.4f}")
print(f"cond(F): {diag.condition_number:.1f}")
print(f"Solver: {diag.solver_method}")
```

## Reproducing Paper Results

> Run experiments as modules from the repository root (after `pip install -e .`)
> so that both the installed `gsa_cpd` package and the local `experiments`
> package resolve correctly.

### Monte Carlo Simulations (Section 4)

```bash
python -m experiments.monte_carlo.run_all --quick    # fast (~5 min)
python -m experiments.monte_carlo.run_all            # full (~2 hours)

# Or a single experiment, e.g. the Gaussian-limit check (Theorem 1):
python -m experiments.monte_carlo.exp01_gaussian_limit --quick
```

Results are written to timestamped directories under `results/`.

### Real-Data Benchmarks (Section 5)

The real-data results in Section 5 (NASA IMS Bearing, NSL-KDD, SKAB, TCPD,
FRED macro series, PhysioNet 2019, etc.) rely on large external datasets,
several of which require registration or per-dataset download. **Those
benchmark scripts are not bundled in this supplement.** This repository
focuses on the fully self-contained, immediately reproducible parts of the
paper: the Monte Carlo study (Section 4), the formal proofs (Section 2), and
the `gsa_cpd` package with its test suite. See the paper's Section 5 and
[docs/DATASETS.md](docs/DATASETS.md) for dataset sources and the reported
real-data results.

### Formal Proofs (Section 2)

```bash
cd Lean && lake build GSA
```

Key files: `InfoFunctional.lean` (Theorem 2), `Convergence.lean` (Theorem 4), `FAR_ADD.lean` (Theorem 6).

### Running Tests

```bash
pytest tests/ -v
```

## Key Results

| Scenario | GSA Advantage | Reproduced here |
|---|---|---|
| Gaussian limit | S=1 poly = classical CUSUM (exact match, validated) | Yes (MC + test) |
| Non-Gaussian (gamma_3 >= 8) | 30--36% ADD reduction vs. classical CUSUM | Yes (MC) |
| Heavy tails (kurtosis > 20) | Only working method (NASA IMS Bearing dataset) | Paper Section 5 |
| Cybersecurity (NSL-KDD) | FAR = 0%, DetRate = 100% | Paper Section 5 |

Rows marked "Paper Section 5" are real-data benchmarks reported in the paper;
the synthetic (Monte Carlo) and Gaussian-limit results are reproducible
directly from this repository (see above).

## Package Structure

```
src/gsa_cpd/
    __init__.py             # Public API: GSADetector, BasisType, ThresholdType
    core/
        detector.py         # GSADetector -- main detector class
        basis.py            # BasisType enum + basis evaluation
        threshold.py        # ThresholdType enum + PE-criterion thresholds
        solver.py           # FK=Y linear system solver (direct / ridge / SVD)
        moments.py          # Moment computation utilities
        diagnostics.py      # GSADiagnostics dataclass
    stopping_rules/
        cusum.py            # CUSUMRule -- minimax (Lorden)
        grsh.py             # GRShRule -- Bayesian (Girshick-Rubin-Shiryaev)
        srp.py              # SRPRule -- quasi-minimax (Shiryaev-Roberts)
    baselines/
        oracle_cusum.py     # OracleCUSUM -- upper bound (known LLR)
        sign_cusum.py       # SignCUSUM -- nonparametric sign-based
        mad_cusum.py        # MADCUSUM -- MAD-normalized robust
        ewma.py             # EWMA -- exponentially weighted moving average
    data/
        schema.py           # TimeSeriesData, DatasetInfo containers
        preprocessing.py    # Data loading and preprocessing
    utils/
        distributions.py    # Distribution factory and sampling
        metrics.py          # FAR, ADD, J(s) evaluation metrics
experiments/
    monte_carlo/            # Section 4: Monte Carlo simulations (run via -m)
    real_data/              # Section 5: notes only (datasets external; see paper)
results/                    # Generated outputs (created on first run)
tests/
    test_detector.py        # GSADetector unit tests
    test_solver.py          # FK=Y solver tests
    test_threshold.py       # Threshold computation tests
    test_stopping_rules.py  # Stopping rule tests
    test_gaussian_limit.py  # Gaussian limit property test (S=1 poly = CUSUM)
Lean/                       # Lean 4 formal proofs (Mathlib)
```

## Citation

```bibtex
@article{zabolotnii2025gsa,
  title   = {Generalized Stochastic Approximation of Log-Likelihood Ratio
             for Robust Sequential Change-Point Detection},
  author  = {Zabolotnii, Serhii},
  year    = {2025},
  note    = {Software available at https://github.com/SZabolotnii/KuYuPe-Change_Point-code-supplement}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
