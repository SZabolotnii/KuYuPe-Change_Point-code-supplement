# GSA Change-Point Detection

Companion code for the paper:

> **"Generalized Stochastic Approximation of the Log-Likelihood Ratio for Robust Sequential Change-Point Detection"**
>
> Serhii Zabolotnii

## Highlights

- **Moment-based LLR approximation** -- no PDF knowledge required, only moments up to order 2s
- **Three basis types:** polynomial, fractional, logarithmic (plus Hermite)
- **Kunchenko PE-criterion** for per-step FAR bounding via Chebyshev / Vysochansky-Petunin / Cantelli inequalities
- **No shortest-delay claim.** On the benchmark corpora the detector is never the fastest, and on
  a change of shape at an unchanged mean it does not fire at all; see
  [`experiments/real_data/README.md`](experiments/real_data/README.md)
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

### Monte Carlo Simulations (Section 7 of the article)

```bash
python -m experiments.monte_carlo.run_all --quick    # fast (~5 min)
python -m experiments.monte_carlo.run_all            # full (~2 hours)

# Or a single experiment, e.g. the Gaussian-limit check (Theorem 1):
python -m experiments.monte_carlo.exp01_gaussian_limit --quick
```

Results are written to timestamped directories under `results/`.

### Benchmark study on measured streams

**This study is not part of the article.** It was withdrawn from the manuscript
on 2026-09-05 after a scoring defect was found in the harness that produced it:
a run in which the detector never raised an alarm was recorded as a detection at
the last sample of the test segment. Baselines were unaffected, so the defect
inflated the GSA rows alone — on NSL-KDD the detection rate falls from 1.00 to
**0.00** once the study's own protocol is applied.

The corpora themselves (NASA IMS Bearing, NSL-KDD, SKAB, TCPD, FRED macro
series, PhysioNet 2019, US RealInt) are large, and several need registration, so
the loaders are still not bundled. The **result files are**, together with the
script that re-scores them, so every corrected number is checkable offline:

```bash
python experiments/real_data/score.py             # both scoring rules, six corpora
python experiments/real_data/score.py --set refit # after the standardisation fix
```

What the study establishes is a scope result — where the detector holds its
false-alarm level and where it cannot see the change at all — not a superiority
claim. Read [`experiments/real_data/README.md`](experiments/real_data/README.md)
before citing any number from it; the earlier public version of the paper
(`arXiv:2605.23419`) reports the uncorrected values.

### Formal Proofs

```bash
cd Lean && lake build GSA
```

Key files: `InfoFunctional.lean` (Parseval partial sums over an **orthonormal** basis — see the
warning in its header: it does NOT certify Theorem 2 of the preprint), `BridgeGap.lean` (the
finite-dimensional bridge `Yᵀ F⁻¹ Y = max_K (KᵀY)²/(KᵀFK)`, plus the explicit `NOT FORMALISED`
inventory), `Convergence.lean` (Theorem 4), `FAR_ADD.lean` (Theorem 6).

> **Correction 2026-08-23.** Theorem 2(a) and 2(c) of the preprint are false; only part (b),
> monotonicity, holds. The complete-basis limit is `(1/c)·2Δ/(2−Δ)` in the triangular
> discrimination, not the Jeffreys divergence. Do not describe Theorem 2 as "Lean-verified".

### Running Tests

```bash
pytest tests/ -v
```

## Key Results

| Scenario | Result | Reproduced here |
|---|---|---|
| Gaussian limit | S=1 poly = classical CUSUM (exact match, validated) | Yes (MC + test) |
| Non-Gaussian (gamma_3 >= 8) | 30--36% ADD reduction vs. classical CUSUM | Yes (MC) |
| Measured streams | Scope result only: the level is held where the change is large relative to the statistic's in-control excursion (NAB EC2, FAR 0.10 at DR 0.43); a change of shape at an unchanged mean on a short calibration segment is not seen at all (NSL-KDD, SKAB) | Yes, from the stored result files |

> **Correction 2026-09-05.** Two rows previously claimed "only working method
> (NASA IMS Bearing)" and "FAR = 0%, DetRate = 100% (NSL-KDD)". Both were
> artefacts of the scoring defect described in
> [`experiments/real_data/README.md`](experiments/real_data/README.md): a run in
> which no alarm fired was counted as a detection. Under the protocol the study
> states, the NSL-KDD detection rate is **0.00**, not 1.00, and the detector
> raises no alarm on any of the 30 trials. Do not cite the old figures.

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
    monte_carlo/            # Section 7 of the article (run via -m)
    real_data/              # benchmark study: result files + re-scoring script
                            #   (not part of the article; corpora external)
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
