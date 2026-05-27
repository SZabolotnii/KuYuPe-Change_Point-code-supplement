# API Reference

Brief reference for the public classes in `gsa_cpd`.

## Core

### GSADetector

```python
from gsa_cpd import GSADetector
```

Main detector class. Approximates the log-likelihood ratio using basis functions and moment-based optimization, then applies CUSUM stopping rule.

**Constructor:**

```python
GSADetector(
    basis: BasisType = BasisType.POLY,
    degree: int = 2,
    epsilon: float = 0.01,
    threshold_type: ThresholdType = ThresholdType.CHEBYSHEV,
    threshold_scale: float = 1.0,
    ridge_lambda: float = 1e-6,
    phi_max: float = 10.0,
    winsor_pct: float = 0.05,
)
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `basis` | `BasisType` | `POLY` | Basis function family |
| `degree` | `int` | `2` | Approximation order s (number of basis functions) |
| `epsilon` | `float` | `0.01` | Target false alarm rate |
| `threshold_type` | `ThresholdType` | `CHEBYSHEV` | Threshold computation method |
| `threshold_scale` | `float` | `1.0` | Safety multiplier for threshold |
| `ridge_lambda` | `float` | `1e-6` | Tikhonov regularization for FK=Y solver |
| `phi_max` | `float` | `10.0` | Clip basis function values to [-phi_max, phi_max] |
| `winsor_pct` | `float` | `0.05` | Winsorization percentage for each tail (0 to disable) |

**Methods:**

| Method | Signature | Description |
|---|---|---|
| `fit` | `fit(calibration_data: np.ndarray, delta: float = 0.2) -> GSADetector` | Calibrate on H0 data. `delta` is the MDE coefficient (expected relative change in moments). Returns self for chaining. |
| `predict` | `predict(x: float) -> bool` | Process one observation. Returns `True` if change detected. |
| `reset` | `reset() -> None` | Reset runtime state for a new monitoring session. Does not clear fitted parameters. |

**Properties:**

| Property | Type | Description |
|---|---|---|
| `diagnostics` | `GSADiagnostics` | Diagnostic information from calibration |
| `coefficients` | `np.ndarray` or `None` | Coefficient vector K from FK=Y |
| `alarm_time` | `int` or `None` | Time index of the alarm |

---

### BasisType

```python
from gsa_cpd import BasisType
```

Enum of available basis function families.

| Value | Basis Functions | Best For |
|---|---|---|
| `BasisType.POLY` | {x, x^2, ..., x^s} | General purpose, Gaussian limit (s=1) |
| `BasisType.LOG` | {x, ln\|x\|, x ln\|x\|, (ln\|x\|)^2, ...} | Log-scale variations |
| `BasisType.FRAC` | {sgn(x)\|x\|^a1, sgn(x)\|x\|^a2, ...} | Heavy-tailed data (kurtosis > 6) |
| `BasisType.HERMITE` | {He_1(x), He_2(x), ..., He_s(x)} | Gaussian or near-Gaussian data |

**Notes:**
- `POLY` with s=1 exactly reproduces classical Gaussian CUSUM (Theorem 1).
- `FRAC` uses predefined exponents [1.0, 0.5, 1/3, 2/3, 0.25, 0.75] and is generally the most robust choice for non-Gaussian data.
- `HERMITE` uses probabilist's Hermite polynomials orthogonal w.r.t. N(0,1).

---

### ThresholdType

```python
from gsa_cpd import ThresholdType
```

Enum of threshold computation methods. All analytic methods use the PE (Probability Error) criterion.

| Value | Formula | Assumption |
|---|---|---|
| `ThresholdType.CHEBYSHEV` | h = E + sigma / sqrt(epsilon) | Universal (any distribution) |
| `ThresholdType.VP` | h = E + (2/3) sigma / sqrt(epsilon) | Unimodal distributions |
| `ThresholdType.CANTELLI` | h = E + sigma sqrt(1/epsilon - 1) | One-sided bound |
| `ThresholdType.SIMULATION` | MC-calibrated via bootstrap | No assumptions; slowest |

**Notes:**
- `CHEBYSHEV` is the default and provides a per-step FAR bound for any distribution via Chebyshev's inequality (a proven per-step bound at `s=1`; for higher orders the in-control run length is set by Monte Carlo calibration).
- `VP` (Vysochansky-Petunin) gives a tighter threshold when the unimodal assumption holds, reducing ADD.
- `SIMULATION` uses bootstrap resampling of calibration data and binary search on ARL.

---

### GSADiagnostics

```python
from gsa_cpd import GSADiagnostics
```

Dataclass with diagnostic information collected during `GSADetector.fit()`.

| Field | Type | Description |
|---|---|---|
| `condition_number` | `float` | Condition number of system matrix F. Values > 1e6 indicate ill-conditioning. |
| `J_s` | `float` | Information functional J(s) = K^T Y (Kunchenko definition) |
| `J_s_lean` | `float` | J(s) = sum(K_i^2) (Lean/Parseval definition) |
| `E_L_H0` | `float` | Expected value of Lambda under H0 |
| `E_L_H1` | `float` | Expected value of Lambda under H1 |
| `Var_L_H0` | `float` | Variance of Lambda under H0 |
| `threshold` | `float` | Detection threshold h |
| `coeffs` | `np.ndarray` | Coefficient vector K |
| `k0` | `float` | Bias term |
| `eta` | `float` | Efficiency coefficient (E[L\|H1] - E[L\|H0]) / sqrt(Var[L\|H0]) |
| `solver_method` | `str` | Method used to solve FK=Y: "direct", "ridge", or "svd" |

---

## Stopping Rules

Modular stopping rules that can be used independently of `GSADetector`. Each takes LLR increments as input.

### CUSUMRule

```python
from gsa_cpd.stopping_rules import CUSUMRule
```

Cumulative Sum stopping rule (Lorden minimax).

g_t = max(0, g_{t-1} + lambda_t); alarm when g_t > threshold.

```python
rule = CUSUMRule(threshold=5.0)
rule.update(lambda_t)  # -> bool (True = alarm)
rule.statistic         # current g_t
rule.reset()
```

### GRShRule

```python
from gsa_cpd.stopping_rules import GRShRule
```

Additive Bayesian stopping rule (Girshick-Rubin-Shiryaev).

S_t = S_{t-1} + lambda_t (no reset to zero); alarm when S_t > threshold.

Bayesian-optimal for geometric prior on change time. More sensitive to drift but less robust to past false spikes.

```python
rule = GRShRule(threshold=10.0)
rule.update(lambda_t)  # -> bool
rule.statistic         # current S_t
rule.reset()
```

### SRPRule

```python
from gsa_cpd.stopping_rules import SRPRule
```

Shiryaev-Roberts Procedure. Quasi-minimax optimal.

R_t = (1 + R_{t-1}) * exp(lambda_t); alarm when R_t > H.

Includes a static method for threshold calibration via ARL matching:

```python
rule = SRPRule(threshold=100.0)
rule.update(lambda_t)  # -> bool
rule.statistic         # current R_t
rule.reset()

# Calibrate threshold via Monte Carlo ARL estimation
H = SRPRule.calibrate_arl(
    compute_llr=lambda x: detector._compute_llr(x),
    sample_h0=lambda n: np.random.normal(0, 1, n),
    target_arl=500,
    n_runs=200,
)
```

**`calibrate_arl` parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compute_llr` | `Callable[[float], float]` | required | LLR function |
| `sample_h0` | `Callable[[int], np.ndarray]` | required | H0 sampler |
| `target_arl` | `float` | required | Target Average Run Length under H0 |
| `max_iter` | `int` | `12` | Binary search iterations |
| `tol` | `float` | `0.05` | Relative tolerance for ARL match |
| `n_runs` | `int` | `200` | Number of MC runs |
| `max_run_length` | `int` | `5000` | Max observations per run |

---

## Baselines

Baseline detectors for comparison. All follow the same `fit` / `predict` / `reset` interface as `GSADetector`.

### SignCUSUM

```python
from gsa_cpd.baselines import SignCUSUM
```

Nonparametric sign-based CUSUM. Uses sign(x - median) as the test statistic. No distributional assumptions whatsoever.

```python
det = SignCUSUM(epsilon=0.01)
det.fit(calibration_data)
det.predict(x)       # -> bool
det.alarm_time       # int or None
```

### MADCUSUM

```python
from gsa_cpd.baselines import MADCUSUM
```

MAD-normalized robust CUSUM. Normalizes observations using Median Absolute Deviation, then applies CUSUM on squared z-scores. Robust to outliers in calibration data.

```python
det = MADCUSUM(epsilon=0.01)
det.fit(calibration_data)
det.predict(x)       # -> bool
det.alarm_time       # int or None
```

### EWMA

```python
from gsa_cpd.baselines import EWMA
```

Exponentially Weighted Moving Average control chart. Signals when the EWMA deviates from the calibration mean by more than L * sigma.

```python
det = EWMA(lam=0.1, L=3.5)
det.fit(calibration_data)
det.predict(x)       # -> bool
det.alarm_time       # int or None
```

| Parameter | Type | Default | Description |
|---|---|---|---|
| `lam` | `float` | `0.1` | Smoothing parameter (0 < lam <= 1) |
| `L` | `float` | `3.5` | Control limit multiplier (number of sigmas) |

### OracleCUSUM

```python
from gsa_cpd.baselines import OracleCUSUM
```

Oracle CUSUM using the true analytical LLR. Requires full knowledge of both f0 and f1. Used as an upper-bound baseline to measure GSA efficiency loss.

```python
from scipy.stats import norm

true_llr = lambda x: norm.logpdf(x, loc=1) - norm.logpdf(x, loc=0)
det = OracleCUSUM(true_llr_func=true_llr, epsilon=0.01)
det.fit(calibration_data)
det.predict(x)       # -> bool
det.alarm_time       # int or None
```

---

## Data Schema

### TimeSeriesData

```python
from gsa_cpd.data.schema import TimeSeriesData, DatasetInfo, AnnotatedChangepoint
```

Container for a time series with annotated change points.

```python
ts = TimeSeriesData(
    values=np.array([...]),
    changepoints=[
        AnnotatedChangepoint(index=100, label="mean shift", confidence=1.0),
        AnnotatedChangepoint(index=500, label="variance change", confidence=0.9),
    ],
    info=DatasetInfo(
        name="my_dataset",
        domain="finance",
        source="Yahoo Finance",
        n_observations=1000,
        kurtosis=5.2,
        skewness=0.8,
    ),
)
ts.n  # -> 1000
```
