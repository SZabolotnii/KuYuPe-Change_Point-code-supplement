# Benchmark study on measured streams

**This study is not part of the article.** It was withdrawn from the manuscript on
2026-09-05, after a scoring defect was found in the harness that produced it. What it
establishes is a *scope* result — where the detector holds its false-alarm level and where it
does not see a change at all — not a superiority claim. It is kept here in full, with the defect,
the corrected numbers and the code that recomputes them, because the earlier public version of the
paper (`arXiv:2605.23419`) reports the uncorrected numbers and a reader deserves to be able to
check them.

Everything below is reproducible offline from the stored result files:

```bash
python experiments/real_data/score.py            # the original six corpora, both scoring rules
python experiments/real_data/score.py --set refit  # after the standardisation fix
```

---

## 1. The defect

The detector reported `detection_time = len(test_segment)` when no alarm ever fired, and the
harness scored any stop at or after the change point as a detection. A run in which the statistic
never crossed its threshold was therefore recorded as a **detection at the last sample of the
segment**, with a delay equal to the distance from the change point to the end of the segment.

The baseline wrappers report `None` in that situation and were scored correctly. The defect is
confined to the GSA rows.

The protocol the study itself states is unambiguous:

> A trial ends in a false alarm if the detector stops before the change point; in a detection if
> it stops after it; or in a missed target if no alarm is raised before the end of the test
> segment.

The arithmetic tell: the published ADD on NSL-KDD is 91.9, and the mean distance from the change
point to the end of the segment over those 30 trials is 91.93. Every "detection" was the
end-of-segment fallback. The same signature appears as `DR = 1 - FAR` **exactly**, in every GSA
row of the published table.

## 2. The corrected numbers

Re-scored under the protocol as stated (`score.py`, one representative GSA configuration per
corpus):

| Corpus | trials | published rule (ADD, FAR, DR) | protocol rule (ADD, FAR, DR) | silent runs |
|---|---:|---|---|---:|
| US RealInt | 3 | (20.3, 0.000, **1.00**) | (18.0, 0.000, **0.67**) | 1 |
| SKAB | 128 | (347.6, 0.000, **1.00**) | (261.7, 0.000, **0.02**) | 125 |
| NASA IMS | 15 | (160.2, 0.000, **1.00**) | (61.6, 0.000, **0.60**) | 6 |
| TCPD | 89 | (152.5, 0.202, **0.80**) | (29.1, 0.202, **0.11**) | 61 |
| NAB EC2 | 49 | (636.5, 0.102, **0.90**) | (473.9, 0.102, **0.65**) | 12 |
| NSL-KDD | 30 | (91.9, 0.000, **1.00**) | (∞, 0.000, **0.00**) | 30 |

No baseline in any corpus has a single silent run. The false-alarm rates are unaffected: a false
alarm was determined by a stop *before* the change point, which the defect does not touch.

## 3. Why the detector was silent

Two separate causes, both now measured.

**A flattened dictionary.** `evaluate_basis_*` clips every basis value to `[-phi_max, phi_max]`
with `phi_max = 10`. TCPD series are supplied unnormalised — `brent_spot` lives on [16.9, 36.0],
`businv` on [8.0·10⁵, 1.1·10⁶] — so every basis value is clipped to the same constant 10. The
dictionary then carries no information, the orthonormalisation is singular, the approximated
log-likelihood ratio is a constant, and the CUSUM statistic can never leave zero: the detector is
**inert at every threshold**, not conservative. This flattens the basis on **53 of 89** TCPD
trials, and on 0 of 49 NAB and 0 of 128 SKAB trials, both of which are normalised upstream.

`GSADetector` now takes `standardize="robust"` (median and 1.4826·MAD, fitted on the calibration
sample alone), warns when the clip flattens the dictionary, and records
`diagnostics.basis_degenerate`. The default remains `None`, so the Monte-Carlo study of the paper
reproduces unchanged.

**A change too small to move the statistic.** Standardising removes the inertness but does not
make the detector good. Comparing the peak of the statistic on the test segment against its own
threshold, after the fix:

| Corpus | median `max g / h` | best trial | trials crossing |
|---|---:|---:|---:|
| NSL-KDD | 0.251 | 0.287 | **0 of 30** |
| SKAB | 0.282 | 2.22 | **1 of 128** |
| NAB EC2 | 1.25 | 438 | 26 of 49 |
| TCPD | 4.82 | 1757 | 70 of 89 |

On NSL-KDD the statistic reaches a quarter of its threshold and never crosses it in any trial. No
threshold rescues that; the change simply does not drive the statistic far enough over the test
segment. NSL-KDD and SKAB are exactly the corpora whose change is a change of **shape at a nearly
unchanged mean** on a short calibration segment (NSL-KDD: mean 3.22, standard deviation 0.147),
which is the empirical face of the parity theorem — and the measurement says the surrogate is
blind there too, not only the linear detectors.

## 4. What the study does and does not show

After the standardisation fix (`results/refit/`, quantile threshold calibration with 2000
in-control runs):

| Corpus | (ADD, FAR, DR) |
|---|---|
| TCPD | (76.1, 0.506, 0.28) |
| NAB EC2 | (521.0, 0.102, 0.43) |
| SKAB | (236.0, 0.000, 0.01) |
| NSL-KDD | (∞, 0.000, 0.00) |

**Shows.** The detector holds its false-alarm level where the change is large relative to the
in-control excursion of the statistic (NAB EC2). Its per-sample cost is `O(s)`.

**Does not show.** Any superiority over the retrospective and Bayesian baselines. On TCPD the
kernel change-point detector reaches DR 0.58 against 0.28 here; on NSL-KDD only the sign and MAD
CUSUM variants detect anything at all, and they do so at false-alarm rates above 0.66.

**Open.** Holding the level on TCPD needs a threshold multiplier of eight to sixteen. That is not
a result and no such constant is recommended anywhere: it is a *measurement* of how badly an
i.i.d. bootstrap of a short, non-stationary calibration segment understates the real excursions of
the statistic. A moving-block bootstrap that preserves short-range dependence is the natural next
step, and it is a change of assumption rather than a tuning choice.

**Not investigated.** On SKAB the `GSA-frac-S2` configuration calibrates to a threshold of
2.03·10¹⁸ — a numerical blow-up — where the other five configurations calibrate to about 16.3.

## 5. What is in this directory

```
score.py                       recompute every number above from the stored files
results/published/*.json       the six corpora as originally run (six GSA configurations
                               plus sign/MAD CUSUM, EWMA, BOCPD, kernel CPD and PELT)
results/rescored/*.json        both scoring rules applied to those files
results/refit/*.json           four corpora re-run with standardize="robust" and
                               quantile threshold calibration, carrying per-trial
                               diagnostics: alarm_raised, max_g_stat, basis_degenerate,
                               the calibrated threshold and the standardising scale
```

The corpus loaders are **not** bundled: NASA IMS, NSL-KDD, SKAB, TCPD, PhysioNet 2019, FTSE 100
and US RealInt each need a separate download, several behind registration. Sources, versions and
citations are in [`../../docs/DATASETS.md`](../../docs/DATASETS.md). The `gsa_cpd` package is
dataset-agnostic: once a univariate series is a NumPy array, the `GSADetector` API from the
top-level README applies directly — and on a series that is not already on a unit scale, pass
`standardize="robust"`.

## 6. Reproducible parts of the paper

- Monte-Carlo study (Section 7 of the article) — [`../monte_carlo/`](../monte_carlo/)
- Formal proofs — [`../../Lean/`](../../Lean/)
- The `gsa_cpd` package and its test suite — [`../../tests/`](../../tests/)
