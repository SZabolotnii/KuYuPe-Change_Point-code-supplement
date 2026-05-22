# Real-Data Benchmark Experiments (Section 5)

The real-data results reported in Section 5 of the paper were produced on
large external datasets (NASA IMS Bearing, NSL-KDD, SKAB, TCPD, FRED macro
series, PhysioNet 2019, FTSE 100, US RealInt). Several of these require
registration or per-dataset manual download.

**These benchmark scripts are not bundled in this supplement.** This
repository is scoped to the parts of the paper that are fully self-contained
and immediately reproducible without external data:

- Monte Carlo study (Section 4) — see [`../monte_carlo/`](../monte_carlo/)
- Formal proofs (Section 2) — see [`../../Lean/`](../../Lean/)
- The `gsa_cpd` package and its test suite — see [`../../tests/`](../../tests/)

## Where to find the real-data details

- **Dataset sources, versions, and citations:** [`../../docs/DATASETS.md`](../../docs/DATASETS.md)
- **Reported results and discussion:** Section 5 of the paper

The `gsa_cpd` package itself is dataset-agnostic: once a univariate series is
loaded into a NumPy array, the same `GSADetector` API used in the Monte Carlo
experiments applies directly. See the Quick Start in the top-level README.
