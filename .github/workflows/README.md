# CI workflows

This directory holds GitHub Actions workflows for the public release of
`gsa_cpd`.

## `test.yml` — main CI

Runs on every push and pull request to `main`, plus `workflow_dispatch`
for manual triggers.

Two jobs:

1. **pytest** — Python test suite under `tests/` against Python 3.10,
   3.11, 3.12 (matrix). Installs the package via `pip install -e .[dev]`
   so that the dev-dependency block in `pyproject.toml` covers
   `pytest`, `hypothesis`, etc.

2. **lean** — Lean 4 + mathlib build of `Lean/GSA/`. Installs `elan`,
   caches the `.lake` build directory keyed on `lean-toolchain` +
   `lake-manifest.json`, then runs `lake build GSA`. A grep guard fails
   the job if any `sorry` is detected in committed `*.lean` sources
   (the project asserts 0 `sorry` / 1 documented `axiom`).

## Adding more checks

* **Type-checking**: add a `mypy --strict src/gsa_cpd/` step to the
  pytest job after dependencies install.
* **Coverage**: add `pip install pytest-cov` and run
  `pytest --cov=src/gsa_cpd --cov-report=xml`, then upload via
  `codecov/codecov-action`.
* **Reproducibility audit**: a third job that re-runs
  `python -m experiments.monte_carlo.exp01_gaussian_limit` with a small
  `--n_trials` and asserts the printed K₁ ≈ 0.5 ± 0.05 — guards against
  silent regressions in the GSA detector.
