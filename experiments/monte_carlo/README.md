# Monte Carlo Experiments

Reproducible Monte Carlo experiments for validating the GSA change-point detection framework.

## Prerequisites

```bash
pip install -r requirements.txt
```

## Running all experiments

From the repository root:

```bash
# Quick run (reduced trial counts, ~5 min)
python -m experiments.monte_carlo.run_all --quick

# Full run (~1-2 hours)
python -m experiments.monte_carlo.run_all
```

## Individual experiments

Each experiment can be run independently:

```bash
# Exp01: Gaussian limit verification (Theorem 1)
python -m experiments.monte_carlo.exp01_gaussian_limit --quick

# Exp02: ADD vs approximation degree (Pearson III)
python -m experiments.monte_carlo.exp02_add_vs_degree --quick

# Exp03: Basis function comparison (Student-t, Pearson III)
python -m experiments.monte_carlo.exp03_basis_comparison --quick

# Exp04: FAR control verification (all configurations)
python -m experiments.monte_carlo.exp04_far_control --quick

# Exp05: Stopping rule comparison (CUSUM/GRSh/SRP)
python -m experiments.monte_carlo.exp05_stopping_rules --quick
```

## Experiment descriptions

| Experiment | What it tests | Key parameters |
|---|---|---|
| exp01 | S=1 poly GSA matches classical CUSUM on Gaussian data | N_cal=1000, 100 trials |
| exp02 | ADD decreases with degree s=1..4 on Pearson III | gamma_3={0,2,10}, 500 trials |
| exp03 | Poly vs frac vs log bases at s=2 | Student-t(5), Pearson III(10) |
| exp04 | Empirical FAR <= epsilon across all configs | 1000 H0-only trials |
| exp05 | CUSUM vs GRSh vs SRP with GSA-LLR | ARL-matched thresholds |

## Output

Results are saved as JSON files in timestamped directories under `results/`.
Each experiment also prints a formatted summary table to stdout.

## Customization

All scripts accept `--help` for available options. Common flags:

- `--quick` -- fewer trials for fast iteration
- `--n_trials N` -- override trial count
- `--seed N` -- set random seed (default: 42)
- `--results_dir DIR` -- override output directory
