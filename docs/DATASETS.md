# Benchmark Datasets

Descriptions and download instructions for the 8 benchmark datasets used in Section 5 of the paper. All datasets are publicly available.

> **Note.** The real-data benchmark scripts (including the
> `experiments/real_data/download_datasets.py` helper referenced in the
> commands below) are **not bundled in this supplement** — see the README.
> This document is provided so that the datasets, versions, and citations
> behind the paper's Section 5 results are fully traceable. Obtain each
> dataset directly from the **Source** / URL listed for it.

## Overview

| # | Dataset | Domain | N obs | Kurtosis | Change Types | Source |
|---|---|---|---|---|---|---|
| 1 | US RealInt | Economics | 228 | ~3.2 | Mean shift (structural breaks) | Bai & Perron, 2003 |
| 2 | NASA IMS Bearing | Engineering | ~20k | >20 | Variance + distribution shift | Qiu et al., 2006 |
| 3 | NSL-KDD | Cybersecurity | ~150k | varies | Mean + distribution shift | Tavallaee et al., 2009 |
| 4 | SKAB | Industrial IoT | ~35k | ~4--8 | Mean + variance shift | Skoltech, 2020 |
| 5 | TCPD | Multi-domain | varies | varies | Mixed (42 time series) | van den Burg & Williams, 2020 |
| 6 | FTSE 100 | Finance | ~10k | ~8--15 | Volatility regime change | Yahoo Finance |
| 7 | FEDFUNDS | Economics | ~780 | ~3.5 | Mean shift (policy regimes) | FRED, St. Louis Fed |
| 8 | PhysioNet 2019 | Medical | varies | varies | Distribution shift (sepsis) | Reyna et al., 2019 |

---

## 1. US Real Interest Rate (US RealInt)

**Description:** Quarterly U.S. ex-post real interest rate, 1961Q1--2017Q4. Contains 3 known structural breaks identified by Bai-Perron test methodology. A classical benchmark for offline change-point detection.

**Citation:**
> Bai, J., & Perron, P. (2003). Computation and analysis of multiple structural change models. *Journal of Applied Econometrics*, 18(1), 1--22.

**Download:**
```bash
# Automatic (via experiment script)
python experiments/real_data/download_datasets.py --dataset us_realint

# Manual
# Available from the Journal of Applied Econometrics data archive:
# https://onlinelibrary.wiley.com/journal/10991255
# Also bundled with the R package 'strucchange'
```

**Expected file:** `experiments/real_data/data/us_realint.csv`

---

## 2. NASA IMS Bearing Dataset

**Description:** Vibration signals from 4 bearings on a loaded shaft running to failure. Collected by the Center for Intelligent Maintenance Systems (IMS), University of Cincinnati. Features extreme heavy tails (kurtosis > 20) in the degradation phase, making it an ideal test for heavy-tailed robustness.

**Citation:**
> Qiu, H., Lee, J., Lin, J., & Yu, G. (2006). Wavelet filter-based weak signature detection method and its application on rolling element bearing prognostics. *Journal of Sound and Vibration*, 289(4--5), 1066--1090.

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset nasa_ims

# Manual
# NASA Prognostics Data Repository:
# https://www.nasa.gov/content/prognostics-center-of-excellence-data-set-repository
# Direct: https://ti.arc.nasa.gov/tech/dash/groups/pcoe/prognostic-data-repository/
# Download "IMS Bearing Data" (6 GB compressed)
```

**Expected files:** `experiments/real_data/data/nasa_ims/` (directory with test 1, 2, 3 subdirectories)

**Preprocessing notes:**
- Use channel 1 (bearing 1, accelerometer 1) from Test 2
- Compute RMS of each 20,480-point snapshot to get a scalar time series
- Known failure: outer race defect at approximately observation 984

---

## 3. NSL-KDD

**Description:** An improved version of KDD Cup 1999 dataset for network intrusion detection. Contains labeled normal and attack traffic records. Used to evaluate GSA as a one-class anomaly detector where calibration uses only normal traffic.

**Citation:**
> Tavallaee, M., Bagheri, E., Lu, W., & Ghorbani, A. A. (2009). A detailed analysis of the KDD CUP 99 data set. *Proceedings of the 2nd IEEE Symposium on Computational Intelligence for Security and Defense Applications (CISDA)*.

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset nsl_kdd

# Manual
# University of New Brunswick:
# https://www.unb.ca/cic/datasets/nsl.html
# Files needed: KDDTrain+.txt, KDDTest+.txt
```

**Expected files:** `experiments/real_data/data/nsl_kdd/KDDTrain+.txt`, `KDDTest+.txt`

**Preprocessing notes:**
- Extract numeric features (columns 1--34)
- Calibrate on "normal" class from training set
- Test on full test set (mixed normal + attack)

---

## 4. SKAB (Skoltech Anomaly Benchmark)

**Description:** Sensor data from a testbed for water circulation system with flow, pressure, and temperature sensors. Contains labeled anomalies from valve manipulation and pump failures. Published by Skolkovo Institute of Science and Technology.

**Citation:**
> Katser, I., & Kozitsin, V. (2020). Skoltech Anomaly Benchmark (SKAB). https://github.com/waico/SKAB

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset skab

# Manual
# GitHub repository:
# https://github.com/waico/SKAB
# Clone and copy the 'data' directory
```

**Expected files:** `experiments/real_data/data/skab/` (directory with CSV files per anomaly scenario)

**Preprocessing notes:**
- Use the `anomaly-free/` subdirectory for calibration
- Evaluate on individual scenario files (valve1, valve2, etc.)
- Primary signal: "Volume Flow RateRMS" or combined feature vector

---

## 5. TCPD (Turing Change Point Dataset)

**Description:** A curated collection of 42 annotated time series across multiple domains (climate, finance, biology, etc.) created for benchmarking change-point detection algorithms. Maintained by The Alan Turing Institute.

**Citation:**
> van den Burg, G. J. J., & Williams, C. K. I. (2020). An evaluation of change point detection algorithms. *arXiv preprint arXiv:2003.06222*.

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset tcpd

# Manual
# GitHub repository:
# https://github.com/alan-turing-institute/TCPD
# Requires: pip install tcpd
# Then: python -c "from tcpd import datasets; datasets.download()"
```

**Expected files:** `experiments/real_data/data/tcpd/` (directory with JSON files per series)

**Preprocessing notes:**
- Each series has expert-annotated change points with confidence scores
- We report median F1 score across all 42 series
- Default tolerance: 5% of series length for matching detected vs. true change points

---

## 6. FTSE 100

**Description:** Daily closing prices of the FTSE 100 index (London Stock Exchange). Used for detecting volatility regime changes, including the 2008 financial crisis and COVID-19 market crash.

**Citation:**
> Data sourced from Yahoo Finance. FTSE 100 ticker: ^FTSE.

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset ftse100

# Manual (Python)
# pip install yfinance
# import yfinance as yf
# data = yf.download("^FTSE", start="2000-01-01", end="2024-12-31")
# data["Close"].to_csv("experiments/real_data/data/ftse100.csv")
```

**Expected file:** `experiments/real_data/data/ftse100.csv`

**Preprocessing notes:**
- Compute log-returns: r_t = ln(P_t / P_{t-1})
- Use squared log-returns for volatility change detection
- Known regime changes: 2008-09 (financial crisis), 2016-06 (Brexit vote), 2020-03 (COVID-19)

---

## 7. FEDFUNDS (Federal Funds Rate)

**Description:** Monthly effective federal funds rate from the Federal Reserve Economic Data (FRED) database. Contains well-documented policy regime shifts (Volcker tightening, Great Moderation, zero-lower-bound era).

**Citation:**
> Board of Governors of the Federal Reserve System. Effective Federal Funds Rate [FEDFUNDS]. Retrieved from FRED, Federal Reserve Bank of St. Louis. https://fred.stlouisfed.org/series/FEDFUNDS

**Download:**
```bash
# Automatic
python experiments/real_data/download_datasets.py --dataset fedfunds

# Manual
# FRED direct download:
# https://fred.stlouisfed.org/series/FEDFUNDS
# Click "Download" -> CSV
```

**Expected file:** `experiments/real_data/data/fedfunds.csv`

**Preprocessing notes:**
- Monthly frequency, 1954-07 to present
- Use the rate values directly (no differencing needed for level-shift detection)
- Known regime changes: 1979-10 (Volcker), 1982, 2001, 2008-12, 2015-12, 2022-03

---

## 8. PhysioNet Computing in Cardiology Challenge 2019

**Description:** Clinical data from ICU patients for early prediction of sepsis. Each patient record is a multivariate time series of vital signs and lab values. Used to evaluate GSA for distributional shift detection in medical monitoring.

**Citation:**
> Reyna, M. A., Josef, C. S., Jeter, R., Shashikumar, S. P., Westover, M. B., Nemati, S., Clifford, G. D., & Sharma, A. (2019). Early prediction of sepsis from clinical data: the PhysioNet/Computing in Cardiology Challenge 2019. *Critical Care Medicine*, 48(2), 210--217.

**Download:**
```bash
# Automatic (downloads a subset)
python experiments/real_data/download_datasets.py --dataset physionet2019

# Manual
# PhysioNet:
# https://physionet.org/content/challenge-2019/1.0.0/
# Requires PhysioNet credentialed access (free registration)
# Download training set A (hospital system A)
```

**Expected files:** `experiments/real_data/data/physionet2019/` (directory with per-patient PSV files)

**Preprocessing notes:**
- Extract heart rate (HR), mean arterial pressure (MAP), and temperature signals
- Each patient is a separate time series (variable length, 8--336 hours)
- Sepsis onset label column: SepsisLabel (0/1)
- Use pre-sepsis window (6 hours before onset) as the change region

---

## Automated Download

To download all datasets at once:

```bash
python experiments/real_data/download_datasets.py --all
```

Some datasets (NASA IMS, PhysioNet 2019) are large and may require several minutes to download. The script will skip datasets that are already present in the `data/` directory.

## Data Directory Structure

After downloading, the expected layout is:

```
experiments/real_data/data/
    us_realint.csv
    nasa_ims/
        test1/
        test2/
        test3/
    nsl_kdd/
        KDDTrain+.txt
        KDDTest+.txt
    skab/
        anomaly-free/
        valve1/
        valve2/
        ...
    tcpd/
        *.json
    ftse100.csv
    fedfunds.csv
    physionet2019/
        training_setA/
            p*.psv
```
