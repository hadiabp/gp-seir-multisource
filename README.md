# gp-seir-multisource

**Calibrating a Multi-Source Epidemic Model via Generalized Profiling:
Data Stream Valuation Through Cooperative Game Theory**

Hadi A. Peivasti & Martijn Schoot Uiterkamp — Tilburg University

---

## Overview

This repository provides the complete analysis pipeline for fitting a continuous-time SEIR epidemic model with time-varying transmission rate β(t) to four simultaneous Dutch COVID-19 surveillance streams, and for quantifying each stream's informational contribution using Shapley values from cooperative game theory.

**The core question:** How much does each surveillance data stream actually contribute to model accuracy? The answer depends heavily on the epidemic phase — an insight only visible through the rolling-window Shapley analysis.

**Data streams used:**

| Stream | Observation model | Estimated parameter |
|--------|-------------------|---------------------|
| Daily reported cases | μ_case = N · ρ · σ · E(t) | ρ — detection fraction |
| Hospital admissions | μ_hosp = N · p_H · γ · I(t) | p_H — hospitalisation fraction |
| ICU admissions | μ_ICU = N · p_ICU · γ · I(t) | p_ICU — ICU fraction |
| COVID RADAR app (I-fraction) | μ_radar = α_R · I(t) | α_R — app scaling |

**Validation targets** (held out — never used in fitting):
- RIVM R_t estimates (Wallinga–Lipsitch method)
- RIVM prevalence estimates (hospitalisation + seroprevalence)

---

## Repository Structure

```
gp-seir-multisource/
├── README.md
├── LICENSE
├── gp_input_final_v2.csv             Main input data (repository root)
├── R/
│   ├── gp_seir_complete_v2.R         Full-period analysis (primary entry point)
│   ├── gp_seir_wave_analysis_v2.R    Per-wave fits + Shapley/Banzhaf values
│   └── rolling_shapley_v2.R          Rolling-window Shapley + Banzhaf (48 windows)
└── Python/
    └── compute_banzhaf.py            Post-hoc Banzhaf from ablation CSVs (stdlib only)
```

---

## Requirements

### R (≥ 4.1.0)

```r
install.packages(c("data.table", "fda", "minpack.lm", "zoo"))
install.packages("EpiEstim")   # optional — for Cori Rt comparison
```

### Python (≥ 3.9)

Standard library only — no external packages required.

---

## Data

Place `gp_input_final_v2.csv` in the **repository root** before running any script. All commands below assume you run them from the repository root.

The file contains daily Dutch COVID-19 surveillance data for 2020-02-27 to 2023-04-01. Key columns used in the analysis:

| Column | Description |
|--------|-------------|
| `date` | ISO date (YYYY-MM-DD) |
| `cases_raw` | Daily confirmed cases |
| `hosp_nice` | Daily hospital admissions (NICE registry) |
| `icu_daily` | Daily ICU admissions |
| `radar_I_frac_sm7` | COVID RADAR infectious fraction (7-day smoothed) |
| `rt_rivm` | RIVM R_t point estimate (**validation only**) |
| `prev_json_avg` | RIVM prevalence estimate (**validation only**) |

The analysis window is **2020-03-01 to 2021-03-01** (366 days).

---

## Usage

### 1. Full-period analysis

Runs all 15 source-subset ablation fits, computes Shapley and Banzhaf values, profile likelihoods, and sensitivity analyses. Exports 26+ TikZ-ready CSV files.

```bash
# Without delay kernels (output prefix: FALSE_)
Rscript R/gp_seir_complete_v2.R

# With log-normal reporting delay kernels (output prefix: TRUE_) — used in paper
Rscript R/gp_seir_complete_v2.R --delay
```

**Typical runtime:** 2–4 hours (15 subsets × 8 starts each).

### 2. Per-wave analysis

Independent fits for three epidemic periods: Wave 1 (Mar–Jun 2020), Inter-wave (Jun–Oct 2020), Wave 2 (Oct 2020–Mar 2021).

```bash
Rscript R/gp_seir_wave_analysis_v2.R               # all three waves
Rscript R/gp_seir_wave_analysis_v2.R --wave=wave1  # single wave
```

`USE_DELAY_KERNEL = TRUE` is hardcoded — per-wave script always uses delay kernels.

> **Known issue in wave-1 start-date sensitivity:** The loop at line ~1073 (`for (ds in starts)`) coerces R `Date` objects to numeric. Fix: replace with `for (i in seq_along(starts)) { ds <- starts[i]; ... }`. This only affects the wave-1 start-date sensitivity section, not the main fits or Shapley values.

### 3. Rolling-window analysis

48 sliding windows (W = 5 weeks, S = 1 week step) across the full epidemic year.

```bash
Rscript R/rolling_shapley_v2.R              # default settings
Rscript R/rolling_shapley_v2.R --resume     # resume from checkpoint
```

**Typical runtime:** 4–8 hours (~7,200 model fits total).

### 4. Post-hoc Banzhaf (Python)

Recompute Banzhaf and Shapley values from an existing ablation CSV:

```bash
python Python/compute_banzhaf.py
python Python/compute_banzhaf.py my_ablation.csv my_output.csv
```

---

## Method Summary

### Generalized Profiling

SEIR state variables and β(t) are represented as cubic B-splines (14-day and 28-day knots). A two-level cascade optimization avoids explicit ODE solving:

- **Inner loop** (Levenberg–Marquardt): optimize state spline coefficients for fixed β(t)
- **Outer loop** (gradient descent + line search): optimize β(t) and {ρ, p_H, p_ICU, α_R}

ODE fidelity is enforced by a quadratic penalty with weight η = 10⁶ (~1700× ODE/data ratio).

### Data Source Valuation

The characteristic function is the **Pearson correlation between the model's infectious fraction and the held-out RIVM prevalence series**, after fitting with each of the 15 non-empty source subsets. Shapley values give each source's weighted average marginal contribution. Banzhaf values (equal weights) serve as a robustness check.

---

## Model Parameters

### Fixed

| Parameter | Value | Interpretation |
|-----------|-------|----------------|
| σ | 1/5.5 day⁻¹ | Latent period ≈ 5.5 d |
| γ | 1/9.5 day⁻¹ | Infectious period ≈ 9.5 d |
| N | 17,400,000 | Dutch population 2020 |
| η | 10⁶ | ODE penalty weight |
| State knots | 14 days | B-spline resolution |
| β(t) knots | 28 days | Transmission B-spline |
| Delay (cases) | mean 3.0 d, SD 1.5 d | Symptom-onset to report |
| Delay (hosp) | mean 2.0 d, SD 1.0 d | Removal-flow to admission |
| Delay (ICU) | mean 3.5 d, SD 1.5 d | Removal-flow to ICU record |

### Estimated

| Parameter | Bounds | Interpretation |
|-----------|--------|----------------|
| ρ | [0.05, 0.60] | Case detection fraction |
| p_H | [0.002, 0.05] | Hospitalisation fraction |
| p_ICU | [0.0003, 0.02] | ICU admission fraction |
| α_R | [0.10, 1.00] | RADAR app-to-prevalence scaling |

---

## Key Results

### Full-period (366 days, delay kernel, η = 10⁶)

| Metric | Value |
|--------|-------|
| R_t Pearson (GP vs RIVM) | 0.728 |
| R_t Pearson (Cori vs RIVM) | 0.332 |
| Prevalence Pearson | 0.885 |
| ODE RMSE | 2.5 × 10⁻⁴ |
| ρ (detection) | 0.217 |
| α_R (RADAR scaling) | 0.456 |

### Per-wave Shapley values — prevalence Pearson characteristic function

| Source | Wave 1 | Inter-wave | Wave 2 |
|--------|--------|------------|--------|
| Cases | −0.143 | +0.237 | +0.207 |
| Hospital | **+0.409** | +0.247 | +0.243 |
| **ICU** | **+0.410** | +0.245 | **+0.266** |
| RADAR | +0.138 | +0.243 | +0.096 |
| v(all 4) | 0.812 | 0.972 | 0.812 |

### Per-wave Shapley values — R_t Pearson characteristic function

| Source | Wave 1 | Inter-wave | Wave 2 |
|--------|--------|------------|--------|
| Cases | −0.177 | +0.072 | +0.159 |
| **Hospital** | **+0.405** | +0.132 | +0.189 |
| ICU | +0.362 | +0.112 | +0.178 |
| RADAR | +0.278 | **+0.163** | **+0.222** |
| v(all 4) | 0.869 | 0.480 | 0.748 |

Shapley and normalised Banzhaf rankings are identical in all 6 wave × metric combinations.

### Rolling Shapley (48 windows, prevalence)

- RADAR rank-1 in **40%** of windows (19/48); Cases 27% (13/48); Hosp 23% (11/48); ICU 10% (5/48)
- Rank-1 source changes **23 times** across 48 windows
- RADAR arc: φ = **+0.518** (Oct 14 — wave-2 onset) → φ = **−0.622** (Feb 10) at unchanged ~10k/day engagement, detecting representativeness erosion invisible to conventional monitoring

---

## Citation

> Peivasti, H.A. & Schoot Uiterkamp, M. (2026).
> *Calibrating a Multi-Source Epidemic Model via Generalized Profiling:
> Data Stream Valuation Through Cooperative Game Theory.*
> Tilburg University.

---

## License

MIT — see [LICENSE](LICENSE).
