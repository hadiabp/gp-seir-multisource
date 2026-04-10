# gp-seir-multisource

**Calibration of Epidemic Models via Multi-Source Generalized Profiling:
Data Source Valuation Through Cooperative Game Theory**

Hadi A. Peivasti & Martijn Schoot Uiterkamp — Tilburg University

---

## Overview

This repository provides the full analysis pipeline for fitting a continuous-time
SEIR epidemic model with time-varying transmission rate β(t) to four simultaneous
Dutch COVID-19 surveillance streams, and for quantifying each stream's informational
contribution using Shapley and Banzhaf values from cooperative game theory.

The four data streams are:

| Stream | Observation model | Fitted parameter |
|--------|------------------|-----------------|
| Daily reported cases | μ_case = N · ρ · σ · E(t) | ρ (detection fraction) |
| Hospital admissions | μ_hosp = N · p_H · γ · I(t) | p_H (hosp. fraction) |
| ICU admissions | μ_icu = N · p_ICU · γ · I(t) | p_ICU (ICU fraction) |
| COVID RADAR app (I-fraction) | μ_radar = α_R · I(t) | α_R (app scaling) |

Validation targets (held-out — never used in fitting):
- RIVM R_t estimates (Wallinga–Lipsitch method)
- RIVM prevalence estimates (hospitalisation + seroprevalence)

---

## Repository Structure

```
gp-seir-multisource/
├── README.md
├── LICENSE
├── .gitignore
├── R/
│   ├── gp_seir_complete_v2.R         Full-period analysis (primary entry point)
│   ├── gp_seir_wave_analysis_v2.R    Per-wave independent fits + game theory
│   └── rolling_shapley_v2.R          Rolling-window Shapley + Banzhaf
└── python/
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

Place `gp_input_final_v2.csv` in the **repository root** before running any script.

All commands in this README assume you run them from the repository root, so the R scripts can find `gp_input_final_v2.csv` via a relative path.

The file contains daily Dutch COVID-19 surveillance data for the period
2020-02-27 to 2023-04-01 with the following columns (among others):

| Column | Description |
|--------|-------------|
| `date` | ISO date (YYYY-MM-DD) |
| `cases_raw` | Daily confirmed cases |
| `hosp_nice` | Daily hospital admissions |
| `icu_daily` | Daily ICU admissions |
| `radar_I_frac_sm7` | COVID RADAR app infectious fraction (7-day smoothed) |
| `rt_rivm` | RIVM R_t point estimate (**validation only**) |
| `rt_low`, `rt_up` | RIVM R_t 95% interval bounds |
| `prev_json_avg` | RIVM prevalence estimate (**validation only**) |
| `prev_json_low`, `prev_json_up` | Prevalence interval bounds |
| `tests_total`, `tests_positive` | Daily testing volumes |
| `positivity_rate` | Test positivity rate |

The analysis window used in the paper is 2020-03-01 to 2021-03-01 (366 days).

---

## Usage

### 1. Full-period analysis

The primary script. Runs all 15 source-subset ablation fits, computes Shapley and
Banzhaf values for 16 characteristic functions, profile likelihoods, sensitivity
analyses, Cori comparison, and exports 26+ TikZ-ready CSV files.

```bash
# Run from the repository root

# Without delay kernels (output prefix: FALSE_)
Rscript R/gp_seir_complete_v2.R

# With log-normal reporting delay kernels (output prefix: TRUE_)
Rscript R/gp_seir_complete_v2.R --delay
```

Typical runtime: 2–4 hours on a modern laptop (15 subsets × 8 starts each).

---

### 2. Per-wave analysis

Independent fits for each of three epidemic periods
(Wave 1: Mar–Jun 2020; Inter-wave: Jun–Oct 2020; Wave 2: Oct 2020–Mar 2021).
Produces per-wave Shapley/Banzhaf values, sensitivity analyses, and Cori comparison.

```bash
# Run from the repository root
Rscript R/gp_seir_wave_analysis_v2.R                         # all three waves
Rscript R/gp_seir_wave_analysis_v2.R --wave=wave1            # single wave
Rscript R/gp_seir_wave_analysis_v2.R --waves=wave1,wave2     # two waves
```

Note: `USE_DELAY_KERNEL = TRUE` is hardcoded — the per-wave script always uses delay kernels.

---

### 3. Rolling-window analysis

Sliding-window (W = 5 weeks, S = 1 week step) Shapley and Banzhaf values
across 48 windows spanning the full year.

```bash
# Run from the repository root
Rscript R/rolling_shapley_v2.R                               # default settings
Rscript R/rolling_shapley_v2.R --window=7 --step=2           # wider window
Rscript R/rolling_shapley_v2.R --resume                      # resume from checkpoint
```

Typical runtime: 4–8 hours (48 windows × 15 subsets × 10 starts each).

---

### 4. Post-hoc Banzhaf (Python)

Recompute Banzhaf and Shapley values from an existing ablation CSV without re-running R:

```bash
# Run from the repository root
python python/compute_banzhaf.py                                  # processes FALSE_ and TRUE_ defaults
python python/compute_banzhaf.py my_ablation.csv my_output.csv    # custom paths
```

---

## Method

### Generalized Profiling

State variables S(t), E(t), I(t), R(t) and transmission rate β(t) are represented as
cubic B-splines (14-day and 28-day knots respectively). A two-level cascade optimization
avoids explicit ODE solving during the fit:

- **Inner loop** (Levenberg–Marquardt): fix β(t) and scaling parameters; optimize spline
  coefficients for the state variables
- **Outer loop** (gradient descent): optimize β(t) and {ρ, p_H, p_ICU, α_R}

ODE fidelity is enforced by a quadratic penalty with weight η = 10⁶ (~1700× ODE/data ratio).

### Data Source Valuation

The characteristic function is Pearson correlation between the model R_t trajectory and
the held-out RIVM R_t series, evaluated after fitting with only the sources in coalition S.

**Why Pearson, not RMSE:** Pearson is shift/scale-invariant and measures directional R_t
tracking only. RMSE conflates tracking quality with a systematic level offset (bias up to
±1.1 R_t units) arising from the constant detection-fraction assumption. A source that
corrects this calibration artifact — rather than improving epidemic tracking — would receive
spuriously high RMSE-Shapley value.

**Alternative:** Skill score v(S) = 1 − MSE(S)/Var(R_t^RIVM) — natural zero, avoids negative
singletons. Both Pearson and skill score are computed by all v2 scripts.

### Shapley Values

```
φ_i = Σ_{S ⊆ N\{i}} [|S|!(n−|S|−1)!/n!] · [v(S∪{i}) − v(S)]
```

Requires all 2⁴ − 1 = 15 non-empty subset fits. Satisfies efficiency (Σφ_i = v(N)),
symmetry, linearity, and null-player axioms.

### Banzhaf Values

```
β_i = (1/2^{n−1}) · Σ_{S ⊆ N\{i}} [v(S∪{i}) − v(S)]
```

Same 15 fits — no additional model runs. Normalised Banzhaf:
β̃_i = β_i · v(N) / Σ_j β_j (requires |Σβ_j| ≥ 0.05 and sign(Σβ_j) = sign(v(N))).

---

## Model Parameters

### Fixed (working assumptions — same across all scripts)

| Parameter | Value | Interpretation |
|-----------|-------|---------------|
| σ | 1/5.5 day⁻¹ | Latent-to-infectious transition rate (latent period ≈ 5.5 d) |
| γ | 1/9.5 day⁻¹ | Infectious-to-removed transition rate (inf. period ≈ 9.5 d) |
| N | 17,400,000 | Dutch population 2020 |
| η | 10⁶ | ODE penalty weight (≈ 1700× ODE/data ratio) |
| κ_β | 10³ | β(t) roughness penalty |
| κ_mass | 10⁶ | Mass conservation penalty |
| β_max | 0.6 | Soft upper cap on transmission rate |
| State knots | 14 days | B-spline resolution for S, E, I, R |
| β(t) knots | 28 days | B-spline resolution for transmission rate |
| Delay (cases) | mean 3.0 d, SD 1.5 d | Symptom-onset to case-report delay (log-normal) |
| Delay (hosp) | mean 2.0 d, SD 1.0 d | Removal-flow to admission-record delay |
| Delay (ICU) | mean 3.5 d, SD 1.5 d | Removal-flow to ICU-record delay |

### Estimated per-analysis

| Parameter | Bounds | Interpretation |
|-----------|--------|---------------|
| ρ | [0.05, 0.60] | Case detection fraction |
| p_H | [0.002, 0.05] | Hospitalisation fraction |
| p_ICU | [0.0003, 0.02] | ICU admission fraction |
| α_R | [0.10, 1.00] | RADAR app–to–prevalence scaling |

### Rolling-window settings

| Parameter | Value | Notes |
|-----------|-------|-------|
| Window width W | 5 weeks (35 days) | Default |
| Step size S | 1 week (7 days) | 48 windows over 366 days |
| β initialisation | constant (β₀ = γ everywhere) | Avoids hospital-init bias in rolling context |
| Warm starts per subset | 4 | |
| Cold starts per subset | 3 | From default init, wider perturbation |
| Level starts per subset | 3 | Uninformed flat β |
| Total starts | 10 per subset | 150 fits per window |
| Stability threshold | 0.05 | Flag if \|Σβ_i\| < 0.05 or sign mismatch |

### Cori serial interval (EpiEstim)

| Parameter | Value | Source |
|-----------|-------|--------|
| Mean SI | 5.1 days | Bi et al. (2020), COVID-19 serial interval |
| SD SI | 2.8 days | Same |
| Window τ | 7 days | |

---

## Output Files

### gp_seir_complete_v2.R (prefix: `FALSE_` or `TRUE_`)

**Model fit and game theory:**
- `*_gp_seir_results.csv` — Full time series: S, E, I, R, β, R_t, model predictions, residuals, flows
- `*_ablation_results.csv` — All 15 source subsets × all metrics (Pearson, Spearman, skill, RMSE, bias, per-wave)
- `*_shapley_banzhaf_values.csv` — φ, β, β̃ for 16 characteristic functions × 4 sources
- `*_shapley_values.csv` — Shapley only (backward-compatible wide format)
- `*_shapley_interactions_all.csv` — Pairwise Shapley interaction indices I_ij for all metrics
- `*_leave_one_out.csv` — Drop-one-source analysis
- `*_superadditivity.csv` — v(N) vs Σv({i}), redundancy %, best pair

**Validation:**
- `*_rt_comprehensive.csv` — R_t: Pearson/Spearman/skill/R²/RMSE/MAE/bias/threshold, pre-June, per-wave
- `*_prev_comprehensive.csv` — Prevalence: same metrics
- `*_rt_method_comparison.csv` — GP-SEIR vs Cori vs RIVM: Pearson/Spearman/RMSE/bias/skill/threshold
- `*_wave_params.csv` — Per-wave sliced metrics from main fit
- `*_wave_refit_params.csv` — Independent per-wave refits with rho, pH, alphaR
- `*_wave_summary.csv` — Attack rates, peak timing, mean R_t, %time above 1

**Sensitivity:**
- `*_wdata_sensitivity.csv` — 6 ODE-weight values: rt_corr, skill, ode_rmse, cases_rmse
- `*_beta_init_sensitivity.csv` — 3 initialisation strategies: hospital, constant, random
- `*_profile_likelihood.csv` — Approximate profile likelihood for ρ, p_H, p_ICU, α_R
- `*_sigma_gamma_sensitivity.csv` — 4 alternate σ/γ configurations

**Diagnostics:**
- `*_residual_diagnostics.csv` — ACF lag-1, Shapiro–Wilk, Ljung–Box per stream
- `*_stream_fit_metrics.csv` — Per-stream Pearson, RMSE, MAE, bias
- `*_cross_correlations.csv` — Pairwise CCF and best lag between streams
- `*_rt_lag_analysis.csv` — Model–RIVM Pearson/RMSE at lags −10..+20 days
- `*_ode_compliance.csv` — ODE residual RMSE and max per compartment
- `*_derived_timeseries.csv` — New infections, onset/removal flows, growth rates, doubling times
- `*_multistart_params.csv` — Per-start J, ρ, pH, pICU, α_R, rt_corr
- `*_summary.csv` — All key numbers in a single long-format file

**TikZ-ready (for LaTeX figures):**
- `*_tikz_seir.csv`, `*_tikz_rt.csv`, `*_tikz_prev.csv` — Main trajectories
- `*_tikz_cases.csv`, `*_tikz_hosp.csv`, `*_tikz_icu.csv`, `*_tikz_radar.csv` — Stream fits
- `*_tikz_residuals.csv`, `*_tikz_acf.csv`, `*_tikz_qq.csv` — Diagnostics
- `*_tikz_detection.csv`, `*_tikz_cumulative.csv` — Detection and R(t)
- `*_tikz_ablation.csv`, `*_tikz_beta_ablation.csv` — Ablation
- `*_tikz_shapley_banzhaf.csv` — Shapley + normBanzhaf for bar charts
- `*_tikz_rt_scatter.csv`, `*_tikz_bland_altman.csv` — R_t validation plots
- `*_tikz_epiestim.csv`, `*_tikz_rt_threeway.csv` — GP vs Cori vs RIVM
- `*_tikz_profile.csv`, `*_tikz_rt_wdata.csv`, `*_tikz_sigma_gamma.csv`

**Workspace:**
- `*_gp_seir_plots.pdf` — 9-page diagnostic PDF
- `*_gp_seir_workspace.RData` — Complete R workspace (load to resume)
- `*_gp_seir_checkpoint_*.RData` — Checkpoints after each expensive step

### rolling_shapley_v2.R (prefix: `rolling_shapley_w{W}_s{S}_{init}`)

- `*_subset_results.csv` — All 15 subsets × 48 windows × all metrics (primary output)
- `*_shapley_all.csv`, `*_banzhaf_all.csv` — Long-format values with stability flags
- `*_tikz_rt.csv`, `*_tikz_prev.csv` — Wide-format Shapley for pgfplots
- `*_tikz_banzhaf_rt.csv`, `*_tikz_banzhaf_prev.csv` — Wide-format normBanzhaf (NaN = unstable)
- `*_rank1_agreement.csv` — Shapley/Banzhaf agreement at 3 restriction levels
- `*_windows.csv` — All-4-source summary: R_t corr, params, peak, attack rate
- `*_redundancy.csv`, `*_interactions.csv` — Redundancy and pairwise interactions over time
- `*_cori_comparison.csv` — GP vs Cori per window (requires EpiEstim)

---

## Key Results Summary

Per-wave Shapley values φ_i under R_t Pearson characteristic function:

| Source | Wave 1 | Inter-wave | Wave 2 |
|--------|--------|------------|--------|
| Cases | +0.077 | +0.072 | +0.159 |
| Hospital | +0.423 | +0.132 | +0.189 |
| **ICU** | **+0.513** | +0.112 | +0.178 |
| RADAR | −0.193 | **+0.163** | **+0.222** |

All-4-source R_t correlation: Wave 1 = 0.820, Inter-wave = 0.480, Wave 2 = 0.748.

Rankings are identical under Shapley and normalised Banzhaf in all 6 period × metric combinations.
Rolling-window rank-1 agreement: **86%** (R_t) and **95%** (prevalence) across stable signal windows.

---

## Citation

If you use this code, please cite:

> Peivasti, H.A. & Schoot Uiterkamp, M. (2026).
> *Calibration of Epidemic Models via Multi-Source Generalized Profiling:
> Data Source Valuation Through Cooperative Game Theory.*
> Tilburg University.

---

## License

MIT — see [LICENSE](LICENSE).
