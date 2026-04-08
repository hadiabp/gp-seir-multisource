# gp-seir-multisource

**Calibration of Epidemic Models via Multi-Source Generalized Profiling: Data Source Valuation Through Cooperative Game Theory**

Authors: Hadi A. Peivasti & Martijn Schoot Uiterkamp (Tilburg University)

This repository contains the R and Python code accompanying the paper. We fit an SEIR model with time-varying transmission rate to four Dutch COVID-19 surveillance streams (cases, hospital admissions, ICU admissions, and the COVID RADAR contact-monitoring app) using generalized profiling, and then quantify the information value of each data source via Shapley and Banzhaf values from cooperative game theory.

## What the code does

1. **Multi-source GP-SEIR fit.** Fits an SEIR model with time-varying β(t) to all four data streams simultaneously using B-spline state representations and parameter cascading. No numerical ODE solver inside the optimization loop.
2. **Source ablation.** Re-fits the model on each of the 2⁴−1 = 15 non-empty subsets of the four data sources by setting binary source indicators z_m ∈ {0,1}.
3. **Cooperative game-theoretic analysis.** Computes exact Shapley values, Banzhaf values, Shapley interaction indices, leave-one-out impacts, and superadditivity ratios from the 15 ablation results, treating each data source as a player and the external-validation correlation against independent RIVM Rt or prevalence estimates as the characteristic function.
4. **Sensitivity analyses.** Multi-start optimization (5 warm + 3 cold), profile likelihoods for the four scaling parameters, σ/γ sensitivity, ODE penalty weight sensitivity, three β(t) initialization strategies, optional log-normal delay kernels.

## Repository structure

```
gp-seir-multisource/
├── README.md                       This file
├── LICENSE                         MIT
├── R/
│   ├── gp_seir_complete.R          Full-period (366-day) analysis (primary entry point)
│   ├── gp_seir_wave_analysis.R     Per-wave independent fits (alternative analysis)
│   └── helpers/
│       ├── seir_rhs.R              SEIR vector field on log scale
│       ├── delay_kernels.R         Log-normal delay convolution
│       ├── shapley.R               Exact Shapley value computation
│       └── interactions.R          Shapley interaction index
├── python/
│   └── compute_banzhaf.py          Banzhaf values from existing ablation CSVs
├── docs/
│   ├── paper.tex                   The paper
│   ├── references.bib              Bibliography
│   └── PROJECT_REFERENCE.md        Comprehensive project reference
└── data/
    └── (input data not committed; see Data section below)
```

## Quick start

### Requirements

R packages: `data.table`, `fda`, `minpack.lm`, `zoo`
Python (only for the Banzhaf post-hoc script): standard library, no external dependencies

### Run the full analysis

```r
# In R, working directory containing the data CSV
source("R/gp_seir_complete.R")
# Produces FALSE_*.csv files in the working directory
```

To re-run with explicit log-normal delay kernels (robustness check):

```r
USE_DELAY_KERNEL <- TRUE
source("R/gp_seir_complete.R")
# Produces TRUE_*.csv files
```

To run the per-wave independent-fit analysis (different aggregation strategy):

```r
source("R/gp_seir_wave_analysis.R")
```

### Compute Banzhaf values from existing ablation CSVs

The ablation CSVs already contain the per-coalition validation values for both Rt and prevalence in 12 different time windows (full period, pre/post-June, three waves × two metrics). The Python script reads these and computes both Shapley (sanity check; output should match `*_shapley_values.csv` exactly) and Banzhaf values for all metrics.

```bash
python3 python/compute_banzhaf.py
# Reads FALSE_ablation_results.csv and TRUE_ablation_results.csv
# Outputs FALSE_banzhaf_values.csv and TRUE_banzhaf_values.csv
```

## Two configurations

The code supports two configurations selected by a global flag:

| Flag | Behavior | Used in paper as |
|---|---|---|
| `USE_DELAY_KERNEL=FALSE` | Instantaneous observation models. | **Primary** results |
| `USE_DELAY_KERNEL=TRUE` | Each count-based observation is convolved with a log-normal delay kernel (case=3.0d mean, hosp=2.0d, ICU=3.5d). | Robustness check |

We use FALSE as primary because at the standard B-spline knot spacing (14 days), explicit delay kernels concentrating within 1–3 days are absorbed by the spline flexibility. The TRUE configuration produces virtually identical parameter estimates and external-validation correlations (to three decimal places), with only the cost function value differing by ~10%.

## Two analyses

The code contains two related but distinct analyses:

| Script | Approach | Used by |
|---|---|---|
| `gp_seir_complete.R` | **Full-period (366-day) fit** with sliced Shapley evaluation per wave. Same fit used everywhere; only the validation window changes. | The paper |
| `gp_seir_wave_analysis.R` | **Independent per-wave fits**, separate B-spline bases, separate scaling parameters per wave. | Technical report (separate document) |

The per-wave refits reveal substantial parameter heterogeneity (e.g., the detection fraction ρ̂ varies 2.8-fold from 0.114 in wave 1 to 0.319 in wave 2), exposing the constant-parameter assumption of the full-period fit as a source of misspecification.

## Data

The analysis uses public Dutch COVID-19 data:

- **Cases, hospital admissions, ICU admissions** — RIVM Dutch national COVID-19 database (Geubbels et al., 2023, *Scientific Data*)
- **COVID RADAR contact data** — DANS dataset (van Dijk, 2020; doi:10.17026/dans-zcd-m9dh)
- **RIVM Rt and prevalence estimates** — RIVM open data, used exclusively for external validation (zero weight in the fitting objective)

Population: 17,400,000 (Statistics Netherlands StatLine table 03759ned, 1 January 2020).

The processed input file `gp_input_final_v2.csv` (pre-aligned daily series) is not committed to the repository; the script `prepare_data.R` (not included) was used to download and align the original sources.

## Headline results

| Quantity | Value |
|---|---|
| Detection fraction ρ̂ | 0.217 |
| Hospital fraction p̂_H | 0.0108 |
| ICU fraction p̂_ICU | 0.00267 |
| RADAR scaling α̂_R | 0.457 |
| Rt Pearson correlation with RIVM | 0.728 |
| Prevalence Pearson correlation with RIVM | 0.885 |
| Best single source for prevalence | **hospital alone** (r = 0.907) |
| All-four-source prevalence correlation | 0.904 |
| Shapley redundancy (Rt) | 75.5% |
| Shapley redundancy (Prev) | 72.5% |
| Wave-2 case Shapley value (Prev) | **−0.093** |
| Wave-2 RADAR Shapley value (Prev) | **+0.212** |

All twelve Shapley interaction indices are negative, indicating universal substitutability between data sources. Banzhaf values agree with Shapley to within ±0.04 after normalization, confirming the rankings are robust to the choice of cooperative solution concept.

## Citation

If you use this code or find the methodology useful, please cite:

```bibtex
@article{Peivasti2026GPSEIR,
  author  = {Peivasti, Hadi A. and Schoot Uiterkamp, Martijn},
  title   = {Calibration of Epidemic Models via Multi-Source Generalized Profiling:
             Data Source Valuation Through Cooperative Game Theory},
  year    = {2026},
  journal = {arXiv preprint}
}
```

## License

MIT License — see `LICENSE` for details.
