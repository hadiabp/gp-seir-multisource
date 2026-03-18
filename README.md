# GP-SEIR: Multi-Source Generalized Profiling for Epidemic Model Calibration

Code and data for: **Calibration of Epidemic Models via Multi-Source Generalized Profiling** (Peivasti & Schoot Uiterkamp, 2026, Tilburg University).

## Key Results

- Rt correlation with RIVM: **0.82** (not used in fitting)
- Prevalence correlation: **0.85** (fully independent validation)
- RADAR contact data provides **highest marginal information** among all streams
- Hospital+ICU+RADAR without cases achieves **0.90** prevalence correlation

## Structure

```
R/fit_gp_seir_analysis.R   — Complete analysis (modes: full/manual/ablation_only)
data/gp_input_final_v2.csv — Processed input data
results/                    — Output CSVs (ablation, weights, sensitivity, wave refit)
tikz/                       — CSV data for pgfplots/LaTeX figures
plots/gp_seir_plots.pdf     — Generated figures (10 pages)
```

## Data Sources

| Source | URL |
|--------|-----|
| RIVM COVID-19 | https://data.rivm.nl/covid-19/ |
| COVID RADAR (DANS) | https://doi.org/10.17026/dans-zcd-m9dh |
| Schoot Uiterkamp et al. (2025) | https://github.com/mhhschootuiterkamp/DCM-data-integration-in-SEIR-models |

## Requirements

```r
install.packages(c("data.table", "fda", "minpack.lm", "zoo"))
```

## Usage

```r
# Full analysis (several hours): set RUN_MODE <- "full"
# Quick fit with custom weights: set RUN_MODE <- "manual"
# Ablation study only: set RUN_MODE <- "ablation_only"
# Change ranking: set SORT_CRITERION <- "prev_corr" or "rt_corr"
Rscript R/fit_gp_seir_analysis.R
```

## License

MIT
