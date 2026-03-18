#!/usr/bin/env Rscript
##############################################################################
# supplement_rt_analysis.R
#
# Run AFTER fit_gp_seir_analysis.R to produce additional Rt outputs.
# Reads: gp_seir_results.csv, ablation_results.csv, weight_optimization.csv,
#        wdata_sensitivity.csv, wave_refit_params.csv, wave_params.csv
# Produces: tikz_rt_scatter.csv, tikz_rt_waves.csv, tikz_rt_prev_tradeoff.csv,
#           rt_comprehensive.csv, gp_seir_plots_supplement.pdf
##############################################################################

suppressPackageStartupMessages(library(data.table))

cat("Supplementary Rt Analysis\n\n")

gp   <- fread("gp_seir_results.csv")
abl  <- fread("ablation_results.csv")
wopt <- fread("weight_optimization.csv")
ws   <- fread("wdata_sensitivity.csv")

dates <- as.Date(gp$date)
N <- 17400000

# ══════════════════════════════════════════════════════════════════
# 1. Comprehensive Rt metrics for main fit
# ══════════════════════════════════════════════════════════════════

ok <- is.finite(gp$rt_rivm) & is.finite(gp$Rt_model)

rt_metrics <- data.frame(
  metric = c("pearson", "spearman", "R2", "RMSE", "MAE", "bias",
             "threshold_agreement",
             "pearson_pre_jun", "spearman_pre_jun", "RMSE_pre_jun", "bias_pre_jun",
             "pearson_post_jun", "spearman_post_jun", "RMSE_post_jun", "bias_post_jun",
             "pearson_wave1", "pearson_interwave", "pearson_wave2",
             "RMSE_wave1", "RMSE_interwave", "RMSE_wave2"),
  value = NA_real_
)

r_rivm <- gp$rt_rivm[ok]; r_mod <- gp$Rt_model[ok]
rt_metrics$value[1] <- cor(r_rivm, r_mod)
rt_metrics$value[2] <- cor(r_rivm, r_mod, method = "spearman")
ss_res <- sum((r_rivm - r_mod)^2); ss_tot <- sum((r_rivm - mean(r_rivm))^2)
rt_metrics$value[3] <- 1 - ss_res / ss_tot
rt_metrics$value[4] <- sqrt(mean((r_rivm - r_mod)^2))
rt_metrics$value[5] <- mean(abs(r_rivm - r_mod))
rt_metrics$value[6] <- mean(r_mod - r_rivm)
rt_metrics$value[7] <- mean((r_mod > 1) == (r_rivm > 1))

# Pre/post June 12
cutoff <- as.Date("2020-06-12")
pre <- ok & dates < cutoff; post <- ok & dates >= cutoff
rt_metrics$value[8]  <- cor(gp$rt_rivm[pre], gp$Rt_model[pre])
rt_metrics$value[9]  <- cor(gp$rt_rivm[pre], gp$Rt_model[pre], method = "spearman")
rt_metrics$value[10] <- sqrt(mean((gp$rt_rivm[pre] - gp$Rt_model[pre])^2))
rt_metrics$value[11] <- mean(gp$Rt_model[pre] - gp$rt_rivm[pre])
rt_metrics$value[12] <- cor(gp$rt_rivm[post], gp$Rt_model[post])
rt_metrics$value[13] <- cor(gp$rt_rivm[post], gp$Rt_model[post], method = "spearman")
rt_metrics$value[14] <- sqrt(mean((gp$rt_rivm[post] - gp$Rt_model[post])^2))
rt_metrics$value[15] <- mean(gp$Rt_model[post] - gp$rt_rivm[post])

# Per-wave
w1 <- ok & dates < "2020-06-01"
iw <- ok & dates >= "2020-06-01" & dates < "2020-10-01"
w2 <- ok & dates >= "2020-10-01"
rt_metrics$value[16] <- cor(gp$rt_rivm[w1], gp$Rt_model[w1])
rt_metrics$value[17] <- cor(gp$rt_rivm[iw], gp$Rt_model[iw])
rt_metrics$value[18] <- cor(gp$rt_rivm[w2], gp$Rt_model[w2])
rt_metrics$value[19] <- sqrt(mean((gp$rt_rivm[w1] - gp$Rt_model[w1])^2))
rt_metrics$value[20] <- sqrt(mean((gp$rt_rivm[iw] - gp$Rt_model[iw])^2))
rt_metrics$value[21] <- sqrt(mean((gp$rt_rivm[w2] - gp$Rt_model[w2])^2))

fwrite(rt_metrics, "rt_comprehensive.csv")
cat("  Saved rt_comprehensive.csv\n")

# ══════════════════════════════════════════════════════════════════
# 2. TikZ: Rt scatter plot data
# ══════════════════════════════════════════════════════════════════

# Add wave labels for coloring
wave_label <- ifelse(dates < "2020-06-01", "wave1",
              ifelse(dates < "2020-10-01", "inter", "wave2"))
fwrite(data.table(
  rt_rivm = gp$rt_rivm, rt_model = gp$Rt_model,
  wave = wave_label, day = gp$day
), "tikz_rt_scatter.csv")
cat("  Saved tikz_rt_scatter.csv\n")

# ══════════════════════════════════════════════════════════════════
# 3. TikZ: Per-wave Rt bar chart data
# ══════════════════════════════════════════════════════════════════

wave_rt <- data.frame(
  period = c("Wave 1\n(Mar-Jun)", "Inter-wave\n(Jun-Oct)", "Wave 2\n(Oct-Mar)", "Full"),
  pearson = c(rt_metrics$value[16], rt_metrics$value[17], rt_metrics$value[18], rt_metrics$value[1]),
  rmse = c(rt_metrics$value[19], rt_metrics$value[20], rt_metrics$value[21], rt_metrics$value[4])
)
fwrite(wave_rt, "tikz_rt_waves.csv")
cat("  Saved tikz_rt_waves.csv\n")

# ══════════════════════════════════════════════════════════════════
# 4. TikZ: Rt vs Prev tradeoff across all configs
# ══════════════════════════════════════════════════════════════════

# Ablation configs
abl_trade <- data.table(
  config = abl$sources, type = "ablation",
  rt_corr = abl$rt_corr, prev_corr = abl$prev_corr
)

# Weight optimization - tag Pareto-optimal
wopt_trade <- data.table(
  config = paste0("c=", wopt$w_cases, " h=", wopt$w_hosp,
                  " i=", wopt$w_icu, " r=", wopt$w_radar),
  type = "weight_opt",
  rt_corr = wopt$rt_corr, prev_corr = wopt$prev_corr
)

# Mark Pareto-optimal points in weight optimization
wopt_trade[, pareto := {
  p <- rep(FALSE, .N)
  for (i in seq_len(.N)) {
    dominated <- FALSE
    for (j in seq_len(.N)) {
      if (j != i && rt_corr[j] >= rt_corr[i] && prev_corr[j] >= prev_corr[i] &&
          (rt_corr[j] > rt_corr[i] || prev_corr[j] > prev_corr[i])) {
        dominated <- TRUE; break
      }
    }
    p[i] <- !dominated
  }
  p
}]

trade <- rbind(abl_trade[, pareto := TRUE], wopt_trade, fill = TRUE)
fwrite(trade, "tikz_rt_prev_tradeoff.csv")
cat("  Saved tikz_rt_prev_tradeoff.csv\n")

# ══════════════════════════════════════════════════════════════════
# 5. Ablation sorted by Rt (alternative ranking)
# ══════════════════════════════════════════════════════════════════

abl_rt <- abl[order(-abl$rt_corr), ]
fwrite(abl_rt[, .(sources, rt_corr, rt_rmse, prev_corr)], "tikz_ablation_rt.csv")
cat("  Saved tikz_ablation_rt.csv\n")

# ══════════════════════════════════════════════════════════════════
# 6. PDF: Supplementary Rt figures
# ══════════════════════════════════════════════════════════════════

pdf("gp_seir_plots_supplement.pdf", width = 14, height = 10)

# Page 1: Rt scatter plot (colored by wave)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
cols <- ifelse(wave_label == "wave1", "red",
        ifelse(wave_label == "inter", "blue", "darkgreen"))
plot(gp$rt_rivm[ok], gp$Rt_model[ok], pch = 16, cex = 0.5, col = cols[ok],
     xlab = expression(R[t]^RIVM), ylab = expression(R[t]^model),
     main = sprintf("A. Rt scatter (r=%.3f, R²=%.3f)", rt_metrics$value[1], rt_metrics$value[3]),
     xlim = c(0.3, 2.3), ylim = c(0.3, 2.3))
abline(0, 1, col = "gray40", lty = 2)
abline(h = 1, v = 1, col = "gray80", lty = 3)
legend("topleft", c("Wave 1", "Inter-wave", "Wave 2"),
       col = c("red", "blue", "darkgreen"), pch = 16, cex = 0.7)

# Bland-Altman plot
mean_rt <- (gp$rt_rivm[ok] + gp$Rt_model[ok]) / 2
diff_rt <- gp$Rt_model[ok] - gp$rt_rivm[ok]
plot(mean_rt, diff_rt, pch = 16, cex = 0.4, col = cols[ok],
     xlab = expression("Mean of " * R[t]^RIVM * " and " * R[t]^model),
     ylab = expression(R[t]^model - R[t]^RIVM),
     main = sprintf("B. Bland-Altman (bias=%.3f, LOA=[%.2f,%.2f])",
                    mean(diff_rt), mean(diff_rt)-1.96*sd(diff_rt),
                    mean(diff_rt)+1.96*sd(diff_rt)))
abline(h = mean(diff_rt), col = "red")
abline(h = mean(diff_rt) + c(-1.96, 1.96) * sd(diff_rt), col = "red", lty = 2)
abline(h = 0, col = "gray60", lty = 3)

# Page 2: Per-wave Rt comparison
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (wl in c("wave1", "inter", "wave2")) {
  idx <- ok & wave_label == wl
  lbl <- switch(wl, wave1="Wave 1 (Mar-Jun)", inter="Inter-wave (Jun-Oct)", wave2="Wave 2 (Oct-Mar)")
  pc <- cor(gp$rt_rivm[idx], gp$Rt_model[idx])
  rm <- sqrt(mean((gp$rt_rivm[idx] - gp$Rt_model[idx])^2))
  plot(gp$rt_rivm[idx], gp$Rt_model[idx], pch = 16, cex = 0.6,
       col = switch(wl, wave1="red", inter="blue", wave2="darkgreen"),
       xlab = expression(R[t]^RIVM), ylab = expression(R[t]^model),
       main = sprintf("%s\nr=%.3f, RMSE=%.3f", lbl, pc, rm),
       xlim = c(0.3, 2.3), ylim = c(0.3, 2.3))
  abline(0, 1, col = "gray40", lty = 2)
  abline(h = 1, v = 1, col = "gray80", lty = 3)
}

# Page 3: Rt vs Prevalence tradeoff
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(wopt$rt_corr, wopt$prev_corr, pch = 16, cex = 0.7, col = "gray60",
     xlab = expression(R[t] ~ "correlation"), ylab = "Prevalence correlation",
     main = "C. Rt vs Prevalence tradeoff (64 weight configs)",
     xlim = c(0.67, 0.84), ylim = c(0.68, 0.92))
# Pareto frontier
pareto_idx <- wopt_trade$pareto
points(wopt$rt_corr[pareto_idx], wopt$prev_corr[pareto_idx],
       pch = 17, cex = 1.2, col = "red")
# Ablation points
points(abl$rt_corr, abl$prev_corr, pch = 15, cex = 1.0, col = "blue")
text(abl$rt_corr + 0.003, abl$prev_corr, abl$sources, cex = 0.45, adj = 0)
# Main fit
main_row <- wopt[wopt$w_cases == 1 & wopt$w_hosp == 1 & wopt$w_icu == 1 & wopt$w_radar == 1, ]
if (nrow(main_row) > 0)
  points(main_row$rt_corr, main_row$prev_corr, pch = 8, cex = 2, col = "black", lwd = 2)
legend("bottomleft", c("Weight configs", "Pareto-optimal", "Ablation configs", "Main fit (1,1,1,1)"),
       col = c("gray60", "red", "blue", "black"), pch = c(16, 17, 15, 8), cex = 0.65)

# Rt heatmap (aggregated best Rt over hosp/icu)
wc_vals <- sort(unique(wopt$w_cases))
wr_vals <- sort(unique(wopt$w_radar))
mat_rt <- matrix(NA, nrow = length(wc_vals), ncol = length(wr_vals))
for (i in seq_along(wc_vals)) {
  for (j in seq_along(wr_vals)) {
    rows <- wopt[wopt$w_cases == wc_vals[i] & wopt$w_radar == wr_vals[j], ]
    if (nrow(rows) > 0) mat_rt[i, j] <- max(rows$rt_corr)
  }
}
image(wc_vals, wr_vals, mat_rt,
      xlab = expression(w[cases]), ylab = expression(w[radar]),
      main = "D. Best Rt correlation (over w_hosp, w_icu)",
      col = hcl.colors(20, "Greens3", rev = TRUE))
for (i in seq_along(wc_vals))
  for (j in seq_along(wr_vals))
    if (!is.na(mat_rt[i, j]))
      text(wc_vals[i], wr_vals[j], sprintf("%.3f", mat_rt[i, j]), cex = 0.7)

# Page 4: Ablation bar charts (both metrics side by side)
par(mfrow = c(1, 2), mar = c(4, 10, 3, 1))
abl_prev <- abl[order(abl$prev_corr), ]
bp1 <- barplot(abl_prev$prev_corr, horiz = TRUE, names.arg = abl_prev$sources,
               las = 1, col = "steelblue", xlim = c(0, 1),
               main = "E. Ablation: prevalence corr", xlab = "Correlation")
text(abl_prev$prev_corr - 0.03, bp1, sprintf("%.3f", abl_prev$prev_corr), cex = 0.65)

abl_rt2 <- abl[order(abl$rt_corr), ]
bp2 <- barplot(abl_rt2$rt_corr, horiz = TRUE, names.arg = abl_rt2$sources,
               las = 1, col = "darkgreen", xlim = c(0, 1),
               main = "F. Ablation: Rt corr", xlab = "Correlation")
text(abl_rt2$rt_corr - 0.03, bp2, sprintf("%.3f", abl_rt2$rt_corr), cex = 0.65)

dev.off()
cat("  Saved gp_seir_plots_supplement.pdf (4 pages)\n")

cat("\nDone.\n")
