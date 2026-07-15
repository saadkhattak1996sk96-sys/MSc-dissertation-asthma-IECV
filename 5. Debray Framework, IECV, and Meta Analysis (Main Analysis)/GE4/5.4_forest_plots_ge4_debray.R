# =====================================================================
# forest_plots_ge4_debray.R
#
# GE4 True Debray — forest plots for C-statistic and calibration slope
# by fold, alongside the pooled (Bayesian, half-Student-t) estimate.
#
# Confidence intervals for the C-statistic are built on the logit
# scale and back-transformed, matching the scale used for cross-fold
# pooling. Calibration slope is left untransformed, per the pooling
# method used (Wood 2015).
#
# No O:E forest plot is produced for true Debray: O:E is exactly 1.0
# in every fold by construction (see fit_ge4_debray_all_folds.R), so a
# plot of it would carry no information.
#
# Input:  ge4_debray_v3_all_folds_summary.rds,
#         ge4_debray_v3_pooled_bayesian_results.rds
# Output: GE4_debray_v3_forest_cstat.png, GE4_debray_v3_forest_slope.png
# =====================================================================

library(ggplot2)
library(dplyr)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge4_debray_v3_all_folds_summary.rds")
pooled <- readRDS("/users/hlskhatt/outputs/ge4_debray_v3_pooled_bayesian_results.rds")
output_dir <- "/users/hlskhatt/outputs"

make_forest_plot <- function(data, pooled_est, pooled_lb, pooled_ub, title, xlab) {
  plot_data <- data.frame(
    label = c(as.character(data$fold), "Pooled (Bayesian, half-Student-t)"),
    estimate = c(data$estimate, pooled_est),
    lower = c(data$lower, pooled_lb),
    upper = c(data$upper, pooled_ub),
    type = c(rep("Fold", nrow(data)), "Pooled")
  )
  plot_data$label <- factor(plot_data$label, levels = rev(plot_data$label))

  ggplot(plot_data, aes(x = estimate, y = label, color = type)) +
    geom_point(aes(shape = type, size = type)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    scale_shape_manual(values = c("Fold" = 15, "Pooled" = 18)) +
    scale_size_manual(values = c("Fold" = 3, "Pooled" = 5)) +
    scale_color_manual(values = c("Fold" = "black", "Pooled" = "darkred")) +
    labs(title = title, x = xlab, y = NULL) +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13))
}

# ---------------------------------------------------------------------
# C-statistic: per-fold CIs built on the logit scale, back-transformed
# ---------------------------------------------------------------------
cstat_data <- data.frame(
  fold = fixed_results$fold,
  estimate = qlogis(fixed_results$c_stat),
  lower = qlogis(fixed_results$c_stat) - 1.96 * (fixed_results$c_stat_se / (fixed_results$c_stat * (1 - fixed_results$c_stat))),
  upper = qlogis(fixed_results$c_stat) + 1.96 * (fixed_results$c_stat_se / (fixed_results$c_stat * (1 - fixed_results$c_stat)))
)
cstat_data$estimate <- plogis(cstat_data$estimate)
cstat_data$lower <- plogis(cstat_data$lower)
cstat_data$upper <- plogis(cstat_data$upper)

p_cstat <- make_forest_plot(cstat_data, pooled$fit_cstat$est, pooled$fit_cstat$ci.lb, pooled$fit_cstat$ci.ub,
                             "GE4 True Debray (v3 folds) — C-statistic by Fold (IECV)", "C-statistic")
ggsave(file.path(output_dir, "GE4_debray_v3_forest_cstat.png"), p_cstat, width = 8, height = 5, dpi = 300)

# ---------------------------------------------------------------------
# Calibration slope: untransformed
# ---------------------------------------------------------------------
slope_data <- data.frame(
  fold = fixed_results$fold,
  estimate = fixed_results$slope,
  lower = fixed_results$slope - 1.96 * fixed_results$slope_se,
  upper = fixed_results$slope + 1.96 * fixed_results$slope_se
)
p_slope <- make_forest_plot(slope_data, pooled$slope_summary["mu", "Median"],
                             pooled$slope_summary["mu", "Lower95"], pooled$slope_summary["mu", "Upper95"],
                             "GE4 True Debray (v3 folds) — Calibration Slope by Fold (IECV)", "Calibration Slope")
p_slope <- p_slope + geom_vline(xintercept = 1, linetype = "dashed", color = "grey50")
ggsave(file.path(output_dir, "GE4_debray_v3_forest_slope.png"), p_slope, width = 8, height = 5, dpi = 300)

cat("Saved: GE4_debray_v3_forest_cstat.png, GE4_debray_v3_forest_slope.png\n")
