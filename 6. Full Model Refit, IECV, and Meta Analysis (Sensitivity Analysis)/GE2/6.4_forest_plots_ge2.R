# =============================================================================
# GE2 FULL MODEL REFIT — FOREST PLOTS
# =============================================================================
# Purpose:
#   Produces forest plots for C-statistic, calibration slope, and O:E across
#   the 7 folds plus the pooled estimate.
#
# Note on confidence interval construction:
#   C-statistic and O:E are pooled on transformed scales (logit and log
#   respectively), so their fold-level confidence intervals are built on the
#   same transformed scale via the delta method, then back-transformed for
#   display — keeping the plotted intervals consistent with how pooling was
#   actually performed. Calibration slope is pooled untransformed, so its
#   intervals are built directly on the raw scale.
# =============================================================================

library(ggplot2)
library(dplyr)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_v3_all_folds_summary.rds")
pooled <- readRDS("/users/hlskhatt/outputs/ge2_v3_full_refit_pooled_results.rds")
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

# C-statistic — logit-scale CI, back-transformed
logit_est <- qlogis(fixed_results$c_stat)
logit_se <- fixed_results$c_stat_se / (fixed_results$c_stat * (1 - fixed_results$c_stat))
cstat_data <- data.frame(
  fold = fixed_results$fold,
  estimate = fixed_results$c_stat,
  lower = plogis(logit_est - 1.96 * logit_se),
  upper = plogis(logit_est + 1.96 * logit_se)
)
p_cstat <- make_forest_plot(cstat_data, pooled$fit_cstat$est, pooled$fit_cstat$ci.lb, pooled$fit_cstat$ci.ub,
                             "GE2 Full Refit (v3 folds) — C-statistic by Fold (IECV)", "C-statistic")
ggsave(file.path(output_dir, "GE2_full_refit_forest_cstat.png"), p_cstat, width = 8, height = 5, dpi = 300)

# Calibration slope — untransformed CI (pooled untransformed by design)
slope_data <- data.frame(
  fold = fixed_results$fold,
  estimate = fixed_results$slope,
  lower = fixed_results$slope - 1.96 * fixed_results$slope_se,
  upper = fixed_results$slope + 1.96 * fixed_results$slope_se
)
p_slope <- make_forest_plot(slope_data, pooled$slope_summary["mu","Median"], pooled$slope_summary["mu","Lower95"], pooled$slope_summary["mu","Upper95"],
                             "GE2 Full Refit (v3 folds) — Calibration Slope by Fold (IECV)", "Calibration Slope")
p_slope <- p_slope + geom_vline(xintercept = 1, linetype = "dashed", color = "grey50")
ggsave(file.path(output_dir, "GE2_full_refit_forest_slope.png"), p_slope, width = 8, height = 5, dpi = 300)

# O:E — log-scale CI, back-transformed
log_est <- log(fixed_results$oe_mean)
log_se <- fixed_results$oe_se / fixed_results$oe_mean
oe_data <- data.frame(
  fold = fixed_results$fold,
  estimate = fixed_results$oe_mean,
  lower = exp(log_est - 1.96 * log_se),
  upper = exp(log_est + 1.96 * log_se)
)
p_oe <- make_forest_plot(oe_data, pooled$fit_oe$est, pooled$fit_oe$ci.lb, pooled$fit_oe$ci.ub,
                          "GE2 Full Refit (v3 folds) — O:E Ratio by Fold (IECV, exploratory)", "O:E Ratio")
p_oe <- p_oe + geom_vline(xintercept = 1, linetype = "dashed", color = "grey50")
ggsave(file.path(output_dir, "GE2_full_refit_forest_oe.png"), p_oe, width = 8, height = 5, dpi = 300)
