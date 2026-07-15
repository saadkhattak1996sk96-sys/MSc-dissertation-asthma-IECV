# =====================================================================
# pool_ge2_debray_sensitivity_reml_uniform.R
#
# GE2 True Debray — sensitivity comparison. Tests whether the pooled
# result depends on the specific pooling method chosen, by comparing
# the primary Bayesian (half-Student-t) result against two alternative,
# equally defensible approaches:
#
#   - Uniform(0, 2) prior on tau — Debray et al. (2019)'s own documented
#     alternative, chosen as the "loosest reasonable" bookend.
#   - REML (frequentist random-effects meta-analysis).
#
# Note: the uniform prior's tau.max = 2 is not the same quantity as the
# primary model's tau.max = 10. In the primary (half-Student-t) model,
# tau.max is an outer safety ceiling on a distribution already
# concentrated near zero. In the uniform model, tau.max IS the entire
# prior (every value between 0 and tau.max is equally likely), so it
# must match the specific value in the cited source rather than the
# primary model's ceiling.
#
# Input:  ge2_debray_v3_all_folds_summary.rds,
#         ge2_debray_v3_pooled_bayesian_results.rds
# Output: ge2_debray_v3_sensitivity_comparison.rds
# =====================================================================

library(metafor)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_debray_v3_all_folds_summary.rds")

# ---------------------------------------------------------------------
# REML (frequentist random-effects meta-analysis)
# ---------------------------------------------------------------------
logit_cstat <- qlogis(fixed_results$c_stat)
logit_cstat_se <- fixed_results$c_stat_se / (fixed_results$c_stat * (1 - fixed_results$c_stat))
reml_cstat <- rma(yi = logit_cstat, sei = logit_cstat_se, method = "REML")
pooled_cstat_reml <- plogis(reml_cstat$b[1])

reml_slope <- rma(yi = fixed_results$slope, sei = fixed_results$slope_se, method = "REML")

# ---------------------------------------------------------------------
# Uniform(0, 2) prior Bayesian sensitivity model
# ---------------------------------------------------------------------
generic_model <- "
model {
  for (i in 1:k) {
    theta[i] ~ dnorm(theta_true[i], theta.prec[i])
    theta.prec[i] <- 1 / (theta.se[i] * theta.se[i])
    theta_true[i] ~ dnorm(mu, tau.prec)
  }
  pred ~ dnorm(mu, tau.prec)
  tau.prec <- 1 / (bsTau * bsTau)
  bsTau ~ dunif(0, tau.max)
  mu ~ dnorm(hp.mu.mean, hp.mu.prec)
  hp.mu.prec <- 1 / hp.mu.var
}
"

cstat_data_unif <- list(k = nrow(fixed_results), theta = logit_cstat, theta.se = logit_cstat_se,
                         tau.max = 2, hp.mu.mean = 0, hp.mu.var = 1000)
set.seed(123)
fit_cstat_unif <- run.jags(model = generic_model, data = cstat_data_unif, monitor = c("mu", "bsTau", "pred"),
                            n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
cstat_unif_summary <- summary(fit_cstat_unif)
pooled_cstat_unif <- plogis(cstat_unif_summary["mu", "Median"])

slope_data_unif <- list(k = nrow(fixed_results), theta = fixed_results$slope, theta.se = fixed_results$slope_se,
                         tau.max = 2, hp.mu.mean = 0, hp.mu.var = 1000)
set.seed(123)
fit_slope_unif <- run.jags(model = generic_model, data = slope_data_unif, monitor = c("mu", "bsTau", "pred"),
                            n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
slope_unif_summary <- summary(fit_slope_unif)

# ---------------------------------------------------------------------
# Assemble and save the three-way comparison
# ---------------------------------------------------------------------
primary <- readRDS("/users/hlskhatt/outputs/ge2_debray_v3_pooled_bayesian_results.rds")$final_summary

comparison <- data.frame(
  measure = c("C-statistic", "Calibration slope"),
  half_student_t_primary = primary$pooled_estimate,
  uniform_prior = c(round(pooled_cstat_unif, 4), round(slope_unif_summary["mu", "Median"], 4)),
  reml = c(round(pooled_cstat_reml, 4), round(reml_slope$b[1], 4)),
  tau2_half_student_t = primary$tau2,
  tau2_uniform = c(round(cstat_unif_summary["bsTau", "Median"]^2, 5),
                   round(slope_unif_summary["bsTau", "Median"]^2, 5)),
  tau2_reml = c(round(reml_cstat$tau2, 5), round(reml_slope$tau2, 5))
)

cat("\n=== TRUE DEBRAY GE2 — SENSITIVITY COMPARISON ===\n")
print(comparison)

saveRDS(comparison, "/users/hlskhatt/outputs/ge2_debray_v3_sensitivity_comparison.rds")
