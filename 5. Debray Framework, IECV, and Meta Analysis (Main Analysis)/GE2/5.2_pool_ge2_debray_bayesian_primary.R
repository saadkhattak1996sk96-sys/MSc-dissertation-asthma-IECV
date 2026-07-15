# =====================================================================
# pool_ge2_debray_bayesian_primary.R
#
# GE2 True Debray — cross-fold Bayesian random-effects pooling.
# Primary pooling method per the pre-specified analysis plan: one
# justified, citable method per performance measure.
#
# Prior specification (Debray et al. 2019, Supporting Information
# 4.2-4.3): half-Student-t prior on tau, location 0, 3 degrees of
# freedom, truncated at [0, 10]. Scale 0.5 for the logit-transformed
# C-statistic and calibration slope.
#
# O:E ratio is intentionally excluded from pooling here. Under true
# Debray, each fold's O:E is a mathematical identity (~1.0) rather than
# an independent estimate, since the held-out fold's intercept is
# re-estimated from its own observed outcome before evaluation — there
# is no genuine between-fold variation to pool.
#
# Input:  ge2_debray_v3_all_folds_summary.rds (7-fold results table)
# Output: ge2_debray_v3_pooled_bayesian_results.rds
# =====================================================================

library(metamisc)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_debray_v3_all_folds_summary.rds")

# ---------------------------------------------------------------------
# C-statistic: logit-scale Bayesian random-effects pooling via
# metamisc::valmeta(), half-Student-t prior on heterogeneity
# ---------------------------------------------------------------------
fit_cstat <- valmeta(
  cstat    = fixed_results$c_stat,
  cstat.se = fixed_results$c_stat_se,
  slab     = fixed_results$fold,
  method   = "BAYES",
  pars     = list(hp.tau.dist = "dhalft", hp.tau.sigma = 0.5, hp.tau.df = 3, hp.tau.max = 10),
  verbose  = FALSE
)

# ---------------------------------------------------------------------
# Calibration slope: untransformed Bayesian random-effects pooling,
# same prior specification, via a hand-built JAGS model (valmeta does
# not support calibration slope pooling directly)
# ---------------------------------------------------------------------
slope_model <- "
model {
  for (i in 1:k) {
    theta[i] ~ dnorm(theta_true[i], theta.prec[i])
    theta.prec[i] <- 1 / (theta.se[i] * theta.se[i])
    theta_true[i] ~ dnorm(mu, tau.prec)
  }
  pred ~ dnorm(mu, tau.prec)
  tau.prec <- 1 / (bsTau * bsTau)
  bsTau ~ dt(0, tau.dt.prec, tau.df) T(0, tau.max)
  tau.dt.prec <- 1 / (tau.sigma * tau.sigma)
  mu ~ dnorm(hp.mu.mean, hp.mu.prec)
  hp.mu.prec <- 1 / hp.mu.var
}
"

slope_data <- list(
  k = nrow(fixed_results),
  theta = fixed_results$slope,
  theta.se = fixed_results$slope_se,
  tau.sigma = 0.5, tau.df = 3, tau.max = 10,
  hp.mu.mean = 0, hp.mu.var = 1000
)

set.seed(123)
fit_slope <- run.jags(
  model = slope_model, data = slope_data,
  monitor = c("mu", "bsTau", "pred"), n.chains = 4,
  adapt = 1000, burnin = 4000, sample = 10000, method = "rjags"
)
slope_summary <- summary(fit_slope)

# ---------------------------------------------------------------------
# Assemble and save the pooled summary
# ---------------------------------------------------------------------
final_summary <- data.frame(
  measure = c("C-statistic", "Calibration slope"),
  pooled_estimate = c(round(fit_cstat$est, 4), round(slope_summary["mu", "Median"], 4)),
  ci_lower = c(round(fit_cstat$ci.lb, 4), round(slope_summary["mu", "Lower95"], 4)),
  ci_upper = c(round(fit_cstat$ci.ub, 4), round(slope_summary["mu", "Upper95"], 4)),
  tau2 = c(round(fit_cstat$tau2, 5), round(slope_summary["bsTau", "Median"]^2, 5))
)

cat("\n=== TRUE DEBRAY GE2 — POOLED RESULTS ===\n")
print(final_summary)

saveRDS(list(fit_cstat = fit_cstat, fit_slope = fit_slope, slope_summary = slope_summary,
             final_summary = final_summary),
        "/users/hlskhatt/outputs/ge2_debray_v3_pooled_bayesian_results.rds")
