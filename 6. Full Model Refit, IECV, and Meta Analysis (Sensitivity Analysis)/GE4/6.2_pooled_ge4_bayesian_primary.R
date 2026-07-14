# =====================================================================
# pool_ge4_bayesian_primary.R
# Cross-fold Bayesian random-effects pooling, half-Student-t prior (primary
# method per Laura's constraint: one justified, citable method per measure).
#
# NOTE: O:E is pooled via a hand-rolled JAGS model (log-transformed), NOT
# via valmeta()'s built-in BAYES method for O:E. valmeta()'s internal
# defaults failed to converge for this outcome (psrf = 1.06 on bsTau,
# above the 1.05 threshold) and do not expose adapt/burnin/sample for
# tuning. The hand-rolled model gives full control over MCMC settings.
# =====================================================================

library(metamisc)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge4_v3_all_folds_summary.rds")

# --- C-statistic: valmeta's BAYES method, converged cleanly at defaults ---
fit_cstat <- valmeta(
  cstat    = fixed_results$c_stat,
  cstat.se = fixed_results$c_stat_se,
  slab     = fixed_results$fold,
  method   = "BAYES",
  pars     = list(hp.tau.dist = "dhalft", hp.tau.sigma = 0.5, hp.tau.df = 3, hp.tau.max = 10),
  verbose  = FALSE
)

safe_tau2 <- function(fit) {
  val <- tryCatch(fit$tau2, error = function(e) NA)
  val <- suppressWarnings(as.numeric(val))
  if (length(val) == 1 && !is.na(val)) return(val)

  val <- tryCatch(fit$tau, error = function(e) NA)
  val <- suppressWarnings(as.numeric(val))
  if (length(val) == 1 && !is.na(val)) return(val^2)

  val <- tryCatch(fit$fit$summaries["bsTau", "Median"], error = function(e) NA)
  val <- suppressWarnings(as.numeric(val))
  if (length(val) == 1 && !is.na(val)) return(val^2)

  NA
}

# --- Generic random-effects model, used for both slope (raw scale) and
#     O:E (log scale) — gives explicit control over adapt/burnin/sample ---
generic_model <- "
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

# --- Calibration slope: plain Rubin's-rules-consistent random-effects
#     model, untransformed (Wood 2015 justification) ---
slope_data <- list(
  k = nrow(fixed_results),
  theta = fixed_results$slope,
  theta.se = fixed_results$slope_se,
  tau.sigma = 0.5, tau.df = 3, tau.max = 10,
  hp.mu.mean = 0, hp.mu.var = 1000
)

set.seed(123)
fit_slope <- run.jags(
  model = generic_model, data = slope_data,
  monitor = c("mu", "bsTau", "pred"), n.chains = 4,
  adapt = 1000, burnin = 4000, sample = 10000, method = "rjags"
)
slope_summary <- summary(fit_slope)

# --- O:E: log-transformed, hand-rolled model with increased
#     adapt/burnin/sample specifically to clear the convergence threshold
#     valmeta()'s defaults could not reach ---
log_oe <- log(fixed_results$oe_mean)
log_oe_se <- fixed_results$oe_se / fixed_results$oe_mean

oe_data <- list(
  k = nrow(fixed_results),
  theta = log_oe,
  theta.se = log_oe_se,
  tau.sigma = 1.5, tau.df = 3, tau.max = 10,
  hp.mu.mean = 0, hp.mu.var = 1000
)

set.seed(123)
fit_oe <- run.jags(
  model = generic_model, data = oe_data,
  monitor = c("mu", "bsTau", "pred"), n.chains = 4,
  adapt = 2000, burnin = 10000, sample = 30000, method = "rjags"
)
oe_summary <- summary(fit_oe)

# --- Explicit convergence check, printed rather than left to a warning ---
cat("\n=== O:E convergence check (psrf, want < 1.05 for all rows) ===\n")
print(fit_oe$psrf$psrf)

oe_tau2 <- oe_summary["bsTau", "Median"]^2

final_summary <- data.frame(
  measure = c("C-statistic", "Calibration slope", "O:E ratio (exploratory)"),
  pooled_estimate = c(round(fit_cstat$est, 4),
                      round(slope_summary["mu","Median"], 4),
                      round(exp(oe_summary["mu","Median"]), 4)),
  ci_lower = c(round(fit_cstat$ci.lb, 4),
               round(slope_summary["mu","Lower95"], 4),
               round(exp(oe_summary["mu","Lower95"]), 4)),
  ci_upper = c(round(fit_cstat$ci.ub, 4),
               round(slope_summary["mu","Upper95"], 4),
               round(exp(oe_summary["mu","Upper95"]), 4)),
  tau2 = c(round(safe_tau2(fit_cstat), 5),
           round(slope_summary["bsTau","Median"]^2, 5),
           round(oe_tau2, 5))
)

saveRDS(list(fit_cstat = fit_cstat, fit_slope = fit_slope, slope_summary = slope_summary,
             fit_oe = fit_oe, oe_summary = oe_summary, final_summary = final_summary),
        "/users/hlskhatt/outputs/ge4_v3_pooled_bayesian_results.rds")

print(final_summary)
