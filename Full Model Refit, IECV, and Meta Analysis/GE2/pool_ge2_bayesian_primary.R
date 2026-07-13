library(metafor)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_all_folds_summary_FIXED.rds")

cat("=========== C-STATISTIC: BAYESIAN POOLING (logit-transformed) ===========\n")
fit_cstat <- valmeta(
  cstat    = fixed_results$c_stat,
  cstat.se = fixed_results$c_stat_se,
  slab     = fixed_results$fold,
  method   = "BAYES",
  pars     = list(hp.tau.dist = "dhalft", hp.tau.sigma = 0.5, hp.tau.df = 3, hp.tau.max = 10),
  verbose  = FALSE
)
cat("Pooled C-statistic:", round(fit_cstat$est, 4), "| 95% CrI:", round(fit_cstat$ci.lb, 4), "-", round(fit_cstat$ci.ub, 4), "| Tau2:", round(fit_cstat$tau2, 5), "\n")

cat("\n=========== CALIBRATION SLOPE: BAYESIAN POOLING (hand-built JAGS) ===========\n")
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
  k = nrow(fixed_results), theta = fixed_results$slope, theta.se = fixed_results$slope_se,
  tau.sigma = 0.5, tau.df = 3, tau.max = 10, hp.mu.mean = 0, hp.mu.var = 1000
)
fit_slope <- run.jags(model = slope_model, data = slope_data, monitor = c("mu", "bsTau", "pred"),
                       n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
slope_summary <- summary(fit_slope)
cat("Pooled Slope:", round(slope_summary["mu","Median"], 4), "| 95% CrI:", round(slope_summary["mu","Lower95"], 4), "-", round(slope_summary["mu","Upper95"], 4), "| Tau2:", round(slope_summary["bsTau","Median"]^2, 5), "\n")

cat("\n=========== O:E RATIO: BAYESIAN POOLING (log-transformed JAGS, EXPLORATORY) ===========\n")
log_oe <- log(fixed_results$oe_mean)
log_oe_se <- fixed_results$oe_se / fixed_results$oe_mean
oe_data <- list(
  k = nrow(fixed_results), theta = log_oe, theta.se = log_oe_se,
  tau.sigma = 1.5, tau.df = 3, tau.max = 10, hp.mu.mean = 0, hp.mu.var = 1000
)
fit_oe_jags <- run.jags(model = slope_model, data = oe_data, monitor = c("mu", "bsTau", "pred"),
                         n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
oe_summary <- summary(fit_oe_jags)
cat("Pooled O:E:", round(exp(oe_summary["mu","Median"]), 4), "| 95% CrI:", round(exp(oe_summary["mu","Lower95"]), 4), "-", round(exp(oe_summary["mu","Upper95"]), 4), "| Tau2 (log scale):", round(oe_summary["bsTau","Median"]^2, 5), "\n")

cat("\n\n=========== FINAL SUMMARY, ALL THREE MEASURES ===========\n")
final_summary <- data.frame(
  measure = c("C-statistic", "Calibration slope", "O:E ratio (exploratory)"),
  pooled_estimate = c(round(fit_cstat$est, 4), round(slope_summary["mu","Median"], 4), round(exp(oe_summary["mu","Median"]), 4)),
  ci_lower = c(round(fit_cstat$ci.lb, 4), round(slope_summary["mu","Lower95"], 4), round(exp(oe_summary["mu","Lower95"]), 4)),
  ci_upper = c(round(fit_cstat$ci.ub, 4), round(slope_summary["mu","Upper95"], 4), round(exp(oe_summary["mu","Upper95"]), 4)),
  tau2 = c(round(fit_cstat$tau2, 5), round(slope_summary["bsTau","Median"]^2, 5), round(oe_summary["bsTau","Median"]^2, 5))
)
print(final_summary, row.names = FALSE)

saveRDS(list(fit_cstat = fit_cstat, fit_slope = fit_slope, slope_summary = slope_summary,
             fit_oe_jags = fit_oe_jags, oe_summary = oe_summary, final_summary = final_summary),
        "/users/hlskhatt/outputs/ge2_pooled_bayesian_results_FINAL.rds")
cat("\nSaved: ge2_pooled_bayesian_results_FINAL.rds\n")
