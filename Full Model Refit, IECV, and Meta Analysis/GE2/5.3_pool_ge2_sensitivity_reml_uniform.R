library(metafor)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_all_folds_summary_FIXED.rds")

cat("=========== REML: C-STATISTIC (logit-transformed) ===========\n")
logit_cstat <- qlogis(fixed_results$c_stat)
logit_cstat_se <- fixed_results$c_stat_se / (fixed_results$c_stat * (1 - fixed_results$c_stat))
reml_cstat <- rma(yi = logit_cstat, sei = logit_cstat_se, method = "REML")
pooled_cstat_reml <- plogis(reml_cstat$b[1])
cat("Pooled C-stat (REML):", round(pooled_cstat_reml, 4), "| 95% CI:", round(plogis(reml_cstat$ci.lb), 4), "-", round(plogis(reml_cstat$ci.ub), 4), "| Tau2:", round(reml_cstat$tau2, 5), "\n")

cat("\n=========== REML: CALIBRATION SLOPE (untransformed) ===========\n")
reml_slope <- rma(yi = fixed_results$slope, sei = fixed_results$slope_se, method = "REML")
cat("Pooled Slope (REML):", round(reml_slope$b[1], 4), "| 95% CI:", round(reml_slope$ci.lb, 4), "-", round(reml_slope$ci.ub, 4), "| Tau2:", round(reml_slope$tau2, 5), "\n")

cat("\n=========== REML: O:E RATIO (log-transformed) ===========\n")
log_oe <- log(fixed_results$oe_mean)
log_oe_se <- fixed_results$oe_se / fixed_results$oe_mean
reml_oe <- rma(yi = log_oe, sei = log_oe_se, method = "REML")
pooled_oe_reml <- exp(reml_oe$b[1])
cat("Pooled O:E (REML):", round(pooled_oe_reml, 4), "| 95% CI:", round(exp(reml_oe$ci.lb), 4), "-", round(exp(reml_oe$ci.ub), 4), "| Tau2:", round(reml_oe$tau2, 5), "\n")

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

cat("\n=========== UNIFORM(0,2) PRIOR: C-STATISTIC (logit-transformed) ===========\n")
cstat_data_unif <- list(k = nrow(fixed_results), theta = logit_cstat, theta.se = logit_cstat_se,
                         tau.max = 2, hp.mu.mean = 0, hp.mu.var = 1000)
fit_cstat_unif <- run.jags(model = generic_model, data = cstat_data_unif, monitor = c("mu", "bsTau", "pred"),
                            n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
cstat_unif_summary <- summary(fit_cstat_unif)
pooled_cstat_unif <- plogis(cstat_unif_summary["mu","Median"])
cat("Pooled C-stat (Uniform prior):", round(pooled_cstat_unif, 4),
    "| 95% CrI:", round(plogis(cstat_unif_summary["mu","Lower95"]), 4), "-", round(plogis(cstat_unif_summary["mu","Upper95"]), 4),
    "| Tau2:", round(cstat_unif_summary["bsTau","Median"]^2, 5), "\n")

cat("\n=========== UNIFORM(0,2) PRIOR: CALIBRATION SLOPE ===========\n")
slope_data_unif <- list(k = nrow(fixed_results), theta = fixed_results$slope, theta.se = fixed_results$slope_se,
                         tau.max = 2, hp.mu.mean = 0, hp.mu.var = 1000)
fit_slope_unif <- run.jags(model = generic_model, data = slope_data_unif, monitor = c("mu", "bsTau", "pred"),
                            n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
slope_unif_summary <- summary(fit_slope_unif)
cat("Pooled Slope (Uniform prior):", round(slope_unif_summary["mu","Median"], 4),
    "| 95% CrI:", round(slope_unif_summary["mu","Lower95"], 4), "-", round(slope_unif_summary["mu","Upper95"], 4),
    "| Tau2:", round(slope_unif_summary["bsTau","Median"]^2, 5), "\n")

cat("\n=========== UNIFORM(0,2) PRIOR: O:E RATIO (log-transformed) ===========\n")
oe_data_unif <- list(k = nrow(fixed_results), theta = log_oe, theta.se = log_oe_se,
                      tau.max = 2, hp.mu.mean = 0, hp.mu.var = 1000)
fit_oe_unif <- run.jags(model = generic_model, data = oe_data_unif, monitor = c("mu", "bsTau", "pred"),
                         n.chains = 4, adapt = 1000, burnin = 4000, sample = 10000, method = "rjags")
oe_unif_summary <- summary(fit_oe_unif)
pooled_oe_unif <- exp(oe_unif_summary["mu","Median"])
cat("Pooled O:E (Uniform prior):", round(pooled_oe_unif, 4),
    "| 95% CrI:", round(exp(oe_unif_summary["mu","Lower95"]), 4), "-", round(exp(oe_unif_summary["mu","Upper95"]), 4),
    "| Tau2 (log scale):", round(oe_unif_summary["bsTau","Median"]^2, 5), "\n")

cat("\n\n=========== THREE-WAY COMPARISON: HALF-STUDENT-T (PRIMARY) vs UNIFORM vs REML ===========\n")
primary <- readRDS("/users/hlskhatt/outputs/ge2_pooled_bayesian_results_FINAL.rds")$final_summary
comparison <- data.frame(
  measure = c("C-statistic", "Calibration slope", "O:E ratio (exploratory)"),
  half_student_t_primary = primary$pooled_estimate,
  uniform_prior = c(round(pooled_cstat_unif, 4), round(slope_unif_summary["mu","Median"], 4), round(pooled_oe_unif, 4)),
  reml = c(round(pooled_cstat_reml, 4), round(reml_slope$b[1], 4), round(pooled_oe_reml, 4)),
  tau2_half_student_t = primary$tau2,
  tau2_uniform = c(round(cstat_unif_summary["bsTau","Median"]^2, 5), round(slope_unif_summary["bsTau","Median"]^2, 5), round(oe_unif_summary["bsTau","Median"]^2, 5)),
  tau2_reml = c(round(reml_cstat$tau2, 5), round(reml_slope$tau2, 5), round(reml_oe$tau2, 5))
)
print(comparison, row.names = FALSE)

saveRDS(comparison, "/users/hlskhatt/outputs/ge2_sensitivity_comparison.rds")
cat("\nSaved: ge2_sensitivity_comparison.rds\n")
