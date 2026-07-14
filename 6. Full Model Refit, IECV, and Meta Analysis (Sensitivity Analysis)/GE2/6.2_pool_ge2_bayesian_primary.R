library(metamisc)
library(runjags)

fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_v3_all_folds_summary.rds")

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

fit_slope <- run.jags(
  model = slope_model, data = slope_data,
  monitor = c("mu", "bsTau", "pred"), n.chains = 4,
  adapt = 1000, burnin = 4000, sample = 10000, method = "rjags"
)
slope_summary <- summary(fit_slope)

fit_oe <- valmeta(
  measure = "OE",
  OE      = fixed_results$oe_mean,
  OE.se   = fixed_results$oe_se,
  slab    = fixed_results$fold,
  method  = "BAYES",
  pars    = list(hp.tau.dist = "dhalft", hp.tau.sigma = 1.5, hp.tau.df = 3, hp.tau.max = 10),
  verbose = FALSE
)

oe_tau2 <- safe_tau2(fit_oe)

final_summary <- data.frame(
  measure = c("C-statistic", "Calibration slope", "O:E ratio (exploratory)"),
  pooled_estimate = c(round(fit_cstat$est, 4), round(slope_summary["mu","Median"], 4), round(fit_oe$est, 4)),
  ci_lower = c(round(fit_cstat$ci.lb, 4), round(slope_summary["mu","Lower95"], 4), round(fit_oe$ci.lb, 4)),
  ci_upper = c(round(fit_cstat$ci.ub, 4), round(slope_summary["mu","Upper95"], 4), round(fit_oe$ci.ub, 4)),
  tau2 = c(round(safe_tau2(fit_cstat), 5), round(slope_summary["bsTau","Median"]^2, 5), round(oe_tau2, 5))
)

saveRDS(list(fit_cstat = fit_cstat, fit_slope = fit_slope, slope_summary = slope_summary,
             fit_oe = fit_oe, final_summary = final_summary),
        "/users/hlskhatt/outputs/ge2_v3_pooled_bayesian_results.rds")

print(final_summary)
