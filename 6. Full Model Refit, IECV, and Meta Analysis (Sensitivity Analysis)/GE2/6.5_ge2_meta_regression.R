# =============================================================================
# GE2 FULL MODEL REFIT — META-REGRESSION
# =============================================================================
# Purpose:
#   Explores whether fold-level heterogeneity in discrimination or
#   calibration is associated with two candidate fold-level moderators:
#   event rate and mean patient age. Exploratory only — with 7 folds this
#   analysis is underpowered relative to standard meta-regression guidance.
#
# Note:
#   Each of the six regressions (2 moderators x 3 outcome measures) is
#   wrapped so that a failed fit is reported explicitly (via a `converged`
#   flag and NA values) rather than silently dropped from the results table.
# =============================================================================

library(dplyr)
library(metafor)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")
fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_v3_all_folds_summary.rds")

fold_moderators <- year1 %>%
  filter(!is.na(blakey_outcome_ge2)) %>%
  group_by(iecv_fold_v3) %>%
  summarise(event_rate = mean(blakey_outcome_ge2, na.rm = TRUE),
            mean_age = mean(age, na.rm = TRUE)) %>%
  rename(fold = iecv_fold_v3)

merged <- merge(fixed_results, fold_moderators, by = "fold", sort = FALSE)
merged <- merged[match(fixed_results$fold, merged$fold), ]

logit_cstat <- qlogis(merged$c_stat)
logit_cstat_se <- merged$c_stat_se / (merged$c_stat * (1 - merged$c_stat))
log_oe <- log(merged$oe_mean)
log_oe_se <- merged$oe_se / merged$oe_mean

run_metareg <- function(yi, sei, moderator) {
  tryCatch(rma(yi = yi, sei = sei, mods = ~ moderator, method = "REML"),
           error = function(e) NULL)
}

get_pval <- function(fit) {
  if (is.null(fit)) return(NA)
  tryCatch(fit$pval[2], error = function(e) NA)
}

get_r2 <- function(fit) {
  if (is.null(fit)) return(NA)
  tryCatch(fit$R2, error = function(e) NA)
}

get_converged <- function(fit) {
  !is.null(fit)
}

fit_cstat_event <- run_metareg(logit_cstat, logit_cstat_se, merged$event_rate)
fit_slope_event <- run_metareg(merged$slope, merged$slope_se, merged$event_rate)
fit_oe_event <- run_metareg(log_oe, log_oe_se, merged$event_rate)
fit_cstat_age <- run_metareg(logit_cstat, logit_cstat_se, merged$mean_age)
fit_slope_age <- run_metareg(merged$slope, merged$slope_se, merged$mean_age)
fit_oe_age <- run_metareg(log_oe, log_oe_se, merged$mean_age)

all_fits <- list(fit_cstat_event, fit_slope_event, fit_oe_event,
                  fit_cstat_age, fit_slope_age, fit_oe_age)

summary_table <- data.frame(
  outcome = rep(c("C-statistic", "Calibration slope", "O:E"), 2),
  moderator = c(rep("event_rate", 3), rep("mean_age", 3)),
  p_value = sapply(all_fits, get_pval),
  r_squared_pct = sapply(all_fits, get_r2),
  converged = sapply(all_fits, get_converged)
)

saveRDS(list(fold_moderators = fold_moderators, summary_table = summary_table),
        "/users/hlskhatt/outputs/ge2_v3_full_refit_meta_regression_results.rds")
