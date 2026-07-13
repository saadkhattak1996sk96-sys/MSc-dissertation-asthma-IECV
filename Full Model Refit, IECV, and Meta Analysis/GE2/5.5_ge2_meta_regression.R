library(dplyr)
library(metafor)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean.rds")
fixed_results <- readRDS("/users/hlskhatt/outputs/ge2_all_folds_summary_FIXED.rds")

fold_moderators <- year1 %>%
  filter(!is.na(blakey_outcome_ge2)) %>%
  group_by(iecv_fold_v2) %>%
  summarise(
    event_rate = mean(blakey_outcome_ge2, na.rm = TRUE),
    mean_age = mean(age, na.rm = TRUE)
  ) %>%
  rename(fold = iecv_fold_v2)

merged <- merge(fixed_results, fold_moderators, by = "fold")

logit_cstat <- qlogis(merged$c_stat)
logit_cstat_se <- merged$c_stat_se / (merged$c_stat * (1 - merged$c_stat))
log_oe <- log(merged$oe_mean)
log_oe_se <- merged$oe_se / merged$oe_mean

run_metareg <- function(yi, sei, moderator, mod_name, outcome_name) {
  cat("\n---", outcome_name, "~", mod_name, "---\n")
  fit <- tryCatch(
    rma(yi = yi, sei = sei, mods = ~ moderator, method = "REML"),
    error = function(e) { cat("Model failed to fit:", conditionMessage(e), "\n"); return(NULL) }
  )
  if (!is.null(fit)) {
    cat("Moderator coefficient:", round(fit$b[2], 5), "| p-value:", round(fit$pval[2], 4), "\n")
    cat("R-squared (% heterogeneity explained):", round(fit$R2, 2), "%\n")
    cat("Residual Tau2:", round(fit$tau2, 6), "\n")
  }
  fit
}

fit_cstat_event <- run_metareg(logit_cstat, logit_cstat_se, merged$event_rate, "event_rate", "C-statistic (logit)")
fit_slope_event <- run_metareg(merged$slope, merged$slope_se, merged$event_rate, "event_rate", "Calibration slope")
fit_oe_event <- run_metareg(log_oe, log_oe_se, merged$event_rate, "event_rate", "O:E (log)")

fit_cstat_age <- run_metareg(logit_cstat, logit_cstat_se, merged$mean_age, "mean_age", "C-statistic (logit)")
fit_slope_age <- run_metareg(merged$slope, merged$slope_se, merged$mean_age, "mean_age", "Calibration slope")
fit_oe_age <- run_metareg(log_oe, log_oe_se, merged$mean_age, "mean_age", "O:E (log)")

summary_table <- data.frame(
  outcome = rep(c("C-statistic", "Calibration slope", "O:E"), 2),
  moderator = c(rep("event_rate", 3), rep("mean_age", 3)),
  p_value = c(fit_cstat_event$pval[2], fit_slope_event$pval[2], fit_oe_event$pval[2],
              fit_cstat_age$pval[2], fit_slope_age$pval[2], fit_oe_age$pval[2]),
  r_squared_pct = c(fit_cstat_event$R2, fit_slope_event$R2, fit_oe_event$R2,
                     fit_cstat_age$R2, fit_slope_age$R2, fit_oe_age$R2)
)
print(summary_table, row.names = FALSE)

cat("\nNOTE: K=7, meta-regression is explicitly underpowered per Debray et al. (2019).\n")
cat("Any 'significant' C-statistic moderator here is very likely driven by Fold 1 being\n")
cat("simultaneously the youngest and highest-event-rate subset, not two independent findings.\n")

saveRDS(list(fold_moderators = fold_moderators, summary_table = summary_table),
        "/users/hlskhatt/outputs/ge2_meta_regression_results.rds")
