# =====================================================================
# ge4_debray_meta_regression.R
#
# GE4 True Debray — exploratory meta-regression. Tests whether
# fold-level event rate or mean patient age explain the between-fold
# variation observed in discrimination and calibration.
#
# O:E is not tested here: under true Debray, O:E is exactly 1.0 in
# every fold by construction (a forced identity, not real information),
# so there is no between-fold variation for a moderator to explain.
#
# NOTE: with only K = 7 folds, this analysis is underpowered relative
# to the benchmark in the literature for this kind of test (Debray
# et al. 2019's own worked example required K = 21 studies to detect a
# significant moderator). Four tests are run here (2 measures x 2
# moderators) with no multiple-comparison correction. R-squared values
# should be interpreted cautiously alongside the p-value: a high R²
# accompanying a non-significant p-value can reflect near-zero
# underlying heterogeneity left to explain, rather than a genuinely
# strong moderating relationship, and should not be reported as a
# confirmed effect on its own.
#
# Input:  year1_ge_clean_v3.rds, ge4_debray_v3_all_folds_summary.rds
# Output: ge4_debray_v3_meta_regression_results.rds
# =====================================================================

library(dplyr)
library(metafor)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")
fixed_results <- readRDS("/users/hlskhatt/outputs/ge4_debray_v3_all_folds_summary.rds")

# ---------------------------------------------------------------------
# Fold-level moderators: event rate and mean age within each v3 fold
# ---------------------------------------------------------------------
fold_moderators <- year1 %>%
  filter(!is.na(blakey_outcome_ge4)) %>%
  group_by(iecv_fold_v3) %>%
  summarise(event_rate = mean(blakey_outcome_ge4, na.rm = TRUE),
            mean_age = mean(age, na.rm = TRUE)) %>%
  rename(fold = iecv_fold_v3)

# --- Explicit match()-based reorder after merge(), to guarantee row
#     alignment between fold-level results and fold-level moderators
#     regardless of merge()'s internal sort behaviour ---
merged <- merge(fixed_results, fold_moderators, by = "fold", sort = FALSE)
merged <- merged[match(fixed_results$fold, merged$fold), ]

logit_cstat <- qlogis(merged$c_stat)
logit_cstat_se <- merged$c_stat_se / (merged$c_stat * (1 - merged$c_stat))

# ---------------------------------------------------------------------
# Meta-regression models. Wrapped in tryCatch so that a single failed
# fit returns NA rather than silently shifting subsequent values.
# ---------------------------------------------------------------------
run_metareg <- function(yi, sei, moderator) {
  tryCatch(rma(yi = yi, sei = sei, mods = ~ moderator, method = "REML"),
           error = function(e) NULL)
}

fit_cstat_event <- run_metareg(logit_cstat, logit_cstat_se, merged$event_rate)
fit_slope_event <- run_metareg(merged$slope, merged$slope_se, merged$event_rate)
fit_cstat_age   <- run_metareg(logit_cstat, logit_cstat_se, merged$mean_age)
fit_slope_age   <- run_metareg(merged$slope, merged$slope_se, merged$mean_age)

get_pval <- function(fit) if (is.null(fit)) NA else fit$pval[2]
get_r2   <- function(fit) if (is.null(fit)) NA else fit$R2

summary_table <- data.frame(
  outcome = rep(c("C-statistic", "Calibration slope"), 2),
  moderator = c(rep("event_rate", 2), rep("mean_age", 2)),
  p_value = c(get_pval(fit_cstat_event), get_pval(fit_slope_event),
              get_pval(fit_cstat_age), get_pval(fit_slope_age)),
  r_squared_pct = c(get_r2(fit_cstat_event), get_r2(fit_slope_event),
                     get_r2(fit_cstat_age), get_r2(fit_slope_age)),
  converged = c(!is.null(fit_cstat_event), !is.null(fit_slope_event),
                !is.null(fit_cstat_age), !is.null(fit_slope_age))
)

cat("\n=== TRUE DEBRAY GE4 — META-REGRESSION (K=7, 4 tests, O:E excluded) ===\n")
print(summary_table)

cat("\n--- Fold-level moderators used ---\n")
print(fold_moderators)

saveRDS(list(fold_moderators = fold_moderators, summary_table = summary_table),
        "/users/hlskhatt/outputs/ge4_debray_v3_meta_regression_results.rds")
