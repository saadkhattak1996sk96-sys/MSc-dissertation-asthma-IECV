# =====================================================================
# ge4_meta_regression.R
# Meta-regression: does fold-level event rate or mean age explain
# variation in C-statistic, calibration slope, or O:E across folds?
#
# NOTE: K=7 folds — this analysis is underpowered relative to Debray's
# own K=21 benchmark for this kind of test. Six tests run here (3
# measures x 2 moderators) with no multiple-comparison correction —
# treat any "significant" result with appropriate caution in write-up.
# =====================================================================

library(dplyr)
library(metafor)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")
fixed_results <- readRDS("/users/hlskhatt/outputs/ge4_v3_all_folds_summary.rds")

fold_moderators <- year1 %>%
  filter(!is.na(blakey_outcome_ge4)) %>%
  group_by(iecv_fold_v3) %>%
  summarise(event_rate = mean(blakey_outcome_ge4, na.rm = TRUE),
            mean_age = mean(age, na.rm = TRUE)) %>%
  rename(fold = iecv_fold_v3)

# --- sort = FALSE + explicit match() reorder: avoids relying on merge()'s
#     default alphabetical sort silently happening to match fold order ---
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

fit_cstat_event <- run_metareg(logit_cstat, logit_cstat_se, merged$event_rate)
fit_slope_event <- run_metareg(merged$slope, merged$slope_se, merged$event_rate)
fit_oe_event <- run_metareg(log_oe, log_oe_se, merged$event_rate)
fit_cstat_age <- run_metareg(logit_cstat, logit_cstat_se, merged$mean_age)
fit_slope_age <- run_metareg(merged$slope, merged$slope_se, merged$mean_age)
fit_oe_age <- run_metareg(log_oe, log_oe_se, merged$mean_age)

# --- Defensive helpers: return NA + flag rather than crashing if any
#     rma() call failed (returns NULL from the tryCatch above) ---
get_pval <- function(fit) if (is.null(fit)) NA else fit$pval[2]
get_r2 <- function(fit) if (is.null(fit)) NA else fit$R2

summary_table <- data.frame(
  outcome = rep(c("C-statistic", "Calibration slope", "O:E"), 2),
  moderator = c(rep("event_rate", 3), rep("mean_age", 3)),
  p_value = c(get_pval(fit_cstat_event), get_pval(fit_slope_event), get_pval(fit_oe_event),
              get_pval(fit_cstat_age), get_pval(fit_slope_age), get_pval(fit_oe_age)),
  r_squared_pct = c(get_r2(fit_cstat_event), get_r2(fit_slope_event), get_r2(fit_oe_event),
                     get_r2(fit_cstat_age), get_r2(fit_slope_age), get_r2(fit_oe_age)),
  converged = c(!is.null(fit_cstat_event), !is.null(fit_slope_event), !is.null(fit_oe_event),
                !is.null(fit_cstat_age), !is.null(fit_slope_age), !is.null(fit_oe_age))
)

print(summary_table)

saveRDS(list(fold_moderators = fold_moderators, summary_table = summary_table),
        "/users/hlskhatt/outputs/ge4_v3_meta_regression_results.rds")
