# =====================================================================
# fit_ge4_debray_final_model.R
#
# GE4 True Debray — FINAL COMBINED MODEL.
#
# All 7 IECV rounds in fit_ge4_debray_all_folds.R serve as a stress
# test of the method, confirming Debray's stability assumption holds
# for this predictor set. This script builds the actual model the
# dissertation delivers: one shared coefficient set fitted using every
# patient across all 7 eras simultaneously, with no fold held out,
# retaining one stratified intercept per era.
#
# A stability check is included: the final model's shared coefficients
# are compared against the average of the 7 development-fold
# coefficients already estimated in fit_ge4_debray_all_folds.R. Close
# agreement confirms the final model is a stable summary of what the
# 7-fold IECV already demonstrated, not an artefact of combining all
# the data together.
#
# Full refit does not require an equivalent final-model build (see
# project decision log): it is a sensitivity analysis testing whether
# letting coefficients vary freely by era changes the answer, a
# question already answered by the fold-level comparison. Only true
# Debray, as the primary analysis, needs a hand-over-able final model.
#
# Input:  year1_ge_clean_v3.rds,
#         ge4_debray_v3_shared_coefs_by_fold.rds
# Output: ge4_debray_v3_FINAL_MODEL.rds
# =====================================================================

library(dplyr)
library(mice)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")

base_select <- function(df) {
  df %>% select(blakey_outcome_ge4, iecv_fold_v3, age_cat, Gender_Coded, Smoking_clean,
                 ics_cat, Rhinitis_Dx,
                 ba_Acute_OS_Courses, ba_asthma_AE, ba_all_consultations,
                 ba_SABA_Dosage, ba_LTRA_Rx, ba_LABA_Rx,
                 BMI_clean, ba_PF_Percent_Pred, Eosinophil_clean)
}

meth_spec <- function(df) {
  m <- make.method(df)
  m[c("BMI_clean", "ba_PF_Percent_Pred", "Eosinophil_clean")] <- "pmm"
  m["Smoking_clean"] <- "polyreg"
  m
}

recode_derived_cats <- function(df) {
  df$bmi_cat <- cut(df$BMI_clean, breaks = c(-Inf, 18.5, 25, 30, Inf),
                     labels = c("Underweight", "Normal", "Overweight", "Obese"))
  df$bmi_cat <- relevel(df$bmi_cat, ref = "Normal")
  df$pef_cat <- cut(df$ba_PF_Percent_Pred, breaks = c(-Inf, 60, 79, Inf),
                     labels = c("<=60", "61-79", ">=80"))
  df$pef_cat <- relevel(df$pef_cat, ref = ">=80")
  df$eos_cat <- cut(df$Eosinophil_clean, breaks = c(-Inf, 0.4, Inf),
                     labels = c("<=0.4", ">0.4"))
  df$eos_cat <- relevel(df$eos_cat, ref = "<=0.4")
  df
}

formula_stratified <- as.formula(
  "blakey_outcome_ge4 ~ 0 + iecv_fold_v3 + age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx"
)

# ---------------------------------------------------------------------
# All patients with a known GE4 outcome, across all 7 eras — no fold
# held out
# ---------------------------------------------------------------------
final_data <- year1 %>% filter(!is.na(blakey_outcome_ge4)) %>% base_select()
final_data$iecv_fold_v3 <- factor(final_data$iecv_fold_v3)

cat("--- Folds present in final model data (should be all 7) ---\n")
print(table(final_data$iecv_fold_v3))

# ---------------------------------------------------------------------
# Multiple imputation (m = 10) on the full dataset, single-core
# ---------------------------------------------------------------------
imp_final <- mice(final_data, m = 10, method = meth_spec(final_data),
                   maxit = 5, seed = 123, printFlag = FALSE)

cat("--- mice() logged events, FINAL MODEL data ---\n")
print(imp_final$loggedEvents)

# --- Verification: confirm imputed copies are genuinely diverse
#     across all four imputed variables before proceeding ---
check_vars <- c("Eosinophil_clean", "BMI_clean", "ba_PF_Percent_Pred", "Smoking_clean")
comp1 <- complete(imp_final, 1)
total_diff <- 0
for (v in check_vars) {
  for (k in 2:imp_final$m) {
    compk <- complete(imp_final, k)
    total_diff <- total_diff + sum(comp1[[v]] != compk[[v]], na.rm = TRUE)
  }
}
cat("Imputation diversity check (sum of differences, copy 1 vs all others, all 4 vars):", total_diff, "\n")
if (total_diff == 0) stop("Imputed copies are not sufficiently diverse -- halting.")

# ---------------------------------------------------------------------
# Fit the stratified-intercept model on each imputed copy, pool via
# Rubin's rules
# ---------------------------------------------------------------------
fits_final <- list()
for (i in 1:imp_final$m) {
  train_i <- recode_derived_cats(complete(imp_final, i))
  fits_final[[i]] <- glm(formula_stratified, data = train_i, family = binomial)
}

pooled_final <- pool(as.mira(fits_final))
pooled_coefs_final <- summary(pooled_final)$estimate
names(pooled_coefs_final) <- summary(pooled_final)$term

fold_intercept_names <- grep("^iecv_fold_v3", names(pooled_coefs_final), value = TRUE)
shared_coefs_final <- pooled_coefs_final[setdiff(names(pooled_coefs_final), fold_intercept_names)]

cat("\n--- Total coefficient count (7 fold intercepts + shared predictor terms) ---\n")
cat("Count:", length(pooled_coefs_final), "\n")

# ---------------------------------------------------------------------
# Stability check: final model's shared coefficients vs the average of
# the 7 development-fold coefficients from the IECV loop
# ---------------------------------------------------------------------
shared_coefs_by_fold <- readRDS("/users/hlskhatt/outputs/ge4_debray_v3_shared_coefs_by_fold.rds")
avg_dev_coefs <- Reduce(`+`, lapply(shared_coefs_by_fold, function(x) x[names(shared_coefs_final)])) /
  length(shared_coefs_by_fold)

comparison_coefs <- data.frame(
  term = names(shared_coefs_final),
  avg_7fold_dev = round(avg_dev_coefs, 4),
  final_7fold = round(shared_coefs_final, 4),
  abs_diff = round(abs(avg_dev_coefs - shared_coefs_final), 4)
)

cat("\n--- Stability check: final model vs average of the 7 development fits ---\n")
print(comparison_coefs, row.names = FALSE)

saveRDS(list(pooled_coefs_final = pooled_coefs_final,
             fold_intercepts = pooled_coefs_final[fold_intercept_names],
             shared_coefs_final = shared_coefs_final,
             comparison_coefs = comparison_coefs),
        "/users/hlskhatt/outputs/ge4_debray_v3_FINAL_MODEL.rds")

cat("\n=== FINAL COMBINED MODEL SAVED ===\n")
