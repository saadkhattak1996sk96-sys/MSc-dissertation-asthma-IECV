# =====================================================================
# fit_ge2_debray_all_folds.R
#
# GE2 (>=2 exacerbations) — TRUE DEBRAY IECV PIPELINE, v3 fold structure
#
# This is the primary analysis for the dissertation. Unlike a full model
# refit (run separately as a sensitivity analysis), Debray's framework
# freezes one shared set of predictor coefficients across all folds and
# allows only the intercept (baseline risk) to be re-estimated fresh in
# each held-out fold. This script performs the full internal-external
# cross-validation (IECV) cycle across all 7 v3 folds:
#
#   1. Fit a stratified-intercept development model on the 6 training
#      folds (one intercept per training fold, one shared coefficient
#      set across all of them).
#   2. Freeze the shared coefficients, discard the training intercepts.
#   3. Re-estimate a fresh intercept for the held-out fold using the
#      offset-trick (Debray 2013, Section 2.2.4).
#   4. Evaluate discrimination (C-statistic), calibration (slope), and
#      O:E ratio on the held-out fold.
#
# Missing data (BMI, PEF, eosinophil, smoking) is handled via multiple
# imputation (m = 10), performed separately within the training and
# validation data of each fold to avoid information leakage.
#
# Method notes:
#   - C-statistic: median across the 10 imputed copies (Marshall 2009),
#     with its standard error obtained via a 200-repetition bootstrap.
#   - Calibration slope: pooled across the 10 imputed copies using
#     Rubin's rules (Wood 2015).
#   - O:E ratio: under true Debray, this is a mathematical identity
#     (~1.0 in every fold) rather than an independent estimate, since
#     the held-out fold's intercept is re-estimated from its own
#     observed outcome. It is retained here for completeness and as a
#     confirmation that recalibration behaved as expected, not as an
#     informative performance measure in its own right.
#
# Output: one row per fold (fit statistics) plus each fold's frozen
# shared coefficients, saved to disk for use by the pooling, forest
# plot, meta-regression, and final-model scripts.
# =====================================================================

library(dplyr)
library(mice)
library(pROC)

# ---------------------------------------------------------------------
# Load the cleaned, fold-labelled cohort (v3 fold boundaries)
# ---------------------------------------------------------------------
year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")

# ---------------------------------------------------------------------
# Column selection: outcome, fold label, and all 17 GE2 predictors
# (age, sex, BMI, smoking, PEF, ICS adherence, rhinitis, eosinophil,
# acute OCS courses, A&E attendance, consultations, SABA/LTRA/LABA
# prescribing, eczema, GERD, NSAIDs)
# ---------------------------------------------------------------------
base_select <- function(df) {
  df %>% select(blakey_outcome_ge2, iecv_fold_v3, age_cat, Gender_Coded, Smoking_clean,
                 ics_cat, Rhinitis_Dx,
                 ba_Acute_OS_Courses, ba_asthma_AE, ba_all_consultations,
                 ba_SABA_Dosage, ba_LTRA_Rx, ba_LABA_Rx,
                 Eczema_Dx, GERD_Dx, NSAIDS,
                 BMI_clean, ba_PF_Percent_Pred, Eosinophil_clean)
}

# ---------------------------------------------------------------------
# Imputation method specification: predictive mean matching for the
# three continuous predictors with genuine missingness, polytomous
# regression for smoking status. Categorical derivatives (bmi_cat,
# pef_cat, eos_cat) are deliberately NOT included here — they are
# constructed after imputation from the completed continuous values,
# so that a patient's category can never disagree with their own
# imputed continuous value.
# ---------------------------------------------------------------------
meth_spec <- function(df) {
  m <- make.method(df)
  m[c("BMI_clean", "ba_PF_Percent_Pred", "Eosinophil_clean")] <- "pmm"
  m["Smoking_clean"] <- "polyreg"
  m
}

# ---------------------------------------------------------------------
# Derive categorical predictors from the completed continuous values,
# applied once per imputed copy after imputation
# ---------------------------------------------------------------------
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

# ---------------------------------------------------------------------
# Model formulae
#   - Stratified-intercept development formula: "0 +" removes the
#     single shared intercept; iecv_fold_v3 gives each training fold
#     its own intercept term (Debray 2013, Section 2.1.3)
#   - Shared-only formula (no intercept term): used to build the
#     design matrix for the held-out fold once coefficients are frozen
# ---------------------------------------------------------------------
formula_stratified <- as.formula(
  "blakey_outcome_ge2 ~ 0 + iecv_fold_v3 + age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx +
   Eczema_Dx + GERD_Dx + NSAIDS"
)

formula_shared_normal <- ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
  pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
  ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
  ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx +
  Eczema_Dx + GERD_Dx + NSAIDS

all_folds <- c("Fold1_2005_Apr2008", "Fold2_May08_Mar09", "Fold3_Apr09_Aug10",
               "Fold4_Sep10_Aug11", "Fold5_Sep11_Dec11", "Fold6_2012", "Fold7_2013")

results_dir <- "/users/hlskhatt/outputs/ge2_debray_v3_fold_results"

full_results_ge2_debray <- data.frame()
shared_coefs_by_fold <- list()

# =====================================================================
# MAIN IECV LOOP — one iteration per held-out fold
# =====================================================================
for (fold in all_folds) {

  cat("\n=========== STARTING FOLD:", fold, "===========\n")

  # --- Build training (6 folds) and validation (held-out fold) data ---
  train_raw <- year1 %>% filter(iecv_fold_v3 != fold, !is.na(blakey_outcome_ge2)) %>% base_select()
  valid_raw <- year1 %>% filter(iecv_fold_v3 == fold, !is.na(blakey_outcome_ge2)) %>% base_select()
  train_raw$iecv_fold_v3 <- factor(train_raw$iecv_fold_v3)

  # --- Multiple imputation, training data (m = 10) ---
  imp_train <- mice(train_raw, m = 10, method = meth_spec(train_raw),
                     maxit = 5, seed = 123, printFlag = FALSE)

  # --- Verification: confirm the 10 imputed copies are genuinely
  #     distinct from one another before proceeding to model fitting ---
  comp1 <- complete(imp_train, 1)
  comp2 <- complete(imp_train, 3)
  n_diff_check <- sum(comp1$Eosinophil_clean != comp2$Eosinophil_clean, na.rm = TRUE)
  cat("Imputation diversity check (copy 1 vs copy 3, Eosinophil):", n_diff_check, "\n")
  if (n_diff_check == 0) stop("Imputed copies are not sufficiently diverse in fold ", fold, " -- halting.")

  # --- Fit the stratified-intercept development model on each imputed
  #     copy, then pool via Rubin's rules ---
  fits_strat <- list()
  for (i in 1:imp_train$m) {
    train_i <- recode_derived_cats(complete(imp_train, i))
    fits_strat[[i]] <- glm(formula_stratified, data = train_i, family = binomial)
  }
  pooled_strat <- pool(as.mira(fits_strat))
  pooled_coefs_strat <- summary(pooled_strat)$estimate
  names(pooled_coefs_strat) <- summary(pooled_strat)$term

  # --- Freeze the shared coefficients; discard the training-fold
  #     intercepts entirely (they represent the training pile's baseline
  #     risk, not the held-out fold's) ---
  fold_intercept_names <- grep("^iecv_fold_v3", names(pooled_coefs_strat), value = TRUE)
  shared_coefs <- pooled_coefs_strat[setdiff(names(pooled_coefs_strat), fold_intercept_names)]
  shared_coefs_by_fold[[fold]] <- shared_coefs

  # --- Multiple imputation, validation (held-out) data, m = 10,
  #     performed independently of the training imputation ---
  imp_valid <- mice(valid_raw, m = 10, method = meth_spec(valid_raw),
                     maxit = 5, seed = 123, printFlag = FALSE)

  # --- Build the linear predictor (frozen coefficients, no intercept)
  #     for each imputed copy of the held-out fold ---
  lp_no_int_list <- list()
  actual_list <- list()
  for (i in 1:imp_valid$m) {
    valid_i <- recode_derived_cats(complete(imp_valid, i))
    mm_full <- model.matrix(formula_shared_normal, data = valid_i)
    mm_no_int <- mm_full[, colnames(mm_full) != "(Intercept)", drop = FALSE]

    missing_check <- setdiff(colnames(mm_no_int), names(shared_coefs))
    if (length(missing_check) > 0) {
      stop("Fold ", fold, " missing coefficients for: ", paste(missing_check, collapse = ", "))
    }

    lp_no_int_list[[i]] <- mm_no_int %*% shared_coefs[colnames(mm_no_int)]
    actual_list[[i]] <- valid_i$blakey_outcome_ge2
  }

  n_fold <- length(actual_list[[1]])
  actual_outcome <- actual_list[[1]]

  # --- Recalibrate a fresh intercept for the held-out fold via the
  #     offset trick (Debray 2013, Section 2.2.4), pooled across the
  #     10 imputed copies using Rubin's rules ---
  new_int_fits <- list()
  for (i in 1:imp_valid$m) {
    new_int_fits[[i]] <- glm(actual_list[[i]] ~ offset(lp_no_int_list[[i]]), family = binomial)
  }
  new_intercept <- summary(pool(as.mira(new_int_fits)))$estimate[1]

  # --- Predicted probabilities under the recalibrated model ---
  pred_prob_matrix <- matrix(NA, nrow = n_fold, ncol = imp_valid$m)
  for (i in 1:imp_valid$m) pred_prob_matrix[, i] <- 1 / (1 + exp(-(lp_no_int_list[[i]] + new_intercept)))

  # --- Discrimination: C-statistic per imputed copy, median across
  #     copies (Marshall 2009) ---
  c_stat_per_copy <- sapply(1:imp_valid$m, function(i)
    as.numeric(auc(roc(actual_outcome, pred_prob_matrix[, i], quiet = TRUE))))
  c_stat_median <- median(c_stat_per_copy)

  # --- Calibration: slope per imputed copy, pooled via Rubin's rules
  #     (Wood 2015) ---
  slope_est_per_copy <- numeric(imp_valid$m)
  slope_se_per_copy <- numeric(imp_valid$m)
  for (i in 1:imp_valid$m) {
    slope_fit <- glm(actual_outcome ~ qlogis(pred_prob_matrix[, i]), family = binomial)
    slope_est_per_copy[i] <- coef(slope_fit)[2]
    slope_se_per_copy[i] <- sqrt(vcov(slope_fit)[2, 2])
  }
  slope_pooled <- pool.scalar(Q = slope_est_per_copy, U = slope_se_per_copy^2, n = n_fold)
  slope_estimate <- slope_pooled$qbar
  slope_se <- sqrt(slope_pooled$t)

  # --- O:E ratio: expected to equal ~1.0 under correct recalibration
  #     (Van Calster et al. 2016) ---
  oe_ratio <- sum(actual_outcome) / sum(rowMeans(pred_prob_matrix))

  # --- Bootstrap (200 reps) for the C-statistic's standard error:
  #     resample patients with replacement, recompute discrimination
  #     on the resampled set ---
  set.seed(123)
  B <- 200
  boot_cstat <- numeric(B)
  for (b in 1:B) {
    resample_rows <- sample(1:n_fold, n_fold, replace = TRUE)
    boot_outcome <- actual_outcome[resample_rows]
    boot_preds <- pred_prob_matrix[resample_rows, , drop = FALSE]
    c_stat_copies_b <- sapply(1:imp_valid$m, function(i)
      as.numeric(auc(roc(boot_outcome, boot_preds[, i], quiet = TRUE))))
    boot_cstat[b] <- median(c_stat_copies_b)
  }
  c_stat_se <- sd(boot_cstat)

  # --- Record and save this fold's results ---
  fold_result <- data.frame(
    fold = fold, n_train = nrow(train_raw), n_valid = n_fold,
    c_stat = round(c_stat_median, 4), c_stat_se = round(c_stat_se, 4),
    slope = round(slope_estimate, 4), slope_se = round(slope_se, 4),
    oe_ratio = round(oe_ratio, 4)
  )

  full_results_ge2_debray <- rbind(full_results_ge2_debray, fold_result)

  saveRDS(list(fold = fold, shared_coefs = shared_coefs,
               c_stat = c_stat_median, c_stat_se = c_stat_se,
               slope = slope_estimate, slope_se = slope_se,
               boot_cstat = boot_cstat, oe_ratio = oe_ratio),
          file.path(results_dir, paste0("ge2_debray_v3_", fold, "_result.rds")))

  saveRDS(full_results_ge2_debray, "/users/hlskhatt/outputs/ge2_debray_v3_all_folds_summary.rds")
  saveRDS(shared_coefs_by_fold, "/users/hlskhatt/outputs/ge2_debray_v3_shared_coefs_by_fold.rds")

  cat("COMPLETED FOLD:", fold, "| C-stat:", round(c_stat_median, 4),
      "| Slope:", round(slope_estimate, 4), "| O:E:", round(oe_ratio, 4), "\n")
}

cat("\n\n=== ALL 7 FOLDS COMPLETE ===\n")
print(full_results_ge2_debray, row.names = FALSE)
