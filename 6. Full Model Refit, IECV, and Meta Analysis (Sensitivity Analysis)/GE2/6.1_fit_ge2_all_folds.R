# =============================================================================
# GE2 FULL MODEL REFIT — FOLD FITTING (SENSITIVITY ANALYSIS)
# =============================================================================
# Purpose:
#   Implements the full-refit IECV sensitivity analysis for the GE2 outcome.
#   Unlike the primary (Debray) framework, which shares one fixed coefficient
#   set across all folds and only re-estimates the intercept, this script
#   re-estimates BOTH coefficients and intercept independently in each of the
#   7 temporal folds. This tests whether Blakey's predictor-outcome
#   relationships are stable over time, rather than assuming it.
#
# Method:
#   - Leave-one-fold-out cross-validation across 7 temporal folds (v3 boundaries)
#   - Multiple imputation (m = 10) via single-core mice(), with an explicit
#     post-imputation diversity check to confirm the 10 imputed copies are
#     genuinely distinct before they are used for anything downstream
#   - Rubin's rules pooling of development-fold coefficients
#   - Offset-trick intercept recalibration on the held-out fold
#   - C-statistic and calibration slope evaluated per imputed copy, then pooled
#   - O:E ratio recovered via Steyerberg & Harrell bootstrap recalibration
#     (200 reps): intercept is recalibrated on a bootstrap resample of the
#     held-out fold, then evaluated against that fold's original data
#   - C-statistic bootstrap SE is computed by evaluating discrimination on the
#     resampled data itself, not the full original fold, since AUC is
#     invariant to a constant intercept shift on unresampled data
#
# Output:
#   - Per-fold summary table (n, C-statistic, calibration slope, O:E, SEs)
#   - Per-fold .rds files containing pooled coefficients and bootstrap draws
# =============================================================================

library(dplyr)
library(mice)
library(pROC)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")

# -----------------------------------------------------------------------
# Diversity check: confirms multiple imputation produced genuinely distinct
# copies before they are used for any downstream model fitting or pooling.
# Halts the script immediately if any two copies are found to be identical.
# -----------------------------------------------------------------------
check_imputation_diversity <- function(imp, vars_to_check, label) {
  m <- imp$m
  min_diff <- Inf
  for (a in 1:(m - 1)) {
    for (b in (a + 1):m) {
      copy_a <- complete(imp, a)
      copy_b <- complete(imp, b)
      diff_count <- 0
      for (v in vars_to_check) {
        diff_count <- diff_count + sum(copy_a[[v]] != copy_b[[v]], na.rm = TRUE)
      }
      if (diff_count < min_diff) min_diff <- diff_count
    }
  }
  cat(sprintf("[%s] Smallest difference found between any two imputed copies: %d\n", label, min_diff))
  if (min_diff == 0) {
    stop(sprintf("STOPPING: [%s] Two imputed copies are completely identical. Imputation diversity check failed — do not continue.", label))
  }
  invisible(min_diff)
}

build_split_v3 <- function(year1, held_out_fold, outcome_col) {
  base_select <- function(df) {
    df %>% select(all_of(outcome_col), age_cat, Gender_Coded, Smoking_clean,
                   ics_cat, Rhinitis_Dx,
                   ba_Acute_OS_Courses, ba_asthma_AE, ba_all_consultations,
                   ba_SABA_Dosage, ba_LTRA_Rx, ba_LABA_Rx,
                   Eczema_Dx, GERD_Dx, NSAIDS,
                   BMI_clean, ba_PF_Percent_Pred, Eosinophil_clean)
  }
  train_raw <- year1 %>% filter(iecv_fold_v3 != held_out_fold, !is.na(.data[[outcome_col]])) %>% base_select()
  valid_raw <- year1 %>% filter(iecv_fold_v3 == held_out_fold, !is.na(.data[[outcome_col]])) %>% base_select()
  list(train_raw = train_raw, valid_raw = valid_raw)
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

model_formula <- as.formula(
  "blakey_outcome_ge2 ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx +
   Eczema_Dx + GERD_Dx + NSAIDS"
)

formula_no_outcome <- ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
  pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
  ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
  ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx +
  Eczema_Dx + GERD_Dx + NSAIDS

impute_vars <- c("BMI_clean", "ba_PF_Percent_Pred", "Eosinophil_clean", "Smoking_clean")

all_folds <- c("Fold1_2005_Apr2008", "Fold2_May08_Mar09", "Fold3_Apr09_Aug10",
               "Fold4_Sep10_Aug11", "Fold5_Sep11_Dec11", "Fold6_2012", "Fold7_2013")

results_dir <- "/users/hlskhatt/outputs/ge2_v3_fold_results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

full_results_ge2_v3 <- data.frame()

for (fold in all_folds) {

  split_f <- build_split_v3(year1, fold, "blakey_outcome_ge2")

  # Training-fold imputation: single-core mice(), diversity-checked
  imp_train <- mice(split_f$train_raw, m = 10, method = meth_spec(split_f$train_raw),
                     maxit = 5, seed = 123, printFlag = FALSE)
  check_imputation_diversity(imp_train, impute_vars, paste0(fold, "_train"))

  fits <- list()
  for (i in 1:imp_train$m) {
    train_i <- complete(imp_train, i)
    train_i <- recode_derived_cats(train_i)
    fits[[i]] <- glm(model_formula, data = train_i, family = binomial)
  }
  pooled <- pool(as.mira(fits))
  pooled_summary <- summary(pooled, type = "all")
  pooled_coefs <- pooled_summary$estimate
  names(pooled_coefs) <- pooled_summary$term
  train_intercept <- pooled_coefs["(Intercept)"]

  # Validation-fold imputation: single-core mice(), diversity-checked
  imp_valid <- mice(split_f$valid_raw, m = 10, method = meth_spec(split_f$valid_raw),
                     maxit = 5, seed = 123, printFlag = FALSE)
  check_imputation_diversity(imp_valid, impute_vars, paste0(fold, "_valid"))

  lp_no_int_list <- list()
  actual_list <- list()
  for (i in 1:imp_valid$m) {
    valid_i <- complete(imp_valid, i)
    valid_i <- recode_derived_cats(valid_i)

    complete_rows <- complete.cases(valid_i[, all.vars(formula_no_outcome)])
    if (any(!complete_rows)) valid_i <- valid_i[complete_rows, ]

    mm <- model.matrix(formula_no_outcome, data = valid_i)
    lp_no_int_list[[i]] <- (mm %*% pooled_coefs[colnames(mm)]) - train_intercept
    actual_list[[i]] <- valid_i$blakey_outcome_ge2
  }

  n_fold <- length(actual_list[[1]])
  actual_outcome <- actual_list[[1]]

  # Offset-trick recalibration of the intercept on the held-out fold
  new_int_fits <- list()
  for (i in 1:imp_valid$m) {
    new_int_fits[[i]] <- glm(actual_list[[i]] ~ offset(lp_no_int_list[[i]]), family = binomial)
  }
  new_intercept <- summary(pool(as.mira(new_int_fits)))$estimate[1]

  pred_prob_matrix <- matrix(NA, nrow = n_fold, ncol = imp_valid$m)
  for (i in 1:imp_valid$m) pred_prob_matrix[, i] <- 1 / (1 + exp(-(lp_no_int_list[[i]] + new_intercept)))

  c_stat_per_copy <- sapply(1:imp_valid$m, function(i) as.numeric(auc(roc(actual_outcome, pred_prob_matrix[, i], quiet = TRUE))))
  c_stat_median <- median(c_stat_per_copy)

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

  # Bootstrap: recalibrate intercept on a resample, evaluate against original
  # fold data (Steyerberg & Harrell method) — recovers a genuine per-fold O:E
  set.seed(123)
  B <- 200
  boot_oe <- numeric(B)
  boot_cstat <- numeric(B)

  for (b in 1:B) {
    resample_rows <- sample(1:n_fold, n_fold, replace = TRUE)
    pred_prob_copies <- matrix(NA, nrow = n_fold, ncol = imp_valid$m)
    for (i in 1:imp_valid$m) {
      recal_fit <- glm(actual_list[[i]][resample_rows] ~ offset(lp_no_int_list[[i]][resample_rows]), family = binomial)
      new_int_b <- coef(recal_fit)[1]
      pred_prob_copies[, i] <- 1 / (1 + exp(-(lp_no_int_list[[i]] + new_int_b)))
    }

    boot_outcome <- actual_outcome[resample_rows]
    c_stat_copies_b <- sapply(1:imp_valid$m, function(i)
      as.numeric(auc(roc(boot_outcome, pred_prob_copies[resample_rows, i], quiet = TRUE))))
    boot_cstat[b] <- median(c_stat_copies_b)

    pred_prob_avg <- rowMeans(pred_prob_copies)
    boot_oe[b] <- sum(actual_outcome) / sum(pred_prob_avg)
  }

  c_stat_se <- sd(boot_cstat)
  oe_mean <- mean(boot_oe)
  oe_se <- sd(boot_oe)

  fold_result <- data.frame(
    fold = fold, n_train = nrow(split_f$train_raw), n_valid = n_fold,
    c_stat = round(c_stat_median, 4), c_stat_se = round(c_stat_se, 4),
    slope = round(slope_estimate, 4), slope_se = round(slope_se, 4),
    oe_mean = round(oe_mean, 4), oe_se = round(oe_se, 5)
  )

  full_results_ge2_v3 <- rbind(full_results_ge2_v3, fold_result)

  saveRDS(list(fold = fold, pooled_coefs = pooled_coefs,
               c_stat = c_stat_median, c_stat_se = c_stat_se,
               slope = slope_estimate, slope_se = slope_se,
               boot_oe = boot_oe, boot_cstat = boot_cstat,
               oe_mean = oe_mean, oe_se = oe_se),
          file.path(results_dir, paste0("ge2_v3_", fold, "_result.rds")))

  saveRDS(full_results_ge2_v3, "/users/hlskhatt/outputs/ge2_v3_all_folds_summary.rds")
}
