# =====================================================================
# fit_ge4_all_folds.R
# GE4 (>=4 exacerbations) full-refit IECV pipeline, v3 fold structure
# Mirrors fit_ge2_all_folds.R exactly, with the outcome swapped to
# blakey_outcome_ge4 and Eczema_Dx/GERD_Dx/NSAIDS dropped, since Blakey's
# own GE4 model excludes these three predictors by his stepwise selection.
# =====================================================================

library(dplyr)
library(mice)
library(pROC)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean_v3.rds")

# --- Build train/validation split for a given held-out fold ---
build_split_v3_ge4 <- function(year1, held_out_fold, outcome_col) {
  base_select <- function(df) {
    df %>% select(all_of(outcome_col), age_cat, Gender_Coded, Smoking_clean,
                   ics_cat, Rhinitis_Dx,
                   ba_Acute_OS_Courses, ba_asthma_AE, ba_all_consultations,
                   ba_SABA_Dosage, ba_LTRA_Rx, ba_LABA_Rx,
                   BMI_clean, ba_PF_Percent_Pred, Eosinophil_clean)
  }
  train_raw <- year1 %>% filter(iecv_fold_v3 != held_out_fold, !is.na(.data[[outcome_col]])) %>% base_select()
  valid_raw <- year1 %>% filter(iecv_fold_v3 == held_out_fold, !is.na(.data[[outcome_col]])) %>% base_select()
  list(train_raw = train_raw, valid_raw = valid_raw)
}

# --- Imputation method spec: PMM for continuous predictors with real
#     missingness (BMI, PEF, eosinophil), polytomous regression for smoking ---
meth_spec <- function(df) {
  m <- make.method(df)
  m[c("BMI_clean", "ba_PF_Percent_Pred", "Eosinophil_clean")] <- "pmm"
  m["Smoking_clean"] <- "polyreg"
  m
}

# --- Derive categorical predictors from imputed continuous variables,
#     applied AFTER imputation so bmi_cat/pef_cat/eos_cat never carry
#     independent missingness of their own (fixes the GE2 imputation-
#     consistency bug found during QA — same fix applied here from the start) ---
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

# --- GE4 model formula: identical predictor set to GE2, minus
#     Eczema_Dx/GERD_Dx/NSAIDS (Blakey's GE4 stepwise selection excludes these) ---
model_formula_ge4 <- as.formula(
  "blakey_outcome_ge4 ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx"
)

formula_no_outcome_ge4 <- ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
  pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
  ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
  ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx

all_folds <- c("Fold1_2005_Apr2008", "Fold2_May08_Mar09", "Fold3_Apr09_Aug10",
               "Fold4_Sep10_Aug11", "Fold5_Sep11_Dec11", "Fold6_2012", "Fold7_2013")

n_cores <- min(10, parallel::detectCores())
results_dir <- "/users/hlskhatt/outputs/ge4_v3_fold_results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

full_results_ge4_v3 <- data.frame()

for (fold in all_folds) {

  split_f <- build_split_v3_ge4(year1, fold, "blakey_outcome_ge4")

  # --- Impute training data (m=10), fit GE4 model on each copy, pool
  #     coefficients via Rubin's rules ---
  imp_train <- futuremice(split_f$train_raw, m = 10, method = meth_spec(split_f$train_raw),
                           maxit = 5, seed = 123, parallelseed = 123, n.core = n_cores)

  fits <- list()
  for (i in 1:imp_train$m) {
    train_i <- complete(imp_train, i)
    train_i <- recode_derived_cats(train_i)
    fits[[i]] <- glm(model_formula_ge4, data = train_i, family = binomial)
  }
  pooled <- pool(as.mira(fits))
  pooled_summary <- summary(pooled, type = "all")
  pooled_coefs <- pooled_summary$estimate
  names(pooled_coefs) <- pooled_summary$term
  train_intercept <- pooled_coefs["(Intercept)"]

  # --- Impute validation (held-out) data separately, m=10 ---
  imp_valid <- futuremice(split_f$valid_raw, m = 10, method = meth_spec(split_f$valid_raw),
                           maxit = 5, seed = 123, parallelseed = 123, n.core = n_cores)

  # --- Linear predictor with intercept stripped out, per imputed copy ---
  lp_no_int_list <- list()
  actual_list <- list()
  for (i in 1:imp_valid$m) {
    valid_i <- complete(imp_valid, i)
    valid_i <- recode_derived_cats(valid_i)

    complete_rows <- complete.cases(valid_i[, all.vars(formula_no_outcome_ge4)])
    if (any(!complete_rows)) valid_i <- valid_i[complete_rows, ]

    mm <- model.matrix(formula_no_outcome_ge4, data = valid_i)
    lp_no_int_list[[i]] <- (mm %*% pooled_coefs[colnames(mm)]) - train_intercept
    actual_list[[i]] <- valid_i$blakey_outcome_ge4
  }

  n_fold <- length(actual_list[[1]])
  actual_outcome <- actual_list[[1]]

  # --- Mandatory per-fold intercept recalibration (fixes the "no
  #     recalibration" bug found in the GE2 build — applied correctly here
  #     from the start) ---
  new_int_fits <- list()
  for (i in 1:imp_valid$m) {
    new_int_fits[[i]] <- glm(actual_list[[i]] ~ offset(lp_no_int_list[[i]]), family = binomial)
  }
  new_intercept <- summary(pool(as.mira(new_int_fits)))$estimate[1]

  pred_prob_matrix <- matrix(NA, nrow = n_fold, ncol = imp_valid$m)
  for (i in 1:imp_valid$m) pred_prob_matrix[, i] <- 1 / (1 + exp(-(lp_no_int_list[[i]] + new_intercept)))

  # --- C-statistic: median across MI copies (Marshall/Wahl-style
  #     treatment for C-stat's non-normal sampling distribution) ---
  c_stat_per_copy <- sapply(1:imp_valid$m, function(i) as.numeric(auc(roc(actual_outcome, pred_prob_matrix[, i], quiet = TRUE))))
  c_stat_median <- median(c_stat_per_copy)

  # --- Calibration slope: proper Rubin's rules pooling via pool.scalar ---
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

  # --- Bootstrap (200 reps): C-stat SE via resample-and-evaluate;
  #     O:E via Steyerberg & Harrell method (recalibrate on resample,
  #     evaluate on the ORIGINAL full fold) ---
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

  full_results_ge4_v3 <- rbind(full_results_ge4_v3, fold_result)

  saveRDS(list(fold = fold, pooled_coefs = pooled_coefs,
               c_stat = c_stat_median, c_stat_se = c_stat_se,
               slope = slope_estimate, slope_se = slope_se,
               boot_oe = boot_oe, boot_cstat = boot_cstat,
               oe_mean = oe_mean, oe_se = oe_se),
          file.path(results_dir, paste0("ge4_v3_", fold, "_result.rds")))

  saveRDS(full_results_ge4_v3, "/users/hlskhatt/outputs/ge4_v3_all_folds_summary.rds")
}
