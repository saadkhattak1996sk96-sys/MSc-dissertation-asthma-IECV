library(dplyr)
library(pmvalsampsize)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean.rds")

all_folds <- c("Fold1_2005_2007", "Fold2_Jan08_Mar09", "Fold3_Apr09_Aug10",
               "Fold4_Sep10_Aug11", "Fold5_Sep11_Dec11", "Fold6_2012", "Fold7_2013")

anticipated_cstat <- list(
  GE2 = 0.785,
  GE4 = 0.867
)

formula_ge2 <- as.formula(
  "blakey_outcome_ge2 ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx +
   Eczema_Dx + GERD_Dx + NSAIDS"
)

formula_ge4 <- as.formula(
  "blakey_outcome_ge4 ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx"
)

predictor_cols <- function(f) all.vars(f)[-1]

build_lp_stats <- function(year1, held_out_fold, outcome_col, formula_obj) {
  cols_needed <- c(outcome_col, predictor_cols(formula_obj),
                   "BMI_clean", "ba_PF_Percent_Pred", "Eosinophil_clean")
  cols_needed <- unique(cols_needed)

  df <- year1 %>%
    filter(iecv_fold_v2 != held_out_fold, !is.na(.data[[outcome_col]])) %>%
    select(all_of(cols_needed))

  df$bmi_cat <- cut(df$BMI_clean, breaks = c(-Inf, 18.5, 25, 30, Inf),
                     labels = c("Underweight", "Normal", "Overweight", "Obese"))
  df$bmi_cat <- relevel(df$bmi_cat, ref = "Normal")
  df$pef_cat <- cut(df$ba_PF_Percent_Pred, breaks = c(-Inf, 60, 79, Inf),
                     labels = c("<=60", "61-79", ">=80"))
  df$pef_cat <- relevel(df$pef_cat, ref = ">=80")
  df$eos_cat <- cut(df$Eosinophil_clean, breaks = c(-Inf, 0.4, Inf),
                     labels = c("<=0.4", ">0.4"))
  df$eos_cat <- relevel(df$eos_cat, ref = "<=0.4")

  df <- df[complete.cases(df[, all.vars(formula_obj)]), ]

  fit <- glm(formula_obj, data = df, family = binomial)
  lp <- predict(fit, type = "link")

  valid_df <- year1 %>%
    filter(iecv_fold_v2 == held_out_fold, !is.na(.data[[outcome_col]]))
  prevalence <- mean(valid_df[[outcome_col]], na.rm = TRUE)

  list(lp_mean = mean(lp), lp_sd = sd(lp), prevalence = prevalence,
       n_dev = nrow(df), n_valid = nrow(valid_df))
}

run_adequacy <- function(outcome_label, outcome_col, formula_obj) {
  cstat_anticipated <- anticipated_cstat[[outcome_label]]
  results <- data.frame()

  for (fold in all_folds) {
    stats_f <- build_lp_stats(year1, fold, outcome_col, formula_obj)

    ss <- pmvalsampsize(
      type       = "b",
      cstatistic = cstat_anticipated,
      prevalence = stats_f$prevalence,
      lpnormal   = c(stats_f$lp_mean, stats_f$lp_sd)
    )

    results <- rbind(results, data.frame(
      outcome = outcome_label, fold = fold,
      n_dev = stats_f$n_dev, n_valid = stats_f$n_valid,
      prevalence = round(stats_f$prevalence, 4),
      cstat_anticipated = cstat_anticipated,
      n_required = ss$results_table["Final SS", "Samp_size"],
      adequacy_ratio = round(stats_f$n_valid / ss$results_table["Final SS", "Samp_size"], 2)
    ))
  }
  results
}

ge2_adequacy <- run_adequacy("GE2", "blakey_outcome_ge2", formula_ge2)
ge4_adequacy <- run_adequacy("GE4", "blakey_outcome_ge4", formula_ge4)

adequacy_all <- rbind(ge2_adequacy, ge4_adequacy)

print(adequacy_all)

saveRDS(adequacy_all, "/users/hlskhatt/outputs/pmvalsampsize_ge2_ge4_all_folds.rds")
write.csv(adequacy_all, "/users/hlskhatt/outputs/pmvalsampsize_ge2_ge4_all_folds.csv", row.names = FALSE)
