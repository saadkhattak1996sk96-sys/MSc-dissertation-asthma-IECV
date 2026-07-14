library(dplyr)
library(lubridate)
library(pmvalsampsize)

year1 <- readRDS("/users/hlskhatt/outputs/year1_ge_clean.rds")

if (!"entry_month" %in% names(year1)) {
  year1$entry_month <- floor_date(dmy(year1$Start_Year_1), "month")
}

anticipated_cstat_ge4 <- 0.867

formula_ge4 <- as.formula(
  "blakey_outcome_ge4 ~ age_cat + Gender_Coded + bmi_cat + Smoking_clean +
   pef_cat + ics_cat + Rhinitis_Dx + eos_cat +
   ba_Acute_OS_Courses + ba_asthma_AE + ba_all_consultations +
   ba_SABA_Dosage + ba_LTRA_Rx + ba_LABA_Rx"
)

predictor_cols_ge4 <- all.vars(formula_ge4)[-1]

current_boundary <- as.Date("2007-12-01")
fold2_ceiling     <- as.Date("2009-03-01")
candidate_boundaries <- seq(current_boundary, fold2_ceiling %m-% months(1), by = "1 month")

assign_split_folds <- function(year1, boundary) {
  year1 %>%
    mutate(
      test_fold = case_when(
        entry_month <= boundary       ~ "Fold1_test",
        entry_month <= fold2_ceiling  ~ "Fold2_test",
        TRUE ~ NA_character_
      )
    )
}

build_lp_and_adequacy <- function(df_split, target_fold, all_other_folds_train) {
  train_df <- df_split %>%
    filter(!(test_fold %in% c("Fold1_test", "Fold2_test")) | test_fold != target_fold,
           !is.na(blakey_outcome_ge4)) %>%
    select(blakey_outcome_ge4, all_of(predictor_cols_ge4))

  valid_df <- df_split %>%
    filter(test_fold == target_fold, !is.na(blakey_outcome_ge4))

  train_df <- train_df[complete.cases(train_df), ]

  fit <- glm(formula_ge4, data = train_df, family = binomial)
  lp  <- predict(fit, type = "link")

  prevalence <- mean(valid_df$blakey_outcome_ge4, na.rm = TRUE)
  n_valid    <- nrow(valid_df)

  ss <- pmvalsampsize(
    type       = "b",
    cstatistic = anticipated_cstat_ge4,
    prevalence = prevalence,
    lpnormal   = c(mean(lp), sd(lp))
  )

  n_required <- ss$results_table["Final SS", "Samp_size"]

  list(n_valid = n_valid, prevalence = round(prevalence, 4),
       n_required = n_required, adequacy_ratio = round(n_valid / n_required, 2))
}

results <- data.frame()

for (b in candidate_boundaries) {
  boundary_date <- as.Date(b, origin = "1970-01-01")
  df_split <- assign_split_folds(year1, boundary_date)

  fold1_res <- build_lp_and_adequacy(df_split, "Fold1_test", "Fold2_test")
  fold2_res <- build_lp_and_adequacy(df_split, "Fold2_test", "Fold1_test")

  results <- rbind(results, data.frame(
    boundary = as.character(boundary_date),
    fold1_n_valid = fold1_res$n_valid, fold1_prevalence = fold1_res$prevalence,
    fold1_adequacy = fold1_res$adequacy_ratio,
    fold2_n_valid = fold2_res$n_valid, fold2_prevalence = fold2_res$prevalence,
    fold2_adequacy = fold2_res$adequacy_ratio
  ))
}

print(results)

both_ok <- results[results$fold1_adequacy >= 1.0 & results$fold2_adequacy >= 1.0, ]
cat("\n--- Boundaries where BOTH Fold1 and Fold2 are adequate (>=1.0) ---\n")
print(both_ok)

saveRDS(results, "/users/hlskhatt/outputs/ge4_fold1_boundary_shift_test.rds")
write.csv(results, "/users/hlskhatt/outputs/ge4_fold1_boundary_shift_test.csv", row.names = FALSE)
