# =============================================================================
# 01_data_quality_and_preprocessing.R
#
# PURPOSE
# This script documents the data quality screening, recoding, and imputation
# steps applied to the OPCRD extract prior to model fitting. It reflects the
# preprocessing pipeline actually used in this project's analysis scripts
# (see 01_prep_cohort/ for the full cohort-construction code this integrates
# into). It is included here, separately and explicitly, so that a reader can
# see exactly what was checked and how, without needing to read every
# downstream script to reconstruct it.
#
# GOVERNANCE NOTE
# This script contains no patient-level data, no output, and no results.
# It is written against representative OPCRD column names and is intended
# to run only within the University of Liverpool's Barkla HPC environment,
# where the underlying dataset is held. See DATA_GOVERNANCE.md (repository
# root) for the accompanying methods-level summary.
# =============================================================================

library(dplyr)
library(mice)

# -----------------------------------------------------------------------------
# SECTION 1 — CLINICALLY IMPLAUSIBLE VALUE SCREENING
#
# OPCRD is a routinely collected primary care dataset. It is not curated for
# research use at the point of entry, and several continuous predictors were
# found, on inspection, to contain values that are not physiologically
# possible (e.g. a BMI of over 190, or a peak flow recorded as several
# thousand percent of predicted). These are not the dataset's standard
# missing-data sentinel (-1, handled separately below) — they are genuine,
# non-missing values that fall outside any plausible clinical range.
#
# Left unaddressed, these values would be treated as real data by multiple
# imputation and by the models themselves, silently distorting both. The
# function below recodes implausible values to NA so that they are handled
# by the same imputation process as genuine missingness, rather than treated
# as informative.
#
# Thresholds below reflect plausible clinical ranges for each variable, not
# arbitrary cutoffs: BMI outside 10-100 kg/m^2, PEF % predicted outside
# 10-200%, and eosinophil count above 10 (10^9/L) are not attainable in a
# living patient and were the boundaries actually used in this project.
# -----------------------------------------------------------------------------

screen_implausible_values <- function(df) {

  df <- df %>%
    mutate(
      BMI = if_else(BMI < 10 | BMI > 100, NA_real_, BMI),

      ba_PF_Percent_Pred = if_else(
        ba_PF_Percent_Pred < 10 | ba_PF_Percent_Pred > 200,
        NA_real_,
        ba_PF_Percent_Pred
      ),

      ba_Eosinophil_Count = if_else(
        ba_Eosinophil_Count > 10,
        NA_real_,
        ba_Eosinophil_Count
      )
      # Note: no floor applied to eosinophil count — a value of 0 is
      # clinically plausible and is not treated as implausible here.
    )

  df
}

# -----------------------------------------------------------------------------
# SECTION 2 — ICS ADHERENCE: SENTINEL HANDLING, NOT IMPUTATION
#
# `ba_ICS_adherence_fixed` uses -1 as a sentinel value. Critically, this does
# NOT mean "adherence unknown" — it means the patient was never prescribed
# an inhaled corticosteroid at all (confirmed via prescription count, drug
# name, and device fields all being simultaneously absent for this group).
# Treating -1 as missing-and-to-be-imputed would be wrong: it is not a gap
# in knowledge, it is a real, informative clinical state ("never
# prescribed"), and imputing a plausible-looking adherence percentage for
# these patients would manufacture information that does not exist.
#
# The function below instead builds an explicit "No ICS" category alongside
# the genuine adherence bands, so that never-prescribed status is retained
# as real information rather than discarded or imputed over.
# -----------------------------------------------------------------------------

build_ics_category <- function(df) {

  df <- df %>%
    mutate(
      ics_cat = case_when(
        ba_ICS_adherence_fixed == -1                                  ~ "No ICS",
        ba_ICS_adherence_fixed > 0  & ba_ICS_adherence_fixed < 40     ~ ">0-39.9",
        ba_ICS_adherence_fixed >= 40 & ba_ICS_adherence_fixed < 60    ~ "40-59.9",
        ba_ICS_adherence_fixed >= 60 & ba_ICS_adherence_fixed < 80    ~ "60-79.9",
        ba_ICS_adherence_fixed >= 80 & ba_ICS_adherence_fixed < 100   ~ "80-100",
        ba_ICS_adherence_fixed >= 100                                 ~ ">=100",
        TRUE                                                          ~ NA_character_
      ),
      ics_cat = factor(
        ics_cat,
        levels = c("No ICS", ">0-39.9", "40-59.9", "60-79.9", "80-100", ">=100")
      )
    )

  # Because "No ICS" is now its own explicit category, ICS adherence carries
  # zero genuine missingness among patients who were ever prescribed one —
  # this variable is therefore NOT included in the multiple imputation model
  # below; it is complete by construction, not by assumption.

  df
}

# -----------------------------------------------------------------------------
# SECTION 3 — MULTIPLE IMPUTATION, WITH AN EXPLICIT DIVERSITY SAFEGUARD
#
# Four variables carry genuine missingness after the steps above: eosinophil
# count, PEF % predicted, BMI, and smoking status. These are imputed using
# predictive mean matching (continuous) / polytomous regression (smoking),
# independently within each cross-validation fold, so that no information
# from a validation fold leaks into its own imputation model.
#
# This project intentionally uses single-core mice() rather than a
# multi-core parallel imputation wrapper. A parallel imputation approach was
# tested during development and found, on inspection, to produce imputed
# copies that were not genuinely independent of one another across cores —
# a failure mode that would silently understate imputation uncertainty in
# every downstream pooled estimate. Single-core mice() does not share that
# failure mode. The diversity check below is retained as a permanent
# safeguard: it will halt execution rather than allow a non-diverse set of
# imputations to pass silently into the analysis.
# -----------------------------------------------------------------------------

impute_fold_data <- function(df, m = 10, seed = 123) {

  vars_to_impute <- c("ba_Eosinophil_Count", "ba_PF_Percent_Pred", "BMI", "Smoking_Status")

  method <- make.method(df)
  method[!names(method) %in% vars_to_impute] <- ""
  method["Smoking_Status"] <- "polyreg"
  method[c("ba_Eosinophil_Count", "ba_PF_Percent_Pred", "BMI")] <- "pmm"

  imp <- mice(
    df,
    m = m,
    method = method,
    seed = seed,
    printFlag = FALSE
  )

  check_imputation_diversity(imp, vars_to_impute)

  imp
}

# Compares every pair of imputed copies against each other for the imputed
# variables, and stops execution if any two copies are found to be
# identical. Genuine imputation diversity is required for Rubin's rules to
# correctly represent imputation uncertainty; this check exists specifically
# to catch a failure mode that would otherwise be silent.
check_imputation_diversity <- function(imp, vars_to_impute) {

  completed_sets <- lapply(seq_len(imp$m), function(i) {
    complete(imp, action = i) %>% select(all_of(vars_to_impute))
  })

  n_copies <- length(completed_sets)

  for (i in 1:(n_copies - 1)) {
    for (j in (i + 1):n_copies) {
      identical_copy <- all(
        mapply(
          function(a, b) isTRUE(all.equal(a, b)),
          completed_sets[[i]],
          completed_sets[[j]]
        )
      )
      if (identical_copy) {
        stop(sprintf(
          "Imputation diversity check failed: copies %d and %d are identical. Halting — do not proceed with non-diverse imputed data.",
          i, j
        ))
      }
    }
  }

  invisible(TRUE)
}

# =============================================================================
# END OF FILE
# No results, fitted objects, or patient-level output are generated or
# stored by this script. It is intended to be sourced by the fold-level
# model-fitting scripts, not run standalone against public data.
# =============================================================================
