# MSc Dissertation — Asthma Exacerbation Prediction Models: External Validation and IECV

**Student:** Saad Ishaq Khattak (dentist by background)
**Degree:** MSc Health Data Science, University of Liverpool
**Primary Supervisor:** Dr Laura Bonnett
**Co-Supervisor:** Mr Thomas Spain

## Project summary

This repository holds analysis code for an MSc dissertation identifying and
externally validating clinical prediction models for asthma exacerbation risk.
The project has four components: a systematic review, a classical external
validation (Noble et al., 2021), and Internal-External Cross-Validation (IECV)
of Blakey et al.'s model under two complementary frameworks — a primary
framework in which predictor–outcome relationships are held fixed across time
periods with period-specific baseline risk only (following Debray et al.'s
IECV methodology), and a full-model-refit framework, run as a sensitivity
analysis, in which the model is re-estimated independently within each time
period. The IECV analysis is the dissertation's main analytical contribution.

Data source: Optimum Patient Care Research Database (OPCRD).

## Data governance

**No patient data, derived datasets, or model outputs are stored in this
repository.** All code operates on OPCRD data accessed under an approved data
sharing agreement. Data files are excluded via `.gitignore` and must never be
committed. Anyone running this code needs their own approved OPCRD access —
this repo does not provide or imply access to any dataset.

The one exception is `0. PRISMA/PRISMA_diagram.png` — a diagram of published
literature screening counts from the systematic review, not derived from
OPCRD or any patient-level data, so it is not subject to the restriction above.

See **`1. Cohort_preparation/1.1_DATA_GOVERNANCE.md`** for a methods-level
summary of the data quality checks, missingness handling, and imputation
approach used throughout the IECV analysis, and
`1.2_Data_quality_and_preprocessing.R` in the same folder for the
corresponding implementation. Consistent with a deliberate, logged project
decision, published code reflects the final, corrected pipeline — it does
not include the development-time debugging or discovery trail behind
individual fixes.

## Repository structure

- **`0. PRISMA/`** — PRISMA 2020 flow diagram for the systematic review
  (`PRISMA.Rmd`, `PRISMA_diagram.png`). Runs locally in RStudio; does not
  touch OPCRD data or require Barkla.

- **`1. Cohort_preparation/`** — cohort construction and data quality
  screening (`1.1_DATA_GOVERNANCE.md`, `1.2_Data_quality_and_preprocessing.R`).
  Loads OPCRD, filters to the active-asthma analysis cohort, screens
  continuous predictors for physiologically implausible values, and handles
  the ICS-adherence missing-data sentinel. Feeds into all downstream
  analyses below.

- **`2. baseline_characteristics_table/`** — baseline characteristics table,
  Table 1 (`baseline_table.Rmd`). Exports to Word via a standalone pandoc
  binary (Barkla has no system-wide pandoc). Confirmed running via SLURM
  batch job on Barkla (job 9762366, 12 July 2026).

- **`3. noble_external_validation/`** — external validation of Noble et al.
  (2021) (`noble_validation.Rmd`). Constructs Noble's 14 predictors and a
  substitute outcome (OPCRD lacks Noble's true outcome variables),
  reproduces Noble's own published worked example as a correctness check,
  and reports discrimination and calibration against Noble's reported
  figures. Retained as a worked example demonstrating why classical external
  validation is insufficient for this literature, rather than as the
  dissertation's primary contribution. Confirmed running via SLURM batch job
  on Barkla (job 9763685, 12 July 2026).

- **`4. Sample_size_calculation_IECV/`** — per-fold sample size adequacy for
  the IECV design (`4.1_pmvalsampsize_ge2_ge4_v2.R`,
  `4.2_ge4_fold1_boundary_shift_test.R`,
  `4.3_pmvalsampsize_ge2_ge4_v3_confirmed.R`). Documents the original
  (v2) fold boundary check, the GE4 Fold 1 sample-size failure, the
  boundary-shift fix, and the final (v3) confirmation that all fold/outcome
  combinations are adequate. All later analyses use the v3 fold structure
  confirmed here.

- **`5. Debray Framework, IECV, and Meta Analysis (Main Analysis)/`** — the
  dissertation's primary analytical contribution, true Debray framework,
  both outcomes:
  - **`GE2/`** — `5.1_fit_ge2_debray_all_folds.R`,
    `5.2_pool_ge2_debray_bayesian_primary.R`,
    `5.3_pool_ge2_debray_sensitivity_reml_uniform.R`,
    `5.4_forest_plots_ge2_debray.R`, `5.5_ge2_debray_meta_regression.R`,
    `5.6_fit_ge2_debray_final_model.R` — the last script builds the final
    combined model (all folds, none held out), the actual object handed
    over as the dissertation's GE2 deliverable.
  - **`GE4/`** — mirrors the GE2 structure and naming convention exactly
    (`5.1`–`5.6`).

- **`6. Full Model Refit, IECV, and Meta Analysis (Sensitivity Analysis)/`**
  — the full-refit sensitivity analysis, both outcomes:
  - **`GE2/`** — `6.1_fit_ge2_all_folds.R`, `6.2_pool_ge2_bayesian_primary.R`,
    `6.3_pool_ge2_sensitivity_reml_uniform.R`, `6.4_forest_plots_ge2.R`,
    `6.5_ge2_meta_regression.R`.
  - **`GE4/`** — `6.1_fit_ge4_all_folds.R`,
    `6.2_pooled_ge4_bayesian_primary.R`,
    `6.3_pool_ge4_sensitivity_reml_uniform.R`, `6.4_forest_plots_ge4.R`,
    `6.5_ge4_meta_regression.R`. Full refit deliberately does not include an
    equivalent to `5.6`'s final combined model: full refit's purpose is to
    test whether letting coefficients vary freely by era changes the
    answer, a question already fully answered by the fold-level and pooled
    comparison against the primary (Debray) analysis. A no-fold-held-out
    version would have no principled connection to that question.

All `.Rmd`/`.R` files in this repository document code that has already been
run and confirmed on Barkla — they are not scripts meant to be re-executed
directly from this repository. File paths inside them reference the
author's own Barkla home directory and will not resolve on another machine
without adjustment.

## Status

All four outcome/framework combinations for the IECV analysis — GE2 Debray,
GE4 Debray, GE2 full refit, GE4 full refit — are built, run at full
production scale on Barkla, and cross-checked against one another. This
completes the analysis-build phase of the project. 
