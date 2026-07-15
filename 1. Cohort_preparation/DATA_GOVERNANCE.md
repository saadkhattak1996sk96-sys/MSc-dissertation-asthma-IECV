# Data Quality, Cleaning, and Governance Note

This file documents, at a methods level, how the OPCRD dataset was checked and
processed prior to model fitting. It is intended to sit alongside the analysis
code in this repository so that a reader understands what preprocessing
occurred without requiring access to the underlying patient-level data, which
is not included in this repository under the project's data governance
agreement.

## Scope of this note

This document describes **what was done and why**, at the level of a methods
section. It does not reproduce:
- patient-level data, counts, or identifiers,
- intermediate debugging output, or
- the step-by-step discovery process behind any individual fix.

Those are not appropriate for a public repository and are not included here.
This note exists so that the absence of that material is a deliberate,
documented choice rather than an unexplained gap.

## 1. Data governance

No patient-level data, derived datasets, or model outputs are stored in this
repository. All scripts operate on data held securely on the University of
Liverpool's Barkla HPC system and are not designed to execute outside that
environment. Scripts are included here for methodological transparency and
reproducibility of approach, not as a runnable pipeline against public data.

## 2. Data quality checks performed

Prior to model fitting, all continuous predictor variables (including BMI,
peak expiratory flow, and eosinophil count) were screened for
physiologically implausible values. Where a recorded value fell outside a
clinically plausible range for that variable, it was recoded to missing
rather than retained or silently excluded. This screening step is applied
within the cohort-construction code in this repository (see the `01_prep_cohort`
scripts) and is a standard, expected step in working with large routinely
collected primary care datasets, which are not curated for research use at
the point of entry.

## 3. Missing data handling

A subset of predictor variables had genuine missingness in the source data.
Missingness mechanisms were assessed empirically (comparing outcome rates
between complete and incomplete cases) prior to selecting a handling
strategy. Multiple imputation was adopted as the primary approach for the
IECV/Blakey analysis, with complete-case analysis retained as a sensitivity
comparison. Imputation was performed independently within each
cross-validation fold to avoid information leakage between development and
validation data. Imputation diversity (i.e., that generated imputed datasets
are genuinely distinct from one another) is explicitly checked in code and
will halt execution if a diversity check fails, rather than proceeding
silently.

## 4. Sample size adequacy

Per-fold sample size adequacy for the internal-external cross-validation
design was assessed using established sample size methodology for clinical
prediction model validation (`pmvalsampsize`), based on each outcome's
anticipated discrimination performance as reported in the source
development paper. Fold boundaries were selected to ensure adequate sample
size across all folds for both outcomes evaluated.

## 5. Model validation framework

Two complementary modelling frameworks are implemented and reported: a
primary framework in which predictor–outcome relationships are held fixed
across time periods with period-specific baseline risk only (following
Debray et al.'s IECV methodology), and a sensitivity framework in which the
full model is re-estimated within each time period. Both frameworks are
evaluated using the same internal-external cross-validation structure, and
pooled performance estimates are combined using Bayesian random-effects
meta-analysis, with a restricted maximum likelihood approach reported
separately as a sensitivity comparison.

## 6. A note on transparency

Data cleaning, missingness handling, and validation methodology are
described above at the level a reader needs to assess the soundness of the
approach. The process by which individual issues were identified during
development (e.g., specific software behaviour, intermediate diagnostic
output) is not included, consistent with standard practice for public
research code repositories, which document a finished, correct pipeline
rather than its development history.

---

*This note should be read alongside the dissertation's Methods section,
which provides full detail on cohort definition, outcome definitions, and
statistical methodology.*
