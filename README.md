# MSc Dissertation — Asthma Exacerbation Prediction Models: External Validation and IECV

**Student:** Saad Ishaq Khattak (dentist by background)
**Degree:** MSc Health Data Science, University of Liverpool
**Primary Supervisor:** Dr Laura Bonnett
**Co-Supervisor:** Mr Thomas Spain

## Project summary
This repository holds analysis code for an MSc dissertation identifying and
externally validating clinical prediction models for asthma exacerbation risk.
The project has four components: a systematic review, a classical external
validation (Noble et al., 2021), Blakey intercept recalibration, and
Internal-External Cross-Validation (IECV) with Debray pooling — the main
analytical contribution.

Data source: Optimum Patient Care Research Database (OPCRD).

## Data governance
**No patient data, derived datasets, or model outputs are stored in this
repository.** All code operates on OPCRD data accessed under an approved data
sharing agreement. Data files are excluded via `.gitignore` and must never be
committed. Anyone running this code needs their own approved OPCRD access —
this repo does not provide or imply access to any dataset.

The one exception is `04_prisma/PRISMA_diagram.png` — a diagram of published
literature screening counts from the systematic review, not derived from
OPCRD or any patient-level data, so it is not subject to the restriction above.

## Repository structure

- **`01_prisma/`** — PRISMA 2020 flow diagram for the systematic review
  (`prisma.Rmd`, `PRISMA_diagram.png`). Runs locally in RStudio; does not
  touch OPCRD data or require Barkla.

- **`02_baseline_table/`** — baseline characteristics table, Table 1
  (`baseline_table.Rmd`). Contains the full pipeline as a single combined
  document: Steps 1-7 (cohort construction — loading OPCRD, filtering,
  recoding, labelling, deriving the PWP-based exacerbation stratification)
  followed by Steps 8-9 (building and exporting Table 1). On Barkla, Steps
  1-7 and Steps 8-9 were run as two separate scripts (`prep_cohort.R` and
  `baseline_table.R`), not one — this combined `.Rmd` is a documentation
  convenience, not a literal transcript of two separate executions. Exports
  to Word via a standalone pandoc binary (Barkla has no system-wide pandoc).
  Confirmed running via SLURM batch jobs on Barkla — prep (job 9762235) and
  baseline table (job 9762366) — 12 July 2026.

- **`03_noble_validation/`** — external validation of Noble et al. (2021)
  (`noble_validation.R`). Loads `year1.rds` and `pwp.rds`. Constructs
  Noble's 14 predictors and a substitute outcome (OPCRD lacks Noble's true
  outcome variables), reproduces Noble's own published worked example as a
  correctness check, and reports discrimination and calibration against
  Noble's reported figures. Retained as a worked example demonstrating why
  classical external validation is insufficient for this literature, rather
  than as the dissertation's primary contribution. Confirmed running via
  SLURM batch job on Barkla (job 9763685, 12 July 2026).

All `.Rmd` files in this repository use `eval=FALSE` throughout — they are
records of code that has already been run and confirmed elsewhere (Barkla
or local RStudio), not scripts meant to be re-executed from this repository.
File paths inside them reference the author's own Barkla home directory and
will not resolve on another machine without adjustment.

## Status
Repository created 12 July 2026 as part of testing the Barkla batch
submission pipeline ahead of the full IECV analysis. All three components
above are written, tested, and confirmed working as of 12 July 2026.
Remaining work: Blakey intercept recalibration and IECV with Debray pooling
(the dissertation's primary analytical contribution) — not yet started in
this repository.
