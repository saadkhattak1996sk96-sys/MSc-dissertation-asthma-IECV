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

## Repository structure
Added incrementally as each component is built and confirmed running on the
University of Liverpool's Barkla HPC cluster:
- `/prep/` — shared cohort construction (run once, output not tracked)
- `/prisma/` — PRISMA 2020 flow diagram
- `/baseline_table/` — baseline characteristics table (Table 1)
- `/noble_validation/` — external validation of Noble et al. (2021)

## Status
Repository created [12/07/2026] as part of testing the Barkla batch
submission pipeline ahead of the full IECV analysis. Code is added here only
after being confirmed to run successfully on Barkla.
