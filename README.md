# Tooth Loss and Liver Fibrosis Among US Adults With MASLD (NHANES 2017–March 2020)

Analysis code for the cross-sectional study *"Tooth Loss and Liver Fibrosis Among US Adults With Metabolic Dysfunction–Associated Steatotic Liver Disease: A Cross-Sectional Analysis of NHANES 2017–March 2020."*

**Author:** Kindson Abone, BDS

## Summary

Using nationally representative NHANES 2017–March 2020 data (n = 6,271 adults), this study examines whether clinically examined tooth loss and untreated dental caries are associated with metabolic dysfunction–associated steatotic liver disease (MASLD) and with liver fibrosis measured by vibration-controlled transient elastography.

## Data

All data are publicly available from the National Center for Health Statistics (<https://wwwn.cdc.gov/nchs/nhanes/>). No data use agreement is required. Raw and processed data files are **not** stored in this repository; `scripts/01_setup_and_download.R` downloads them directly from the CDC.

## How to reproduce

1. Install R (version 4.6.1 was used) and RStudio.
2. Open `dental-masld.Rproj` in RStudio.
3. Run the scripts in order:

| Script | What it does | Main outputs |
|---|---|---|
| `01_setup_and_download.R` | Installs packages, creates folders, downloads 22 NHANES files | `data/raw/` |
| `02_build_analytic_dataset.R` | Merges files, derives dental and liver variables, applies exclusions | `data/processed/nhanes_merged.rds`, `output/tables/flowchart_counts.csv` |
| `03_table1.R` | Weighted participant characteristics and missing-data summary | `output/tables/table1_*.csv` |
| `04_main_models.R` | Survey-weighted logistic regression (Models 1–4) and dose–response splines | `output/tables/main_models*.csv`, `output/figures/spline_*.png` |
| `05_sensitivity_aim3.R` | Nine sensitivity analyses, effect modification, discrimination (AUC), main figure | `output/tables/sensitivity.csv`, `aim3_*.csv`, `output/figures/spline_fibrosis_masld.png` |

Script 05 takes about 5–10 minutes because of multiple imputation.

## Key analytic choices

- Survey design: MEC examination weights (`WTMECPRP`), strata and PSUs; subgroups analyzed as domains.
- Steatosis: CAP ≥ 285 dB/m; MASLD: steatosis plus ≥ 1 cardiometabolic risk factor (2023 consensus).
- Significant fibrosis: liver stiffness ≥ 8.0 kPa.
- Tooth loss: 28 teeth, third molars excluded; implants and root fragments counted as missing.
- Untreated caries: carious coronal surface (code Z) or retained root tip (code J); dentate adults only.

## R packages

tidyverse, haven, survey, srvyr, gtsummary, rms, pROC, here, janitor, mice, mitools, splines.

## License

To be decided by the author before the repository is made public.
