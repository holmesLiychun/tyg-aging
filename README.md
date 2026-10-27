# Decoupling Framework for Composite Metabolic Indices and Biological Aging

Analysis code for BMC Medicine submission:  
"Algorithmic collinearity, not biology? A decoupling framework for composite metabolic indices and biological aging clocks"

## Analysis Plan

Lock date: 2026-09-15 (pre-specified; no post-hoc modification)

## Repository Structure

- **code/CHARLS/** — Longitudinal IPCW-weighted mixed models (Steps 0–10)
- **code/NHANES/** — Cross-sectional decoupling & sensitivity (Steps 01–70)
- **code/MR/** — Mendelian randomization (two-sample & multivariable)
- **data/example/** — Synthetic test data (same schema, no real subjects)
- **data/derived/** — Intermediate datasets (reproducibly generated)
- **docs/** — Variable dictionaries (CHARLS & NHANES, with cross-references)
- **output/figures/** — Main and supplementary figures
- **output/tables/** — Regression and mediation tables

## Data Access

Raw CHARLS and NHANES data are restricted-access and **not included**.  
Apply at: https://charls.charlsdata.com (CHARLS) | CDC NHANES portal (NHANES).  
Place authorized files in `data/raw/` (git-ignored). Synthetic data in `data/example/` allows full pipeline execution without raw data.

## Environment

- R >= 4.5.0
- Package versions locked in `renv.lock` — run `renv::restore()`
- Session info: `sessionInfo.txt`

## Run Order
NHANES (cross-sectional)
source("code/NHANES/00_main_NHANES.R")
CHARLS (longitudinal)
source("code/CHARLS/00_run_all_CHARLS.R")
MR
source("code/MR/MR_analysis_v3.R")
source("code/MR/mvmr_fg_tg_phenoage.R")
## Output Map

| Deliverable | Script |
|---|---|
| Participant flowchart | `NHANES/30_flowchart.R` |
| Table 1 (baseline) | `NHANES/40_descriptive.R`, `CHARLS/Step2b_table1.R` |
| Main decoupling models | `NHANES/50_regression.R`, `CHARLS/Step4_BA_mod.R` |
| Sensitivity analyses | `NHANES/60_sensitivity.R` |
| Stratified forest plots | `NHANES/65b_stratified_analysis.R`, `NHANES/65c_stratified_figures.R` |
| NRI/IDI | `NHANES/66_nri_idi.R` |
| MR results | `MR/MR_analysis_v3.R` |

## Key Variable Coding

| Variable | Coding |
|---|---|
| `bl_fasting` | 1 = fasting, 2 = non-fasting |
| `gender` | 0 = female, 1 = male |
| `diabe` | 0 = no, 1 = yes |
| `bl_crea`, `bl_ua` | mg/dL |

## License

Code: MIT. Data products: CC BY 4.0.

## Citation

See `CITATION.cff`.