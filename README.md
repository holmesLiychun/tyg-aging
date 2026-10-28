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

## Output Map — NHANES

| Deliverable | Script | Key Output |
|---|---|---|
| Participant flowchart | `30_flowchart.R` | `flowchart_data.csv`, `supp_P0_exclusion_summary.txt` |
| Table 1 (baseline) | `40_descriptive.R` | `tables/table1_tyg_baseline.docx`, `tyg_tertile_distribution.csv` |
| BioAge construction (6 alg) | `20_bioage.R`, `25_phenoage.R` | KDM-11, PhenoAge-9, etc. residuals merged into `nhanes` |
| Main effects (6 alg × 5 models) | `50_regression.R` | `main_tyg_ba_effects.csv`, `main_effects_model5.csv` |
| Decoupling matrix | `60_sensitivity.R` §60.5 | `decoupling_matrix_standardized.csv` |
| Sensitivity analyses | `60_sensitivity.R` | ±3SD winsor, race-stratified, mediation |
| Stratified forest plot | `65b_stratified_analysis.R`, `65c_stratified_figures.R` | `forest_stratified_PhenoAge9_FINAL_FIX.png` |
| NRI/IDI | `66_nri_idi.R` | `table_S5_NRI_IDI.csv` (unweighted) |
| Correlation matrix | `67_corr.R` | `table_Sx_TyG_correlation.csv` |
| Nonlinearity & heatmap | `70_visualization.R` | `Figure_RCS_6in1.tiff`, `decoupling_heatmap.png` |

## Output Map — CHARLS (Longitudinal)

| Deliverable | Script Section | Key Output Files |
|---|---|---|
| Table 1 Part B (BA outcomes by TyG quartile) | Step 4b | `tables/Table1_PartB_BA_baseline_by_TyGq.csv`, `tables/Table1_complete.csv` |
| Cross-sectional decoupling (KDM_core6) | Step 5 | `tables/CS_TyG_BA_results.csv`, `figures/CS_Spline_KDMcore6_2in1.pdf` |
| Longitudinal LME (diagonal RE) | Step 6 | `tables/LME_main_unweighted.csv` (β for TyG_bl × time) |
| Time-varying TyG: GEE vs LME + within-person centering | Step 6b | `tables/Step6b_TimeVarying_GEE_LME.csv` |
| Sensitivity (paired, trimmed, minimal, sex-stratified) | Step 7 (7.1–7.5) | Results fed into Fig 3 forest plot |
| Cox survival (all-cause, CVD, cancer) | Step 8 | `05_cox_summary.csv`, `figures/Fig4_KM_by_TyG_quartile.pdf` |
| IPCW-weighted LME (3 weight models) | Step 9 | `tables/08_ipcw_results.csv`, `figures/Fig1_IPCW_Estimates.pdf` |
| Bootstrap (500 resamples, cluster-level) | Step 9b | `figures/Fig2_Bootstrap_Distribution.pdf` |
| Sensitivity forest plot | Step 9c | `tables/Fig3_forest_data.csv`, `figures/Fig3_Sensitivity_Forest.pdf` |
| Final export & follow-up summary | Step 10 | `01_dt_unit_long.csv`, `03_surv_final.csv`, `04_ipcw_data.csv`, `09_followup_summary.csv` |

> **Note:** CHARLS pipeline uses `KDM_core6` only (not 6-algorithm sweep as in NHANES). All longitudinal models use `BA_resid_KDM_core6` and `BA_resid_sex_KDM_core6`.

## Key Variable Coding

### CHARLS

| Variable | Coding | Note |
|---|---|---|
| `ID` | Character, subject identifier | Rows with NA/empty ID removed at pipeline start |
| `wave` | 1–4 | `time = wave - 1` for LME |
| `bl_fasting` | 1 = fasting, 2 = non-fasting | Analysis restricted to `== 1` |
| `gender` | 0 = female, 1 = male → `factor()` | |
| `tyg` | Continuous (log transform of TG×FPG) | Winsorized 1%–99% in `dt_unit` |
| `tyg_bl` | Baseline TyG (wave 1) | Standardized: `tyg_bl_sd = tyg_bl / SD` |
| `tyg_q` | Q1–Q4 (quartiles of `tyg`) | Reference = Q1 |
| `BA_resid_KDM_core6` | Residual from age ~ 6 markers | Sex-stratified; winsorized 1%–99% |
| `BA_resid_sex_KDM_core6` | Sex-specific KDM core6 residual | Sensitivity outcome |
| `R_final` | 1 = completed ≥3 waves, 0 = otherwise | IPCW attrition indicator |
| `death` | 1 = dead, 0 = alive/censored | From `radyear` |
| `cvd_death`, `cancer_death` | 1/0 | Cause-specific (ragcod: 2=CVD, 1=cancer) |
| `fu_years` | Follow-up (days/365.25) | From `base_dt` to `event_dt` |

### NHANES (key additions)

| Variable | Coding | Note |
|---|---|---|
| `tyg_tert` | T1/T2/T3 | T1 = reference |
| `kdm_accel`, `pheno_accel`, etc. | Residuals (Mean≈0) | 6 algorithms |

## Methodological Notes

- **CHARLS IPCW**: Three weight specifications (SW1–SW3) via logistic propensity for completing ≥3 waves; capped at 99th percentile.
- **CHARLS LME**: Diagonal random-effects structure `pdDiag(~1)`; `lmeControl(maxIter=1000, opt="nlminb")`.
- **CHARLS Bootstrap**: Cluster-level resampling (subject-level, 500 iterations, `set.seed(2024)`).
- **CHARLS Winsorize**: 1%–99% on all `BA_resid_*` variables in `dt_unit` (post-Step 4b, pre-Step 5).
- **NHANES NRI/IDI**: Unweighted (package limitation); declared as sensitivity.
- **NHANES weights**: `svyglm` with fasting weights + strata/PSU.

## License

Code: MIT. Data products: CC BY 4.0.

## Citation

See `CITATION.cff`.