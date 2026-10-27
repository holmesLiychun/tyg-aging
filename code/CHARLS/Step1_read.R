# ============================================================
# Step 1: 读入 CHARLS 各波数据
# ============================================================
CH$reset()
source("Step0_config.R")

cat("--- Step 1: Data Loading ---\n")
# ---- Step 1: Data Loading ----
ck_1_path <- file.path(PATH_RDATA, "ck_1_read_dt_unit.rds")

if (file.exists(ck_1_path)) {
  cat("  Found Step 1 checkpoint, loading (skip reading DTA)...\n")
  load_ck("dt_unit", "1_read")
  cat(sprintf("  ✓ dt_unit loaded from checkpoint (%d rows, %d unique IDs)\n", 
              nrow(dt_unit), n_distinct(dt_unit$ID)))
  
} else {
  cat("  Reading DTA files...\n")
  dt_raw <- read_dta(PATH_DATA)
  CH$log("OK", sprintf("Loaded raw: %d rows, %d cols", nrow(dt_raw), ncol(dt_raw)))
  
  # 长表确认
  dt_unit <- dt_raw %>%
    mutate(
      wave = as.integer(wave),
      wave_year = as.integer(WAVE_YEAR[as.character(wave)]),
      time_years = as.numeric(TIME_YEARS[as.character(wave)]),
      ID = as.character(ID)
    ) %>%
    filter(age >= 45)
  
  CH$log("OK", sprintf("dt_unit: %d rows, %d unique IDs", nrow(dt_unit), n_distinct(dt_unit$ID)))
  
  # 变量存在性审计
  all_needed <- unique(unlist(lapply(BA_VERSIONS, function(x) x$markers)))
  all_needed <- c(all_needed, "systo", "gender", "age", "bmi", "smokev", 
                  "drinkl", "edu", "exercise", "bl_fasting", "DIED")
  
  orig_vars <- all_needed[!grepl("^ln_|^tg_|^TyG", all_needed)]
  missing_vars <- setdiff(orig_vars, names(dt_unit))
  
  if (length(missing_vars) > 0) {
    CH$log("WARN", sprintf("Missing RAW variables: %s", paste(missing_vars, collapse=", ")))
  } else {
    CH$log("OK", "All raw audit variables present")
  }
  
  # 保存 checkpoint
  save_ck("dt_unit", "1_read")
  cat("  ✓ Checkpoint saved: ck_1_read_dt_unit.rds\n")
}

