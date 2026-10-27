# ============================================================================
# 模块 65b（修正版）：分层分析
# ============================================================================

timer_start("模块65b 分层分析（修正BMI）")

strat_results <- list()

# ---- 1. BMI 分层 ----
for (lvl in levels(main_dat$bmi_cat2)) {
  sub_idx <- main_dat$bmi_cat2 == lvl
  n_sub <- sum(sub_idx, na.rm = TRUE)
  if (n_sub < 100) { log_write(sprintf("  ⚠ BMI %s: n=%d, 跳过", lvl, n_sub)); next }
  
  ds_sub <- subset(main_ds, sub_idx)
  
  form_sub <- as.formula(
    "pheno9_accel ~ tyg_z + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
  )
  
  fit_sub <- tryCatch(survey::svyglm(form_sub, design = ds_sub, family = gaussian()), error = function(e) NULL)
  if (is.null(fit_sub)) next
  
  sm <- summary(fit_sub)$coefficients
  tyg_row <- grep("tyg_z", rownames(sm))
  strat_results[[length(strat_results) + 1]] <- data.frame(
    Stratification = "BMI category",
    Level = lvl,
    Algorithm = "PhenoAge-9",
    N = n_sub,
    Beta = round(sm[tyg_row, "Estimate"], 4),
    SE = round(sm[tyg_row, "Std. Error"], 4),
    P_value = round(sm[tyg_row, "Pr(>|t|)"], 4),
    stringsAsFactors = FALSE
  )
  log_write(sprintf("  ✓ BMI %s: β=%.4f, P=%.4f, n=%d", lvl, sm[tyg_row,"Estimate"], sm[tyg_row,"Pr(>|t|)"], n_sub))
}

# ---- 2. dm_status 分层 ----
for (lvl in levels(main_dat$dm_status)) {
  sub_idx <- main_dat$dm_status == lvl
  n_sub <- sum(sub_idx, na.rm = TRUE)
  if (n_sub < 100) { log_write(sprintf("  ⚠ DM %s: n=%d, 跳过", lvl, n_sub)); next }
  
  ds_sub <- subset(main_ds, sub_idx)
  
  for (algo in c("pheno9_accel", "pheno_accel")) {
    algo_name <- ifelse(algo == "pheno9_accel", "PhenoAge-9", "Modified PhenoAge")
    
    form_sub <- as.formula(
      paste(algo, "~ tyg_z + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified")
      # 注意：分层后dm_status不再需要作为协变量（因为层内同质）
    )
    
    fit_sub <- tryCatch(survey::svyglm(form_sub, design = ds_sub, family = gaussian()), error = function(e) NULL)
    if (is.null(fit_sub)) next
    
    sm <- summary(fit_sub)$coefficients
    tyg_row <- grep("tyg_z", rownames(sm))
    strat_results[[length(strat_results) + 1]] <- data.frame(
      Stratification = "Diabetes status",
      Level = lvl,
      Algorithm = algo_name,
      N = n_sub,
      Beta = round(sm[tyg_row, "Estimate"], 4),
      SE = round(sm[tyg_row, "Std. Error"], 4),
      P_value = round(sm[tyg_row, "Pr(>|t|)"], 4),
      stringsAsFactors = FALSE
    )
    log_write(sprintf("  ✓ DM %s, %s: β=%.4f, P=%.4f, n=%d", lvl, algo_name, sm[tyg_row,"Estimate"], sm[tyg_row,"Pr(>|t|)"], n_sub))
  }
}

# ---- 输出 ----
if (length(strat_results) > 0) {
  strat_df <- do.call(rbind, strat_results)
  write.csv(strat_df, file.path(OUT_DIR, "table_S4_stratified_analysis.csv"), row.names = FALSE)
  log_write("=== 分层分析汇总 ===")
  print(strat_df, row.names = FALSE)
  log_write(sprintf("✓ Table S4 已保存 (%d 行)", nrow(strat_df)))
}

timer_end("模块65b 分层分析")