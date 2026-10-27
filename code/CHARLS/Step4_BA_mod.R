
# ============================================================
# Step 4: BA 投影（四版·分性别·残差统一）
# ============================================================
CH$reset()
if (!exists("dt_unit")) load_ck("dt_unit", "2_audit")
if (!exists("dt_base")) load_ck("dt_base", "2_audit")
source("ba_algorithms.R")

cat("--- Step 4: BA Projection ---\n")

# 合并 BA 容器
dt_ba <- dt_unit %>% dplyr::select(ID, wave, wave_year, age, gender)

# ---- 变量别名（让 NHANES 算法认得你的 bl_* 变量）----
if (!"bun" %in% names(dt_unit)) dt_unit$bun <- dt_unit$bl_bun
if (!"uap" %in% names(dt_unit)) dt_unit$uap <- dt_unit$bl_ua
if (!"wbc" %in% names(dt_unit)) dt_unit$wbc <- dt_unit$bl_wbc
if (!"mcv" %in% names(dt_unit)) dt_unit$mcv <- dt_unit$bl_mcv
if (!"sbp" %in% names(dt_unit)) dt_unit$sbp <- dt_unit$systo   # ← 关键！收缩压
if (!"creat" %in% names(dt_unit)) dt_unit$creat <- dt_unit$bl_crea
if (!"lncreat" %in% names(dt_unit)) dt_unit$lncreat <- ifelse(dt_unit$bl_crea > 0, log(dt_unit$bl_crea), NA_real_)


# ---- 投影 ----
for (v in names(BA_VERSIONS)) {
  cfg <- BA_VERSIONS[[v]]
  CH$log("INFO", sprintf("Projecting %s (%s)...", v, cfg$label), "4_BA")
  
  
  if (cfg$method == "kdm_core6_fixed") {
    # 加载 NHANES 系数
    if (!file.exists(cfg$coef_path)) {
      CH$log("FAIL", sprintf("%s: coefficient file not found: %s", v, cfg$coef_path), "4_BA")
      next
    }
    coef_list <- readRDS(cfg$coef_path)
    
    # 确保 lncreat 存在（如果 Step2 里只做了 ln_crea）
    if (!"lncreat" %in% names(dt_unit)) {
      dt_unit$lncreat <- ifelse(dt_unit$bl_crea > 0, log(dt_unit$bl_crea + 1), NA_real_)
    }
    
    raw <- calc_kdm_core6_fixed(dt_unit, coef_list, v)
    dt_unit[[paste0("BA_raw_", v)]] <- raw
    dt_unit <- add_resid(dt_unit, paste0("BA_raw_", v), v)
    next  # 跳过下面的常规 calc_kdm / calc_spa
  }

  

  
  # R2(age) of raw BA
  # 替换原来的 R2 报告

}

# ---- 互相关（自动找所有 BA_resid_* 列）----
rc <- names(dt_unit)[grepl("^BA_resid_", names(dt_unit))]

if (length(rc) >= 2) {
  cm <- cor(dt_unit[, rc], use = "pairwise.complete.obs")
  CH$log("OK", "Cross-version resid correlation:", "4_BA")
  print(round(cm, 3))
} else if (length(rc) == 1) {
  CH$log("WARN", sprintf("Only 1 residual column (%s), skip cross-correlation", rc), "4_BA")
} else {
  CH$log("FAIL", "No BA_resid_* columns found!", "4_BA")
}

# ---- 残差-年龄正交性诊断（也改成动态）----
for (v in rc) {
  valid <- complete.cases(dt_unit[[v]], dt_unit$age)
  if (sum(valid) > 10) {
    r <- cor(dt_unit[[v]][valid], dt_unit$age[valid])
    m <- mean(dt_unit[[v]][valid])
    cat(sprintf("%s resid-age r=%.4f mean=%.4f\n", v, r, m))
  } else {
    cat(sprintf("%s: insufficient data for correlation\n", v))
  }
}

# 只保留 KDM_core6 相关列
keep_cols <- c("ID", "wave", "wave_year", "age", "gender",
               "BA_raw_KDM_core6", "BA_resid_KDM_core6", "BA_resid_sex_KDM_core6")

# 如果你后续分析只需要这几个，可以精简 dt_unit
dt_ba_final <- dt_unit[, intersect(names(dt_unit), keep_cols)]

# 保存最终结果
#save_ck("4_BA_final", list(
#  dt_ba = dt_ba_final,
#  dt_unit = dt_unit   # 保留完整版以防万一
#))

#save_ck("4_BA", list(dt_unit=dt_unit, BA_VERSIONS=BA_VERSIONS))
