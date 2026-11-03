
# ============================================================================
# 60.5b ★ DMRC 解耦：TyG 残差化（去除 TG + glucose，按性别分层）
# ============================================================================
timer_start("模块60.5b DMRC解耦")
log_write(">>> 60.5b DMRC 解耦：TyG 残差化（去 TG + glucose）")

# ---- 1. 构造解耦基底 ----
dmrc_base <- ana_kdmcore6 %>%
  filter(!is.na(tyg), !is.na(age), !is.na(sex))

log_write(sprintf("  基底 n=%d（去缺 tyg/age/sex 后）", nrow(dmrc_base)))

# ---- 2. 从 nhanes 重新拉 trig / glucose（强制覆盖，不管 ana_kdmcore6 里有没有）----
dmrc_base <- dmrc_base %>%
  dplyr::select(-trig, -glucose) %>%          # 先删掉旧的（如果是全NA）
  left_join(
    nhanes %>% dplyr::select(SEQN, trig, glucose), 
    by = "SEQN"
  )

log_write(sprintf("  join 后: trig 非NA=%d, glucose 非NA=%d / 共 %d 行",
                  sum(!is.na(dmrc_base$trig)),
                  sum(!is.na(dmrc_base$glucose)),
                  nrow(dmrc_base)))

# 过滤缺值
dmrc_base <- dmrc_base %>%
  filter(!is.na(trig), !is.na(glucose))

log_write(sprintf("  去缺后 n=%d", nrow(dmrc_base)))

if (nrow(dmrc_base) < 100) {
  stop("  ❌ dmrc_base 不足100行，请检查 SEQN 匹配情况")
}

# ---- 3. 按性别分层残差化 ----
dmrc_base$tyg_decoupled <- NA_real_

for (sx in c("Male", "Female")) {
  idx <- dmrc_base$sex == sx
  sub_sx <- dmrc_base[idx, ]
  
  if (nrow(sub_sx) < 50) {
    log_write(sprintf("  ⚠ %s: 仅 %d 行，跳过", sx, nrow(sub_sx)))
    next
  }
  
  fit_dmrc <- lm(tyg ~ trig + glucose + age, data = sub_sx)
  dmrc_base$tyg_decoupled[idx] <- residuals(fit_dmrc)
  
  log_write(sprintf("  ✓ %s: n=%d, R²=%.3f, 残差范围 [%.3f, %.3f]",
                    sx, nrow(sub_sx), summary(fit_dmrc)$r.squared,
                    min(residuals(fit_dmrc)), max(residuals(fit_dmrc))))
}

# Z-score 标准化
dmrc_base$tyg_decoupled_z <- as.numeric(
  scale(dmrc_base$tyg_decoupled, center = TRUE, scale = TRUE)
)

# 原始 tyg_z 也放进基底
dmrc_base$tyg_z <- as.numeric(scale(dmrc_base$tyg, center = TRUE, scale = TRUE))

log_write("✓ 残差化完成")

# ---- 3. 构建 available_algos 子集（不依赖全局）----
# 从 ana_kdmcore6 出发，join 其他算法的 accel
dmrc_algos <- list(
  "KDM_core6"         = "kdm_core6_accel",
  "KDM-11"           = "kdm_accel",
  "KDM-9"            = "kdm9_accel",
  "Modified PhenoAge"= "pheno_accel",
  "PhenoAge-9"       = "pheno9_accel",
  "Pheno-noG"        = "pheno_nog_accel",
  "Pheno-noCRP"      = "phenocrp_accel"
)

# join 所有 accel 到 dmrc_base
for (nm in names(dmrc_algos)) {
  accel_col <- dmrc_algos[[nm]]
  if (accel_col %in% names(dmrc_base)) next
  
  # 尝试从对应的 ana_ 对象中取
  ana_obj_name <- switch(nm,
                         "KDM_core6"          = "ana_kdmcore6",
                         "KDM-11"            = "ana_kdm",
                         "KDM-9"             = "ana_kdm9",
                         "Modified PhenoAge" = "ana_pheno",
                         "PhenoAge-9"        = "ana_pheno9",
                         "Pheno-noG"         = "ana_pheno_nog",
                         "Pheno-noCRP"       = "ana_phenocrp"
  )
  
  if (exists(ana_obj_name)) {
    d_src <- get(ana_obj_name)
    if (accel_col %in% names(d_src)) {
      dmrc_base <- dmrc_base %>%
        left_join(d_src %>% dplyr::select(SEQN, all_of(accel_col)), by = "SEQN")
      log_write(sprintf("  ✓ %s: %s 已 join", nm, accel_col))
    } else {
      log_write(sprintf("  ⚠ %s: %s 不在 %s", nm, accel_col, ana_obj_name))
    }
  } else {
    log_write(sprintf("  ⚠ %s: %s 不存在", nm, ana_obj_name))
  }
}

# ---- 4. 跑回归：Before vs After（Model 5 协变量）----
dmrc_covs <- c("age", "sex", "race", "bmi", "smoke", "alcohol", 
               "paq_active_unified", "edu_cat")

dmrc_results <- list()
counter_dmrc <- 0

for (expo_type in c("Before", "After")) {
  expo_var <- ifelse(expo_type == "Before", "tyg_z", "tyg_decoupled_z")
  expo_label <- ifelse(expo_type == "Before", 
                       "Before decoupling (Full TyG)",
                       "After decoupling (TG+glucose removed)")
  
  for (nm in names(dmrc_algos)) {
    accel_col <- dmrc_algos[[nm]]
    
    if (!accel_col %in% names(dmrc_base)) next
    
    # 筛选完整病例
    keep_cols <- c("SEQN", expo_var, accel_col, "weight_fasting", "psu", "strata", dmrc_covs)
    keep_cols <- keep_cols[keep_cols %in% names(dmrc_base)]
    
    sub_data <- dmrc_base %>%
      dplyr::select(all_of(keep_cols)) %>%
      filter(!is.na(.data[[accel_col]]), !is.na(.data[[expo_var]])) %>%
      na.omit()
    
    if (nrow(sub_data) < 100) next
    
    # 构造 survey design
    ds_dmrc <- survey::svydesign(
      id = ~psu, strata = ~strata,
      weights = ~weight_fasting, nest = TRUE, data = sub_data
    )
    
    form_dmrc <- as.formula(
      paste(accel_col, "~", expo_var, "+", paste(dmrc_covs, collapse = " + "))
    )
    
    fit_dmrc <- tryCatch(
      survey::svyglm(form_dmrc, design = ds_dmrc, family = gaussian()),
      error = function(e) NULL
    )
    
    if (is.null(fit_dmrc)) next
    
    sm <- summary(fit_dmrc)$coefficients
    if (!expo_var %in% rownames(sm)) next
    
    counter_dmrc <- counter_dmrc + 1
    dmrc_results[[counter_dmrc]] <- data.frame(
      Exposure = expo_label,
      ExposureVar = expo_var,
      Outcome = nm,
      OutcomeVar = accel_col,
      Model = "Model 5 (DMRC)",
      Beta = round(sm[expo_var, "Estimate"], 4),
      SE = round(sm[expo_var, "Std. Error"], 4),
      CI_low = round(sm[expo_var, "Estimate"] - 1.96 * sm[expo_var, "Std. Error"], 4),
      CI_high = round(sm[expo_var, "Estimate"] + 1.96 * sm[expo_var, "Std. Error"], 4),
      P_value = round(sm[expo_var, "Pr(>|t|)"], 4),
      N_unweighted = nrow(sub_data),
      stringsAsFactors = FALSE
    )
  }
}

# ---- 5. 整理输出 ----
if (length(dmrc_results) > 0) {
  dmrc_combined <- do.call(rbind, dmrc_results)
  
  # 排序因子
  dmrc_combined$Outcome <- factor(dmrc_combined$Outcome, 
                                  levels = c("KDM_core6", "KDM-11", "KDM-9", "Modified PhenoAge", 
                                             "PhenoAge-9", "Pheno-noG", "Pheno-noCRP"))
  dmrc_combined$Exposure <- factor(dmrc_combined$Exposure,
                                   levels = c("Before decoupling (Full TyG)", "After decoupling (TG+glucose removed)"))
  
  write.csv(dmrc_combined, 
            file.path(OUT_DIR, "dmrc_before_after_comparison.csv"), 
            row.names = FALSE)
  
  log_write(sprintf("✓ DMRC 解耦前后对比数据已保存 (%d 行)", nrow(dmrc_combined)))
  
  # 打印 Model 5 汇总
  log_write("--- DMRC Model 5 β 汇总 ---")
  print(dmrc_combined %>% 
          dplyr::select(Outcome, Exposure, Beta, CI_low, CI_high, P_value, N_unweighted))
  
  # 计算衰减百分比
  # ---- 5b. 计算 β 衰减百分比 ----
  decay <- dmrc_combined %>%
    dplyr::select(Outcome, Exposure, Beta) %>%
    tidyr::pivot_wider(names_from = Exposure, values_from = Beta) %>%
    dplyr::rename(
      Beta_Before = `Before decoupling (Full TyG)`,
      Beta_After  = `After decoupling (TG+glucose removed)`
    ) %>%
    mutate(
      Decay_pct = round(100 * (Beta_Before - Beta_After) / abs(Beta_Before), 1),
      Direction_flip = case_when(
        Beta_Before < 0 & Beta_After > 0 ~ "Flip to positive",
        Beta_Before > 0 & Beta_After < 0 ~ "Flip to negative",
        TRUE ~ "Same direction"
      )
    )
  
  log_write("--- β 衰减百分比 + 方向翻转 ---")
  print(decay %>% dplyr::select(Outcome, Beta_Before, Beta_After, Decay_pct, Direction_flip))
  
} else {
  log_write("  ⚠ DMRC 未产生任何结果")
}

#cp_save(dmrc_combined, "60.5b", "dmrc_before_after")
timer_end("模块60.5b DMRC解耦")
log_write("===== 模块 60.5b 结束 =====\n")