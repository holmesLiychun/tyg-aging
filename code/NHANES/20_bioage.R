# ============================================================================
# 模块 20：生物年龄算法计算（KDM + PhenoAge，自训练线性 BA）
# ============================================================================
# 统一创建 uap（uric acid）别名
nhanes <- nhanes %>% mutate(uap = uric_acid)
timer_start("模块20 生物年龄算法计算")
# 清除可能污染的旧列
for (col in c("pheno_ba","pheno_accel","pheno9_ba","pheno9_accel",
              "kdm_core6_ba","kdm_core6_accel")) {
  if (col %in% names(nhanes) && sum(!is.na(nhanes[[col]])) < 100) {
    nhanes[[col]] <- NULL
    log_write(sprintf("  ♻ 清除无效列 %s", col))
  }
}


library(dplyr)
library(tidyr)

# ---- 辅助函数（与模块25共用逻辑）----
train_ba_accel <- function(data, formula_str, ba_name, accel_name) {
  f <- as.formula(formula_str)
  result_list <- list()
  
  for (sx in c("Male", "Female")) {
    sub <- data %>% filter(sex == sx)
    vars_needed <- all.vars(f)
    complete_idx <- complete.cases(sub[, vars_needed])
    sub_clean <- sub[complete_idx, ]
    
    if (nrow(sub_clean) < 50) {
      log_write(sprintf("  ⚠ %s: 仅 %d 行，跳过", sx, nrow(sub_clean)))
      next
    }
    
    fit <- lm(f, data = sub_clean)
    ba_pred <- predict(fit, newdata = sub_clean)
    accel <- residuals(fit)
    
    result_list[[sx]] <- data.frame(
      SEQN = sub_clean$SEQN,
      temp_ba = as.numeric(ba_pred),
      temp_accel = as.numeric(accel),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(result_list) == 0) return(NULL)
  
  bound <- do.call(rbind, result_list)
  colnames(bound)[colnames(bound) == "temp_ba"] <- ba_name
  colnames(bound)[colnames(bound) == "temp_accel"] <- accel_name
  
  log_write(sprintf("  ✓ %s: n=%d, BA mean=%.1f, Accel mean=%.4f, SD=%.2f",
                    ba_name, nrow(bound),
                    mean(bound[[ba_name]], na.rm=TRUE),
                    mean(bound[[accel_name]], na.rm=TRUE),
                    sd(bound[[accel_name]], na.rm=TRUE)))
  bound
}

# ============================================================================
# B.1 KDM-11（11 markers, 不含 CRP, 7 周期）———— 已跑通，保留原逻辑
# ============================================================================
log_write(">>> B.1 KDM-11（11 markers, 不含 CRP, 7 周期）")

# 检查是否已有缓存结果
if (!"kdm_ba" %in% names(nhanes)) {
  # 从模块10的衍生变量计算
  kdm_base <- nhanes %>%
    filter(age >= 20) %>%
    mutate(
      uap = uric_acid,
      lncreat = log(creatinine + 1),
      lncrp = log(crp + 1)
    )
  
  kdm_formula <- "age ~ albumin + alp + lncreat + hba1c + sbp + bun + uap + lymph + mcv + wbc + totchol"
  
  kdm_res <- train_ba_accel(
    kdm_base %>% filter(!is.na(albumin), !is.na(alp), !is.na(creatinine),
                        !is.na(hba1c), !is.na(sbp), !is.na(bun), !is.na(uap),
                        !is.na(lymph), !is.na(mcv), !is.na(wbc), !is.na(totchol)),
    formula_str = kdm_formula,
    ba_name = "kdm_ba",
    accel_name = "kdm_accel"
  )
  
  if (!is.null(kdm_res)) {
    nhanes <- nhanes %>%
      left_join(kdm_res %>% dplyr::select(SEQN, kdm_ba, kdm_accel), by = "SEQN")
  }
} else {
  log_write("  ℹ kdm_ba 已存在，跳过")
}

log_write(sprintf("  KDM-11 完成: n=%d, r(BA,age)=%.3f",
                  sum(!is.na(nhanes$kdm_ba)),
                  ifelse(sum(!is.na(nhanes$kdm_ba)) > 10,
                         cor(nhanes$kdm_ba, nhanes$age, use="complete.obs"), NA)))
# ============================================================================
# B.1b KDM_core6（6 markers, 不含 CRP/FPG/TG, 7 周期）———— 跨队列统一主算法
# ============================================================================
log_write(">>> B.1b KDM_core6（6 markers, 去代谢组分, 7 周期）")

if (!"kdm_core6_ba" %in% names(nhanes)) {
  kdm_core6_base <- nhanes %>%
    filter(age >= 20, exam_stat == 2) %>%
    mutate(
      lncreat = log(creatinine + 1),   # 与 KDM-11 保持一致
      # 其他变量原值：bun, uap(uric_acid), wbc, mcv, sbp
    )
  
  kdm_core6_formula <- "age ~ lncreat + bun + uap + wbc + mcv + sbp"
  
  kdm_core6_res <- train_ba_accel(
    kdm_core6_base %>% filter(!is.na(creatinine), !is.na(bun), !is.na(uap),
                              !is.na(wbc), !is.na(mcv), !is.na(sbp)),
    formula_str = kdm_core6_formula,
    ba_name = "kdm_core6_ba",
    accel_name = "kdm_core6_accel"
  )
  
  if (!is.null(kdm_core6_res)) {
    nhanes <- nhanes %>%
      left_join(kdm_core6_res %>% dplyr::select(SEQN, kdm_core6_ba, kdm_core6_accel), 
                by = "SEQN")
    
    # ★ 提取并保存固定系数（供 CHARLS 使用）
    kdm_core6_coef <- list()
    for (sx in c("Male", "Female")) {
      sub_coef <- kdm_core6_base %>%
        filter(sex == sx, !is.na(creatinine), !is.na(bun), !is.na(uap),
               !is.na(wbc), !is.na(mcv), !is.na(sbp))
      fit_coef <- lm(as.formula(kdm_core6_formula), data = sub_coef)
      kdm_core6_coef[[sx]] <- list(
        intercept   = coef(fit_coef)[1],
        coefficients = coef(fit_coef)[-1],
        markers     = names(coef(fit_coef))[-1],
        n_train     = nrow(sub_coef),
        r_squared   = summary(fit_coef)$r.squared
      )
    }
    saveRDS(kdm_core6_coef, file.path(OUT_DIR, "cache", "kdm_core6_coef_nhanes.rds"))
    
    # 打印系数摘要
    for (sx in names(kdm_core6_coef)) {
      log_write(sprintf("  %s: n=%d, R²=%.3f", 
                        sx, kdm_core6_coef[[sx]]$n_train, 
                        kdm_core6_coef[[sx]]$r_squared))
      log_write(sprintf("    Intercept=%.3f", kdm_core6_coef[[sx]]$intercept))
      log_write(sprintf("    Coefficients: %s", 
                        paste(names(kdm_core6_coef[[sx]]$coefficients), 
                              round(kdm_core6_coef[[sx]]$coefficients, 4), 
                              sep="=", collapse=", ")))
    }
  }
} else {
  log_write("  ℹ kdm_core6_ba 已存在，跳过")
}

log_write(sprintf("  KDM_core6 完成: n=%d, r(BA,age)=%.3f",
                  sum(!is.na(nhanes$kdm_core6_ba)),
                  ifelse(sum(!is.na(nhanes$kdm_core6_ba)) > 10,
                         cor(nhanes$kdm_core6_ba, nhanes$age, use="complete.obs"), NA)))
# ============================================================================
# B.2 Modified PhenoAge（不含 CRP, 7 周期）———— 修复
# ============================================================================
log_write(">>> B.2 Modified PhenoAge（不含 CRP, 7 周期）")

if (!"pheno_ba" %in% names(nhanes)) {
  pheno_mod_base <- nhanes %>%
    filter(age >= 20, exam_stat == 2) %>%
    mutate(
      albumin_gL   = albumin * 10,
      lncreat_umol = log(creatinine * 88.42 + 1)
    )
  
  mod_formula <- "age ~ albumin_gL + alp + lncreat_umol + hba1c + sbp + bun + uap + lymph + mcv + wbc + totchol + trig"
  
  pheno_mod_res <- train_ba_accel(
    pheno_mod_base %>% filter(!is.na(albumin), !is.na(alp), !is.na(creatinine),
                              !is.na(hba1c), !is.na(sbp), !is.na(bun), !is.na(uap),
                              !is.na(lymph), !is.na(mcv), !is.na(wbc),
                              !is.na(totchol), !is.na(trig)),
    formula_str = mod_formula,
    ba_name = "pheno_ba",
    accel_name = "pheno_accel"
  )
  
  if (!is.null(pheno_mod_res)) {
    nhanes <- nhanes %>%
      left_join(pheno_mod_res %>% dplyr::select(SEQN, pheno_ba, pheno_accel), by = "SEQN")
  }
} else {
  log_write("  ℹ pheno_ba 已存在，跳过")
}

n_pheno <- sum(!is.na(nhanes$pheno_ba))
log_write(sprintf("  Modified PhenoAge 完成: n=%d, r(BA,age)=%.3f",
                  n_pheno,
                  ifelse(n_pheno > 10, cor(nhanes$pheno_ba, nhanes$age, use="complete.obs"), NA)))

# ============================================================================
# B.3 PhenoAge-9 原版 Levine（含 CRP, 5 周期）———— 修复
# ============================================================================
log_write(">>> B.3 PhenoAge-9 原版（含 CRP, 5 周期）")

if (!"pheno9_ba" %in% names(nhanes)) {
  pheno9_base <- nhanes %>%
    filter(age >= 20, exam_stat == 2, cycle %in% CRP_CYCLES) %>%
    mutate(
      albumin_gL   = albumin * 10,
      lncreat_umol = log(creatinine * 88.42 + 1),
      lncrp_mgdl   = log(crp / 10 + 1)
    )
  
  pheno9_formula <- "age ~ albumin_gL + alp + lncreat_umol + I(glucose/18.0) + lncrp_mgdl + lymph + mcv + rdw + wbc"
  
  pheno9_res <- train_ba_accel(
    pheno9_base %>% filter(!is.na(albumin), !is.na(alp), !is.na(creatinine),
                           !is.na(glucose), !is.na(crp), !is.na(lymph),
                           !is.na(mcv), !is.na(rdw), !is.na(wbc)),
    formula_str = pheno9_formula,
    ba_name = "pheno9_ba",
    accel_name = "pheno9_accel"
  )
  
  if (!is.null(pheno9_res)) {
    nhanes <- nhanes %>%
      left_join(pheno9_res %>% dplyr::select(SEQN, pheno9_ba, pheno9_accel), by = "SEQN")
  }
} else {
  log_write("  ℹ pheno9_ba 已存在，跳过")
}

n_pheno9 <- sum(!is.na(nhanes$pheno9_ba))
log_write(sprintf("  PhenoAge-9 完成: n=%d (5 周期), r(BA,age)=%.3f",
                  n_pheno9,
                  ifelse(n_pheno9 > 10, cor(nhanes$pheno9_ba, nhanes$age, use="complete.obs"), NA)))

# ============================================================================
# B.4 KDM-9（含 CRP, 5 周期）———— 已跑通，保留
# ============================================================================
log_write(">>> B.4 KDM-9（含 CRP, 5 周期）")

if (!"kdm9_ba" %in% names(nhanes)) {
  kdm9_base <- nhanes %>%
    filter(age >= 20, cycle %in% CRP_CYCLES) %>%
    mutate(
      uap = uric_acid,
      lncreat = log(creatinine + 1),
      lncrp = log(crp + 1)
    )
  
  kdm9_formula <- "age ~ albumin + alp + lncreat + hba1c + sbp + bun + uap + lymph + mcv + wbc + totchol + lncrp"
  
  kdm9_res <- train_ba_accel(
    kdm9_base %>% filter(!is.na(albumin), !is.na(alp), !is.na(creatinine),
                         !is.na(hba1c), !is.na(sbp), !is.na(bun), !is.na(uap),
                         !is.na(lymph), !is.na(mcv), !is.na(wbc), !is.na(totchol),
                         !is.na(crp)),
    formula_str = kdm9_formula,
    ba_name = "kdm9_ba",
    accel_name = "kdm9_accel"
  )
  
  if (!is.null(kdm9_res)) {
    nhanes <- nhanes %>%
      left_join(kdm9_res %>% dplyr::select(SEQN, kdm9_ba, kdm9_accel), by = "SEQN")
  }
} else {
  log_write("  ℹ kdm9_ba 已存在，跳过")
}

n_kdm9 <- sum(!is.na(nhanes$kdm9_ba))
log_write(sprintf("  KDM-9 完成: n=%d (5 周期), r(BA,age)=%.3f",
                  n_kdm9,
                  ifelse(n_kdm9 > 10, cor(nhanes$kdm9_ba, nhanes$age, use="complete.obs"), NA)))

# ============================================================================
# 质量验证
# ============================================================================
log_write("--- BA 质量验证 ---")

ba_summary <- data.frame(
  Metric = c("KDM-11 (7 cycles)", "KDM_core6 (7 cycles)", 
             "Modified PhenoAge (7 cycles)", "PhenoAge-9 (5 cycles)", "KDM-9 (5 cycles)"),
  N = c(sum(!is.na(nhanes$kdm_ba)), sum(!is.na(nhanes$kdm_core6_ba)), 
        sum(!is.na(nhanes$pheno_ba)), sum(!is.na(nhanes$pheno9_ba)), sum(!is.na(nhanes$kdm9_ba))),
  Corr_with_age = c(
    ifelse(sum(!is.na(nhanes$kdm_ba)) > 10, round(cor(nhanes$kdm_ba, nhanes$age, use="complete.obs"), 3), NA),
    ifelse(sum(!is.na(nhanes$kdm_core6_ba)) > 10, round(cor(nhanes$kdm_core6_ba, nhanes$age, use="complete.obs"), 3), NA),
    ifelse(sum(!is.na(nhanes$pheno_ba)) > 10, round(cor(nhanes$pheno_ba, nhanes$age, use="complete.obs"), 3), NA),
    ifelse(sum(!is.na(nhanes$pheno9_ba)) > 10, round(cor(nhanes$pheno9_ba, nhanes$age, use="complete.obs"), 3), NA),
    ifelse(sum(!is.na(nhanes$kdm9_ba)) > 10, round(cor(nhanes$kdm9_ba, nhanes$age, use="complete.obs"), 3), NA)
  ),
  Accel_Mean = c(
    round(mean(nhanes$kdm_accel, na.rm=TRUE), 4),
    round(mean(nhanes$kdm_core6_accel, na.rm=TRUE), 4),
    round(mean(nhanes$pheno_accel, na.rm=TRUE), 4),
    round(mean(nhanes$pheno9_accel, na.rm=TRUE), 4),
    round(mean(nhanes$kdm9_accel, na.rm=TRUE), 4)
  ),
  Accel_SD = c(
    round(sd(nhanes$kdm_accel, na.rm=TRUE), 2),
    round(sd(nhanes$kdm_core6_accel, na.rm=TRUE), 2),
    round(sd(nhanes$pheno_accel, na.rm=TRUE), 2),
    round(sd(nhanes$pheno9_accel, na.rm=TRUE), 2),
    round(sd(nhanes$kdm9_accel, na.rm=TRUE), 2)
  ),
  stringsAsFactors = FALSE
)

print(ba_summary)

# 断言检查
for (accel_var in c("kdm_accel", "kdm_core6_accel", "pheno_accel", "pheno9_accel", "kdm9_accel")) {
  if (accel_var %in% names(nhanes)) {
    m <- mean(nhanes[[accel_var]], na.rm=TRUE)
    if (is.nan(m) || abs(m) > 0.01) {
      log_write(sprintf("  ⚠ 断言警告: %s Mean=%.4f（期望≈0）", accel_var, m))
    } else {
      log_write(sprintf("  ✓ %s Mean=%.4f（居中正常）", accel_var, m))
    }
  }
}

# ============================================================================
# 缓存
# ============================================================================
cp_save(nhanes, "20", "main_with_bioage")
timer_end("模块20 生物年龄算法计算")
log_write("===== 模块 20 结束 =====")