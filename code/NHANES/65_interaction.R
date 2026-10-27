# ============================================================================
# 模块 65：交互分析（TyG × 效应修饰因子 on 生物年龄加速）
# KDM-9 主分析，KDM_core6 敏感性末位
# ============================================================================

timer_start("模块65 交互分析")

library(dplyr)
library(survey)

# ---- 1. 确认主分析 survey design ----
# ds_kdm 是主分析对象 (n=16,018)
if (!exists("ds_kdm")) {
  stop("ds_kdm 不存在，请检查环境")
}

main_ds <- ds_kdm
main_dat <- main_ds$variables
log_write(sprintf("✓ 使用 ds_kdm, n=%d", nrow(main_dat)))

# 确认关键变量
required_vars <- c("tyg_z", "pheno9_accel", "pheno_accel", 
                   "kdm9_accel", "kdm_accel", "kdm_core6_accel",
                   "age", "sex", "bmi", "dm_status")
missing_vars <- setdiff(required_vars, names(main_dat))
if (length(missing_vars) > 0) {
  log_write(sprintf("⚠ 缺失变量: %s", paste(missing_vars, collapse=", ")))
} else {
  log_write("✓ 所有关键变量存在")
}

# ---- 2. 创建分类变量（如果不存在）----
if (!"age_grp" %in% names(main_dat)) {
  main_dat$age_grp <- cut(main_dat$age,
                          breaks = c(0, 40, 50, 60, 100),
                          labels = c("<40", "40-49", "50-59", "≥60"))
  log_write("✓ 创建 age_grp")
}

if (!"bmi_cat" %in% names(main_dat)) {
  main_dat$bmi_cat <- cut(main_dat$bmi,
                          breaks = c(0, 18.5, 25, 30, 100),
                          labels = c("<18.5", "18.5-24.9", "25-29.9", "≥30"))
  log_write("✓ 创建 bmi_cat")
}

# 修正 BMI 分组：<25 为参考（用于稳健性验证）
main_dat$bmi_cat2 <- cut(main_dat$bmi,
                         breaks = c(0, 25, 30, 100),
                         labels = c("<25 (ref)", "25-29.9", "≥30"))
log_write("✓ 创建 bmi_cat2 (<25 ref)")

# 重新绑定
main_ds$variables <- main_dat

# ---- 3. 效应修饰因子 & 结局清单 ----
effect_mods <- c("sex", "age_grp", "bmi_cat", "dm_status")

# ★ 7 算法，KDM-9 主位，KDM_core6 敏感性末位
outcomes <- c("pheno9_accel", "pheno_accel", "kdm9_accel", "kdm_accel",
              "pheno_nog_accel", "phenocrp_accel", "kdm_core6_accel")
outcome_names <- c("PhenoAge-9", "Modified PhenoAge", "KDM-9", "KDM-11",
                   "Pheno-noG", "Pheno-noCRP", "KDM_core6")

# 验证哪些结局可用
available_outcomes <- list()
for (j in seq_along(outcomes)) {
  if (outcomes[j] %in% names(main_dat)) {
    available_outcomes[[outcome_names[j]]] <- outcomes[j]
  } else {
    log_write(sprintf("  ⚠ %s (%s) 不在 main_dat，跳过", outcome_names[j], outcomes[j]))
  }
}
log_write(sprintf("  可用结局: %d 个", length(available_outcomes)))

# ---- 4. 跑交互 ----
interaction_results <- list()

for (i in seq_along(effect_mods)) {
  mod_var <- effect_mods[i]
  log_write(sprintf("--- 效应修饰因子: %s ---", mod_var))
  
  for (algo_name in names(available_outcomes)) {
    algo <- available_outcomes[[algo_name]]
    
    # 构造公式（pir 缺失时自动降级为不含 pir 的版本）
    cov_str <- "age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
    
    # 检查 pir 是否存在，不存在则去掉
    if (!"pir" %in% names(main_dat)) {
      cov_str <- "age + sex + race + edu_cat + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
    }
    
    form_int <- as.formula(
      paste(algo, "~ tyg_z *", mod_var, "+", cov_str)
    )
    
    fit_int <- tryCatch(
      survey::svyglm(form_int, design = main_ds, family = gaussian()),
      error = function(e) {
        log_write(sprintf("  ✗ %s × %s 报错: %s", algo_name, mod_var, e$message))
        NULL
      })
    
    if (is.null(fit_int)) next
    
    sm <- summary(fit_int)$coefficients
    int_rows <- grep(":", rownames(sm))
    
    if (length(int_rows) == 0) {
      log_write(sprintf("  ⚠ %s × %s: 未找到交互项", algo_name, mod_var))
      next
    }
    
    for (r in int_rows) {
      term_name <- rownames(sm)[r]
      p_val <- sm[r, "Pr(>|t|)"]
      
      interaction_results[[length(interaction_results) + 1]] <- data.frame(
        Effect_Modifier = mod_var,
        Algorithm = algo_name,
        Interaction_Term = term_name,
        Beta = round(sm[r, "Estimate"], 4),
        SE = round(sm[r, "Std. Error"], 4),
        P_interaction = round(p_val, 4),
        stringsAsFactors = FALSE
      )
    }
    log_write(sprintf("  ✓ %s × %s 完成", algo_name, mod_var))
  }
}

# ---- 5. 输出 ----
if (length(interaction_results) > 0) {
  int_df <- do.call(rbind, interaction_results)
  
  write.csv(int_df, file.path(OUT_DIR, "table_S3_interaction_pvalues.csv"),
            row.names = FALSE)
  
  log_write("=== 交互分析汇总 ===")
  print(int_df[, c("Effect_Modifier", "Algorithm", "Interaction_Term", "P_interaction")],
        row.names = FALSE)
  log_write(sprintf("✓ Table S3 已保存 (%d 行)", nrow(int_df)))
  
  # ★ 单独导出 KDM_core6 交互结果
  int_kdmcore6 <- int_df %>% filter(Algorithm == "KDM_core6")
  if (nrow(int_kdmcore6) > 0) {
    write.csv(int_kdmcore6,
              file.path(OUT_DIR, "table_S3_interaction_KDMcore6_only.csv"),
              row.names = FALSE)
    log_write(sprintf("✓ KDM_core6 交互单独导出 (%d 行)", nrow(int_kdmcore6)))
  }
  
  # 标记显著的
  sig_int <- int_df[int_df$P_interaction < 0.05, ]
  if (nrow(sig_int) > 0) {
    log_write(sprintf("⚠ %d 个交互项 P < 0.05:", nrow(sig_int)))
    print(sig_int[, c("Effect_Modifier", "Algorithm", "P_interaction")], row.names = FALSE)
    log_write("→ 建议对这些显著交互做分层分析（森林图）")
  } else {
    log_write("✓ 所有交互 P ≥ 0.05，无显著效应修饰")
    log_write("→ 正文可写：No significant effect modification was observed.")
  }
} else {
  log_write("⚠ 未生成任何交互结果")
}

# ============================================================================
# 6. BMI 分组稳健性验证（主算法 PhenoAge-9 + KDM-9）
# ============================================================================
log_write(">>> BMI 分组稳健性验证（bmi_cat2）")

# PhenoAge-9
form_bmi_p9 <- as.formula(
  "pheno9_accel ~ tyg_z * bmi_cat2 + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
)
fit_bmi_p9 <- tryCatch(survey::svyglm(form_bmi_p9, design = main_ds, family = gaussian()),
                       error = function(e) NULL)
if (!is.null(fit_bmi_p9)) {
  sm_p9 <- summary(fit_bmi_p9)$coefficients
  bmi_int_p9 <- sm_p9[grep("bmi_cat2|tyg_z:bmi", rownames(sm_p9)),
                      c("Estimate", "Std. Error", "Pr(>|t|)")]
  log_write("  PhenoAge-9 × bmi_cat2 交互:")
  print(round(bmi_int_p9, 4))
}

# KDM-9（主算法验证）
form_bmi_kdm9 <- as.formula(
  "kdm9_accel ~ tyg_z * bmi_cat2 + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
)
fit_bmi_kdm9 <- tryCatch(survey::svyglm(form_bmi_kdm9, design = main_ds, family = gaussian()),
                         error = function(e) NULL)
if (!is.null(fit_bmi_kdm9)) {
  sm_k9 <- summary(fit_bmi_kdm9)$coefficients
  bmi_int_k9 <- sm_k9[grep("bmi_cat2|tyg_z:bmi", rownames(sm_k9)),
                      c("Estimate", "Std. Error", "Pr(>|t|)")]
  log_write("  KDM-9 × bmi_cat2 交互:")
  print(round(bmi_int_k9, 4))
}

# KDM_core6（敏感性验证）
form_bmi_kdmcore6 <- as.formula(
  "kdm_core6_accel ~ tyg_z * bmi_cat2 + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
)
fit_bmi_kdmcore6 <- tryCatch(survey::svyglm(form_bmi_kdmcore6, design = main_ds, family = gaussian()),
                             error = function(e) NULL)
if (!is.null(fit_bmi_kdmcore6)) {
  sm_kc <- summary(fit_bmi_kdmcore6)$coefficients
  bmi_int_kc <- sm_kc[grep("bmi_cat2|tyg_z:bmi", rownames(sm_kc)),
                      c("Estimate", "Std. Error", "Pr(>|t|)")]
  log_write("  KDM_core6 × bmi_cat2 交互:")
  print(round(bmi_int_kc, 4))
}

timer_end("模块65 交互分析")
log_write("===== 模块 65 结束 =====\n")