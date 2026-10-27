# ============================================================================
# 模块 60：敏感性分析（全算法版，KDM-9 为主 + KDM_core6 敏感性）
# ----------------------------------------------------------------------------
# 输入 ：ana_kdm、ana_pheno、ana_pheno9、ana_kdm9、ana_pheno_nog、ana_phenocrp
#        ana_kdmcore6（新增，KDM_core6 敏感性）
# 输出 ：sensitivity_winsor3sd.csv        —— 截尾（±3SD）
#         sensitivity_race_stratified.csv  —— 种族分层 + 交互
#         six_algorithms_tyg_association.csv—— 7 算法 TyG 关联（含 KDM_core6）
#         mediation_tyg_inflammation_bioage.csv—— 中介分析
#         decoupling_matrix_standardized.csv —— 4暴露 × 7结局 解耦矩阵
# ============================================================================

timer_start("模块60 敏感性分析")

# ---- 算法清单（KDM-9 为主，KDM_core6 放最后作为敏感性）-----------------------
algo_list <- list(
  "KDM-11"            = list(data = "ana_kdm",       accel = "kdm_accel"),
  "KDM-9"             = list(data = "ana_kdm9",      accel = "kdm9_accel"),
  "PhenoAge-9"        = list(data = "ana_pheno9",    accel = "pheno9_accel"),
  "Modified PhenoAge" = list(data = "ana_pheno",     accel = "pheno_accel"),
  "Pheno-noG"         = list(data = "ana_pheno_nog", accel = "pheno_nog_accel"),
  "Pheno-noCRP"       = list(data = "ana_phenocrp",  accel = "phenocrp_accel"),
  "KDM_core6"         = list(data = "ana_kdmcore6",  accel = "kdm_core6_accel")   # ★ 敏感性
)

# 验证可用算法
available_algos <- list()
for (nm in names(algo_list)) {
  dname <- algo_list[[nm]]$data
  if (!exists(dname)) {
    log_write(sprintf("  ⚠ %s: %s 不存在，跳过", nm, dname))
    next
  }
  d <- get(dname)
  if (!algo_list[[nm]]$accel %in% names(d)) {
    log_write(sprintf("  ⚠ %s: %s 中无 %s 列，跳过", nm, dname, algo_list[[nm]]$accel))
    next
  }
  available_algos[[nm]] <- algo_list[[nm]]
}
log_write(sprintf("  可用算法: %d 套（KDM-9 主分析，KDM_core6 敏感性）", length(available_algos)))

# ============================================================================
# 60.1 截尾敏感性分析（BA accel ±3SD）
# ============================================================================
log_write(">>> 60.1 截尾敏感性分析 (BA accel ±3SD)")

winsor_effects <- data.frame(Algorithm = character(), Endpoint = character(),
                             Model = character(), Contrast = character(),
                             Beta = numeric(), SE = numeric(), P = numeric(),
                             stringsAsFactors = FALSE)

for (nm in names(available_algos)) {
  dname <- available_algos[[nm]]$data
  d <- get(dname)
  ep <- available_algos[[nm]]$accel
  
  # 算截尾变量
  m <- mean(d[[ep]], na.rm = TRUE)
  s <- sd(d[[ep]], na.rm = TRUE)
  ep_winsor <- paste0(ep, "_winsor")
  d[[ep_winsor]] <- pmin(pmax(d[[ep]], m - 3 * s), m + 3 * s)
  
  ds_w <- get_design(d, w = "weight_fasting")
  if (!is.null(ds_w$error)) {
    log_write(sprintf("  %s: 设计对象无效，跳过", nm))
    next
  }
  
  for (mod_name in names(covariate_sets)) {
    f <- as.formula(paste(ep_winsor, "~", covariate_sets[[mod_name]]))
    fit_w <- tryCatch(survey::svyglm(f, design = ds_w), error = function(e) NULL)
    if (is.null(fit_w)) next
    
    sm_w <- summary(fit_w)$coefficients
    t3_idx <- grep("tyg_tertT3", rownames(sm_w))
    if (length(t3_idx) == 0) next
    
    t3 <- sm_w[t3_idx[1], c("Estimate", "Std. Error", "Pr(>|t|)")]
    winsor_effects <- rbind(winsor_effects, data.frame(
      Algorithm = nm, Endpoint = ep_winsor,
      Model = mod_name,
      Contrast = "T3 vs T1",
      Beta = round(t3[1], 4), SE = round(t3[2], 4), P = round(t3[3], 4),
      stringsAsFactors = FALSE
    ))
  }
  log_write(sprintf("  ✓ %s 截尾完成", nm))
}

log_write("--- 截尾分析结果（Model 5 T3 vs T1）---")
print(winsor_effects %>% filter(Model == "Model 5"))
write.csv(winsor_effects, file.path(OUT_DIR, "sensitivity_winsor3sd.csv"), row.names = FALSE)

# ============================================================================
# 60.2 种族分层（NHW vs NHB）+ 交互
# ============================================================================
log_write(">>> 60.2 种族分层")

# 构造种族分层分析人群
race_ana <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg)) %>%
  mutate(race_binary = factor(case_when(
    race == "Non-Hispanic Black" ~ "NHB",
    race == "Non-Hispanic White" ~ "NHW",
    TRUE ~ NA_character_
  ))) %>%
  filter(!is.na(race_binary))

# 把各 accel 合并进来
for (nm in names(available_algos)) {
  d <- get(available_algos[[nm]]$data)
  accel_col <- available_algos[[nm]]$accel
  if (accel_col %in% names(race_ana)) next
  d_sub <- d %>% dplyr::select(SEQN, all_of(accel_col))
  race_ana <- race_ana %>% left_join(d_sub, by = "SEQN")
}

race_results <- data.frame(Algorithm = character(), Race = character(),
                           Model = character(),
                           Beta_T3 = numeric(), SE_T3 = numeric(), P_T3 = numeric(),
                           stringsAsFactors = FALSE)

for (nm in names(available_algos)) {
  ep <- available_algos[[nm]]$accel
  if (!ep %in% names(race_ana)) {
    log_write(sprintf("  %s: %s 不在 race_ana 中，跳过", nm, ep))
    next
  }
  
  for (r in c("NHW", "NHB")) {
    sub_r <- race_ana %>% filter(race_binary == r, !is.na(.data[[ep]]))
    if (nrow(sub_r) < 100) {
      log_write(sprintf("  %s %s: 样本量不足（n=%d），跳过", nm, r, nrow(sub_r)))
      next
    }
    
    ds_r <- get_design(sub_r, w = "weight_fasting")
    if (!is.null(ds_r$error)) {
      log_write(sprintf("  %s %s: 设计对象无效", nm, r))
      next
    }
    
    fit_r <- tryCatch(
      survey::svyglm(as.formula(paste(ep, "~ tyg_tert + age + sex + bmi")), design = ds_r),
      error = function(e) NULL)
    if (is.null(fit_r)) next
    
    sm_r <- summary(fit_r)$coefficients
    t3_idx <- grep("tyg_tertT3", rownames(sm_r))
    if (length(t3_idx) == 0) next
    
    t3_row <- sm_r[t3_idx[1], c("Estimate", "Std. Error", "Pr(>|t|)")]
    race_results <- rbind(race_results, data.frame(
      Algorithm = nm, Race = r,
      Model = "tyg_tert + age + sex + bmi",
      Beta_T3 = round(t3_row[1], 4), SE_T3 = round(t3_row[2], 4), P_T3 = round(t3_row[3], 4),
      stringsAsFactors = FALSE
    ))
  }
  log_write(sprintf("  ✓ %s 种族分层完成", nm))
}

log_write("--- 种族分层结果 ---")
print(race_results)
write.csv(race_results, file.path(OUT_DIR, "sensitivity_race_stratified.csv"), row.names = FALSE)

# 交互作用（全人群）
ds_race_all <- get_design(race_ana %>% filter(!is.na(pheno_accel)), w = "weight_fasting")
if (!is.null(ds_race_all$error)) {
  log_write("  种族交互：设计对象无效，跳过")
} else {
  for (nm in names(available_algos)) {
    ep <- available_algos[[nm]]$accel
    if (!ep %in% names(race_ana)) next
    sub_int <- race_ana %>% filter(!is.na(.data[[ep]]))
    ds_int <- get_design(sub_int, w = "weight_fasting")
    if (!is.null(ds_int$error)) next
    
    fit_int <- tryCatch(
      survey::svyglm(as.formula(paste(ep, "~ tyg_tert * race_binary + age + sex + bmi")),
                     design = ds_int),
      error = function(e) NULL)
    if (is.null(fit_int)) next
    
    sm_int <- summary(fit_int)$coefficients
    int_row <- grep("tyg_tertT3:race_binaryNHB", rownames(sm_int))
    if (length(int_row) > 0) {
      int_p <- sm_int[int_row[1], "Pr(>|t|)"]
      log_write(sprintf("  交互 p 值 (%s, T3×NHB): %.4f", nm, int_p))
    }
  }
}

# ============================================================================
# 60.3 七算法 TyG 关联快速检查（KDM-9 主 + KDM_core6 敏感性）
# ============================================================================
log_write(">>> 60.3 七算法 TyG 关联快速检查（含 KDM_core6 敏感性）")

tyg_assoc <- data.frame(Algorithm = character(), Contrast = character(),
                        Beta = numeric(), SE = numeric(), P = numeric(),
                        N = integer(), stringsAsFactors = FALSE)

for (nm in names(available_algos)) {
  dname <- available_algos[[nm]]$data
  d <- get(dname)
  ep <- available_algos[[nm]]$accel
  
  if (!"tyg_tert" %in% names(d)) {
    log_write(sprintf("  %s: 无 tyg_tert，跳过", nm))
    next
  }
  
  ds_tmp <- get_design(d %>% filter(!is.na(.data[[ep]]), !is.na(tyg_tert)),
                       w = "weight_fasting")
  if (!is.null(ds_tmp$error)) {
    log_write(sprintf("  %s: 设计对象无效，跳过", nm))
    next
  }
  
  fit <- tryCatch(
    survey::svyglm(as.formula(paste(ep, "~ tyg_tert + age + sex + race + bmi + smoke + alcohol + paq_active_unified + edu_cat")), design = ds_tmp),
    error = function(e) NULL)
  if (is.null(fit)) {
    log_write(sprintf("  %s: 模型拟合失败", nm))
    next
  }
  
  sm <- summary(fit)$coefficients
  for (tt in c("tyg_tertT2", "tyg_tertT3")) {
    idx <- grep(tt, rownames(sm))
    if (length(idx) > 0) {
      tyg_assoc <- rbind(tyg_assoc, data.frame(
        Algorithm = nm,
        Contrast = ifelse(grepl("T2", tt), "T2 vs T1", "T3 vs T1"),
        Beta = round(sm[idx[1], 1], 4), SE = round(sm[idx[1], 2], 4),
        P = round(sm[idx[1], 4], 4),
        N = nrow(d),
        stringsAsFactors = FALSE
      ))
    }
  }
  log_write(sprintf("  ✓ %s 关联分析完成", nm))
}

log_write("--- 七算法 TyG 关联（KDM-9 主，KDM_core6 敏感性末位） ---")
print(tyg_assoc)
write.csv(tyg_assoc, file.path(OUT_DIR, "six_algorithms_tyg_association.csv"), row.names = FALSE)
# ↑ 文件名保持原名以兼容下游模块70，内容已含7算法

# ============================================================================
# 60.4 中介分析（TyG T3 vs T1 → 炎症 → 生物年龄加速）
# ============================================================================
log_write(">>> 60.4 中介分析")

mediators <- c("lncrp_mgdl", "wbc", "lymph")
outcomes_m <- c()
for (nm in names(available_algos)) {
  outcomes_m <- c(outcomes_m, available_algos[[nm]]$accel)
}
outcomes_m <- unique(outcomes_m)

covariates_m <- c("age", "sex", "race", "bmi", "smoke", "alcohol", "paq_active_unified", "edu_cat")

med_all_results <- list()
counter_med <- 0

# 用 ana_pheno 作为主数据源
med_base <- ana_pheno
for (med in mediators) {
  if (!med %in% names(med_base) && med %in% names(nhanes)) {
    med_base <- med_base %>% left_join(nhanes %>% dplyr::select(SEQN, all_of(med)), by = "SEQN")
  }
}

for (outcome in outcomes_m) {
  if (!outcome %in% names(med_base)) {
    log_write(sprintf("  %s 不在 med_base 中，跳过", outcome))
    next
  }
  
  for (med in mediators) {
    if (!med %in% names(med_base)) {
      log_write(sprintf("  %s 不在 med_base 中，跳过", med))
      next
    }
    
    log_write(sprintf("  >>> TyG_tert -> %s -> %s", med, outcome))
    
    sub_data <- med_base %>%
      filter(!is.na(.data[[outcome]]), !is.na(.data[[med]]), !is.na(tyg_tert)) %>%
      dplyr::select(all_of(c("SEQN", "tyg_tert", outcome, med,
                             "weight_fasting", "psu", "strata", covariates_m))) %>%
      na.omit()
    
    log_write(sprintf("     有效样本量: %d", nrow(sub_data)))
    if (nrow(sub_data) < 100) {
      log_write("     ⚠ 样本量过小，跳过")
      next
    }
    
    design_sub <- survey::svydesign(id = ~psu, strata = ~strata,
                                    weights = ~weight_fasting,
                                    nest = TRUE, data = sub_data)
    
    form_m <- as.formula(paste(med, "~ tyg_tert +", paste(covariates_m, collapse = " + ")))
    form_y <- as.formula(paste(outcome, "~ tyg_tert +", med, "+",
                               paste(covariates_m, collapse = " + ")))
    
    model_m <- tryCatch(survey::svyglm(form_m, design = design_sub, family = gaussian()),
                        error = function(e) NULL)
    model_y <- tryCatch(survey::svyglm(form_y, design = design_sub, family = gaussian()),
                        error = function(e) NULL)
    if (is.null(model_m) || is.null(model_y)) {
      log_write("     ⚠ 中介模型拟合失败，跳过")
      next
    }
    
    med_out <- tryCatch({
      set.seed(123)
      mediation::mediate(model.m = model_m, model.y = model_y,
                         treat = "tyg_tert", mediator = med,
                         control.value = "T1 (low)", treat.value = "T3 (high)",
                         boot = TRUE, sims = 1000)
    }, error = function(e) {
      log_write(sprintf("     ⚠ mediate() 失败: %s", e$message))
      NULL
    })
    if (is.null(med_out)) next
    
    counter_med <- counter_med + 1
    med_all_results[[counter_med]] <- data.frame(
      Outcome = outcome, Mediator = med,
      ACME = round(med_out$d0, 4),
      ACME_CI_low = round(med_out$d0.ci[1], 4),
      ACME_CI_high = round(med_out$d0.ci[2], 4),
      ADE = round(med_out$z0, 4),
      ADE_CI_low = round(med_out$z0.ci[1], 4),
      ADE_CI_high = round(med_out$z0.ci[2], 4),
      Total_Effect = round(med_out$tau.coef, 4),
      Prop_Mediated = round(med_out$n0, 4),
      N = nrow(sub_data), stringsAsFactors = FALSE
    )
    
    log_write(sprintf("     ACME: %.3f (95%% CI: %.3f, %.3f), Prop. Mediated: %.1f%%",
                      med_out$d0, med_out$d0.ci[1], med_out$d0.ci[2],
                      med_out$n0 * 100))
  }
}

if (length(med_all_results) > 0) {
  med_final <- do.call(rbind, med_all_results)
  log_write("--- 中介分析结果 ---")
  print(med_final)
  write.csv(med_final,
            file.path(OUT_DIR, "mediation_tyg_inflammation_bioage.csv"),
            row.names = FALSE)
  log_write(sprintf("✓ 中介分析完成（共 %d 行）", nrow(med_final)))
} else {
  log_write("⚠ 中介分析未产生任何结果")
}

# ============================================================================
# 60.5 解耦矩阵（4 暴露 × 7 结局，标准化 β）
# ============================================================================
log_write(">>> 60.5 解耦矩阵（含 KDM_core6 敏感性）")

# 构造分析基底：用 ana_pheno，left_join 其他 accel
matrix_base <- ana_pheno

for (nm in names(available_algos)) {
  dname <- available_algos[[nm]]$data
  if (dname == "ana_pheno") next
  d <- get(dname)
  accel_col <- available_algos[[nm]]$accel
  if (accel_col %in% names(matrix_base)) next
  d_sub <- d %>% dplyr::select(SEQN, all_of(accel_col))
  matrix_base <- matrix_base %>% left_join(d_sub, by = "SEQN")
}

# 在基底内算 Z-score
matrix_base <- matrix_base %>%
  mutate(
    tyg_z = as.numeric(scale(tyg, center = TRUE, scale = TRUE)),
    tyg_bmi_z = as.numeric(scale(tyg * bmi, center = TRUE, scale = TRUE)),
    tg_hdl_z = ifelse(!is.na(trig) & !is.na(hdl),
                      as.numeric(scale(trig / hdl, center = TRUE, scale = TRUE)), NA_real_),
    fpg_z = as.numeric(scale(glucose, center = TRUE, scale = TRUE))
  ) %>%
  filter(!is.na(tg_hdl_z))

exposures_z <- c("tyg_z", "tyg_bmi_z", "tg_hdl_z", "fpg_z")
expo_label_z <- c("tyg_z" = "TyG (per SD)", "tyg_bmi_z" = "TyG-BMI (per SD)",
                  "tg_hdl_z" = "TG/HDL-C (per SD)", "fpg_z" = "FPG-only (per SD)")

outcome_col_map <- list()
for (nm in names(available_algos)) {
  outcome_col_map[[nm]] <- available_algos[[nm]]$accel
}

models_cov <- list(
  "Model 1" = character(0),                           # 仅暴露，对应主分析 Model 1（仅tyg_tert）
  "Model 2" = c("age", "sex", "race"),                # 对应主分析 Model 2
  "Model 3" = c("age", "sex", "race", "bmi"),         # 对应主分析 Model 3
  "Model 4" = c("age", "sex", "race", "bmi", "smoke"), # 对应主分析 Model 4
  "Model 5" = c("age", "sex", "race", "bmi", "smoke", "alcohol", "paq_active_unified", "edu_cat")  # 全调整
)

matrix_results <- list()
counter_mat <- 0
skipped_mat <- 0

for (expo in exposures_z) {
  for (outcome_name in names(outcome_col_map)) {
    outcome_col <- outcome_col_map[[outcome_name]]
    
    if (!outcome_col %in% names(matrix_base)) {
      skipped_mat <- skipped_mat + 1
      next
    }
    
    for (m_name in names(models_cov)) {
      covs <- models_cov[[m_name]]
      keep_cols <- c("SEQN", expo, outcome_col, "weight_fasting", "psu", "strata", covs)
      
      missing_cols <- setdiff(keep_cols, names(matrix_base))
      if (length(missing_cols) > 0) {
        skipped_mat <- skipped_mat + 1
        next
      }
      
      sub_data <- matrix_base %>%
        dplyr::select(all_of(keep_cols)) %>%
        na.omit()
      
      if (nrow(sub_data) < 100) {
        skipped_mat <- skipped_mat + 1
        next
      }
      
      sub_data <- droplevels(sub_data)
      if (length(unique(sub_data$strata)) < 2 || length(unique(sub_data$psu)) < 2) {
        skipped_mat <- skipped_mat + 1
        next
      }
      
      design_sub <- survey::svydesign(id = ~psu, strata = ~strata,
                                      weights = ~weight_fasting,
                                      nest = TRUE, data = sub_data)
      # 改成：
      if (length(covs) == 0) {
        form <- as.formula(paste(outcome_col, "~", expo))
      } else {
        form <- as.formula(paste(outcome_col, "~", expo, "+", paste(covs, collapse = " + ")))
      }
      
      fit <- tryCatch(survey::svyglm(form, design = design_sub, family = gaussian()),
                      error = function(e) {
                        skipped_mat <<- skipped_mat + 1
                        NULL
                      })
      if (is.null(fit)) next
      
      sm <- summary(fit)$coefficients
      expo_row <- which(rownames(sm) == expo)
      if (length(expo_row) == 0) {
        skipped_mat <- skipped_mat + 1
        next
      }
      
      counter_mat <- counter_mat + 1
      matrix_results[[counter_mat]] <- data.frame(
        Exposure = expo_label_z[expo], ExposureVar = expo,
        Outcome = outcome_name, OutcomeVar = outcome_col,
        Model = m_name,
        Beta = round(sm[expo_row[1], 1], 4),
        SE = round(sm[expo_row[1], 2], 4),
        CI_low = round(sm[expo_row[1], 1] - 1.96 * sm[expo_row[1], 2], 4),
        CI_high = round(sm[expo_row[1], 1] + 1.96 * sm[expo_row[1], 2], 4),
        P_value = round(sm[expo_row[1], 4], 4),
        N_unweighted = nrow(sub_data),
        stringsAsFactors = FALSE
      )
    }
  }
}

log_write(sprintf("✓ 解耦矩阵完成 %d 个模型（4暴露×7结局×5模型），跳过 %d 个", counter_mat, skipped_mat))

if (counter_mat > 0) {
  mat_final <- do.call(rbind, matrix_results)
  log_write("--- Model 5 标准化 β（28 格 = 4暴露×7结局）---")
  print(mat_final %>% filter(Model == "Model 5") %>% 
          dplyr::select(Exposure, Outcome, Beta, CI_low, CI_high, P_value, N_unweighted))
  write.csv(mat_final, file.path(OUT_DIR, "decoupling_matrix_standardized.csv"),
            row.names = FALSE)
  log_write("✓ 已保存: decoupling_matrix_standardized.csv")
}

timer_end("模块60 敏感性分析")
log_write("===== 模块 60 结束 =====\n")