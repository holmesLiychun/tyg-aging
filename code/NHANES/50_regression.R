# ============================================================================
# 模块 50：主回归分析（TyG 三分位 → 生物年龄加速残差）
# ----------------------------------------------------------------------------
# 输入 ：ana_kdm、ana_pheno、ana_pheno9、ana_kdm9、ana_pheno_nog、ana_phenocrp
# 输出 ：main_effects —— 6 算法 × 5 模型的 β/SE/P 值长表
#                            → main_tyg_ba_effects.csv
#        Model 5 汇总表 → main_effects_model5.csv
# ============================================================================

timer_start("模块50 主回归分析")

# ---- 50.1 协变量集与建模函数 -----------------------------------------------
covariate_sets <- list(
  "Model 1" = "tyg_tert",
  "Model 2" = "tyg_tert + age + sex + race",
  "Model 3" = "tyg_tert + age + sex + race + bmi",
  "Model 4" = "tyg_tert + age + sex + race + bmi + smoke",
  "Model 5" = paste("tyg_tert + age + sex + race + bmi + smoke + alcohol +",
                    "paq_active_unified + edu_cat")
)

run_sequential_models <- function(endpoint, design_obj) {
  models <- list()
  for (nm in names(covariate_sets)) {
    f <- as.formula(paste(endpoint, "~", covariate_sets[[nm]]))
    models[[nm]] <- survey::svyglm(f, design = design_obj)
  }
  models
}

extract_effects <- function(models, endpoint_name) {
  rows <- lapply(seq_along(models), function(i) {
    sm <- summary(models[[i]])$coefficients
    kp <- grep("tyg_tert", rownames(sm))
    if (length(kp) == 0) return(NULL)
    sub <- sm[kp, c("Estimate", "Std. Error", "Pr(>|t|)"), drop = FALSE]
    data.frame(Endpoint = endpoint_name, Model = paste0("Model ", i),
               Contrast = rownames(sub),
               Beta = round(sub[, "Estimate"], 4),
               SE = round(sub[, "Std. Error"], 4),
               P_value = round(sub[, "Pr(>|t|)"], 4),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ---- 50.2 抽样设计对象（供模块 60/70 复用）---------------------------------
ds_kdm       <- get_design(ana_kdm,    w = "weight_fasting")
ds_pheno     <- get_design(ana_pheno,  w = "weight_fasting")
ds_pheno9    <- get_design(ana_pheno9, w = "weight_fasting")
ds_kdm9      <- get_design(ana_kdm9,   w = "weight_fasting")
ds_pheno_nog <- if (exists("ana_pheno_nog")) get_design(ana_pheno_nog, w = "weight_fasting") else NULL
ds_phenocrp  <- if (exists("ana_phenocrp"))  get_design(ana_phenocrp,  w = "weight_fasting") else NULL

# 有效性检查
for (d_name in c("ds_kdm", "ds_pheno", "ds_pheno9", "ds_kdm9")) {
  d <- get(d_name)
  if (!is.null(d$error)) {
    log_write(sprintf("  ⚠ %s 构造失败: %s", d_name, d$error))
  } else {
    log_write(sprintf("  ✓ %s: n=%d, degf=%d", d_name, nrow(d), survey::degf(d)))
  }
}
if (!is.null(ds_pheno_nog) && is.null(ds_pheno_nog$error)) {
  log_write(sprintf("  ✓ ds_pheno_nog: n=%d, degf=%d", nrow(ds_pheno_nog), survey::degf(ds_pheno_nog)))
} else {
  log_write("  ⚠ ds_pheno_nog 不可用")
}
if (!is.null(ds_phenocrp) && is.null(ds_phenocrp$error)) {
  log_write(sprintf("  ✓ ds_phenocrp: n=%d, degf=%d", nrow(ds_phenocrp), survey::degf(ds_phenocrp)))
} else {
  log_write("  ⚠ ds_phenocrp 不可用")
}

# ---- 50.3 KDM-11 主分析（7 周期）-------------------------------------------
kdm_models   <- run_sequential_models("kdm_accel", ds_kdm)
kdm_results  <- extract_effects(kdm_models, "KDM-11 Acceleration (7 cycles)")
# ---- 50.3b KDM_core6 主分析（7 周期）----
ds_kdmcore6 <- get_design(ana_kdmcore6, w = "weight_fasting")
kdmcore6_models  <- run_sequential_models("kdm_core6_accel", ds_kdmcore6)
kdmcore6_results <- extract_effects(kdmcore6_models, "KDM_core6 Acceleration (7 cycles)")


# ---- 50.4 Modified PhenoAge（7 周期，主分析）-------------------------------
pheno_models  <- run_sequential_models("pheno_accel", ds_pheno)
pheno_results <- extract_effects(pheno_models, "Modified PhenoAge Acceleration (7 cycles)")

# ---- 50.5 PhenoAge-9 原版（5 周期，敏感性）---------------------------------
pheno9_results <- data.frame()
if (!is.null(ds_pheno9$error)) {
  log_write("  PhenoAge-9 设计对象无效，跳过")
} else if (nrow(ana_pheno9) > 100) {
  pheno9_models  <- run_sequential_models("pheno9_accel", ds_pheno9)
  pheno9_results <- extract_effects(pheno9_models, "PhenoAge-9 Original Acceleration (5 cycles)")
} else {
  log_write(sprintf("  PhenoAge-9 样本量不足（n=%d），跳过", nrow(ana_pheno9)))
}

# ---- 50.6 KDM-9（5 周期，敏感性）-------------------------------------------
kdm9_results <- data.frame()
if (!is.null(ds_kdm9$error)) {
  log_write("  KDM-9 设计对象无效，跳过")
} else if (nrow(ana_kdm9) > 100) {
  kdm9_models  <- run_sequential_models("kdm9_accel", ds_kdm9)
  kdm9_results <- extract_effects(kdm9_models, "KDM-9 Acceleration (with CRP, 5 cycles)")
} else {
  log_write(sprintf("  KDM-9 样本量不足（n=%d），跳过", nrow(ana_kdm9)))
}

# ---- 50.7 Pheno-noG（7 周期，解耦核心敏感性）-------------------------------
phenog_results <- data.frame()
if (is.null(ds_pheno_nog)) {
  log_write("  Pheno-noG: ds_pheno_nog 不存在，跳过")
} else if (!is.null(ds_pheno_nog$error)) {
  log_write(sprintf("  Pheno-noG 设计对象无效: %s，跳过", ds_pheno_nog$error))
} else if (nrow(ana_pheno_nog) > 100) {
  phenog_models  <- run_sequential_models("pheno_nog_accel", ds_pheno_nog)
  phenog_results <- extract_effects(phenog_models, "Pheno-noG Acceleration (7 cycles)")
  log_write(sprintf("  ✓ Pheno-noG 回归完成: %d 行", nrow(phenog_results)))
} else {
  log_write(sprintf("  Pheno-noG 样本量不足（n=%d），跳过", nrow(ana_pheno_nog)))
}

# ---- 50.8 Pheno-noCRP（5 周期，敏感性）-------------------------------------
phenocrp_results <- data.frame()
if (is.null(ds_phenocrp)) {
  log_write("  Pheno-noCRP: ds_phenocrp 不存在，跳过（模块30未生成ana_phenocrp）")
} else if (!is.null(ds_phenocrp$error)) {
  log_write(sprintf("  Pheno-noCRP 设计对象无效: %s，跳过", ds_phenocrp$error))
} else if (nrow(ana_phenocrp) > 100) {
  phenocrp_models  <- run_sequential_models("phenocrp_accel", ds_phenocrp)
  phenocrp_results <- extract_effects(phenocrp_models, "Pheno-noCRP Acceleration (5 cycles)")
  log_write(sprintf("  ✓ Pheno-noCRP 回归完成: %d 行", nrow(phenocrp_results)))
} else {
  log_write(sprintf("  Pheno-noCRP 样本量不足（n=%d），跳过", nrow(ana_phenocrp)))
}

# ---- 50.9 合并所有结果 ------------------------------------------------------
main_effects <- rbind(kdm_results, pheno_results)
if (nrow(pheno9_results) > 0)   main_effects <- rbind(main_effects, pheno9_results)
if (nrow(kdm9_results) > 0)     main_effects <- rbind(main_effects, kdm9_results)
if (nrow(kdmcore6_results) > 0) main_effects <- rbind(main_effects, kdmcore6_results)
if (nrow(phenog_results) > 0)   main_effects <- rbind(main_effects, phenog_results)
if (nrow(phenocrp_results) > 0) main_effects <- rbind(main_effects, phenocrp_results)



log_write("--- 主回归效应估计 ---")
print(main_effects, row.names = FALSE)
write.csv(main_effects, file.path(OUT_DIR, "main_tyg_ba_effects.csv"), row.names = FALSE)

# ---- 50.10 Model 5 汇总表 --------------------------------------------------
if (nrow(main_effects) > 0) {
  m5 <- main_effects %>% filter(grepl("Model 5", Model))
  log_write("--- Model 5 β 汇总 ---")
  print(m5)
  write.csv(m5, file.path(OUT_DIR, "main_effects_model5.csv"), row.names = FALSE)
}

cp_save(main_effects, "50", "main_effects")
cp_save(kdmcore6_results, "50", "kdmcore6_effects")  # 方便模块60/70单独调用
timer_end("模块50 主回归分析")
log_write("===== 模块 50 结束 =====\n")