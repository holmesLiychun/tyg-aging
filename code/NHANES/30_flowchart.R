# ============================================================================
# 模块 30：分析人群流程图
# ============================================================================

timer_start("模块30 分析人群流程图")

# ---- 30.1 抽样设计构造函数 ----
get_design <- function(d, w = "weight_fasting", psu = "psu", stra = "strata") {
  if (!is.data.frame(d) || nrow(d) == 0) {
    return(list(error = "数据为空，无法构造设计对象"))
  }
  d_sub <- d %>% dplyr::filter(!is.na(.data[[w]]), .data[[w]] > 0)
  if (nrow(d_sub) == 0) {
    return(list(error = sprintf("有效权重观测为 0（权重列 %s 全 NA 或非正）", w)))
  }
  if (length(unique(d_sub[[stra]])) < 2 || length(unique(d_sub[[psu]])) < 2) {
    return(list(error = sprintf("strata 或 psu 水平数不足（stra=%d, psu=%d），无法估计方差",
                                length(unique(d_sub[[stra]])), length(unique(d_sub[[psu]])))))
  }
  ds <- tryCatch(
    survey::svydesign(ids = stats::reformulate(psu),
                      strata = stats::reformulate(stra),
                      weights = stats::reformulate(w),
                      nest = TRUE, data = d_sub),
    error = function(e) list(error = e$message)
  )
  ds
}

weighted_n <- function(d, w = "weight_fasting") {
  if (!is.data.frame(d) || nrow(d) == 0) return(0)
  d_sub <- d %>% dplyr::filter(!is.na(.data[[w]]), .data[[w]] > 0)
  if (nrow(d_sub) == 0) return(0)
  ds <- get_design(d_sub, w = w)
  if (!is.null(ds$error)) return(0)
  tryCatch(as.numeric(survey::svytotal(rep(1, nrow(d_sub)), ds)), error = function(e) 0)
}

# ---- 30.2 逐阶段样本量统计 ----
n0_total   <- nrow(nhanes)
n1_adult   <- nhanes %>% filter(age >= 20) %>% nrow()
n1_mec     <- nhanes %>% filter(age >= 20, exam_stat == 2) %>% nrow()
n2_fast    <- nhanes %>% filter(age >= 20, exam_stat == 2, fasting == 1) %>% nrow()
n3_tyg     <- nhanes %>% filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg)) %>% nrow()

n4_kdm <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg)) %>%
  filter_at(vars(all_of(KDM11_MARKERS)), all_vars(!is.na(.))) %>% nrow()

n4_kdmcore6 <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg)) %>%
  filter_at(vars(creatinine, bun, uap, wbc, mcv, sbp), all_vars(!is.na(.))) %>% 
  nrow()

n4_pheno9 <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg), cycle %in% CRP_CYCLES) %>%
  filter_at(vars(all_of(PHENO9_MARKERS)), all_vars(!is.na(.))) %>% nrow()

n4_phenocrp <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg), cycle %in% CRP_CYCLES) %>%
  filter_at(vars(albumin, creatinine, glucose, lymph, mcv, rdw, alp, wbc), 
            all_vars(!is.na(.))) %>% nrow()

n5_kdm    <- n4_kdm
n5_pheno9 <- n4_pheno9

# ============================================================================
# ★ 30.3 先构造所有分析人群数据框（必须在 30.4/30.5/30.7 之前！）
# ============================================================================
add_tyg_tert <- function(df) {
  df %>% mutate(tyg_tert = factor(ntile(tyg, 3), levels = 1:3,
                                  labels = c("T1 (low)", "T2 (mid)", "T3 (high)")))
}

ana_kdm <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(kdm_accel)) %>%
  add_tyg_tert()

ana_kdmcore6 <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(kdm_core6_accel)) %>%
  add_tyg_tert()

ana_pheno <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(pheno_accel)) %>%
  add_tyg_tert()

ana_pheno9 <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(pheno9_accel), cycle %in% CRP_CYCLES) %>%
  add_tyg_tert()

ana_kdm9 <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(kdm9_accel), cycle %in% CRP_CYCLES) %>%
  add_tyg_tert()

ana_phenocrp <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(phenocrp_accel), cycle %in% CRP_CYCLES) %>%
  add_tyg_tert()

ana_pheno_nog <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1,
         !is.na(weight_fasting), weight_fasting > 0,
         !is.na(tyg), !is.na(pheno_nog_accel)) %>%
  add_tyg_tert()

# 保存检查点（★ 补上 ana_kdmcore6 ★）
cp_save(ana_kdm, "30", "ana_kdm")
cp_save(ana_kdmcore6, "30", "ana_kdmcore6")
cp_save(ana_pheno, "30", "ana_pheno")
cp_save(ana_pheno9, "30", "ana_pheno9")
cp_save(ana_kdm9, "30", "ana_kdm9")
cp_save(ana_phenocrp, "30", "ana_phenocrp")
cp_save(ana_pheno_nog, "30", "ana_pheno_nog")

# ---- 30.4 日志输出（重排编号，加入 KDM_core6）----
log_write("\n========== 分析人群筛选流程（原始未加权人数）==========")
log_write(sprintf("Step 0 : 原始合并样本（全部年龄）:                    n = %d", n0_total))
log_write(sprintf("Step 1 : 年龄 >= 20 岁:                               n = %d (排除 %d)", n1_adult, n0_total - n1_adult))
log_write(sprintf("Step 1': 年龄 >= 20 岁 & MEC 检查合格:               n = %d (排除 %d)", n1_mec, n1_adult - n1_mec))
log_write(sprintf("Step 2 : 空腹合格 (fasting == 1):                     n = %d (排除 %d)", n2_fast, n1_mec - n2_fast))
log_write(sprintf("Step 3 : TyG 非缺失:                                  n = %d (排除 %d)", n3_tyg, n2_fast - n3_tyg))
log_write(sprintf("Step 4a: KDM-11 marker 完整（不含CRP, 7周期）:        n = %d (排除 %d)", n4_kdm, n3_tyg - n4_kdm))
log_write(sprintf("Step 4b: KDM_core6 marker 完整（6项, 7周期）:         n = %d (排除 %d)", n4_kdmcore6, n3_tyg - n4_kdmcore6))
log_write(sprintf("Step 4c: PhenoAge-9 marker 完整（含CRP, 5周期）:      n = %d (排除 %d)", n4_pheno9, n3_tyg - n4_pheno9))
log_write(sprintf("Step 4d: Pheno-noCRP marker 完整（去CRP, 5周期）:     n = %d (排除 %d)", n4_phenocrp, n3_tyg - n4_phenocrp))
log_write(sprintf("Step 5a: KDM-11 最终分析人群（敏感性, 7周期）:       n = %d", nrow(ana_kdm)))
log_write(sprintf("Step 5b: KDM_core6 最终分析人群（跨队列主算法）:     n = %d", nrow(ana_kdmcore6)))
log_write(sprintf("Step 5c: PhenoAge-9 最终分析人群（敏感性, 5周期）:    n = %d", nrow(ana_pheno9)))
log_write(sprintf("Step 5d: Pheno-noCRP 最终分析人群（敏感性, 5周期）:   n = %d", nrow(ana_phenocrp)))
log_write(sprintf("Step 5e: Pheno-noG 最终分析人群（敏感性, 7周期）:    n = %d", nrow(ana_pheno_nog)))

# ---- 30.5 加权人口（★ 加入 KDM_core6 ★）----
log_write(sprintf("\n========== 加权人口 =========="))
log_write(sprintf("MEC 加权人口（全部成人）: %.0f",
                  weighted_n(nhanes %>% filter(age >= 20, exam_stat == 2), w = "weight_mec")))
log_write(sprintf("空腹加权人口: %.0f",
                  weighted_n(nhanes %>% filter(age >= 20, exam_stat == 2, fasting == 1),
                             w = "weight_fasting")))
log_write(sprintf("KDM-11 分析人群空腹加权: %.0f", weighted_n(ana_kdm)))
log_write(sprintf("KDM_core6 分析人群空腹加权（跨队列主算法）: %.0f", weighted_n(ana_kdmcore6)))
log_write(sprintf("Modified PhenoAge 分析人群空腹加权: %.0f", weighted_n(ana_pheno)))
log_write(sprintf("PhenoAge-9 分析人群空腹加权（5周期）: %.0f", weighted_n(ana_pheno9)))
log_write(sprintf("KDM-9 分析人群空腹加权（5周期）: %.0f", weighted_n(ana_kdm9)))
log_write(sprintf("Pheno-noCRP 分析人群空腹加权（5周期）: %.0f", weighted_n(ana_phenocrp)))
log_write(sprintf("Pheno-noG 分析人群空腹加权（7周期）: %.0f", weighted_n(ana_pheno_nog)))

# ---- 30.6 各周期样本量（★ 加入 KDMcore6_complete 列 ★）----
log_write(sprintf("\n========== 各周期样本量（按筛选阶段）=========="))
n_by_cycle <- nhanes %>%
  filter(age >= 20, exam_stat == 2) %>%
  group_by(cycle) %>% summarise(
    MEC = n(),
    Fasting = sum(fasting == 1, na.rm = TRUE),
    TyG = sum(!is.na(tyg), na.rm = TRUE),
    CRP_avail = sum(!is.na(crp), na.rm = TRUE),
    KDM_complete = sum(!is.na(albumin) & !is.na(alp) & !is.na(lncreat) & !is.na(hba1c) &
                         !is.na(sbp) & !is.na(bun) & !is.na(uap) & !is.na(lymph) &
                         !is.na(mcv) & !is.na(wbc) & !is.na(totchol), na.rm = TRUE),
    KDMcore6_complete = sum(!is.na(creatinine) & !is.na(bun) & !is.na(uap) & 
                              !is.na(wbc) & !is.na(mcv) & !is.na(sbp), na.rm = TRUE),
    Pheno_complete = sum(!is.na(albumin) & !is.na(lncreat) & !is.na(glucose) &
                           !is.na(hba1c) & !is.na(sbp) & !is.na(bmi) & !is.na(totchol) &
                           !is.na(trig) & !is.na(crp), na.rm = TRUE),
    .groups = "drop"
  )
print(n_by_cycle)

# ---- 30.7 保存流程图数据（★ 加入 KDM_core6 两行 ★）----
flow_data <- data.frame(
  Step = c("原始合并", "成人(>=20岁)", "MEC合格", "空腹合格", "TyG非缺失",
           "KDM-11 marker完整", "KDM_core6 marker完整(6项,7周期)", 
           "PhenoAge-9 marker完整(5周期)", "Pheno-noCRP marker完整(5周期)",
           "KDM-11最终人群(敏感性)", "KDM_core6最终人群(跨队列主算法)",
           "PhenoAge-9最终人群(敏感性)", "Pheno-noCRP最终人群(敏感性)",
           "Pheno-noG最终人群(敏感性)"),
  N_raw = c(n0_total, n1_adult, n1_mec, n2_fast, n3_tyg,
            n4_kdm, n4_kdmcore6, n4_pheno9, n4_phenocrp,
            nrow(ana_kdm), nrow(ana_kdmcore6), 
            nrow(ana_pheno9), nrow(ana_phenocrp), nrow(ana_pheno_nog))
)
write.csv(flow_data, file.path(OUT_DIR, "flowchart_data.csv"), row.names = FALSE)
log_write("✓ 流程图数据已保存: flowchart_data.csv")

# ---- 30.8 论文摘要文本生成（主分析切换为 KDM_core6，保留 KDM-11 并列）----
{
  # 1. 以空腹合格成人为基数
  base_n <- n2_fast
  
  # 2. KDM_core6 排除人数
  excl_tyg_core6      <- n2_fast - n3_tyg
  excl_biomarker_core6<- n3_tyg - n4_kdmcore6
  excl_covar_core6    <- n4_kdmcore6 - nrow(ana_kdmcore6)
  
  # 3. KDM-11 排除人数（原有）
  excl_biomarker_kdm  <- n3_tyg - n4_kdm
  excl_covar_kdm      <- n4_kdm - nrow(ana_kdm)
  
  # 4. 最终样本量与加权
  final_core6_n      <- nrow(ana_kdmcore6)
  final_core6_weight <- weighted_n(ana_kdmcore6)
  final_kdm_n        <- nrow(ana_kdm)
  final_kdm_weight   <- weighted_n(ana_kdm)
  final_pheno9_n     <- nrow(ana_pheno9)
  
  weighted_core6_fmt <- format(round(final_core6_weight, 0), big.mark = ",")
  weighted_kdm_fmt   <- format(round(final_kdm_weight, 0), big.mark = ",")
  
  # 5a. KDM_core6 主分析摘要
  p0_text_core6 <- sprintf(
    "从原始N=%d名受访者中，排除TyG组分缺失n=%d人、生物年龄标志物不完整n=%d人、协变量缺失n=%d人，最终主要分析(KDM_core6, 跨队列统一算法)样本为%d人，加权代表%s名美国成人; PhenoAge-9分析样本为%d人。所有分析均考虑复杂多阶段概率抽样设计，采用调查权重、初级抽样单位和分层变量。",
    base_n, excl_tyg_core6, excl_biomarker_core6, excl_covar_core6,
    final_core6_n, weighted_core6_fmt, final_pheno9_n
  )
  
  # 5b. KDM-11 敏感性分析摘要（原有逻辑保留）
  p0_text_kdm <- sprintf(
    "从原始N=%d名受访者中，排除TyG组分缺失n=%d人、生物年龄标志物不完整n=%d人、协变量缺失n=%d人，最终敏感性分析(KDM-11)样本为%d人，加权代表%s名美国成人; PhenoAge-9分析样本为%d人。所有分析均考虑复杂多阶段概率抽样设计，采用调查权重、初级抽样单位和分层变量。",
    base_n, excl_tyg_core6, excl_biomarker_kdm, excl_covar_kdm,
    final_kdm_n, weighted_kdm_fmt, final_pheno9_n
  )
  
  # 6. 输出到日志
  log_write("\n========== 补-P0 摘要文本（可直接用于论文/补充材料）==========")
  log_write("--- [主分析] KDM_core6 ---")
  log_write(p0_text_core6)
  log_write("--- [敏感性分析] KDM-11 ---")
  log_write(p0_text_kdm)
  
  # 7. 保存为独立 txt 文件
  supp_p0_path <- file.path(OUT_DIR, "supp_P0_exclusion_summary.txt")
  cat("=== 主分析: KDM_core6 ===\n", p0_text_core6, "\n\n",
      "=== 敏感性分析: KDM-11 ===\n", p0_text_kdm, "\n",
      file = supp_p0_path)
  log_write(sprintf("✓ 补-P0 摘要文本已保存至: %s", supp_p0_path))
}

# ---- 清理重复代码，只保留一次 timer_end ----
timer_end("模块30 分析人群流程图")
log_write("===== 模块 30 结束 =====")