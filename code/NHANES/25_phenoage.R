# ============================================================================
# 模块 25：四套 PhenoAge（自训练线性 BA + 性别分层残差加速）
# ============================================================================

timer_start("模块25 四套PhenoAge")

library(dplyr)
library(tidyr)

# ---- 准备基础数据 ----
pheno_base <- nhanes %>%
  filter(age >= 20, exam_stat == 2) %>%
  mutate(
    albumin_gL   = albumin * 10,
    lncreat_umol = log(creatinine * 88.42 + 1),
    lncrp_mgdl   = log(crp / 10 + 1),
    sex_male     = ifelse(sex == "Male", 1, 0)
  )

# ============================================================================
# 辅助函数
# ============================================================================
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
# 1) Modified PhenoAge（去 CRP，7 周期）
# ============================================================================
log_write(">>> 1) Modified PhenoAge（去 CRP, 7 周期）")

pheno_mod_base <- pheno_base %>%
  filter(!is.na(albumin), !is.na(alp), !is.na(creatinine),
         !is.na(hba1c), !is.na(sbp), !is.na(bun), !is.na(uap),
         !is.na(lymph), !is.na(mcv), !is.na(wbc),
         !is.na(totchol), !is.na(trig))

pheno_mod_res <- train_ba_accel(
  pheno_mod_base,
  formula_str = "age ~ albumin_gL + alp + lncreat_umol + hba1c + sbp + bun + uap + lymph + mcv + wbc + totchol + trig",
  ba_name = "pheno_ba",
  accel_name = "pheno_accel"
)

# ============================================================================
# 2) PhenoAge-9 原版 Levine（含 CRP，5 周期）
# ============================================================================
log_write(">>> 2) PhenoAge-9 原版（含 CRP, 5 周期）")

pheno9_base <- pheno_base %>%
  filter(cycle %in% CRP_CYCLES,
         !is.na(crp), !is.na(albumin), !is.na(creatinine), !is.na(glucose),
         !is.na(lymph), !is.na(mcv), !is.na(rdw), !is.na(alp), !is.na(wbc))

pheno9_res <- train_ba_accel(
  pheno9_base,
  formula_str = "age ~ albumin_gL + alp + lncreat_umol + I(glucose/18.0) + lncrp_mgdl + lymph + mcv + rdw + wbc",
  ba_name = "pheno9_ba",
  accel_name = "pheno9_accel"
)

# ============================================================================
# 3) Pheno-noG（去 glucose+trig+totchol，7 周期）
# ============================================================================
log_write(">>> 3) Pheno-noG（去 glucose+trig+totchol, 7 周期）")

nog_base <- pheno_base %>%
  filter(!is.na(albumin), !is.na(alp), !is.na(creatinine), !is.na(hba1c),
         !is.na(sbp), !is.na(bun), !is.na(uap), !is.na(lymph), !is.na(mcv),
         !is.na(wbc), !is.na(rdw))

nog_res <- train_ba_accel(
  nog_base,
  formula_str = "age ~ albumin_gL + alp + lncreat_umol + hba1c + sbp + bun + uap + lymph + mcv + wbc + rdw + lncrp_mgdl",
  ba_name = "pheno_nog_ba",
  accel_name = "pheno_nog_accel"
)

# ============================================================================
# 4) Pheno-noCRP（原版 9 标志物去 CRP，5 周期）
# ============================================================================
log_write(">>> 4) Pheno-noCRP（原版去 CRP, 5 周期）")

phenocrp_base <- pheno_base %>%
  filter(cycle %in% CRP_CYCLES,
         !is.na(albumin), !is.na(creatinine), !is.na(glucose),
         !is.na(lymph), !is.na(mcv), !is.na(rdw), !is.na(alp), !is.na(wbc))
# 注意：故意不包含 crp

phenocrp_res <- train_ba_accel(
  phenocrp_base,
  formula_str = "age ~ albumin_gL + alp + lncreat_umol + I(glucose/18.0) + lymph + mcv + rdw + wbc",
  ba_name = "phenocrp_ba",
  accel_name = "phenocrp_accel"
)

# ============================================================================
# 合并回 nhanes（用 dplyr::select 明确指定，避免 MASS 遮蔽）
# ============================================================================
log_write(">>> 合并四套 PhenoAge 回主数据框...")

# 先删掉已有的旧列（避免重复 join 产生 .x/.y）
# 先删掉已有的旧列（避免重复 join 产生 .x/.y）
nhanes$pheno_ba <- NULL
nhanes$pheno_accel <- NULL
nhanes$pheno9_ba <- NULL
nhanes$pheno9_accel <- NULL
nhanes$pheno_nog_ba <- NULL
nhanes$pheno_nog_accel <- NULL
nhanes$phenocrp_ba <- NULL          # ← 新增
nhanes$phenocrp_accel <- NULL       # ← 新增

if (!is.null(pheno_mod_res)) {
  nhanes <- nhanes %>%
    left_join(dplyr::select(pheno_mod_res, SEQN, pheno_ba, pheno_accel), by = "SEQN")
}
if (!is.null(pheno9_res)) {
  nhanes <- nhanes %>%
    left_join(dplyr::select(pheno9_res, SEQN, pheno9_ba, pheno9_accel), by = "SEQN")
}
if (!is.null(nog_res)) {
  nhanes <- nhanes %>%
    left_join(dplyr::select(nog_res, SEQN, pheno_nog_ba, pheno_nog_accel), by = "SEQN")
}
if (!is.null(phenocrp_res)) {       # ← 新增
  nhanes <- nhanes %>%
    left_join(dplyr::select(phenocrp_res, SEQN, phenocrp_ba, phenocrp_accel), by = "SEQN")
}

# ---- 诊断 ----
log_write(sprintf(
  "PhenoMod 非缺:%d, Pheno9 非缺:%d, Pheno-noG 非缺:%d, Pheno-noCRP 非缺:%d",
  sum(!is.na(nhanes$pheno_ba)),
  sum(!is.na(nhanes$pheno9_ba)),
  sum(!is.na(nhanes$pheno_nog_ba)),
  sum(!is.na(nhanes$phenocrp_ba))    # ← 新增
))

log_write("--- 加速残差诊断 ---")
for (v in c("pheno_accel", "pheno9_accel", "pheno_nog_accel", "phenocrp_accel")) {  # ← 加了 phenocrp_accel
  if (v %in% names(nhanes)) {
    vals <- nhanes[[v]]
    log_write(sprintf(
      "  %s: n=%d, mean=%.4f, sd=%.2f",
      v, sum(!is.na(vals)), mean(vals, na.rm=TRUE), sd(vals, na.rm=TRUE)
    ))
  }
}

cp_save(nhanes, "25", "phenoage_three")
timer_end("模块25 四套PhenoAge")
log_write("===== 模块 25 结束 =====")