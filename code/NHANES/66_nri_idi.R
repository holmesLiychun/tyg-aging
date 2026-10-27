# ============================================================================
# 模块 66：NRI / IDI 分析（TyG 增量判别改善）【完整版 v4】
# 文件：66_nri_idi.R
# ============================================================================
# 
# ⚠⚠⚠ 重要方法学声明 ⚠⚠⚠
# ----------------------------------------------------------------------------
# nricens 包不支持 survey weights（复杂抽样设计权重）。
# 因此本模块的分析为 **未加权（unweighted）敏感性分析**。
# 
# 论文方法部分必须注明：
# "NRI and IDI were calculated from unweighted logistic regression models 
#  as sensitivity analyses, without applying NHANES survey weights, due to 
#  the limitation of the nricens package."
# 
# 结果解释时需说明：未加权分析可能引入选择偏倚，结论需谨慎外推。
# ----------------------------------------------------------------------------
# 
# 方法学设计：
# - 结局：PhenoAge-9 加速 > 中位数（二分类）
# - 基础模型：age + sex + race + bmi + smoke + paq_active_unified + hdl + hba1c
# - 新模型：基础模型 + tyg_z
# - 脂质变量选择 HDL-C（避免与 TyG 中的 TG 信息重叠）
# - 输出：分类NRI(0.5) / 三分类NRI(0.3,0.7) / 连续NRI / IDI
# - 完整病例分析，n=3868
# ============================================================================

timer_start("模块66 NRI/IDI分析")

library(nricens)
library(rms)
library(dplyr)

# ---- 1. 准备数据 ----
main_dat <- ds_kdm$variables

# 非糖尿病亚组
non_dm_idx <- main_dat$dm_status == "Non-diabetes"
nri_dat <- main_dat[non_dm_idx, ]
log_write(sprintf("✓ 非糖尿病亚组数据（含缺失）: n=%d", nrow(nri_dat)))

# ---- 2. 变量确认（已验证全部存在，无需重命名）----
vars_needed <- c("age", "sex", "race", "bmi", "smoke", "paq_active_unified",
                 "hdl", "hba1c", "tyg_z", "pheno9_accel")

missing_vars <- setdiff(vars_needed, names(nri_dat))
if (length(missing_vars) > 0) {
  stop(paste("缺失变量：", paste(missing_vars, collapse=", ")))
} else {
  log_write("✓ 所有必需变量已就位")
}

# ---- 3. 筛选完整病例 ----
nri_dat <- nri_dat %>% 
  filter(complete.cases(!!!syms(vars_needed)))

log_write(sprintf("✓ 完整病例数据: n=%d", nrow(nri_dat)))

# ---- 4. 创建二分类结局 ----
median_pheno <- median(nri_dat$pheno9_accel, na.rm = TRUE)
log_write(sprintf("  PhenoAge-9 加速中位数 = %.4f", median_pheno))

nri_dat$pheno_binary <- ifelse(nri_dat$pheno9_accel > median_pheno, 1, 0)
table(nri_dat$pheno_binary, useNA = "ifany")

# ---- 5. 拟合 glm 模型（nribin 需要 glm + x=TRUE）----
formula_base <- pheno_binary ~ age + sex + race + bmi + smoke + 
  paq_active_unified + hdl + hba1c

formula_new <- pheno_binary ~ age + sex + race + bmi + smoke + 
  paq_active_unified + hdl + hba1c + tyg_z

glm_base <- glm(formula_base, data = nri_dat, family = binomial, x = TRUE)
glm_new  <- glm(formula_new,  data = nri_dat, family = binomial, x = TRUE)

y <- nri_dat$pheno_binary
log_write("✓ GLM 模型拟合完成")

# ---- 6. NRI：二分分类（cut=0.5）----
set.seed(2026)
log_write("计算分类 NRI (cut=0.5)...")

nri_05 <- nricens::nribin(
  event   = y,
  mdl.std = glm_base,
  mdl.new = glm_new,
  updown  = "category",
  cut     = 0.5,
  niter   = 1000,
  alpha   = 0.05,
  msg     = FALSE
)
print(nri_05)

# ---- 7. NRI：三分类（cut=c(0.3, 0.7)）----
log_write("计算三分类 NRI (cut=0.3, 0.7)...")

nri_3cat <- nricens::nribin(
  event   = y,
  mdl.std = glm_base,
  mdl.new = glm_new,
  updown  = "category",
  cut     = c(0.3, 0.7),
  niter   = 1000,
  alpha   = 0.05,
  msg     = FALSE
)
print(nri_3cat)

# ---- 8. NRI：连续（category-free）----
log_write("计算连续 NRI (category-free)...")

nri_cont <- nricens::nribin(
  event   = y,
  mdl.std = glm_base,
  mdl.new = glm_new,
  updown  = "diff",
  cut     = 0,
  niter   = 1000,
  alpha   = 0.05,
  msg     = FALSE
)
print(nri_cont)

# ---- 9. IDI：手算 + bootstrap ----
log_write("计算 IDI (bootstrap 1000次)...")

p_base <- plogis(predict(glm_base, type = "link"))
p_new  <- plogis(predict(glm_new,  type = "link"))

idi_point <- function(pb, pn, yy) {
  mean(pn[yy == 1]) - mean(pb[yy == 1]) -
    (mean(pn[yy == 0]) - mean(pb[yy == 0]))
}

idi0 <- idi_point(p_base, p_new, y)

# Bootstrap
B <- 1000
idi_b <- numeric(B)
n <- length(y)

for (i in seq_len(B)) {
  idx <- sample.int(n, n, replace = TRUE)
  idi_b[i] <- idi_point(p_base[idx], p_new[idx], y[idx])
}

idi_ci <- quantile(idi_b, c(0.025, 0.975))
idi_p  <- t.test(idi_b)$p.value

log_write(sprintf("IDI = %.6f (95%% CI: %.6f - %.6f), P = %.4f",
                  idi0, idi_ci[1], idi_ci[2], idi_p))

# ---- 10. 提取 NRI 结果（兼容不同 nricens 版本）----
# 尝试多种方式提取，确保不报错
extract_nri <- function(obj, label) {
  if (!is.null(obj$nri)) {
    if (is.matrix(obj$nri) || is.data.frame(obj$nri)) {
      return(data.frame(
        Metric   = label,
        Estimate = ifelse("Estimate" %in% colnames(obj$nri), obj$nri["Estimate", "NRI"], obj$nri[1,1]),
        Lower_CI = ifelse("Lower CI" %in% colnames(obj$nri), obj$nri["Lower CI", "NRI"], NA),
        Upper_CI = ifelse("Upper CI" %in% colnames(obj$nri), obj$nri["Upper CI", "NRI"], NA),
        P_value  = ifelse("P-value" %in% colnames(obj$nri), obj$nri["P-value", "NRI"], NA)
      ))
    }
  }
  # 如果结构不同，返回基本值
  return(data.frame(
    Metric   = label,
    Estimate = NA_real_,
    Lower_CI = NA_real_,
    Upper_CI = NA_real_,
    P_value  = NA_real_
  ))
}

# ---- 11. 保存汇总结果 ----
nri_idi_summary <- data.frame(
  Metric   = c("NRI_binary(0.5)", "NRI_3cat(0.3,0.7)", "NRI_continuous", "IDI"),
  Estimate = c(
    ifelse(is.null(nri_05$nri), NA, nri_05$nri[1,1]),
    ifelse(is.null(nri_3cat$nri), NA, nri_3cat$nri[1,1]),
    ifelse(is.null(nri_cont$nri), NA, nri_cont$nri[1,1]),
    idi0
  ),
  Lower_CI = c(
    ifelse(is.null(nri_05$nri), NA, nri_05$nri[1,2]),
    ifelse(is.null(nri_3cat$nri), NA, nri_3cat$nri[1,2]),
    ifelse(is.null(nri_cont$nri), NA, nri_cont$nri[1,2]),
    idi_ci[1]
  ),
  Upper_CI = c(
    ifelse(is.null(nri_05$nri), NA, nri_05$nri[1,3]),
    ifelse(is.null(nri_3cat$nri), NA, nri_3cat$nri[1,3]),
    ifelse(is.null(nri_cont$nri), NA, nri_cont$nri[1,3]),
    idi_ci[2]
  ),
  P_value  = c(
    ifelse(is.null(nri_05$nri), NA, nri_05$nri[1,4]),
    ifelse(is.null(nri_3cat$nri), NA, nri_3cat$nri[1,4]),
    ifelse(is.null(nri_cont$nri), NA, nri_cont$nri[1,4]),
    idi_p
  ),
  stringsAsFactors = FALSE
)

write.csv(nri_idi_summary, file.path(OUT_DIR, "table_S5_NRI_IDI.csv"), row.names = FALSE)
log_write("✓ Table S5 (NRI/IDI 汇总) 已保存")

# ---- 12. 重分类表 ----
pred_base <- plogis(predict(glm_base, type = "link"))
pred_new  <- plogis(predict(glm_new,  type = "link"))

# 二分重分类
reclass_2cat <- table(
  "Base_model" = ifelse(pred_base > 0.5, "High risk", "Low risk"),
  "New_model"  = ifelse(pred_new > 0.5, "High risk", "Low risk")
)
print("=== 二分重分类表 (0.5) ===")
print(reclass_2cat)
write.csv(as.data.frame(reclass_2cat), 
          file.path(OUT_DIR, "table_S5b_reclass_2cat.csv"), row.names = TRUE)

# 三分类重分类
reclass_3cat <- table(
  "Base_model" = cut(pred_base, breaks = c(0, 0.3, 0.7, 1), 
                     labels = c("Low", "Intermediate", "High")),
  "New_model"  = cut(pred_new, breaks = c(0, 0.3, 0.7, 1), 
                     labels = c("Low", "Intermediate", "High"))
)
print("=== 三分类重分类表 (0.3, 0.7) ===")
print(reclass_3cat)
write.csv(as.data.frame(reclass_3cat), 
          file.path(OUT_DIR, "table_S5c_reclass_3cat.csv"), row.names = TRUE)

# ---- 13. 打印结果摘要 ----
log_write("=== 最终结果摘要 ===")
log_write(sprintf("NRI (0.5)        = %.4f", 
                  ifelse(is.null(nri_05$nri), NA, nri_05$nri[1,1])))
log_write(sprintf("NRI (3-category) = %.4f", 
                  ifelse(is.null(nri_3cat$nri), NA, nri_3cat$nri[1,1])))
log_write(sprintf("NRI (continuous) = %.4f", 
                  ifelse(is.null(nri_cont$nri), NA, nri_cont$nri[1,1])))
log_write(sprintf("IDI              = %.6f (P=%.4f)", idi0, idi_p))

timer_end("模块66 NRI/IDI分析")

# ============================================================================
# 模块 66.5：从 bootstrapsample 计算 NRI bootstrap 精确 P 值
# ============================================================================

# ---- 13.5 二分类 NRI bootstrap P 值（正态近似） ----

bs <- nri_05$bootstrapsample

cat("=== bootstrapsample 结构 ===\n")
str(bs)
cat("维度:", dim(bs), "\n")
if (!is.null(colnames(bs))) cat("列名:", colnames(bs), "\n")
if (!is.null(rownames(nri_05$nri))) cat("$nri 行名:", rownames(nri_05$nri), "\n")

bs_mat <- as.matrix(bs)
n_boot <- nrow(bs_mat)
n_stat <- ncol(bs_mat)

cat(sprintf("\nBootstrap 次数 = %d, 统计量数 = %d\n", n_boot, n_stat))

bs_p_values <- data.frame(
  Statistic      = colnames(bs_mat),
  Estimate       = nri_05$nri$Estimate,
  Bootstrap_SE   = NA_real_,
  Bootstrap_P    = NA_real_,
  stringsAsFactors = FALSE
)

for (j in 1:n_stat) {
  bs_vals <- bs_mat[, j]
  bs_vals <- bs_vals[!is.na(bs_vals)]
  
  if (length(bs_vals) >= 10) {
    bs_se <- sd(bs_vals)
    obs   <- nri_05$nri$Estimate[j]
    
    # 正态近似（最稳健，不会溢出）
    z <- abs(obs) / bs_se
    p_norm <- 2 * (1 - pnorm(z))
    
    # ★ 关键修复：写回数据框
    bs_p_values$Bootstrap_SE[j] <- bs_se
    bs_p_values$Bootstrap_P[j]  <- p_norm
    
    cat(sprintf("  %s: Est=%.6f, SE=%.6f, P=%.4f\n",
                colnames(bs_mat)[j], obs, bs_se, p_norm))
  }
}

cat("\n=== NRI Binary Bootstrap P 值汇总 ===\n")
print(bs_p_values, digits = 6)

write.csv(bs_p_values, file.path(OUT_DIR, "table_S5d_NRI_bootstrap_P.csv"), 
          row.names = FALSE)
log_write("✓ NRI binary bootstrap P 值已计算并保存")


# ---- 13.6 三分类 NRI bootstrap P 值（正态近似） ----

bs_3cat <- nri_3cat$bootstrapsample
bs_3cat_mat <- as.matrix(bs_3cat)
cat("\n=== 三分类 NRI bootstrap 结构 ===\n")
cat("维度:", dim(bs_3cat_mat), "\n")

bs_p_3cat <- data.frame(
  Statistic      = colnames(bs_3cat_mat),
  Estimate       = nri_3cat$nri$Estimate,
  Bootstrap_SE   = NA_real_,
  Bootstrap_P    = NA_real_,
  stringsAsFactors = FALSE
)

for (j in 1:ncol(bs_3cat_mat)) {
  bs_vals <- bs_3cat_mat[, j]
  bs_vals <- bs_vals[!is.na(bs_vals)]
  if (length(bs_vals) >= 10) {
    bs_se <- sd(bs_vals)
    obs   <- nri_3cat$nri$Estimate[j]
    z <- abs(obs) / bs_se
    p_norm <- 2 * (1 - pnorm(z))
    
    bs_p_3cat$Bootstrap_SE[j] <- bs_se
    bs_p_3cat$Bootstrap_P[j]  <- p_norm
    
    cat(sprintf("  %s: Est=%.6f, SE=%.6f, P=%.4f\n",
                colnames(bs_3cat_mat)[j], obs, bs_se, p_norm))
  }
}

cat("\n=== NRI 3cat Bootstrap P 值汇总 ===\n")
print(bs_p_3cat, digits = 6)
write.csv(bs_p_3cat, file.path(OUT_DIR, "table_S5e_NRI3cat_bootstrap_P.csv"), 
          row.names = FALSE)


# ---- 13.7 连续 NRI bootstrap P 值（正态近似） ----

bs_cont <- nri_cont$bootstrapsample
bs_cont_mat <- as.matrix(bs_cont)
cat("\n=== 连续 NRI bootstrap 结构 ===\n")
cat("维度:", dim(bs_cont_mat), "\n")

bs_p_cont <- data.frame(
  Statistic      = colnames(bs_cont_mat),
  Estimate       = nri_cont$nri$Estimate,
  Bootstrap_SE   = NA_real_,
  Bootstrap_P    = NA_real_,
  stringsAsFactors = FALSE
)

for (j in 1:ncol(bs_cont_mat)) {
  bs_vals <- bs_cont_mat[, j]
  bs_vals <- bs_vals[!is.na(bs_vals)]
  if (length(bs_vals) >= 10) {
    bs_se <- sd(bs_vals)
    obs   <- nri_cont$nri$Estimate[j]
    z <- abs(obs) / bs_se
    p_norm <- 2 * (1 - pnorm(z))
    
    bs_p_cont$Bootstrap_SE[j] <- bs_se
    bs_p_cont$Bootstrap_P[j]  <- p_norm
    
    cat(sprintf("  %s: Est=%.6f, SE=%.6f, P=%.4f\n",
                colnames(bs_cont_mat)[j], obs, bs_se, p_norm))
  }
}

cat("\n=== NRI Continuous Bootstrap P 值汇总 ===\n")
print(bs_p_cont, digits = 6)
write.csv(bs_p_cont, file.path(OUT_DIR, "table_S5f_NRIcont_bootstrap_P.csv"), 
          row.names = FALSE)

log_write("✓ 所有 NRI bootstrap P 值计算完成")


# ---- 13.8 修复 nri_idi_summary 的 P_value 列 ----

nri_idi_summary_fixed <- nri_idi_summary

# 现在 bs_p_values$Bootstrap_P 有值了，不再返回 NA
nri_idi_summary_fixed$P_value[1] <- bs_p_values$Bootstrap_P[1]   # NRI binary(0.5)
nri_idi_summary_fixed$P_value[2] <- bs_p_3cat$Bootstrap_P[1]       # NRI 3cat
nri_idi_summary_fixed$P_value[3] <- bs_p_cont$Bootstrap_P[1]       # NRI continuous
nri_idi_summary_fixed$P_value[4] <- idi_p                          # IDI

cat("\n=== 修复后的 nri_idi_summary_fixed ===\n")
print(nri_idi_summary_fixed, digits = 7)

write.csv(nri_idi_summary_fixed, 
          file.path(OUT_DIR, "table_S5_NRI_IDI_with_P.csv"), 
          row.names = FALSE)
log_write("✓ table_S5_NRI_IDI_with_P.csv 已保存（含 bootstrap P 值）")

# ============================================================================
# 模块 66.9：高阶内部验证（重写版 v2）
# 文件：66_nri_idi.R  【追加到末尾】
# ============================================================================
#
# 修复记录：
#   v2 (2026-10-14):
#     - [G] calibration 改用 glm 逻辑校准回归（修复 loess coef 长度为0的报错）
#     - [F] 修复 ΔAUC 日志显示错误
#     - [E] IDI P 值改用正态近似，与 NRI 口径统一
# ============================================================================

# ---- 依赖检查 ----
suppressMessages({
  library(pROC)
  library(rms)
  library(Hmisc)
})
log_write("=== 进入模块 66.9：高阶内部验证 ===")

# ============================================================================
# A. 准备：统一数据 + 两模型
# ============================================================================
dat <- nri_dat   # 非糖尿病、完整病例、n=3868
yy  <- dat$pheno_binary

m0 <- lrm(formula_base, data = dat, x = TRUE, y = TRUE)   # 基础
m1 <- lrm(formula_new,  data = dat, x = TRUE, y = TRUE)   # 新

p0 <- plogis(predict(m0, type = "lp"))
p1 <- plogis(predict(m1, type = "lp"))

# 事件率
event_rate <- mean(yy)
log_write(sprintf("事件率 (PhenoAge-9加速>中位数) = %.3f, n_event=%d",
                  event_rate, sum(yy)))

# ============================================================================
# B. NRI 事件/非事件拆分表（在 0.5 阈值上）
# ============================================================================
reclass_events <- function(p_std, p_new, y, cut = 0.5) {
  # 事件组上移 / 非事件组下移
  up_e   <- mean(p_new[y == 1] > cut & p_std[y == 1] <= cut)
  down_e <- mean(p_new[y == 1] <= cut & p_std[y == 1] > cut)
  up_ne  <- mean(p_new[y == 0] > cut & p_std[y == 0] <= cut)
  down_ne<- mean(p_new[y == 0] <= cut & p_std[y == 0] > cut)
  nri_events <- up_e - down_e
  nri_nonev  <- down_ne - up_ne
  nri_total  <- nri_events + nri_nonev
  data.frame(
    subgroup   = c("Events (acceleration present)", "Non-events", "Total"),
    N_up       = c(round(up_e*sum(yy)),   round(up_ne*(length(yy)-sum(yy))), NA),
    N_down     = c(round(down_e*sum(yy)), round(down_ne*(length(yy)-sum(yy))), NA),
    NRI        = c(nri_events, nri_nonev, nri_total),
    stringsAsFactors = FALSE
  )
}

nri_split <- reclass_events(p0, p1, yy, cut = 0.5)
print("=== NRI 事件/非事件拆分表 (0.5) ===")
print(nri_split, digits = 5)
write.csv(nri_split, file.path(OUT_DIR, "table_S5g_NRI_split.csv"), row.names = FALSE)


# ============================================================================
# C. 单/双侧 NRI 检验（H0: NRI=0  vs  H0: NRI<=0）
# ============================================================================
nri_onesided <- function(bs_mat, stat_col, obs, alt = "greater") {
  # alt: "greater" = H0: NRI<=0（关心是否有改善）
  v <- bs_mat[, stat_col]
  v <- v[!is.na(v)]
  se <- sd(v)
  z_two <- obs / se
  p_two <- 2 * (1 - pnorm(abs(z_two)))
  p_one <- 1 - pnorm(z_two)      # P(NRI <= 0)
  data.frame(SE = se, Z = z_two, P_two_sided = p_two, P_one_sided_greater = p_one)
}

# NRI binary 的 "NRI" 统计量列
nri_test_binary <- nri_onesided(nri_05$bootstrapsample, "NRI", nri_05$nri["NRI","Estimate"])
print("=== NRI(0.5) 单/双侧检验 ===")
print(nri_test_binary)

# 连续 NRI 同样
nri_test_cont <- nri_onesided(nri_cont$bootstrapsample, "NRI", nri_cont$nri["NRI","Estimate"])
print("=== NRI(continuous) 单/双侧检验 ===")
print(nri_test_cont)


# ============================================================================
# D. IDI 检验改用正态近似（与 NRI 口径统一）
# ============================================================================
idi_p_norm <- 2 * (1 - pnorm(abs(idi0) / sd(idi_b)))
idi_p_oneside <- 1 - pnorm(idi0 / sd(idi_b))
log_write(sprintf("IDI P值(正态近似,双侧)=%.4f, 单侧(>0)=%.4f", idi_p_norm, idi_p_oneside))


# ============================================================================
# E. 内部验证：optimism-corrected NRI / IDI / AUC（bootstrap）
# ============================================================================
B_int <- 200   # 内部验证 200 次足以；资源允许可提到 500
set.seed(2026)

run_validate <- function(dat, B) {
  nri_opt <- idi_opt <- auc_opt <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    idx <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    db  <- dat[idx, ]
    oob <- setdiff(seq_len(nrow(dat)), unique(idx))
    if (length(oob) < 20) { next }
    g0b <- glm(formula_base, data = db, family = binomial, x = TRUE)
    g1b <- glm(formula_new,  data = db, family = binomial, x = TRUE)
    
    pa0 <- plogis(predict(g0b, db, type = "link"))
    pa1 <- plogis(predict(g1b, db, type = "link"))
    
    # apparent（拟合样本内）
    ap_nri <- mean(pa1[db$pheno_binary==1] > 0.5 & pa0[db$pheno_binary==1] <= 0.5) -
      mean(pa1[db$pheno_binary==1] <= 0.5 & pa0[db$pheno_binary==1] > 0.5) +
      mean(pa1[db$pheno_binary==0] <= 0.5 & pa0[db$pheno_binary==0] > 0.5) -
      mean(pa1[db$pheno_binary==0] > 0.5 & pa0[db$pheno_binary==0] <= 0.5)
    
    ap_idi <- mean(pa1[db$pheno_binary==1]) - mean(pa0[db$pheno_binary==1]) -
      (mean(pa1[db$pheno_binary==0]) - mean(pa0[db$pheno_binary==0]))
    
    ap_auc <- as.numeric(roc(db$pheno_binary, pa1, quiet = TRUE)$auc) -
      as.numeric(roc(db$pheno_binary, pa0, quiet = TRUE)$auc)
    
    # optimism（OOB 样本上：用 bootstrap 模型预测）
    pr0 <- plogis(predict(g0b, dat[oob, ], type = "link"))
    pr1 <- plogis(predict(g1b, dat[oob, ], type = "link"))
    
    op_nri <- mean(pr1[dat$pheno_binary[oob]==1] > 0.5 & pr0[dat$pheno_binary[oob]==1] <= 0.5) -
      mean(pr1[dat$pheno_binary[oob]==1] <= 0.5 & pr0[dat$pheno_binary[oob]==1] > 0.5) +
      mean(pr1[dat$pheno_binary[oob]==0] <= 0.5 & pr0[dat$pheno_binary[oob]==0] > 0.5) -
      mean(pr1[dat$pheno_binary[oob]==0] > 0.5 & pr0[dat$pheno_binary[oob]==0] <= 0.5)
    
    op_idi <- mean(pr1[dat$pheno_binary[oob]==1]) - mean(pr0[dat$pheno_binary[oob]==1]) -
      (mean(pr1[dat$pheno_binary[oob]==0]) - mean(pr0[dat$pheno_binary[oob]==0]))
    
    op_auc <- as.numeric(roc(dat$pheno_binary[oob], pr1, quiet = TRUE)$auc) -
      as.numeric(roc(dat$pheno_binary[oob], pr0, quiet = TRUE)$auc)
    
    nri_opt[i] <- ap_nri - op_nri
    idi_opt[i] <- ap_idi - op_idi
    auc_opt[i] <- ap_auc - op_auc
  }
  list(
    optimism_nri = median(nri_opt, na.rm = TRUE),
    optimism_idi = median(idi_opt, na.rm = TRUE),
    optimism_auc = median(auc_opt, na.rm = TRUE)
  )
}

val_res <- run_validate(dat, B_int)

# optimism-corrected 估计
apparent_auc_new <- as.numeric(roc(yy, p1, quiet = TRUE)$auc)
apparent_auc_old <- as.numeric(roc(yy, p0, quiet = TRUE)$auc)
apparent_auc_diff <- apparent_auc_new - apparent_auc_old

internal_val <- data.frame(
  Metric             = c("NRI(0.5)", "IDI", "ΔAUC"),
  Apparent           = c(nri_05$nri["NRI","Estimate"], idi0, apparent_auc_diff),
  Optimism           = c(val_res$optimism_nri, val_res$optimism_idi, val_res$optimism_auc),
  Optimism_corrected = c(nri_05$nri["NRI","Estimate"] - val_res$optimism_nri,
                         idi0 - val_res$optimism_idi,
                         apparent_auc_diff - val_res$optimism_auc),
  stringsAsFactors = FALSE
)

print("=== 内部验证（optimism-corrected）===")
print(internal_val, digits = 6)
write.csv(internal_val, file.path(OUT_DIR, "table_S5h_internal_validation.csv"), row.names = FALSE)


# ============================================================================
# F. AUC 比较：基础 vs 新模型（DeLong）
# ============================================================================
roc0 <- roc(yy, p0, quiet = TRUE, ci = TRUE)
roc1 <- roc(yy, p1, quiet = TRUE, ci = TRUE)
roc_diff <- roc.test(roc0, roc1, method = "delong")

auc_compare <- data.frame(
  Model   = c("Base (age+sex+race+bmi+smoke+pa+hdl+hba1c)", "Base + TyG (per SD)"),
  AUC     = c(as.numeric(roc0$auc), as.numeric(roc1$auc)),
  AUC_low = c(roc0$ci[1], roc1$ci[1]),
  AUC_high= c(roc0$ci[3], roc1$ci[3]),
  stringsAsFactors = FALSE
)

delta_auc <- apparent_auc_new - apparent_auc_old
log_write(sprintf("AUC base=%.4f (%.4f-%.4f) | AUC new=%.4f (%.4f-%.4f) | Δ=%.5f, P(DeLong)=%.4f",
                  auc_compare$AUC[1], auc_compare$AUC_low[1], auc_compare$AUC_high[1],
                  auc_compare$AUC[2], auc_compare$AUC_low[2], auc_compare$AUC_high[2],
                  delta_auc, roc_diff$p.value))
write.csv(auc_compare, file.path(OUT_DIR, "table_S5i_AUC_compare.csv"), row.names = FALSE)


# ============================================================================
# G. Calibration：slope & intercept（基础 vs 新模型）【重写 v2】
# ============================================================================
# 方法：logistic 校准回归 logit(P_observed) = α + β × logit(P_predicted)
# 理想校准：β = 1, α = 0
#   β < 1 提示过度自信（overconfident）
#   α ≠ 0 提示系统偏差（calibration-in-the-large）
# ============================================================================

# 将预测概率转换为 logit 尺度
logit <- function(p) log(p / (1 - p))

lp0 <- logit(pmin(pmax(p0, 0.001), 0.999))
lp1 <- logit(pmin(pmax(p1, 0.001), 0.999))

# 事件组上的稳健校准回归
cal_glm0 <- glm(yy ~ lp0, family = binomial)
cal_glm1 <- glm(yy ~ lp1, family = binomial)

cal_meta <- data.frame(
  Model                 = c("Base", "Base + TyG"),
  Calibration_slope     = c(coef(cal_glm0)[2], coef(cal_glm1)[2]),
  Calibration_intercept = c(coef(cal_glm0)[1], coef(cal_glm1)[1]),
  stringsAsFactors = FALSE
)

print("=== Calibration slope / intercept ===")
print(cal_meta, digits = 4)
write.csv(cal_meta, file.path(OUT_DIR, "table_S5j_calibration.csv"), row.names = FALSE)

# ---- G2. 分组校准细节（按预测概率五分位）----
calibrate_detail <- function(pred, y, nk = 5, model_name) {
  cuts <- unique(quantile(pred, probs = seq(0, 1, 1/nk)))
  grp  <- cut(pred, breaks = cuts, include.lowest = TRUE)
  obs_rate  <- tapply(y, grp, mean)
  mean_pred <- tapply(pred, grp, mean)
  n_grp     <- tapply(y, grp, length)
  data.frame(
    Model        = model_name,
    Risk_group   = seq_along(obs_rate),
    N            = as.numeric(n_grp),
    Mean_predicted = as.numeric(mean_pred),
    Observed_rate  = as.numeric(obs_rate),
    stringsAsFactors = FALSE
  )
}

cal_detail <- rbind(
  calibrate_detail(p0, yy, model_name = "Base"),
  calibrate_detail(p1, yy, model_name = "Base + TyG")
)
write.csv(cal_detail, file.path(OUT_DIR, "table_S5j_cal_detail.csv"), row.names = FALSE)

log_write(sprintf("校准斜率: Base=%.4f, Base+TyG=%.4f (理想=1.0)",
                  coef(cal_glm0)[2], coef(cal_glm1)[2]))
log_write(sprintf("校准截距: Base=%.4f, Base+TyG=%.4f (理想=0.0)",
                  coef(cal_glm0)[1], coef(cal_glm1)[1]))

# ============================================================================
# H. Brier score（基础 vs 新）
# ============================================================================
brier <- function(p, y) mean((p - y)^2)
bs0 <- brier(p0, yy)
bs1 <- brier(p1, yy)
# 相对改善比例（scaled Brier）
null_prob <- mean(yy)
bs_null <- brier(rep(null_prob, length(yy)), yy)
bs_improve0 <- (bs_null - bs0) / bs_null
bs_improve1 <- (bs_null - bs1) / bs_null

brier_tab <- data.frame(
  Model        = c("Base", "Base + TyG", "Null model"),
  Brier_score  = c(bs0, bs1, bs_null),
  Scaled_Brier = c(bs_improve0, bs_improve1, 0),
  stringsAsFactors = FALSE
)
log_write(sprintf("Brier: base=%.5f (scaled=%.3f) | new=%.5f (scaled=%.3f) | null=%.5f",
                  bs0, bs_improve0, bs1, bs_improve1, bs_null))
write.csv(brier_tab, file.path(OUT_DIR, "table_S5k_brier.csv"), row.names = FALSE)


# ============================================================================
# I. 亚组 NRI（糖尿病 vs 非糖尿病）
# ============================================================================
dat_all <- main_dat %>% filter(complete.cases(!!!syms(vars_needed)))
yy_all  <- ifelse(dat_all$pheno9_accel > median(dat_all$pheno9_accel, na.rm = TRUE), 1, 0)
dat_all$pheno_binary <- yy_all  # ★ 关键：让 glm(formula_base) 能找到因变量

subgroup_nri <- function(sub_idx, label) {
  dsub <- dat_all[sub_idx, ]
  ysub <- yy_all[sub_idx]
  
  # ---- 安全检查 ----
  n_sub <- length(ysub)
  n_ev  <- sum(ysub)
  
  if (n_ev < 20 || n_ev > n_sub - 20 || n_sub < 50) {
    log_write(sprintf("  [%s] 跳过：n=%d, events=%d（不足）", label, n_sub, n_ev))
    return(data.frame(
      Subgroup   = label,
      n          = n_sub,
      n_event    = n_ev,
      NRI        = NA_real_,
      NRI_low    = NA_real_,
      NRI_high   = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # ---- 拟合模型 ----
  g0 <- glm(formula_base, data = dsub, family = binomial, x = TRUE)
  g1 <- glm(formula_new,  data = dsub, family = binomial, x = TRUE)
  
  # ---- 计算 NRI（带错误捕获）----
  res <- tryCatch({
    nricens::nribin(
      event   = ysub,
      mdl.std = g0,
      mdl.new = g1,
      updown  = "category",
      cut     = 0.5,
      niter   = 500,
      msg     = FALSE
    )
  }, error = function(e) {
    log_write(sprintf("  [%s] nribin 报错: %s", label, e$message))
    return(NULL)
  })
  
  if (is.null(res) || is.null(res$nri)) {
    return(data.frame(
      Subgroup = label, n = n_sub, n_event = n_ev,
      NRI = NA_real_, NRI_low = NA_real_, NRI_high = NA_real_
    ))
  }
  
  # ---- 安全提取结果（兼容不同列名/行名）----
  nri_mat <- as.matrix(res$nri)
  
  # 找行名（兼容 "NRI" 或带空格的变体）
  nri_row <- grep("^NRI$", rownames(nri_mat), ignore.case = TRUE, value = FALSE)
  if (length(nri_row) == 0) {
    nri_row <- 1  # 取第一行作为 NRI 总计
    log_write(sprintf("  [%s] 警告：未找到 'NRI' 行名，使用第1行", label))
  } else {
    nri_row <- nri_row[1]
  }
  
  # 找列名（兼容 Estimate / Lower / Upper 或 Lower CI / Upper CI）
  est_col <- grep("Estimate", colnames(nri_mat), ignore.case = TRUE, value = FALSE)[1]
  low_col <- grep("^Lower",  colnames(nri_mat), ignore.case = TRUE, value = FALSE)[1]
  up_col  <- grep("^Upper",  colnames(nri_mat), ignore.case = TRUE, value = FALSE)[1]
  
  nri_val <- if (!is.na(est_col)) nri_mat[nri_row, est_col] else NA_real_
  nri_lo  <- if (!is.na(low_col)) nri_mat[nri_row, low_col] else NA_real_
  nri_hi  <- if (!is.na(up_col))  nri_mat[nri_row, up_col]  else NA_real_
  
  log_write(sprintf("  [%s] NRI=%.4f (%.4f-%.4f), n=%d, events=%d",
                    label, nri_val, nri_lo, nri_hi, n_sub, n_ev))
  
  data.frame(
    Subgroup = label,
    n        = n_sub,
    n_event  = n_ev,
    NRI      = as.numeric(nri_val),
    NRI_low  = as.numeric(nri_lo),
    NRI_high = as.numeric(nri_hi),
    stringsAsFactors = FALSE
  )
}

# ---- 执行 ----
nri_dm  <- subgroup_nri(dat_all$dm_status == "Diabetes", "Diabetes")
nri_ndm <- subgroup_nri(dat_all$dm_status == "Non-diabetes", "Non-diabetes")
nri_subgroup <- rbind(nri_dm, nri_ndm)

print("=== 亚组 NRI(0.5) ===")
print(nri_subgroup, digits = 4)
write.csv(nri_subgroup, file.path(OUT_DIR, "table_S5l_NRI_subgroup.csv"), row.names = FALSE)

# ============================================================================
# J. 综合汇总表（修复版 v5 —— 函数先定义，逐行构建）
# ============================================================================

# ---- 1. 辅助函数（必须先定义！）----
safe_get <- function(obj, row_pat = "NRI", col_pat = "Estimate", default = NA_real_) {
  if (is.null(obj)) return(default)
  mat <- tryCatch(as.matrix(obj$nri), error = function(e) NULL)
  if (is.null(mat) || length(mat) == 0) return(default)
  r <- grep(row_pat, rownames(mat), ignore.case = TRUE, value = FALSE)
  c <- grep(col_pat, colnames(mat), ignore.case = TRUE, value = FALSE)
  if (length(r) == 0 || length(c) == 0) return(default)
  val <- mat[r[1], c[1]]
  if (length(val) == 0 || is.na(val)) return(default)
  return(as.numeric(val))
}

safe_get_nri_p <- function() {
  if (!exists("bs_p_values")) return(NA_real_)
  if (!is.data.frame(bs_p_values) || nrow(bs_p_values) == 0) return(NA_real_)
  r <- grep("^NRI$", rownames(nri_05$nri), ignore.case = TRUE)
  if (length(r) > 0 && r[1] <= nrow(bs_p_values)) {
    return(as.numeric(bs_p_values$Bootstrap_P[r[1]]))
  }
  return(as.numeric(bs_p_values$Bootstrap_P[1]))
}

# ---- 2. 逐行构建（每行5个标量，绝不会长度不匹配）----
final_summary <- do.call(rbind, list(
  
  data.frame(Metric="AUC: Base model",
             Estimate = as.numeric(apparent_auc_old),
             CI_low   = as.numeric(roc0$ci[1]),
             CI_high  = as.numeric(roc0$ci[3]),
             P_value  = NA_real_),
  
  data.frame(Metric="AUC: Base + TyG (per SD)",
             Estimate = as.numeric(apparent_auc_new),
             CI_low   = as.numeric(roc1$ci[1]),
             CI_high  = as.numeric(roc1$ci[3]),
             P_value  = NA_real_),
  
  data.frame(Metric="ΔAUC (DeLong P)",
             Estimate = as.numeric(apparent_auc_diff),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = ifelse(exists("roc_diff"), as.numeric(roc_diff$p.value), NA_real_)),
  
  data.frame(Metric="NRI, categorical (0.5)",
             Estimate = safe_get(nri_05, "NRI", "Estimate"),
             CI_low   = safe_get(nri_05, "NRI", "Lower"),
             CI_high  = safe_get(nri_05, "NRI", "Upper"),
             P_value  = safe_get_nri_p()),
  
  data.frame(Metric="NRI, events (0.5)",
             Estimate = as.numeric(nri_split$NRI[1]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = ifelse(exists("nri_test_binary") && "P_one" %in% names(nri_test_binary),
                               as.numeric(nri_test_binary$P_one[1]), NA_real_)),
  
  data.frame(Metric="NRI, non-events (0.5)",
             Estimate = as.numeric(nri_split$NRI[2]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = ifelse(exists("nri_test_binary") && "P_one" %in% names(nri_test_binary),
                               as.numeric(nri_test_binary$P_one[1]), NA_real_)),
  
  data.frame(Metric="NRI, three-category (0.3, 0.7)",
             Estimate = safe_get(nri_3cat, "NRI", "Estimate"),
             CI_low   = safe_get(nri_3cat, "NRI", "Lower"),
             CI_high  = safe_get(nri_3cat, "NRI", "Upper"),
             P_value  = ifelse(exists("bs_p_3cat") && nrow(bs_p_3cat) > 0,
                               as.numeric(bs_p_3cat$Bootstrap_P[1]), NA_real_)),
  
  data.frame(Metric="NRI, continuous (category-free)",
             Estimate = safe_get(nri_cont, "NRI", "Estimate"),
             CI_low   = safe_get(nri_cont, "NRI", "Lower"),
             CI_high  = safe_get(nri_cont, "NRI", "Upper"),
             P_value  = ifelse(exists("bs_p_cont") && nrow(bs_p_cont) > 0,
                               as.numeric(bs_p_cont$Bootstrap_P[1]), NA_real_)),
  
  data.frame(Metric="IDI",
             Estimate = as.numeric(idi0),
             CI_low   = as.numeric(idi_ci[1]),
             CI_high  = as.numeric(idi_ci[2]),
             P_value  = as.numeric(idi_p_norm)),
  
  data.frame(Metric="Brier score (base)",
             Estimate = as.numeric(bs0),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_),
  
  data.frame(Metric="Brier score (base + TyG)",
             Estimate = as.numeric(bs1),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_),
  
  data.frame(Metric="Calibration slope (base)",
             Estimate = as.numeric(coef(cal_glm0)[2]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_),
  
  data.frame(Metric="Calibration slope (base + TyG)",
             Estimate = as.numeric(coef(cal_glm1)[2]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_),
  
  data.frame(Metric="Calibration intercept (base)",
             Estimate = as.numeric(coef(cal_glm0)[1]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_),
  
  data.frame(Metric="Calibration intercept (base + TyG)",
             Estimate = as.numeric(coef(cal_glm1)[1]),
             CI_low   = NA_real_,
             CI_high  = NA_real_,
             P_value  = NA_real_)
))

# ---- 3. 验证 ----
cat(sprintf("汇总表: %d 行 × %d 列\n", nrow(final_summary), ncol(final_summary)))

print(final_summary, digits = 6)
write.csv(final_summary, file.path(OUT_DIR, "table_S5m_summary_final.csv"), row.names = FALSE)
log_write("✓ 综合汇总表 table_S5m_summary_final.csv 已保存")

# ============================================================================
# 模块 66.10：C-index 内部验证（修复版）
# ============================================================================

library(rms)
library(survival)

# ---- 设置 datadist（只传入建模相关变量，避免无关列报错）----
vars_for_datadist <- c("pheno_binary", "age", "sex", "race", "bmi", "smoke",
                       "paq_active_unified", "hdl", "hba1c", "tyg_z")
dd <- datadist(dat[, vars_for_datadist, drop = FALSE])
options(datadist = "dd")

# ============================================================================
# A. 二分类结局的 C-index 内部验证（修复版）
# ============================================================================
log_write("=== A. 二分类结局 C-index 内部验证 ===")

# 确认 m0, m1 存在
if (!exists("m0") || !exists("m1")) {
  log_write("重新拟合 lrm 模型...")
  m0 <- lrm(formula_base, data = dat, x = TRUE, y = TRUE)
  m1 <- lrm(formula_new,  data = dat, x = TRUE, y = TRUE)
}

# ---- 关键修复：去掉 dxy=TRUE ----
set.seed(2026)
val_m0 <- validate(m0, method = "boot", B = 200)  # ← 去掉 dxy=TRUE
val_m1 <- validate(m1, method = "boot", B = 200)  # ← 去掉 dxy=TRUE

# ---- 调试：查看 validate 输出结构 ----
log_write("validate(m0) 输出结构：")
log_write(capture.output(print(val_m0)))

# ---- 安全提取函数 ----
extract_cindex <- function(val_obj) {
  # 确认 Dxy 行存在
  if (!"Dxy" %in% rownames(val_obj)) {
    log_write("警告：validate 输出中无 'Dxy' 行")
    return(data.frame(
      C_index_Apparent  = NA_real_,
      C_index_Corrected = NA_real_,
      C_index_Optimism  = NA_real_,
      Dxy_Apparent      = NA_real_,
      Dxy_Corrected     = NA_real_,
      Dxy_Optimism      = NA_real_
    ))
  }
  
  # 提取各列（兼容不同列名）
  dxy_orig <- val_obj["Dxy", "index.orig"]
  
  # index.corrected 可能列名不完全匹配
  corr_col <- grep("correct", colnames(val_obj), ignore.case = TRUE, value = FALSE)
  opt_col  <- grep("optimism", colnames(val_obj), ignore.case = TRUE, value = FALSE)
  
  dxy_corr <- if (length(corr_col) > 0) val_obj["Dxy", corr_col[1]] else NA_real_
  dxy_opt  <- if (length(opt_col) > 0)  val_obj["Dxy", opt_col[1]]  else NA_real_
  
  data.frame(
    C_index_Apparent  = 0.5 + as.numeric(dxy_orig) / 2,
    C_index_Corrected = 0.5 + as.numeric(dxy_corr) / 2,
    C_index_Optimism  = as.numeric(dxy_opt) / 2,
    Dxy_Apparent      = as.numeric(dxy_orig),
    Dxy_Corrected     = as.numeric(dxy_corr),
    Dxy_Optimism      = as.numeric(dxy_opt)
  )
}

cidx_m0 <- extract_cindex(val_m0)
cidx_m1 <- extract_cindex(val_m1)

cindex_compare <- data.frame(
  Model             = c("Base", "Base + TyG (per SD)"),
  C_index_Apparent  = c(cidx_m0$C_index_Apparent,  cidx_m1$C_index_Apparent),
  C_index_Corrected = c(cidx_m0$C_index_Corrected, cidx_m1$C_index_Corrected),
  C_index_Optimism  = c(cidx_m0$C_index_Optimism,  cidx_m1$C_index_Optimism),
  Dxy_Apparent      = c(cidx_m0$Dxy_Apparent,      cidx_m1$Dxy_Apparent),
  Dxy_Corrected     = c(cidx_m0$Dxy_Corrected,     cidx_m1$Dxy_Corrected),
  stringsAsFactors  = FALSE
)

print("=== C-index 内部验证（二分类结局）===")
print(cindex_compare, digits = 5)
write.csv(cindex_compare, file.path(OUT_DIR, "table_S5n_cindex_validation_binary.csv"), row.names = FALSE)

log_write(sprintf("Base C-index: apparent=%.4f, corrected=%.4f", 
                  cidx_m0$C_index_Apparent, cidx_m0$C_index_Corrected))
log_write(sprintf("New C-index:  apparent=%.4f, corrected=%.4f", 
                  cidx_m1$C_index_Apparent, cidx_m1$C_index_Corrected))
# ============================================================================
# B. 生存结局 C-index 内部验证（框架，按需启用）
# ============================================================================
log_write("=== B. 生存结局验证（框架，需替换为实际数据）===")
# 取消注释并替换数据对象名即可运行：
# dat_surv <- dat_charls_complete   # 你的 CHARLS 完整病例
# dd_surv  <- datadist(dat_surv[, c("time", "status", vars_covariates)])
# options(datadist = "dd_surv")
# fit_cox0 <- cph(Surv(time, status) ~ age + sex + ..., data = dat_surv, x=TRUE, y=TRUE, surv=TRUE)
# val_cox0 <- validate(fit_cox0, method="boot", B=200, dxy=TRUE)
# ...

# ============================================================================
# C. 连续结局替代方案（框架）
# ============================================================================
log_write("=== C. 连续结局 R-squared 验证（框架）===")
# ols1 <- ols(pheno9_accel ~ age + sex + race + bmi + smoke + paq_active_unified + hdl + hba1c + tyg_z, data = dat)
# val_ols <- validate(ols1, method="boot", B=200)
# R2_corr <- val_ols["R2", "index.corrected"]

log_write("=== 模块 66.10 完成 ===")


# ============================================================================
# 模块 67：相关系数矩阵（Table S12 各衰老算法内变量与 TyG 的相关性）
# ----------------------------------------------------------------------------
# 输出：table_S12_tyg_correlation_matrix.csv
#       table_S12_tyg_correlation_wide.csv
# ============================================================================
timer_start("模块67 相关系数矩阵(Table S12)")

# ---- 67.1 准备数据：空腹成人，TyG 非缺失 ----
cor_dat <- nhanes %>%
  filter(age >= 20, exam_stat == 2, fasting == 1, !is.na(tyg)) %>%
  dplyr::select(
    SEQN, tyg, glucose, trig, crp, lncrp, wbc, lymph,
    albumin, alp, creatinine, lncreat, hba1c, sbp, bun, uap, mcv, totchol,
    albumin_gL, lncreat_umol, lncrp_mgdl, rdw, glucose_mmol
  ) %>%
  mutate(
    FPG = glucose,
    TG  = trig,
    CRP = crp,
    lnCRP = lncrp,
    WBC_count = wbc,
    Lymphocyte = lymph
  )

# 目标变量列表（按你论文要求的顺序）
target_vars <- c(
  "FPG", "TG", "CRP", "lnCRP", "WBC_count", "Lymphocyte",
  "albumin", "alp", "lncreat", "hba1c", "sbp", "bun", "uap", "mcv", "totchol",
  "albumin_gL", "lncreat_umol", "lncrp_mgdl", "rdw", "glucose_mmol"
)

# ---- 67.2 计算 Pearson + Spearman 相关及 P 值的函数 ----
cor_test <- function(data, var1, var2) {
  complete_dat <- data %>% dplyr::select(all_of(c(var1, var2))) %>% na.omit()
  n <- nrow(complete_dat)
  if (n < 10) {
    return(data.frame(
      N = n, Pearson_r = NA_real_, Pearson_p = NA_real_,
      Spearman_rho = NA_real_, Spearman_p = NA_real_
    ))
  }
  pear  <- cor.test(complete_dat[[var1]], complete_dat[[var2]], method = "pearson")
  spear <- cor.test(complete_dat[[var1]], complete_dat[[var2]], method = "spearman")
  
  data.frame(
    N = n,
    Pearson_r    = round(pear$estimate, 3),
    Pearson_p    = round(pear$p.value, 4),
    Spearman_rho = round(spear$estimate, 3),
    Spearman_p   = round(spear$p.value, 4)
  )
}

# ---- 67.3 循环计算 TyG 与每个变量的相关性 ----
cor_results <- list()
for (v in target_vars) {
  if (v %in% names(cor_dat)) {
    res <- cor_test(cor_dat, "tyg", v)
    cor_results[[v]] <- cbind(Variable = v, res, stringsAsFactors = FALSE)
  } else {
    log_write(sprintf("  ⚠ 变量 %s 不在 cor_dat 中，跳过", v))
  }
}
cor_final <- do.call(rbind, cor_results)

# ---- 67.4 按算法分组标记（方便你写论文时分类）----
cor_final <- as.data.frame(cor_final) %>%
  mutate(
    Algorithm_Component = case_when(
      Variable %in% c("FPG", "TG", "CRP", "lnCRP", "WBC_count", "Lymphocyte") 
      ~ "核心指标/炎症",
      Variable %in% c("albumin", "alp", "lncreat", "hba1c", "sbp", "bun", "uap", 
                      "lymph", "mcv", "wbc", "totchol") 
      ~ "KDM-11 / Modified PhenoAge 组分",
      Variable %in% c("albumin_gL", "lncreat_umol", "lncrp_mgdl", "rdw", "glucose_mmol") 
      ~ "PhenoAge-9 / KDM-9 组分",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::select(Algorithm_Component, Variable, N, Pearson_r, Pearson_p, Spearman_rho, Spearman_p)

# ---- 67.5 保存为 CSV（供你直接贴进论文/补充材料）----
write.csv(cor_final,
          file.path(OUT_DIR, "table_S12_tyg_correlation_matrix.csv"),
          row.names = FALSE)

# 同时输出一个宽表（类似你给的占位表格式，方便直接看）
cor_wide <- cor_final %>%
  dplyr::select(Variable, Pearson_r, Spearman_rho, Pearson_p, N) %>%
  rename(
    `TyG (Pearson r)` = Pearson_r,
    `TyG (Spearman ρ)` = Spearman_rho,
    `P-value` = Pearson_p,
    `N` = N
  )
write.csv(cor_wide,
          file.path(OUT_DIR, "table_S12_tyg_correlation_wide.csv"),
          row.names = FALSE)

# ---- 67.6 日志与打印 ----
log_write("--- Table S12 TyG 与各衰老算法组分相关系数 ---")
print(cor_final, row.names = FALSE)
log_write(sprintf("✓ Table S12 已保存: %s",
                  file.path(OUT_DIR, "table_S12_tyg_correlation_matrix.csv")))
log_write(sprintf("✓ Table S12 宽表已保存: %s",
                  file.path(OUT_DIR, "table_S12_tyg_correlation_wide.csv")))

timer_end("模块67 相关系数矩阵(Table S12)")
log_write("===== 模块 67 结束 =====\n")