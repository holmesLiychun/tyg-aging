# ============================================================================
# 模块 40：描述性统计与 Table 1（修复版 + KDM_core6 整合）
# ============================================================================

library(tableone)
library(survey)
library(flextable)
library(officer)

timer_start("模块40 描述性统计与Table1")

nhanes_adult <- nhanes %>% filter(age >= 20, exam_stat == 2)
ds_adult_mec <- get_design(nhanes_adult, w = "weight_mec")

# ---- 40.1 连续变量加权分布（★ 加入 kdm_core6_ba ★）----
log_write("--- 连续变量加权分布（MEC 权重, 全部成人）---")
cont_vars <- c("age", "bmi", "tyg", "homa_ir", "glucose", "trig", "totchol",
               "sbp", "dbp", "albumin", "crp", "wbc", "lymph", "mcv",
               "kdm_ba", "kdm_core6_ba", "pheno_ba", "pheno9_ba")   # ← 加了 kdm_core6_ba
for (v in cont_vars) {
  if (v %in% names(nhanes_adult)) {
    ds_v <- get_design(nhanes_adult %>% filter(!is.na(.data[[v]])), w = "weight_mec")
    if (!is.null(ds_v$error)) next
    m <- tryCatch(survey::svymean(as.formula(paste0("~", v)), ds_v)[1],
                  error = function(e) NA_real_)
    log_write(sprintf("  %-20s 加权 Mean=%.3f", v, m))
  }
}

# ---- 40.2 TyG 三分位分布 ----
tyg_dist <- nhanes %>% filter(!is.na(tyg)) %>%
  mutate(tyg_tert = factor(ntile(tyg, 3), levels = 1:3,
                           labels = c("T1 (low)", "T2 (mid)", "T3 (high)"))) %>%
  group_by(tyg_tert) %>% summarise(
    N = n(),
    TyG_Mean = round(mean(tyg, na.rm = TRUE), 3),
    TyG_Median = round(median(tyg, na.rm = TRUE), 3),
    TyG_SD = round(sd(tyg, na.rm = TRUE), 3),
    TyG_Min = round(min(tyg, na.rm = TRUE), 3),
    TyG_Max = round(max(tyg, na.rm = TRUE), 3),
    Glucose_Mean = round(mean(glucose, na.rm = TRUE), 2),
    Trig_Mean = round(mean(trig, na.rm = TRUE), 2),
    HOMA_IR_Mean = round(mean(homa_ir, na.rm = TRUE), 3),
    .groups = "drop")
log_write("--- TyG 三分位分布 ---")
print(tyg_dist)
write.csv(tyg_dist, file.path(OUT_DIR, "tyg_tertile_distribution.csv"), row.names = FALSE)

# ---- 40.3 BA 分布（按性别, KDM-11 & Modified PhenoAge & KDM_core6）----
ba_dist <- nhanes %>%
  filter(!is.na(kdm_ba)) %>%
  group_by(sex) %>% summarise(
    N = n(),
    KDM_BA_Mean = round(mean(kdm_ba, na.rm = TRUE), 2),
    KDM_BA_SD = round(sd(kdm_ba, na.rm = TRUE), 2),
    KDM_Accel_Mean = round(mean(kdm_accel, na.rm = TRUE), 4),
    KDM_Accel_SD = round(sd(kdm_accel, na.rm = TRUE), 4),
    KDMcore6_BA_Mean = round(mean(kdm_core6_ba, na.rm = TRUE), 2),
    KDMcore6_BA_SD = round(sd(kdm_core6_ba, na.rm = TRUE), 2),
    KDMcore6_Accel_Mean = round(mean(kdm_core6_accel, na.rm = TRUE), 4),
    KDMcore6_Accel_SD = round(sd(kdm_core6_accel, na.rm = TRUE), 4),
    Pheno_BA_Mean = round(mean(pheno_ba, na.rm = TRUE), 2),
    Pheno_BA_SD = round(sd(pheno_ba, na.rm = TRUE), 2),
    Pheno_Accel_Mean = round(mean(pheno_accel, na.rm = TRUE), 4),
    Pheno_Accel_SD = round(sd(pheno_accel, na.rm = TRUE), 4),
    .groups = "drop")
log_write("--- BA 分布（按性别, KDM-11 & KDM_core6 & Modified PhenoAge）---")
print(ba_dist)
write.csv(ba_dist, file.path(OUT_DIR, "ba_distribution_by_sex.csv"), row.names = FALSE)

# ---- 40.4 PhenoAge-9 / KDM-9 / KDM_core6 单独分布 ----
ba_dist9 <- nhanes %>%
  filter(!is.na(pheno9_ba)) %>%
  group_by(sex) %>% summarise(
    N = n(),
    Pheno9_BA_Mean = round(mean(pheno9_ba, na.rm = TRUE), 2),
    Pheno9_BA_SD = round(sd(pheno9_ba, na.rm = TRUE), 2),
    Pheno9_Accel_Mean = round(mean(pheno9_accel, na.rm = TRUE), 4),
    Pheno9_Accel_SD = round(sd(pheno9_accel, na.rm = TRUE), 4),
    .groups = "drop")
write.csv(ba_dist9, file.path(OUT_DIR, "pheno9_distribution_by_sex.csv"), row.names = FALSE)

ba_dist_kdm9 <- nhanes %>%
  filter(!is.na(kdm9_ba)) %>%
  group_by(sex) %>% summarise(
    N = n(),
    KDM9_BA_Mean = round(mean(kdm9_ba, na.rm = TRUE), 2),
    KDM9_BA_SD = round(sd(kdm9_ba, na.rm = TRUE), 2),
    KDM9_Accel_Mean = round(mean(kdm9_accel, na.rm = TRUE), 4),
    KDM9_Accel_SD = round(sd(kdm9_accel, na.rm = TRUE), 4),
    .groups = "drop")
write.csv(ba_dist_kdm9, file.path(OUT_DIR, "kdm9_distribution_by_sex.csv"), row.names = FALSE)

# ---- 40.5 CRP 分布（按周期）----
crp_dist <- nhanes %>%
  filter(!is.na(crp)) %>%
  group_by(cycle) %>% summarise(
    N = n(),
    CRP_Mean = round(mean(crp, na.rm = TRUE), 3),
    CRP_Median = round(median(crp, na.rm = TRUE), 3),
    CRP_SD = round(sd(crp, na.rm = TRUE), 3),
    CRP_P25 = round(quantile(crp, 0.25, na.rm = TRUE), 3),
    CRP_P75 = round(quantile(crp, 0.75, na.rm = TRUE), 3),
    .groups = "drop")
log_write("--- CRP 分布（按周期）---")
print(crp_dist)
write.csv(crp_dist, file.path(OUT_DIR, "crp_distribution_by_cycle.csv"), row.names = FALSE)

# ============================================================================
# 40.6 Table 1 基线特征（按 TyG 三分位，加权，含 SMD）
# ★ 主分析切换为 KDM_core6，KDM-11 降为敏感性 ★
# ============================================================================
log_write(">>> 生成 Table 1（主分析: KDM_core6）")

# 主分析：用 ana_kdmcore6
ana_table <- ana_kdmcore6 %>% filter(!is.na(tyg_tert))
ds_table  <- get_design(ana_table, w = "weight_fasting")

vars    <- c("age", "sex", "race", "bmi", "smoke", "alcohol", "paq_active_unified", "edu_cat")
catVars <- c("sex", "race", "smoke", "alcohol", "paq_active_unified", "edu_cat")

t1 <- svyCreateTableOne(
  vars = vars, strata = "tyg_tert",
  data = ds_table, factorVars = catVars,
  addOverall = FALSE
)
print(t1, smd = TRUE, test = TRUE)

# 转置为论文格式
t1_mat <- print(t1, smd = TRUE, test = TRUE, showAllLevels = TRUE)
t1_raw <- cbind(Variable = rownames(t1_mat), as.data.frame(t1_mat))
colnames(t1_raw) <- c("Variable", "Level", "T1", "T2", "T3", "P", "test", "SMD")

result <- data.frame(Characteristic = character(), T1 = character(), T2 = character(),
                     T3 = character(), P = character(), SMD = character(),
                     stringsAsFactors = FALSE)

for (i in 1:nrow(t1_raw)) {
  row <- t1_raw[i, ]
  if (row$Variable != "") {
    result <- rbind(result, data.frame(Characteristic = row$Variable,
                                       T1 = "", T2 = "", T3 = "", P = "", SMD = "",
                                       stringsAsFactors = FALSE))
    result <- rbind(result, data.frame(Characteristic = paste0("  ", row$Level),
                                       T1 = as.character(row$T1), T2 = as.character(row$T2),
                                       T3 = as.character(row$T3), P = as.character(row$P),
                                       SMD = as.character(row$SMD),
                                       stringsAsFactors = FALSE))
  } else {
    result <- rbind(result, data.frame(Characteristic = paste0("  ", row$Level),
                                       T1 = as.character(row$T1), T2 = as.character(row$T2),
                                       T3 = as.character(row$T3), P = as.character(row$P),
                                       SMD = as.character(row$SMD),
                                       stringsAsFactors = FALSE))
  }
}

result$Characteristic <- gsub("\\(mean\\(SD\\)\\)", "(mean ± SD)", result$Characteristic, fixed = TRUE)
result$Characteristic <- ifelse(result$Characteristic == "n", "N (weighted)", result$Characteristic)
result$Characteristic <- ifelse(result$Characteristic == "age (mean ± SD)", "Age, years (mean ± SD)", result$Characteristic)
result$Characteristic <- ifelse(result$Characteristic == "bmi (mean ± SD)", "BMI (mean ± SD)", result$Characteristic)
result$Characteristic <- ifelse(grepl("^sex", result$Characteristic), "Sex", result$Characteristic)
result$Characteristic <- ifelse(grepl("^race", result$Characteristic), "Race/ethnicity", result$Characteristic)
result$Characteristic <- ifelse(grepl("^smoke", result$Characteristic), "Smoking status", result$Characteristic)
result$Characteristic <- ifelse(grepl("^alcohol", result$Characteristic), "Alcohol use", result$Characteristic)
result$Characteristic <- ifelse(grepl("^paq_active", result$Characteristic), "Physical activity", result$Characteristic)
result$Characteristic <- ifelse(grepl("^edu_cat", result$Characteristic), "Education", result$Characteristic)

t1_flex <- flextable(result) %>%
  set_header_labels(Characteristic = "Characteristic", T1 = "T1 (low)",
                    T2 = "T2 (mid)", T3 = "T3 (high)",
                    P = "P value", SMD = "SMD") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(align = "left", j = 1, part = "all") %>%
  align(align = "center", j = 2:6, part = "all") %>%
  bold(part = "header") %>%
  border(border.top = fp_border(color = "black", width = 2), part = "header") %>%
  border(border.bottom = fp_border(color = "black", width = 1.5), part = "header") %>%
  border(border.bottom = fp_border(color = "black", width = 2),
         part = "body", i = nrow(result)) %>%
  set_table_properties(layout = "autofit", width = 1)

# 确保 tables/ 目录存在
if (!dir.exists(file.path(OUT_DIR, "tables"))) {
  dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE)
}

read_docx() %>%
  body_add_par("Table 1. Baseline Characteristics by TyG Tertile (KDM_core6, 7 cycles)",
               style = "heading 1") %>%
  body_add_flextable(t1_flex) %>%
  print(target = file.path(OUT_DIR, "tables", "table1_tyg_baseline.docx"))

log_write("✓ Table 1 已导出: tables/table1_tyg_baseline.docx (KDM_core6)")

# ---- 40.7 补充 Table 1（KDM-11 敏感性，可选）----
log_write(">>> 生成 Table 1（敏感性分析: KDM-11）")

ana_table_kdm11 <- ana_kdm %>% filter(!is.na(tyg_tert))
ds_table_kdm11  <- get_design(ana_table_kdm11, w = "weight_fasting")

t1_kdm11 <- svyCreateTableOne(
  vars = vars, strata = "tyg_tert",
  data = ds_table_kdm11, factorVars = catVars,
  addOverall = FALSE
)

t1_mat_kdm11 <- print(t1_kdm11, smd = TRUE, test = TRUE, showAllLevels = TRUE)
t1_raw_kdm11 <- cbind(Variable = rownames(t1_mat_kdm11), as.data.frame(t1_mat_kdm11))
colnames(t1_raw_kdm11) <- c("Variable", "Level", "T1", "T2", "T3", "P", "test", "SMD")

# （同样的转置逻辑，复用上面的 result 构建代码）
result_kdm11 <- result  # 结构相同，重新填充
result_kdm11 <- data.frame(Characteristic = character(), T1 = character(), T2 = character(),
                           T3 = character(), P = character(), SMD = character(),
                           stringsAsFactors = FALSE)
for (i in 1:nrow(t1_raw_kdm11)) {
  row <- t1_raw_kdm11[i, ]
  if (row$Variable != "") {
    result_kdm11 <- rbind(result_kdm11, data.frame(Characteristic = row$Variable,
                                                   T1 = "", T2 = "", T3 = "", P = "", SMD = "",
                                                   stringsAsFactors = FALSE))
    result_kdm11 <- rbind(result_kdm11, data.frame(Characteristic = paste0("  ", row$Level),
                                                   T1 = as.character(row$T1), T2 = as.character(row$T2),
                                                   T3 = as.character(row$T3), P = as.character(row$P),
                                                   SMD = as.character(row$SMD),
                                                   stringsAsFactors = FALSE))
  } else {
    result_kdm11 <- rbind(result_kdm11, data.frame(Characteristic = paste0("  ", row$Level),
                                                   T1 = as.character(row$T1), T2 = as.character(row$T2),
                                                   T3 = as.character(row$T3), P = as.character(row$P),
                                                   SMD = as.character(row$SMD),
                                                   stringsAsFactors = FALSE))
  }
}
# 标签替换（同上）
result_kdm11$Characteristic <- gsub("\\(mean\\(SD\\)\\)", "(mean ± SD)", result_kdm11$Characteristic, fixed = TRUE)
result_kdm11$Characteristic <- ifelse(result_kdm11$Characteristic == "n", "N (weighted)", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(result_kdm11$Characteristic == "age (mean ± SD)", "Age, years (mean ± SD)", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(result_kdm11$Characteristic == "bmi (mean ± SD)", "BMI (mean ± SD)", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^sex", result_kdm11$Characteristic), "Sex", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^race", result_kdm11$Characteristic), "Race/ethnicity", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^smoke", result_kdm11$Characteristic), "Smoking status", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^alcohol", result_kdm11$Characteristic), "Alcohol use", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^paq_active", result_kdm11$Characteristic), "Physical activity", result_kdm11$Characteristic)
result_kdm11$Characteristic <- ifelse(grepl("^edu_cat", result_kdm11$Characteristic), "Education", result_kdm11$Characteristic)

t1_flex_kdm11 <- flextable(result_kdm11) %>%
  set_header_labels(Characteristic = "Characteristic", T1 = "T1 (low)",
                    T2 = "T2 (mid)", T3 = "T3 (high)",
                    P = "P value", SMD = "SMD") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  align(align = "left", j = 1, part = "all") %>%
  align(align = "center", j = 2:6, part = "all") %>%
  bold(part = "header") %>%
  border(border.top = fp_border(color = "black", width = 2), part = "header") %>%
  border(border.bottom = fp_border(color = "black", width = 1.5), part = "header") %>%
  border(border.bottom = fp_border(color = "black", width = 2),
         part = "body", i = nrow(result_kdm11)) %>%
  set_table_properties(layout = "autofit", width = 1)

read_docx() %>%
  body_add_par("Table S1. Baseline Characteristics by TyG Tertile (KDM-11, Sensitivity Analysis)",
               style = "heading 1") %>%
  body_add_flextable(t1_flex_kdm11) %>%
  print(target = file.path(OUT_DIR, "tables", "tableS1_tyg_baseline_kdm11.docx"))

log_write("✓ Table S1 (KDM-11 sensitivity) 已导出: tables/tableS1_tyg_baseline_kdm11.docx")

cp_save(t1_raw, "40", "table1_raw_core6")
cp_save(t1_raw_kdm11, "40", "table1_raw_kdm11")
timer_end("模块40 描述性统计与Table1")
log_write("===== 模块 40 结束 =====")