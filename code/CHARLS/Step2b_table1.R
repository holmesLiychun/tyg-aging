# ============================================================
# Step 2b: Sample Selection Flow + Table 1 Part A
# 基线特征表（分析开始前）
# ============================================================
CH$reset()
if (!exists("dt_raw")) {
  # 如果 dt_raw 已被清掉，从 checkpoint 或重新读
  if (!exists("dt_unit")) load_ck("dt_unit", "2_audit")
  if (!exists("dt_base")) load_ck("dt_base", "2_audit")
}

cat("--- Step 2b: Sample Selection Flow + Table 1 ---\n")

# ============================================================
# PART 0: 样本筛选流程图
# ============================================================
cat("\n========== SAMPLE SELECTION FLOW ==========\n")

flow_df <- data.frame(
  Step = character(),
  n_total = integer(),
  n_excluded = integer(),
  Reason = character(),
  n_remaining = integer(),
  stringsAsFactors = FALSE
)

# ---- Step 0: 原始数据 ----
n0 <- if(exists("dt_raw")) nrow(dt_raw) else nrow(dt_unit)
flow_df <- rbind(flow_df, data.frame(
  Step = "0: Raw CHARLS data (all waves, all respondents)",
  n_total = n0, n_excluded = 0, Reason = "-",
  n_remaining = n0, stringsAsFactors = FALSE
))
cat(sprintf("Step 0: Raw data = %d observations\n", n0))

# ---- Step 1: 年龄筛选 (age ≥ 45) ----
# dt_unit 已经过滤了 age >= 45
if (exists("dt_unit")) {
  n1 <- nrow(dt_unit)
  n_excl1 <- n0 - n1
  flow_df <- rbind(flow_df, data.frame(
    Step = "1: Exclude age < 45",
    n_total = n0, n_excluded = n_excl1,
    Reason = "Age < 45 at interview",
    n_remaining = n1, stringsAsFactors = FALSE
  ))
  cat(sprintf("Step 1: Exclude age < 45 → excluded %d, remaining %d\n", n_excl1, n1))
}

# ---- Step 2: 2011 基线波次 ----
if (exists("dt_base")) {
  n2 <- nrow(dt_base)
  n_excl2 <- n1 - n2
  flow_df <- rbind(flow_df, data.frame(
    Step = "2: Restrict to 2011 baseline wave",
    n_total = n1, n_excluded = n_excl2,
    Reason = "Not in wave 2011",
    n_remaining = n2, stringsAsFactors = FALSE
  ))
  cat(sprintf("Step 2: Restrict to 2011 wave → excluded %d, remaining %d\n", n_excl2, n2))
}

# ---- Step 3: 空腹亚组 ----
if (exists("dt_base")) {
  n3 <- sum(dt_base$bl_fasting == 1, na.rm = TRUE)
  n_excl3 <- n2 - n3
  flow_df <- rbind(flow_df, data.frame(
    Step = "3: Restrict to fasting subsample",
    n_total = n2, n_excluded = n_excl3,
    Reason = "Non-fasting (bl_fasting != 1)",
    n_remaining = n3, stringsAsFactors = FALSE
  ))
  cat(sprintf("Step 3: Fasting only → excluded %d, remaining %d\n", n_excl3, n3))
}

# ---- Step 4: TyG 非缺失 ----
if (exists("dt_base")) {
  n4 <- sum(dt_base$bl_fasting == 1 & !is.na(dt_base$TyG_baseline), na.rm = TRUE)
  n_excl4 <- n3 - n4
  flow_df <- rbind(flow_df, data.frame(
    Step = "4: Exclude missing TyG",
    n_total = n3, n_excluded = n_excl4,
    Reason = "Missing TG or glucose (TyG undefined)",
    n_remaining = n4, stringsAsFactors = FALSE
  ))
  cat(sprintf("Step 4: Exclude missing TyG → excluded %d, remaining %d\n", n_excl4, n4))
}

# ---- Step 5: 核心协变量非缺失 ----
core_covars <- c("age", "gender", "bmi", "smokev", "drinkl", "edu", "exercise")
if (exists("dt_base")) {
  dt_base_tmp <- dt_base %>% filter(bl_fasting == 1 & !is.na(TyG_baseline))
  complete_idx <- complete.cases(dt_base_tmp[, core_covars[core_covars %in% names(dt_base_tmp)]])
  n5 <- sum(complete_idx)
  n_excl5 <- n4 - n5
  flow_df <- rbind(flow_df, data.frame(
    Step = "5: Exclude missing core covariates",
    n_total = n4, n_excluded = n_excl5,
    Reason = paste("Missing any of:", paste(core_covars, collapse=", ")),
    n_remaining = n5, stringsAsFactors = FALSE
  ))
  cat(sprintf("Step 5: Exclude missing covariates → excluded %d, remaining %d\n", n_excl5, n5))
}

# ---- 保存筛选流程 ----
flow_path <- file.path(PATH_OUT, "tables", "Sample_Selection_Flow.csv")
write.csv(flow_df, flow_path, row.names = FALSE)
CH$log("OK", sprintf("Sample selection flow saved: %s", flow_path), "2b_table1")

cat("\n--- SAMPLE SELECTION SUMMARY ---\n")
print(flow_df, row.names = FALSE)

# ============================================================
# PART 1: Table 1 Part A — 基线特征（按 TyG 四分位分层）
# ============================================================
cat("\n========== TABLE 1 PART A ==========\n")

# ---- 准备 Table 1 数据框 ----
tbl_df <- dt_base %>%
  filter(bl_fasting == 1) %>%        # 空腹亚组
  filter(!is.na(TyG_baseline)) %>%   # TyG 非缺失
  filter(complete.cases(.[, c("age", "gender", "bmi", "smokev", "drinkl", "edu", "exercise")])) %>%
  mutate(
    # TyG 四分位
    tyg_q = cut(TyG_baseline,
                quantile(TyG_baseline, c(0, .25, .5, .75, 1), na.rm = TRUE),
                labels = c("Q1", "Q2", "Q3", "Q4"),
                include.lowest = TRUE),
    
    # ★ 糖尿病定义（已确认编码：diabe 0=无, 1=有）
    diabetes = ifelse(
      (!is.na(diabe) & diabe == 1) |
        (!is.na(bl_glu) & bl_glu >= 126) |
        (!is.na(bl_hbalc) & bl_hbalc >= 6.5),
      1, 0),
    diabetes = factor(diabetes, levels = c(0, 1), labels = c("No", "Yes")),
    
    # 分类变量因子化
    gender   = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    smokev   = factor(smokev, levels = c(0, 1), labels = c("Never", "Ever")),
    drinkl   = factor(drinkl, levels = c(0, 1), labels = c("No", "Yes")),
    edu      = factor(edu),
    exercise = factor(exercise, levels = c(0, 1), labels = c("No", "Yes"))
  )

CH$log("OK", sprintf("Table1 Part A: %d rows (final analytic sample)", nrow(tbl_df)), "2b_table1")

# ---- 加载 tableone ----
if (!require(tableone, quietly = TRUE)) {
  install.packages("tableone", repos = "https://cloud.r-project.org")
  library(tableone)
}

# ---- 变量列表（覆盖你清单里的所有 🔴 必须项）----
vars_cont <- c(
  "age",              # 年龄
  "bmi",              # BMI
  "TyG_baseline",     # TyG 指数
  "bl_glu",           # 空腹血糖 (mg/dL) 🔴
  "bl_tg",            # 甘油三酯 (mg/dL) 🔴
  "bl_hbalc",         # HbA1c (%) 🟠
  "bl_wbc",           # WBC (×10⁹/L) 🟠
  "bl_cho",           # 总胆固醇 (mg/dL)
  "bl_crea",          # 肌酐
  "bl_ua",            # 尿酸
  "systo"             # 收缩压
)

vars_cat <- c(
  "gender",           # 性别
  "diabetes",         # 糖尿病患病率 🔴
  "smokev",           # 吸烟
  "drinkl",           # 饮酒
  "edu",              # 教育
  "exercise"          # 运动
)

# ---- 生成 Table 1 ----
tab1 <- CreateTableOne(
  vars       = c(vars_cont, vars_cat),
  strata     = "tyg_q",
  data       = tbl_df,
  factorVars = vars_cat,
  addOverall = TRUE
)

# ---- 格式化输出 ----
tab1_print <- print(tab1,
                    showAllLevels = TRUE,
                    printToggle  = FALSE,
                    noSpaces     = TRUE)

# ---- 保存 ----
tab1_path <- file.path(PATH_OUT, "tables", "Table1_PartA_baseline_by_TyGq.csv")
write.csv(tab1_print, tab1_path, row.names = TRUE)
CH$log("OK", sprintf("Table 1 Part A saved: %s", tab1_path), "2b_table1")

# ---- 控制台预览 ----
cat("\n========== TABLE 1 PART A PREVIEW ==========\n")
print(tab1, showAllLevels = TRUE)

# ============================================================
# PART 2: TyG 四分位切点（用于正文/表注）
# ============================================================
tyg_cuts <- quantile(tbl_df$TyG_baseline, c(.25, .5, .75), na.rm = TRUE)
cat(sprintf("\n--- TyG Quartile Cut-points (n=%d) ---\n", nrow(tbl_df)))
cat(sprintf("Q1: < %.3f\n", tyg_cuts[1]))
cat(sprintf("Q2: %.3f – %.3f\n", tyg_cuts[1], tyg_cuts[2]))
cat(sprintf("Q3: %.3f – %.3f\n", tyg_cuts[2], tyg_cuts[3]))
cat(sprintf("Q4: ≥ %.3f\n", tyg_cuts[3]))

cuts_path <- file.path(PATH_OUT, "tables", "TyG_quartile_cutpoints.txt")
writeLines(
  c(sprintf("TyG Quartile Cut-points (n=%d, fasting, 2011):", nrow(tbl_df)),
    sprintf("Q1: < %.3f", tyg_cuts[1]),
    sprintf("Q2: %.3f – %.3f", tyg_cuts[1], tyg_cuts[2]),
    sprintf("Q3: %.3f – %.3f", tyg_cuts[2], tyg_cuts[3]),
    sprintf("Q4: ≥ %.3f", tyg_cuts[3])),
  cuts_path
)
CH$log("OK", sprintf("TyG cut-points saved: %s", cuts_path), "2b_table1")

# ============================================================
# PART 3: 糖尿病患病率明细（补充表/方法学报告）
# ============================================================
dm_detail <- data.frame(
  Category = c("Self-reported (diabe=1)", "FBG≥126 mg/dL", "HbA1c≥6.5%", "Any of above (combined)"),
  n_cases = c(
    sum(tbl_df$diabe == 1, na.rm = TRUE),
    sum(tbl_df$bl_glu >= 126, na.rm = TRUE),
    sum(tbl_df$bl_hbalc >= 6.5, na.rm = TRUE),
    sum(tbl_df$diabetes == "Yes", na.rm = TRUE)
  ),
  n_total = c(
    sum(!is.na(tbl_df$diabe)),
    sum(!is.na(tbl_df$bl_glu)),
    sum(!is.na(tbl_df$bl_hbalc)),
    nrow(tbl_df)
  )
)
dm_detail$prevalence_pct <- round(100 * dm_detail$n_cases / dm_detail$n_total, 1)

cat("\n--- Diabetes Definition Breakdown ---\n")
print(dm_detail, row.names = FALSE)

dm_path <- file.path(PATH_OUT, "tables", "Diabetes_prevalence_detail.csv")
write.csv(dm_detail, dm_path, row.names = FALSE)
CH$log("OK", sprintf("Diabetes detail saved: %s", dm_path), "2b_table1")

# ============================================================
# PART 4: 空腹 vs 非空腹对比（补充材料）
# ============================================================
tbl_df2 <- dt_base %>%
  mutate(
    fasting_grp = factor(bl_fasting, levels = c(0, 1),
                         labels = c("Non-fasting", "Fasting")),
    gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male")),
    smokev = factor(smokev, levels = c(0, 1), labels = c("Never", "Ever")),
    exercise = factor(exercise, levels = c(0, 1), labels = c("No", "Yes"))
  )

tab1_supp <- CreateTableOne(
  vars       = c("age", "bmi", "TyG_baseline", "bl_glu", "bl_tg",
                 "gender", "smokev", "exercise"),
  strata     = "fasting_grp",
  data       = tbl_df2,
  factorVars = c("gender", "smokev", "exercise"),
  addOverall = TRUE
)

tab1_supp_print <- print(tab1_supp, showAllLevels = TRUE,
                         printToggle = FALSE, noSpaces = TRUE)
tab1_supp_path <- file.path(PATH_OUT, "tables", "Table1_fast_vs_nonfast.csv")
write.csv(tab1_supp_print, tab1_supp_path, row.names = TRUE)
CH$log("OK", sprintf("Fasting vs non-fasting saved: %s", tab1_supp_path), "2b_table1")

cat("\n✓ Step 2b complete. All outputs in:", file.path(PATH_OUT, "tables"), "\n")
CH$report("2b_table1")