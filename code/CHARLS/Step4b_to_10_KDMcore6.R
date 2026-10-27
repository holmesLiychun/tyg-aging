# ============================================================
# Step 4b ~ Step 10 统一更新版（干净版 · PDF 出图）
# ★ KDM_core6 only | IPCW 三权重 | 统一 Winsorize | NO-GO → 数据质量检查
# ★ 改造：开头 load_ck 加载 dt_unit，关键节点 save_ck
# ★ 2026-10-27 更新：Step 8 前置 load_ck、Step 9b 前置 load_ck、4b 末存 checkpoint
# ============================================================

# ---- 前置检查 ----
if (!exists("PATH_OUT")) {
  PATH_OUT <- "your_dir"
}
if (!exists("PATH_DATA")) {
  PATH_DATA <- "your_dir"
}
dir.create(PATH_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH_OUT, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH_OUT, "figures"), showWarnings = FALSE, recursive = TRUE)

# 加载必要包
pkg_list <- c("tableone", "nlme", "lme4", "lmerTest", "survival", "splines", 
              "ggplot2", "patchwork", "dplyr", "tidyr", "haven", "mice", "aod", "cmprsk", "geepack")
for (p in pkg_list) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
}

# Winsorize 函数（全局定义，后续多处调用）
winsor <- function(x, lo = 0.01, hi = 0.99) {
  x <- as.numeric(x)
  q <- quantile(x, c(lo, hi), na.rm = TRUE)
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  x
}

# ============================================================
# ★ 前置加载：从 checkpoint 读取 dt_unit（三级回退）
# ============================================================
if (!exists("dt_unit")) {
  cat("--- Loading dt_unit from checkpoint (3-level fallback) ---\n")
  success <- FALSE
  
  # 优先 Step 4b winsorized
  if (file.exists(file.path(PATH_RDATA, "ck_4b_winsor_dt_unit.rds"))) {
    success <- load_ck("dt_unit", "4b_winsor")
    if (success) cat("  ✓ Loaded dt_unit from Step 4b (winsorized)\n")
  }
  # 其次 Step 2 audit
  if (!success && file.exists(file.path(PATH_RDATA, "ck_2_audit_dt_unit.rds"))) {
    success <- load_ck("dt_unit", "2_audit")
    if (success) cat("  ✓ Loaded dt_unit from Step 2 (audit)\n")
  }
  # 最后 Step 1 原始读入
  if (!success && file.exists(file.path(PATH_RDATA, "ck_1_read_dt_unit.rds"))) {
    success <- load_ck("dt_unit", "1_read")
    if (success) cat("  ✓ Loaded dt_unit from Step 1 (read)\n")
  }
  
  if (!success) {
    stop("ERROR: dt_unit not in memory and no checkpoint found! Run Step 1 first.\n  Expected: ck_1_read_dt_unit.rds / ck_2_audit_dt_unit.rds / ck_4b_winsor_dt_unit.rds in results/checkpoints/")
  }
  cat(sprintf("  dt_unit: %d rows, %d unique IDs\n", nrow(dt_unit), n_distinct(dt_unit$ID)))
}

# ========== ID 清洗 ==========
cat("--- Pre-clean: removing rows with NA ID ---\n")
n_before <- nrow(dt_unit)
dt_unit <- dt_unit[!is.na(dt_unit$ID), ]
n_removed <- n_before - nrow(dt_unit)
cat(sprintf("  Removed %d rows with ID = NA (kept %d rows)\n", n_removed, nrow(dt_unit)))
dt_unit$ID <- as.character(dt_unit$ID)

# ============================================================
# Step 4b: Table 1 Part B — BA 结局变量基线分布
# ============================================================
CH$reset()
cat("--- Step 4b: Generating Table 1 Part B (BA outcomes) ---\n")

ba_bl <- dt_unit %>%
  filter(wave == 1 & bl_fasting == 1) %>%
  mutate(
    tyg_q = cut(tyg,
                quantile(tyg, c(0, .25, .5, .75, 1), na.rm = TRUE),
                labels = c("Q1", "Q2", "Q3", "Q4"),
                include.lowest = TRUE)
  )

ba_vars <- c("BA_resid_KDM_core6", "BA_resid_sex_KDM_core6")
ba_vars <- ba_vars[ba_vars %in% names(ba_bl)]

tab_ba <- CreateTableOne(
  vars   = ba_vars,
  strata = "tyg_q",
  data   = ba_bl,
  addOverall = TRUE
)
tab_ba_print <- print(tab_ba, printToggle = FALSE, noSpaces = TRUE)
tab_ba_path <- file.path(PATH_OUT, "tables", "Table1_PartB_BA_baseline_by_TyGq.csv")
write.csv(tab_ba_print, tab_ba_path, row.names = TRUE)
CH$log("OK", sprintf("Table 1 Part B saved: %s", tab_ba_path), "4b_table1")
cat("\n========== TABLE 1 PART B PREVIEW ==========\n")
print(tab_ba, showAllLevels = TRUE)

# ---- 合并 Part A + Part B ----
if (exists("tab1_path") && file.exists(tab1_path)) {
  partA_raw <- read.csv(tab1_path)
  rownames(partA_raw) <- make.unique(as.character(partA_raw[, 1]))
  partA <- partA_raw[, -1, drop = FALSE]
  
  partB_raw <- read.csv(tab_ba_path)
  rownames(partB_raw) <- make.unique(as.character(partB_raw[, 1]))
  partB <- partB_raw[, -1, drop = FALSE]
  
  if (ncol(partA) == ncol(partB)) {
    colnames(partB) <- colnames(partA)
    full_tab <- rbind(partA, partB)
    write.csv(full_tab, file.path(PATH_OUT, "tables", "Table1_complete.csv"))
    cat(sprintf("✓ Table 1 complete: %d rows × %d cols\n", nrow(full_tab), ncol(full_tab)))
  }
}
CH$report("4b_table1")

# ============================================================
# ★ 统一 Winsorize（放在 4b 末尾、Step 5 之前，dt_unit 层面）
# ============================================================
cat("\n--- Winsorizing BA_resid variables in dt_unit (1%-99%) ---\n")
resid_vars <- grep("^BA_resid_", names(dt_unit), value = TRUE)
for (v in resid_vars) {
  dt_unit[[v]] <- winsor(dt_unit[[v]])
}
cat(sprintf("✓ Winsorized %d BA_resid variables in dt_unit\n", length(resid_vars)))

# ★ 保存 Step 4b checkpoint（供后续步骤三级回退使用）
save_ck("dt_unit", "4b_winsor")
cat("  ✓ Checkpoint saved: ck_4b_winsor_dt_unit.rds\n")


# ============================================================
# Step 5: 横断面解耦（★ KDM_core6 版）
# ============================================================
CH$log("OK", "========== Running 5_CrossSectional ==========", "5_CS")

dt_cs <- dt_unit[dt_unit$wave == 1 & dt_unit$bl_fasting == 1, ]
dt_cs <- dt_cs[!is.na(dt_cs$BA_resid_KDM_core6) & !is.na(dt_cs$tyg), ]
cat(sprintf("Cross-sectional sample: n=%d\n", nrow(dt_cs)))

# 因子化
for (v in c("gender", "smokev", "drinkl", "edu", "exercise")) {
  if (v %in% names(dt_cs)) dt_cs[[v]] <- factor(dt_cs[[v]])
}

covar_vec <- c()
for (v in c("age", "gender", "bmi", "smokev", "drinkl", "edu", "exercise")) {
  if (v %in% names(dt_cs)) covar_vec <- c(covar_vec, v)
}
covars <- paste(covar_vec, collapse = " + ")
cat(sprintf("Covariates used: %s\n", covars))

# ---- 主模型 ----
f_kdm <- as.formula(paste("BA_resid_KDM_core6 ~ tyg +", covars))
cat("\nFitting KDM_core6 model...\n")
m_kdm <- lm(f_kdm, data = dt_cs)
s_kdm <- tryCatch(summary(m_kdm)$coefficients["tyg", ], error = function(e) NA)
if (is.null(dim(s_kdm)) && length(s_kdm) == 4) {
  # 正常
} else {
  cat("⚠️ 'tyg' not found in model coefficients! Row names are:\n")
  print(rownames(summary(m_kdm)$coefficients))
  stop("Fix variable name in formula")
}

# ---- 敏感性：sex_KDM_core6 ----
dt_cs_sex <- dt_cs[!is.na(dt_cs$BA_resid_sex_KDM_core6), ]
f_sex <- as.formula(paste("BA_resid_sex_KDM_core6 ~ tyg +", covars))
m_sex <- lm(f_sex, data = dt_cs_sex)
s_sex <- summary(m_sex)$coefficients["tyg", ]

# ---- 结果输出 ----
cat("\n========== CROSS-SECTIONAL RESULTS ==========\n")
cat(sprintf("KDM_core6:     beta=%.4f, se=%.4f, p=%.4g, R2=%.3f, n=%d\n",
            s_kdm[1], s_kdm[2], s_kdm[4], summary(m_kdm)$r.squared, nobs(m_kdm)))
cat(sprintf("sex_KDM_core6: beta=%.4f, se=%.4f, p=%.4g, R2=%.3f, n=%d\n",
            s_sex[1], s_sex[2], s_sex[4], summary(m_sex)$r.squared, nobs(m_sex)))

# ---- ★ 改为数据质量检查（不再基于方向做 GO/NO-GO）----
if (is.na(s_kdm[1]) || abs(s_kdm[1]) > 50 || s_kdm[2] > 10) {
  cat(sprintf("\n⚠️ DATA QUALITY FLAG: beta=%.4f, SE=%.4f (possible convergence issue)\n", s_kdm[1], s_kdm[2]))
} else {
  direction <- ifelse(s_kdm[1] > 0, "positive", "negative/inverse")
  sig_mark <- ifelse(s_kdm[4] < 0.05, "significant", "non-significant")
  cat(sprintf("\n✅ PROCEED: β=%.4f (%s, %s, p=%.3g) → valid for longitudinal analysis\n",
              s_kdm[1], direction, sig_mark, s_kdm[4]))
}

# ---- 保存（★ 修复：强制标量提取）----
cs_out <- data.frame(
  Model = c("KDM_core6", "sex_KDM_core6"),
  beta  = c(as.numeric(s_kdm["Estimate"]), as.numeric(s_sex["Estimate"])),
  SE    = c(as.numeric(s_kdm["Std. Error"]), as.numeric(s_sex["Std. Error"])),
  P     = c(as.numeric(s_kdm["Pr(>|t|)"]), as.numeric(s_sex["Pr(>|t|)"])),
  R2    = c(summary(m_kdm)$r.squared, summary(m_sex)$r.squared),
  n     = c(nobs(m_kdm), nobs(m_sex))
)
write.csv(cs_out, file.path(PATH_OUT, "tables", "CS_TyG_BA_results.csv"), row.names = FALSE)
CH$log("INFO", "Cross-sectional KDM_core6 beta=%.4f p=%.4g -> %s", s_kdm[1], s_kdm[4],
       ifelse(is.na(s_kdm[1]) || abs(s_kdm[1]) > 50, "FLAG", "OK"))


# ============================================================
# Step 5b: 横断面 RCS 绘图（PDF）
# ============================================================
plot_rcs_cs_fixed <- function(model, data, var = "tyg", ref_val = NULL, title = "") {
  if (is.null(ref_val)) ref_val <- median(data[[var]], na.rm = TRUE)
  seq_var <- seq(
    quantile(data[[var]], 0.01, na.rm = TRUE),
    quantile(data[[var]], 0.99, na.rm = TRUE),
    length.out = 100
  )
  newdata <- data %>% slice(1) %>% slice(rep(1, 100))
  newdata[[var]] <- seq_var
  
  num_vars <- names(data)[sapply(data, is.numeric)]
  for (v in setdiff(num_vars, var)) {
    if (v %in% names(newdata)) newdata[[v]] <- median(data[[v]], na.rm = TRUE)
  }
  fac_vars <- names(data)[sapply(data, is.factor)]
  for (v in fac_vars) {
    if (v %in% names(newdata)) newdata[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]]))
  }
  
  ref_data <- newdata[1, ]; ref_data[[var]] <- ref_val
  pred <- predict(model, newdata = newdata, se.fit = TRUE)
  pred_ref <- predict(model, newdata = ref_data, se.fit = TRUE)
  fit_centered <- pred$fit - pred_ref$fit
  
  plot_df <- data.frame(
    x = seq_var, fit = fit_centered,
    Lower = fit_centered - 1.96 * pred$se.fit,
    Upper = fit_centered + 1.96 * pred$se.fit
  )
  
  ggplot(plot_df, aes(x = x, y = fit)) +
    geom_line(size = 1.2, color = "#E63946") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#E63946") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.8) +
    geom_vline(xintercept = ref_val, linetype = "dotted", color = "gray40", size = 0.6) +
    labs(title = title, x = "TyG Index", y = "BA Residual (95% CI)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())
}

ref_tyg_cs <- median(dt_cs$tyg, na.rm = TRUE)
covars_cs <- "age + gender + bmi + smokev + drinkl + edu + exercise"

p1_cs <- plot_rcs_cs_fixed(
  lm(as.formula(paste("BA_resid_KDM_core6 ~ ns(tyg, df=3) +", covars_cs)), data = dt_cs),
  dt_cs, "tyg", ref_tyg_cs, title = "KDM_core6 (Main)"
)
p2_cs <- plot_rcs_cs_fixed(
  lm(as.formula(paste("BA_resid_sex_KDM_core6 ~ ns(tyg, df=3) +", covars_cs)), data = dt_cs_sex),
  dt_cs_sex, "tyg", ref_tyg_cs, title = "sex_KDM_core6 (Sensitivity)"
)

pdf(file.path(PATH_OUT, "figures", "CS_Spline_KDMcore6_2in1.pdf"), width = 7, height = 8)
print(p1_cs / p2_cs)
dev.off()
cat("✓ CS spline PDF saved\n")

# ============================================================
# Step 6: 纵向 LME（对角随机效应）
# ============================================================
cat("========== STEP 6: Longitudinal LME ==========\n")

dt_unit <- dt_unit[!is.na(dt_unit$ID) & dt_unit$ID != "", ]
dt_unit$ID <- as.character(dt_unit$ID)

dt_long <- dt_unit[dt_unit$bl_fasting == 1 & !is.na(dt_unit$BA_resid_KDM_core6) & !is.na(dt_unit$tyg), ] %>%
  mutate(ID = as.character(ID)) %>%
  group_by(ID) %>%
  mutate(
    time   = wave - 1,
    w1_row = which(wave == min(wave, na.rm = TRUE))[1],
    tyg_bl = tyg[w1_row],
    age_bl = age[w1_row]
  ) %>%
  ungroup()

tyg_sd_val <- sd(dt_long$tyg_bl, na.rm = TRUE)
dt_long$tyg_bl_sd <- dt_long$tyg_bl / tyg_sd_val

cat(sprintf("Longitudinal: %d rows, %d unique IDs\n", nrow(dt_long), length(unique(dt_long$ID))))

dt_long$gender <- factor(dt_long$gender)
for (v in c("smokev", "drinkl", "edu", "exercise")) {
  if (v %in% names(dt_long)) dt_long[[v]] <- factor(dt_long[[v]])
}

dt_long <- dt_long %>%
  dplyr::select(ID, time, BA_resid_KDM_core6, tyg_bl_sd, age_bl, gender, bmi, smokev, drinkl, edu, exercise) %>%
  na.omit() %>%
  mutate(ID = droplevels(as.factor(ID)))

cat(sprintf("After NA removal: %d rows, %d IDs\n", nrow(dt_long), nlevels(dt_long$ID)))

ctrl <- lmeControl(maxIter = 1000, msMaxIter = 1000, opt = "nlminb", returnObject = TRUE)

m_lme_main <- lme(
  BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi + smokev + drinkl + edu,
  random = list(ID = pdDiag(~ 1)),
  data   = dt_long,
  na.action = na.omit,
  control = ctrl
)

cat("\n--- Main LME (unweighted) ---\n")
print(summary(m_lme_main)$tTable)

tt <- summary(m_lme_main)$tTable
key <- "tyg_bl_sd:time"
cat(sprintf("\n*** TyG_baseline × time: β=%.4f, SE=%.4f, p=%.4f ***\n",
            tt[key, "Value"], tt[key, "Std.Error"], tt[key, "p-value"]))

lme_main_int <- data.frame(
  Model = "Unweighted LME",
  Term  = rownames(tt)[grep("tyg_bl_sd", rownames(tt))],
  beta  = tt[grep("tyg_bl_sd", rownames(tt)), "Value"],
  SE    = tt[grep("tyg_bl_sd", rownames(tt)), "Std.Error"],
  p     = tt[grep("tyg_bl_sd", rownames(tt)), "p-value"],
  N_rows = nrow(dt_long), N_ids = nlevels(dt_long$ID), AIC = AIC(m_lme_main)
)
write.csv(lme_main_int, file.path(PATH_OUT, "tables", "LME_main_unweighted.csv"), row.names = FALSE)

# ★ 保存 Step 6 checkpoint
save_ck("dt_long", "6_lme")
cat("  ✓ Checkpoint saved: ck_6_lme_dt_long.rds\n")


# ============================================================
# Step 6b: 时变 TyG + 时变 BA — GEE vs LME
# ============================================================
cat("========== STEP 6b: Time-varying TyG + BA (GEE vs LME) ==========\n")

if (!require("geepack", quietly = TRUE)) {
  install.packages("geepack", repos = "https://cloud.r-project.org")
  library(geepack)
}

dt_tv <- dt_unit %>%
  filter(bl_fasting == 1 & !is.na(BA_resid_KDM_core6) & !is.na(tyg)) %>%
  mutate(ID = as.character(ID)) %>%
  group_by(ID) %>%
  mutate(
    time   = wave - 1,
    tyg_bl = tyg[which(wave == min(wave, na.rm = TRUE))[1]],
    age_bl = age[which(wave == min(wave, na.rm = TRUE))[1]]
  ) %>%
  ungroup() %>%
  mutate(tyg_tv_sd = as.numeric(scale(tyg))) %>%
  group_by(ID) %>%
  mutate(
    tyg_mean_wp = mean(tyg, na.rm = TRUE),
    tyg_wp_dev  = tyg - tyg_mean_wp
  ) %>%
  ungroup()

dt_tv$gender <- as.factor(dt_tv$gender)
for (v in c("smokev", "drinkl", "edu", "exercise")) {
  if (v %in% names(dt_tv)) dt_tv[[v]] <- as.factor(dt_tv[[v]])
}

dt_tv <- dt_tv %>%
  dplyr::select(ID, time, BA_resid_KDM_core6, tyg_tv_sd, tyg_wp_dev, tyg_mean_wp, 
                age_bl, gender, bmi, smokev, drinkl, edu) %>%
  na.omit() %>%
  mutate(ID = droplevels(as.factor(ID)))

cat(sprintf("Time-varying dataset: %d rows, %d IDs\n", nrow(dt_tv), nlevels(dt_tv$ID)))

cat("\n--- GEE 1: Exchangeable ---\n")
gee_exch <- geeglm(
  BA_resid_KDM_core6 ~ tyg_tv_sd * time + age_bl + gender + bmi + smokev + drinkl + edu,
  id = ID, data = dt_tv, family = gaussian, corstr = "exchangeable"
)
print(summary(gee_exch)$coefficients)

cat("\n--- GEE 2: AR(1) ---\n")
gee_ar1 <- geeglm(
  BA_resid_KDM_core6 ~ tyg_tv_sd * time + age_bl + gender + bmi + smokev + drinkl + edu,
  id = ID, data = dt_tv, family = gaussian, corstr = "ar1"
)
print(summary(gee_ar1)$coefficients)

cat("\n--- LME (Time-varying TyG) ---\n")
lme_tv <- lme(
  BA_resid_KDM_core6 ~ tyg_tv_sd * time + age_bl + gender + bmi + smokev + drinkl + edu,
  random = list(ID = pdDiag(~ 1)), data = dt_tv, na.action = na.omit, control = ctrl
)
lme_tv_tt <- summary(lme_tv)$tTable
print(lme_tv_tt)

cat("\n--- DIAGNOSTIC: Within-Person Centering ---\n")
gee_wp <- geeglm(
  BA_resid_KDM_core6 ~ tyg_mean_wp + tyg_wp_dev * time + age_bl + gender + bmi,
  id = ID, data = dt_tv, family = gaussian, corstr = "exchangeable"
)
gee_wp_sum <- summary(gee_wp)
wp_coef <- as.data.frame(gee_wp_sum$coefficients)
print(wp_coef)

if ("tyg_wp_dev:time" %in% rownames(wp_coef)) {
  p_wp    <- wp_coef["tyg_wp_dev:time", "Pr(>|W|)"]
  beta_wp <- wp_coef["tyg_wp_dev:time", "Estimate"]
  cat(sprintf("\n*** Within-person TyG × time: β=%.4f, p=%.4g ***\n", beta_wp, p_wp))
  if (!is.na(p_wp) && p_wp < 0.05) {
    cat("  ✅ CONCLUSION: 个体内交互显著！原模型交互 null 是【模型误设】。\n")
  } else {
    cat("  ⚠️ 个体内交互也不显著，继续检查非线性...\n")
  }
}

cat("\n--- DIAGNOSTIC: Non-linear (ns spline) ---\n")
if (require("splines")) {
  gee_spline <- geeglm(
    BA_resid_KDM_core6 ~ ns(tyg_tv_sd, df=3) * time + age_bl + gender + bmi,
    id = ID, data = dt_tv, family = gaussian, corstr = "exchangeable"
  )
  spl_coef <- as.data.frame(summary(gee_spline)$coefficients)
  print(spl_coef[grep("ns|time", rownames(spl_coef)), ])
}

tv_summary <- data.frame(
  Model = c("GEE_Exch", "GEE_AR1", "LME_tv", "GEE_WithinPerson"),
  Term  = c("tyg_tv_sd:time", "tyg_tv_sd:time", "tyg_tv_sd:time", "tyg_wp_dev:time"),
  Beta  = c(
    ifelse("tyg_tv_sd:time" %in% rownames(as.data.frame(summary(gee_exch)$coefficients)), 
           as.data.frame(summary(gee_exch)$coefficients)["tyg_tv_sd:time", "Estimate"], NA),
    ifelse("tyg_tv_sd:time" %in% rownames(as.data.frame(summary(gee_ar1)$coefficients)), 
           as.data.frame(summary(gee_ar1)$coefficients)["tyg_tv_sd:time", "Estimate"], NA),
    ifelse("tyg_tv_sd:time" %in% rownames(lme_tv_tt), lme_tv_tt["tyg_tv_sd:time", "Value"], NA),
    ifelse("tyg_wp_dev:time" %in% rownames(wp_coef), wp_coef["tyg_wp_dev:time", "Estimate"], NA)
  ),
  P_value = c(
    ifelse("tyg_tv_sd:time" %in% rownames(as.data.frame(summary(gee_exch)$coefficients)), 
           as.data.frame(summary(gee_exch)$coefficients)["tyg_tv_sd:time", "Pr(>|W|)"], NA),
    ifelse("tyg_tv_sd:time" %in% rownames(as.data.frame(summary(gee_ar1)$coefficients)), 
           as.data.frame(summary(gee_ar1)$coefficients)["tyg_tv_sd:time", "Pr(>|W|)"], NA),
    ifelse("tyg_tv_sd:time" %in% rownames(lme_tv_tt), lme_tv_tt["tyg_tv_sd:time", "p-value"], NA),
    ifelse("tyg_wp_dev:time" %in% rownames(wp_coef), wp_coef["tyg_wp_dev:time", "Pr(>|W|)"], NA)
  )
)
write.csv(tv_summary, file.path(PATH_OUT, "tables", "Step6b_TimeVarying_GEE_LME.csv"), row.names = FALSE)
cat(sprintf("\n✓ Step 6b complete.\n"))


# ============================================================
# Step 7: Sensitivity Analyses
# ============================================================
cat("========== Running 7_Sensitivity ==========\n")

dt_paired <- dt_long %>% group_by(ID) %>% filter(n() >= 2 & all(c(0, 2) %in% time)) %>% ungroup()
cat(sprintf("  Paired sample: %d rows (%d IDs)\n", nrow(dt_paired), n_distinct(dt_paired$ID)))
if (nrow(dt_paired) > 10) {
  m_paired <- lmer(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + bmi + (1 | ID), data = dt_paired)
  cat("\n--- 7.1 Paired ---\n"); print(summary(m_paired), corr = FALSE)
  paired_coef <- as.data.frame(summary(m_paired)$coefficients)
}

m_blstrat <- lmer(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi + (1 | ID), data = dt_long)
cat("\n--- 7.2 Baseline stratified ---\n"); print(summary(m_blstrat), corr = FALSE)
blstrat_coef <- as.data.frame(summary(m_blstrat)$coefficients)

if ("BA_resid_sex_KDM_core6" %in% names(dt_long)) {
  m_sex_lme <- lmer(BA_resid_sex_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi + (1 | ID), data = dt_long)
  cat("\n--- 7.3 sex_KDM_core6 ---\n"); print(summary(m_sex_lme), corr = FALSE)
}

dt_trim <- dt_long %>% filter(tyg_bl_sd >= quantile(tyg_bl_sd, 0.01, na.rm=TRUE) & tyg_bl_sd <= quantile(tyg_bl_sd, 0.99, na.rm=TRUE))
m_trim <- lmer(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi + (1 | ID), data = dt_trim)
cat("\n--- 7.4 Trimmed ---\n"); print(summary(m_trim), corr = FALSE)

m_min <- lmer(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + bmi + (1 | ID), data = dt_long)
cat("\n--- 7.5 Minimal ---\n"); print(summary(m_min), corr = FALSE)

CH$log("OK", "Step 7 Sensitivity complete")


# ============================================================
# Step 8: Cox 生存分析（★ 加前置 load_ck）
# ============================================================
cat("========== STEP 8: COX SURVIVAL ANALYSIS ==========\n")

# ★ 前置加载：优先从 checkpoint 读取 surv_final
if (!exists("surv_final")) {
  cat("--- Loading surv_final from checkpoint ---\n")
  success <- FALSE
  
  # 优先 Step 8 自己的 checkpoint
  if (file.exists(file.path(PATH_RDATA, "ck_8_cox_surv_final.rds"))) {
    success <- load_ck("surv_final", "8_cox")
    if (success) cat("  ✓ Loaded surv_final from Step 8 checkpoint\n")
  }
  
  if (!success) {
    # ---- 完整构建 surv_final ----
    cat("  Building surv_final from scratch...\n")
    
    dta_d   <- "D:/charls/charlsdone/H_CHARLS_D_Data.dta"
    dta_eol <- "D:/charls/charlsdone/H_CHARLS_EOL_a.dta"
    
    D   <- read_dta(dta_d)
    EOL <- read_dta(dta_eol)
    cat(sprintf("D file: %d rows | EOL file: %d rows\n", nrow(D), nrow(EOL)))
    
    surv <- D %>%
      filter(!is.na(r1iwy)) %>%
      dplyr::select(ID, r1iwy, r1iwm, r2iwy, r2iwm, r3iwy, r3iwm, r4iwy, r4iwm, radyear, radmonth) %>%
      mutate(ID = as.character(ID)) %>%
      mutate(
        w1_dt = ymd(paste0(as.numeric(r1iwy), "-", 
                           if_else(is.na(as.numeric(r1iwm)) | as.numeric(r1iwm) < 1 | as.numeric(r1iwm) > 12, 7L, as.integer(as.numeric(r1iwm))), "-15")),
        w2_dt = ymd(paste0(as.numeric(r2iwy), "-", 
                           if_else(is.na(as.numeric(r2iwm)) | as.numeric(r2iwm) < 1 | as.numeric(r2iwm) > 12, 7L, as.integer(as.numeric(r2iwm))), "-15")),
        w3_dt = ymd(paste0(as.numeric(r3iwy), "-", 
                           if_else(is.na(as.numeric(r3iwm)) | as.numeric(r3iwm) < 1 | as.numeric(r3iwm) > 12, 7L, as.integer(as.numeric(r3iwm))), "-15")),
        w4_dt = ymd(paste0(as.numeric(r4iwy), "-", 
                           if_else(is.na(as.numeric(r4iwm)) | as.numeric(r4iwm) < 1 | as.numeric(r4iwm) > 12, 7L, as.integer(as.numeric(r4iwm))), "-15")),
        base_dt  = w1_dt,
        death_dt = ymd(paste0(as.numeric(radyear), "-",
                              if_else(is.na(as.numeric(radmonth)) | as.numeric(radmonth) < 1 | as.numeric(radmonth) > 12, 7L, as.integer(as.numeric(radmonth))), "-15")),
        death = if_else(!is.na(as.numeric(radyear)), 1L, 0L)
      ) %>%
      mutate(
        last_alive_num = pmax(
          as.numeric(w1_dt), as.numeric(w2_dt), 
          as.numeric(w3_dt), as.numeric(w4_dt), 
          na.rm = TRUE
        ),
        last_alive_dt = as.Date(last_alive_num, origin = "1970-01-01"),
        event_dt = if_else(death == 1, death_dt, last_alive_dt),
        fu_years = as.numeric(event_dt - base_dt, units = "days") / 365.25
      ) %>%
      dplyr::select(ID, base_dt, death, fu_years)
    
    cat(sprintf("surv constructed: %d rows\n", nrow(surv)))
    
    # 确保 dt_unit 在内存
    if (!exists("dt_unit")) {
      if (file.exists(file.path(PATH_RDATA, "ck_1_read_dt_unit.rds"))) {
        load_ck("dt_unit", "1_read")
        cat("  ✓ Loaded dt_unit for surv_final build\n")
      } else {
        stop("ERROR: Need dt_unit to build surv_final! Run Step 1 first.")
      }
    }
    
    base <- dt_unit %>% filter(wave == 1 & bl_fasting == 1) %>% distinct(ID, .keep_all = TRUE) %>% mutate(ID = as.character(ID))
    surv_final <- base %>% left_join(surv, by = "ID") %>%
      mutate(
        fu_years = ifelse(is.na(fu_years), 7.5, fu_years),
        death = ifelse(is.na(death), 0, death),
        tyg_sd = as.numeric(scale(tyg)),
        tyg_q = cut(tyg, quantile(tyg, c(0, .25, .5, .75, 1), na.rm = TRUE), labels = c("Q1","Q2","Q3","Q4"), include.lowest = TRUE)
      ) %>% filter(!is.na(tyg) & fu_years > 0)
    
    resid_vars_surv <- grep("^BA_resid_", names(surv_final), value = TRUE)
    for (v in resid_vars_surv) { surv_final[[v]] <- winsor(surv_final[[v]]) }
    cat(sprintf("✓ Winsorized %d resid vars in surv_final\n", length(resid_vars_surv)))
    
    id_map <- D %>% dplyr::select(ID, ID_w1) %>% mutate(ID_w1 = as.character(ID_w1), ID = as.character(ID)) %>% distinct()
    eol_cod <- EOL %>% dplyr::select(ID_w1, ragcod, racod_c) %>% mutate(ID_w1 = as.character(ID_w1)) %>%
      left_join(id_map, by = "ID_w1") %>% filter(!is.na(ID)) %>% dplyr::select(ID, ragcod, racod_c)
    surv_final <- surv_final %>% left_join(eol_cod, by = "ID") %>%
      mutate(cvd_death = ifelse(death == 1 & ragcod == 2, 1, 0), cancer_death = ifelse(death == 1 & ragcod == 1, 1, 0))
    
    cat(sprintf("surv_final: %d persons, %d deaths\n", nrow(surv_final), sum(surv_final$death, na.rm=TRUE)))
    
    # ★ 保存 checkpoint
    save_ck("surv_final", "8_cox")
    cat("  ✓ Checkpoint saved: ck_8_cox_surv_final.rds\n")
  }
} else {
  cat("  surv_final already in memory, skipping rebuild\n")
}

# ---- Cox 模型拟合 ----
if (sum(surv_final$death, na.rm=TRUE) > 10) {
  cox_all_sd <- coxph(Surv(fu_years, death) ~ tyg_sd + age + gender + bmi + smokev + drinkl + edu, data = surv_final)
  cat("\n--- All-cause (per SD) ---\n"); print(summary(cox_all_sd))
  
  if ("BA_resid_KDM_core6" %in% names(surv_final)) {
    cox_all_kdm <- coxph(Surv(fu_years, death) ~ tyg_sd + BA_resid_KDM_core6 + age + gender + bmi + smokev + drinkl + edu, data = surv_final)
    cat("\n--- + BA_resid_KDM_core6 ---\n"); print(summary(cox_all_kdm))
  }
}
if (sum(surv_final$cvd_death, na.rm=TRUE) >= 10) {
  cox_cvd_sd <- coxph(Surv(fu_years, cvd_death) ~ tyg_sd + age + gender + bmi + smokev + drinkl + edu, data = surv_final %>% filter(!is.na(cvd_death)))
  cat("\n--- CVD ---\n"); print(summary(cox_cvd_sd))
}
if (sum(surv_final$cancer_death, na.rm=TRUE) >= 5) {
  cox_cancer <- coxph(Surv(fu_years, cancer_death) ~ tyg_sd + age + gender + bmi + smokev + drinkl + edu, data = surv_final %>% filter(!is.na(cancer_death)))
  cat("\n--- Cancer ---\n"); print(summary(cox_cancer))
}

# ---- KM 图 ----
if (sum(surv_final$death, na.rm = TRUE) >= 10) {
  km_fit <- survfit(Surv(fu_years, death) ~ tyg_q, data = surv_final)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv, strata = rep(names(km_fit$strata), km_fit$strata))
  
  p_km <- ggplot() +
    geom_step(data = km_data, aes(x = time, y = surv, color = strata), linewidth = 1.2) +
    scale_color_manual(values = c("#2E86AB", "#E63946", "#06D6A0", "#FFD166"), name = "TyG Quartile") +
    labs(title = "Kaplan-Meier Survival by TyG Quartile", x = "Follow-up (years)", y = "Survival Probability") +
    theme_bw(base_size = 13) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  pdf(file.path(PATH_OUT, "figures", "Fig4_KM_by_TyG_quartile.pdf"), width = 7, height = 5.5)
  print(p_km)
  dev.off()
  cat("✓ KM PDF saved\n")
}

# ★ 保存/刷新 Step 8 checkpoint
save_ck("surv_final", "8_cox")
cat("  ✓ Checkpoint refreshed: ck_8_cox_surv_final.rds\n")


# ============================================================
# Step 9: IPCW
# ============================================================
cat("========== STEP 9: IPCW (Fixed) ==========\n")

# 如果 dt_long 不在内存，尝试加载
if (!exists("dt_long")) {
  cat("--- Loading dt_long from checkpoint ---\n")
  if (!load_ck("dt_long", "6_lme")) {
    stop("ERROR: Need dt_long! Run Step 6 first or check ck_6_lme_dt_long.rds")
  }
  cat("  ✓ dt_long loaded\n")
}

ipcw_data <- dt_unit %>%
  mutate(ID = as.character(ID), time = wave - 1) %>%
  group_by(ID) %>%
  mutate(tyg_bl = tyg[wave == 1][1], age_bl = age[wave == 1][1], max_wave = max(wave, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(R_final = ifelse(max_wave >= 3, 1L, 0L)) %>%
  dplyr::select(ID, wave, time, R_final, age_bl, gender, tyg_bl, bmi, smokev, drinkl, edu, exercise, starts_with("BA_resid_"))

tyg_sd_val <- ipcw_data %>% filter(wave == 1) %>% pull(tyg_bl) %>% na.omit() %>% sd()
ipcw_data <- ipcw_data %>% mutate(tyg_bl_sd = tyg_bl / tyg_sd_val)

baseline_rows <- ipcw_data %>% filter(wave == 1) %>%
  filter(!is.na(age_bl), !is.na(gender), !is.na(tyg_bl_sd)) %>%
  distinct(ID, .keep_all = TRUE)

weight_specs <- list(
  list(name = "sw1", formula = R_final ~ age_bl + gender + tyg_bl_sd),
  list(name = "sw2", formula = R_final ~ age_bl + gender + tyg_bl_sd + bmi + smokev + drinkl + edu),
  list(name = "sw3", formula = R_final ~ age_bl + gender + tyg_bl_sd + bmi + smokev + drinkl + edu + exercise)
)

for (ws in weight_specs) {
  fit_glm <- tryCatch(glm(ws$formula, data = baseline_rows, family = binomial()), 
                      error = function(e) { cat(sprintf("  ⚠️ %s failed: %s\n", ws$name, e$message)); NULL })
  if (is.null(fit_glm)) {
    cat(sprintf("  → %s: using uniform weight\n", ws$name))
    baseline_rows[[ws$name]] <- 1
    next
  }
  pred_prob <- pmin(pmax(predict(fit_glm, type = "response"), 1e-4), 1 - 1e-4)
  sw <- ifelse(baseline_rows$R_final == 1, 1 / pred_prob, 1 / (1 - pred_prob))
  sw <- pmin(sw, quantile(sw, 0.99, na.rm = TRUE))
  baseline_rows[[ws$name]] <- sw
  cat(sprintf("  ✓ %s created (range: %.2f - %.2f)\n", ws$name, min(sw), max(sw)))
}

sw_cols <- c("sw1", "sw2", "sw3")
existing_sw <- sw_cols[sw_cols %in% names(baseline_rows)]
baseline_weights <- baseline_rows %>%
  dplyr::select(ID, all_of(existing_sw)) %>%
  rename_with(~ paste0("ipcw_", .), -ID)

if (!"ipcw_sw3" %in% names(baseline_weights)) {
  cat("  ⚠️ sw3 not found, filling with 1\n")
  baseline_weights$ipcw_sw3 <- 1
}

ipcw_data_fixed <- ipcw_data %>% left_join(baseline_weights, by = "ID")
w0 <- median(baseline_weights$ipcw_sw1, na.rm = TRUE)
eps <- w0 * 0.01
for (w in c("ipcw_sw1", "ipcw_sw2", "ipcw_sw3")) {
  ipcw_data_fixed[[w]] <- ifelse(is.na(ipcw_data_fixed[[w]]), w0, ipcw_data_fixed[[w]])
}

fit_base <- ipcw_data_fixed %>%
  dplyr::select(ID, time, tyg_bl_sd, age_bl, gender, bmi, ipcw_sw1, ipcw_sw2, ipcw_sw3, BA_resid_KDM_core6) %>%
  na.omit() %>% mutate(ID = as.factor(ID))

cat(sprintf("fit_base: %d rows, %d IDs\n", nrow(fit_base), nlevels(fit_base$ID)))

m_base <- lme(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi,
              random = list(ID = pdDiag(~ 1)), data = fit_base, weights = varFixed(~ ipcw_sw1), 
              na.action = na.omit, control = ctrl)

compare_df <- data.frame(Weight_Model = character(), TyG_perSD_time_beta = numeric(), 
                         SE = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

for (w in c("ipcw_sw1", "ipcw_sw2", "ipcw_sw3")) {
  fit_i <- fit_base; fit_i$wt_safe <- 1 / pmax(fit_i[[w]], eps)
  m_i <- tryCatch(update(m_base, data = fit_i, weights = varFixed(~ wt_safe)), error = function(e) NULL)
  if (is.null(m_i)) next
  tt <- summary(m_i)$tTable
  key <- "tyg_bl_sd:time"
  if (!key %in% rownames(tt)) next
  compare_df <- rbind(compare_df, data.frame(
    Weight_Model = switch(w, ipcw_sw1 = "SW1: age+gender+TyG", ipcw_sw2 = "SW2: +bmi+smoke+drink+edu", ipcw_sw3 = "SW3: +exercise"),
    TyG_perSD_time_beta = round(tt[key, "Value"], 4), SE = round(tt[key, "Std.Error"], 4),
    p_value = round(tt[key, "p-value"], 4), stringsAsFactors = FALSE
  ))
}

write.csv(compare_df, file.path(PATH_OUT, "tables", "08_ipcw_results.csv"), row.names = FALSE)
cat("✓ IPCW results saved\n"); print(compare_df)

p_ipcw_forest <- ggplot(compare_df, aes(x = TyG_perSD_time_beta, y = reorder(Weight_Model, TyG_perSD_time_beta))) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_errorbarh(aes(xmin = TyG_perSD_time_beta - SE * 1.96, xmax = TyG_perSD_time_beta + SE * 1.96), height = 0.2, color = "#2E86AB", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "IPCW-Weighted TyG × Time Interaction", x = "β (TyG per SD × Time, years)", y = NULL) +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file.path(PATH_OUT, "figures", "Fig1_IPCW_Estimates.pdf"), width = 7, height = 4.5)
print(p_ipcw_forest)
dev.off()
cat("✓ Fig 1 PDF saved\n")

# ★ 保存 Step 9 checkpoint
save_ck("ipcw_data_fixed", "9_ipcw")
cat("  ✓ Checkpoint saved: ck_9_ipcw_ipcw_data_fixed.rds\n")


# ============================================================
# Step 9b (真实 Bootstrap 版)
# ============================================================
cat("========== STEP 9b: Real Bootstrap 500 ==========\n")

# ★ 前置加载 dt_long
if (!exists("dt_long")) {
  cat("--- Loading dt_long from checkpoint ---\n")
  if (!load_ck("dt_long", "6_lme")) {
    stop("ERROR: Need dt_long for bootstrap! Run Step 6 first.")
  }
  cat("  ✓ dt_long loaded\n")
}

set.seed(2024)

n_boot <- 500
boot_betas <- numeric(n_boot)

id_vec <- levels(dt_long$ID)

for (i in 1:n_boot) {
  boot_ids <- sample(id_vec, replace = TRUE)
  boot_data <- dt_long[dt_long$ID %in% boot_ids, ]
  boot_data$ID <- droplevels(as.factor(as.character(boot_data$ID)))
  
  boot_fit <- tryCatch(
    lme(BA_resid_KDM_core6 ~ tyg_bl_sd * time + age_bl + gender + bmi + smokev + drinkl + edu,
        random = list(ID = pdDiag(~ 1)),
        data = boot_data,
        control = lmeControl(maxIter = 200, msMaxIter = 200, returnObject = TRUE)),
    error = function(e) NULL
  )
  
  if (!is.null(boot_fit)) {
    tt_b <- summary(boot_fit)$tTable
    if ("tyg_bl_sd:time" %in% rownames(tt_b)) {
      boot_betas[i] <- tt_b["tyg_bl_sd:time", "Value"]
    }
  }
  
  if (i %% 50 == 0) cat(sprintf("  Bootstrap %d/%d complete\n", i, n_boot))
}

boot_betas <- boot_betas[!is.na(boot_betas)]
boot_ci_lower   <- quantile(boot_betas, 0.025)
boot_ci_upper   <- quantile(boot_betas, 0.975)
boot_beta_point <- mean(boot_betas)
boot_prop_positive <- mean(boot_betas > 0)

cat(sprintf("Bootstrap: β=%.4f, 95%% CI [%.4f, %.4f], P(>0)=%.3f\n",
            boot_beta_point, boot_ci_lower, boot_ci_upper, boot_prop_positive))

boot_plot_data <- data.frame(beta = boot_betas)

p_bootstrap <- ggplot(boot_plot_data, aes(x = beta)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "#2E86AB", alpha = 0.7, color = "white") +
  geom_density(color = "#E63946", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  geom_vline(xintercept = c(boot_ci_lower, boot_ci_upper), linetype = "dotted", color = "#E63946", size = 0.8) +
  labs(title = "Bootstrap Distribution of TyG × Time (IPCW, 500 resamples)", 
       x = "β (TyG × Time)", y = "Density") +
  theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file.path(PATH_OUT, "figures", "Fig2_Bootstrap_Distribution.pdf"), width = 7, height = 5)
print(p_bootstrap)
dev.off()
cat("✓ Fig 2 PDF saved (real bootstrap)\n")


# ============================================================
# Step 9c: Fig 3 — 敏感性森林图
# ============================================================
cat("========== STEP 9c: Sensitivity Forest Plot ==========\n")

sens_rows <- list()

if (exists("paired_coef") && nrow(paired_coef) > 0) {
  cf <- paired_coef[grep("tyg", rownames(paired_coef)), , drop = FALSE]
  if (nrow(cf) > 0) sens_rows$`Paired (W1+W3)` <- cf
}
if (exists("blstrat_coef") && nrow(blstrat_coef) > 0) {
  cf <- blstrat_coef[grep("tyg_bl_sd", rownames(blstrat_coef)), , drop = FALSE]
  if (nrow(cf) > 0) sens_rows$`Baseline stratified` <- cf
}
if (exists("m_sex_lme")) {
  coef_sex <- as.data.frame(summary(m_sex_lme)$coefficients)
  cf <- coef_sex[grep("tyg", rownames(coef_sex)), , drop = FALSE]
  if (nrow(cf) > 0) sens_rows$`sex_KDM_core6` <- cf
}
if (exists("m_trim")) {
  coef_trim <- as.data.frame(summary(m_trim)$coefficients)
  cf <- coef_trim[grep("tyg", rownames(coef_trim)), , drop = FALSE]
  if (nrow(cf) > 0) sens_rows$`Trimmed (1-99%)` <- cf
}
if (exists("m_min")) {
  coef_min <- as.data.frame(summary(m_min)$coefficients)
  cf <- coef_min[grep("tyg", rownames(coef_min)), , drop = FALSE]
  if (nrow(cf) > 0) sens_rows$`Minimal (age+bmi)` <- cf
}
if (exists("m_lme_main")) {
  tt_main <- as.data.frame(summary(m_lme_main)$tTable)
  cf <- tt_main[grep("tyg_bl_sd", rownames(tt_main)), , drop = FALSE]
  if (nrow(cf) > 0) {
    colnames(cf)[colnames(cf) == "Value"] <- "Estimate"
    colnames(cf)[colnames(cf) == "Std.Error"] <- "Std. Error"
    sens_rows$`Main (unweighted)` <- cf
  }
}

if (length(sens_rows) > 0) {
  forest_df <- data.frame(Model = character(), beta = numeric(), SE = numeric(), 
                          lower = numeric(), upper = numeric(), stringsAsFactors = FALSE)
  
  for (nm in names(sens_rows)) {
    cf <- sens_rows[[nm]]
    int_rows <- grepl(":", rownames(cf))
    if (sum(int_rows) > 0) {
      cf_use <- cf[int_rows, , drop = FALSE]
    } else {
      cf_use <- cf[1, , drop = FALSE]
    }
    
    if (!"Estimate" %in% colnames(cf_use) || !"Std. Error" %in% colnames(cf_use)) next
    
    beta_val <- as.numeric(cf_use[1, "Estimate"])
    se_val   <- as.numeric(cf_use[1, "Std. Error"])
    
    if (is.na(beta_val) || is.na(se_val)) next
    
    forest_df <- rbind(forest_df, data.frame(
      Model = nm, beta = beta_val, SE = se_val,
      lower = beta_val - 1.96 * se_val, upper = beta_val + 1.96 * se_val,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  ✓ %s: beta=%.4f, SE=%.4f\n", nm, beta_val, se_val))
  }
  
  if (nrow(forest_df) > 0) {
    forest_df$Model <- factor(forest_df$Model, levels = rev(forest_df$Model))
    
    p_forest <- ggplot(forest_df, aes(x = beta, y = Model)) +
      geom_point(size = 3.5, color = "#E63946") +
      geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25, color = "#E63946", linewidth = 1.1) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
      labs(title = "Sensitivity: TyG × Time Interaction", x = expression(beta ~ "(TyG × Time)"), y = NULL) +
      theme_bw(base_size = 12) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    pdf(file.path(PATH_OUT, "figures", "Fig3_Sensitivity_Forest.pdf"), width = 7, height = 4.5)
    print(p_forest)
    dev.off()
    write.csv(forest_df, file.path(PATH_OUT, "tables", "Fig3_forest_data.csv"), row.names = FALSE)
    cat(sprintf("✓ Fig 3 PDF saved (%d models)\n", nrow(forest_df)))
  }
}


# ============================================================
# Step 10: Export（修复版）
# ============================================================
cat("========== STEP 10: Export ==========\n")

write.csv(dt_unit, file.path(PATH_OUT, "01_dt_unit_long.csv"), row.names = FALSE)
write.csv(surv_final, file.path(PATH_OUT, "03_surv_final.csv"), row.names = FALSE)
write.csv(ipcw_data_fixed, file.path(PATH_OUT, "04_ipcw_data.csv"), row.names = FALSE)

if (exists("cox_all_sd")) {
  safe_extract_cox <- function(model, term = "tyg_sd") {
    if (!inherits(model, "coxph")) return(c(HR = NA_real_, p = NA_real_))
    cf <- summary(model)$coefficients
    if (!term %in% rownames(cf)) {
      cat(sprintf("  ⚠️ '%s' not in cox model coefficients.\n", term))
      return(c(HR = NA_real_, p = NA_real_))
    }
    c(HR = exp(cf[term, "coef"]), p = cf[term, "Pr(>|z|)"])
  }
  
  cox1 <- safe_extract_cox(cox_all_sd, "tyg_sd")
  cox2 <- if (exists("cox_cvd_sd")) safe_extract_cox(cox_cvd_sd, "tyg_sd") else c(HR = NA_real_, p = NA_real_)
  cox3 <- if (exists("cox_cancer")) safe_extract_cox(cox_cancer, "tyg_sd") else c(HR = NA_real_, p = NA_real_)
  
  cox_df <- data.frame(
    Model = c("AllCause_SD", "CVD_SD", "Cancer_SD"),
    HR    = c(cox1["HR"], cox2["HR"], cox3["HR"]),
    p     = c(cox1["p"], cox2["p"], cox3["p"]),
    stringsAsFactors = FALSE
  )
  write.csv(cox_df, file.path(PATH_OUT, "05_cox_summary.csv"), row.names = FALSE)
  cat("✓ Cox summary saved\n")
  print(cox_df)
}

if (exists("surv_final")) {
  fu_sum <- data.frame(
    Metric = c("Median FU (years)", "Total deaths", "CVD deaths", "Cancer deaths", "Total person-years"),
    Value  = c(
      as.numeric(median(surv_final$fu_years, na.rm = TRUE)),
      as.numeric(sum(surv_final$death, na.rm = TRUE)),
      as.numeric(sum(surv_final$cvd_death, na.rm = TRUE)),
      as.numeric(sum(surv_final$cancer_death, na.rm = TRUE)),
      as.numeric(sum(surv_final$fu_years, na.rm = TRUE))
    ),
    stringsAsFactors = FALSE
  )
  write.csv(fu_sum, file.path(PATH_OUT, "09_followup_summary.csv"), row.names = FALSE)
}

cat(sprintf("\n✓✓✓ All done! Output at: %s\n", PATH_OUT))