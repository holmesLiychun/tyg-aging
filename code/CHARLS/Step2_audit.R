# ============================================================
# Step 2: 变量审计 + 单位转换
# ============================================================
CH$reset()
if (!exists("dt_unit")) load_ck("dt_unit", "1_read")

cat("--- Step 2: Audit & Unit Conversion ---\n")

dt_unit <- dt_unit %>%
  mutate(
    # CRP: mg/L → mg/dL
    crp_mgdl = bl_crp / 10,
    ln_crp = ifelse(crp_mgdl > 0, log(crp_mgdl), NA_real_),
    
    # TyG: ln(TG × Glu / 2)
    TyG = ifelse(bl_tg > 0 & bl_glu > 0, log(bl_tg * bl_glu / 2), NA_real_),
    
    # 空腹标记
    fasting = ifelse(bl_fasting == 1, 1, 0),
    
    # ln 转换（用于 KDM 输入）
    ln_tg   = ifelse(bl_tg > 0, log(bl_tg), NA_real_),
    ln_crea = ifelse(bl_crea > 0, log(bl_crea), NA_real_),
    ln_bun  = ifelse(bl_bun > 0, log(bl_bun), NA_real_),
    ln_ua   = ifelse(bl_ua > 0, log(bl_ua), NA_real_),
    
    # TG/HDL ratio（阴性对照，v4 任务12）
    tg_hdl_ratio = ifelse(bl_hdl > 0, bl_tg / bl_hdl, NA_real_)
  )

# 审计：TyG 自算校验
if ("tyg" %in% names(dt_unit)) {
  r_check <- cor(dt_unit$TyG, dt_unit$tyg, use = "complete.obs")
  # ✅ 直接传格式和变量，不嵌套 sprintf
  CH$log(ifelse(r_check > 0.99, "OK", "WARN"),
         "TyG self-calc vs file field r=%.4f", r_check)
}

# 已有 ln_crea，但 KDM_core6 需要 lncreat（无下划线）
dt_unit$lncreat <- ifelse(dt_unit$bl_crea > 0, log(dt_unit$bl_crea + 1), NA_real_)

# 审计：空腹子样本
n_fast <- sum(dt_unit$fasting == 1, na.rm = TRUE)
# ✅ 直接传格式和变量，让 CH$log 内部处理 %%
CH$log("OK", "Fasting subsample: %d (%.1f%%)", n_fast, n_fast/nrow(dt_unit)*100)

# 基线宽表
dt_base <- dt_unit %>%
  filter(wave_year == 2011) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    age_baseline    = age,
    bmi_baseline    = bmi,
    gender_baseline = gender,
    edu_baseline    = edu,
    smokev_baseline = smokev,
    drinkl_baseline = drinkl,
    exercise_baseline = exercise,
    TyG_baseline    = TyG
  )

# ✅ 同样修复
CH$log("OK", "dt_base (2011): %d rows", nrow(dt_base))

save_ck("dt_unit", "2_audit")
save_ck("dt_base", "2_audit")
CH$report("2_audit")