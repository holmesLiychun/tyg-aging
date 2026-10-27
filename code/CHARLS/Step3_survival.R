# ============================================================
# Step 3: 死亡随访库构建
# ============================================================
CH$reset()

# 确保配置已加载
if (!exists("HAS_MORTALITY")) {
  if (file.exists("Step0_config.R")) {
    source("Step0_config.R")
  } else {
    HAS_MORTALITY <- FALSE
    PATH_OUT <- "results"
    PATH_RDATA <- "results/checkpoints"
  }
}

if (!exists("dt_unit")) load_ck("dt_unit", "2_audit")

cat("--- Step 3: Survival ---\n")

# 探测死亡相关列
cand <- c("RwIWSTAT","DIED","INDV_L_DIED","EXB001_1","EXB001_2",
          "radyear","radmonth","r1iwstat","r2iwstat","r3iwstat","r4iwstat","r5iwstat")
have <- intersect(tolower(cand), tolower(names(dt_unit)))
CH$log("OK", if(length(have)) sprintf("death-like cols: %s", paste(have, collapse=",")) else "no death-like cols in dt_unit", "3_survival")

if (!HAS_MORTALITY || length(have) == 0) {
  CH$log("WARN", "HAS_MORTALITY=FALSE / no death vars → building empty surv shell, Cox will be skipped", "3_survival")
  
  # 构造空壳：每人一行，death=NA，follow_up=NA
  dt_surv <- dt_unit %>%
    distinct(ID) %>%
    mutate(
      death = NA_integer_,
      death_year = NA_real_,
      age_entry = NA_real_,
      follow_up_years = NA_real_,
      max_wave = NA_integer_,
      max_year = NA_integer_
    )
  
  # 补 baseline 年份
  dt_base_y <- dt_unit %>%
    filter(wave_year == min(wave_year, na.rm = TRUE)) %>%
    distinct(ID, base_year = wave_year)
  dt_surv <- dt_surv %>% left_join(dt_base_y, by = "ID")
  
  save_ck("dt_surv", "3_survival")
  CH$report("3_survival")
  # 不 quit()，让主控继续
} else {
  # ---- 真有死亡数据时的逻辑（后续有数据再启用）----
  CH$log("OK", "Mortality data detected, building survival library...", "3_survival")
  
  dt_surv <- dt_unit %>% distinct(ID)
  
  if ("RwIWSTAT" %in% names(dt_unit)) {
    d <- dt_unit %>% group_by(ID) %>%
      summarise(
        die = max(ifelse(RwIWSTAT %in% c(5, 6), 1, 0), na.rm = TRUE),
        yr = min(ifelse(RwIWSTAT %in% c(5, 6), wave_year, Inf), na.rm = TRUE),
        .groups = "drop"
      )
    dt_surv <- dt_surv %>% left_join(d, by = "ID") %>%
      mutate(death = ifelse(is.na(die), 0, die),
             death_year = ifelse(is.finite(yr), yr, NA))
  }
  
  # 基线年份
  dt_base_y <- dt_unit %>%
    filter(wave_year == min(wave_year, na.rm = TRUE)) %>%
    distinct(ID, base_year = wave_year)
  dt_surv <- dt_surv %>% left_join(dt_base_y, by = "ID")
  
  # 各人最大 wave_year（用于删失）
  max_wave_info <- dt_unit %>% group_by(ID) %>%
    summarise(max_year = max(wave_year, na.rm = TRUE), .groups = "drop")
  dt_surv <- dt_surv %>% left_join(max_wave_info, by = "ID")
  
  dt_surv <- dt_surv %>%
    mutate(
      follow_up_years = case_when(
        death == 1 & is.finite(death_year) ~ death_year - base_year,
        !is.na(max_year) & !is.na(base_year) ~ max_year - base_year,
        TRUE ~ NA_real_
      ),
      follow_up_years = pmax(follow_up_years, 0, na.rm = TRUE)
    )
  
  CH$log("OK", sprintf("dt_surv: %d ids, %d deaths", nrow(dt_surv), sum(dt_surv$death, na.rm = TRUE)), "3_survival")
  save_ck("dt_surv", "3_survival")
  CH$report("3_survival")
}