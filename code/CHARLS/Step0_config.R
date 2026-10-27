# ============================================================
# Step0_config.R — 唯一需要修改的文件
# ============================================================
# ---- 路径 ----
PATH_DATA   <- "~/charls.dta"
PATH_OUT    <- here::here("data")
PATH_RDATA  <- here::here("data/example")

# ---- 波次年份映射 ----
WAVE_YEAR <- c("1" = 2011, "2" = 2013, "3" = 2015, "4" = 2018, "5" = 2020)
TIME_YEARS <- c("1" = 0, "2" = 2, "3" = 4, "4" = 7, "5" = 9)

# ---- 四版 BA 定义（v4 终审）----
BA_VERSIONS <- list(
  KDMf = list(label="KDM-full (10 markers, no ln_crp)", method="kdm_nhanes",
              markers=c("bl_crea","bl_glu","bl_hbalc","bl_cho","bl_tg","bl_bun","bl_ua","bl_wbc","bl_mcv","systo")),
  KDMnog = list(label="KDM-noG (8 markers)", method="kdm_nhanes",
                markers=c("bl_cho","bl_tg","bl_crea","bl_bun","bl_ua","bl_wbc","bl_mcv","systo")),
  Pheno6 = list(label="PhenoAge-mod-6", method="pheno_mod",
                markers=c("bl_crea","bl_glu","ln_crp","bl_mcv","bl_wbc","age")),
  Pheno5nog = list(label="PhenoAge-mod-noG-5", method="pheno_mod",
                   markers=c("bl_crea","ln_crp","bl_mcv","bl_wbc","age")),
  KDM_core6 = list(
    label = "KDM_core6 (6 markers, fixed NHANES coeff)", 
    method = "kdm_core6_fixed",
    markers = c("lncreat", "bl_bun", "bl_ua", "bl_wbc", "bl_mcv", "systo"),
    coef_path = "C:/Users/user/Desktop/胰岛素抵抗-衰老项目/10.07/output(new)/cache/kdm_core6_coef_nhanes.rds"
  )
)

# ---- 协变量池 ----
COVARS_BASE <- c("age_baseline","gender","bmi_baseline",
                 "smokev_baseline","drinkl_baseline",
                 "edu_baseline","exercise_baseline")

# ---- 决策阈值 ----
DECISION_POSITIVE_P <- 0.05   # KDM-full 横断面 β>0 且 p<此值 → go

HAS_MORTALITY <- FALSE   # 合并好死亡库后改 TRUE
PATH_DEATH     <- NULL    # 死亡/退出dta路径；NULL=不读

# ---- Debug 控制 ----
STOP_ON_FAIL <- FALSE   # TRUE = 任何 FAIL 立即停止；FALSE = 只 WARN
VERBOSE      <- TRUE    # 打印中间结果
# ---- 目录初始化（确保存在）----
for (d in c(PATH_OUT, PATH_RDATA, 
            paste0(PATH_OUT, "/tables"), 
            paste0(PATH_OUT, "/figures"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
cat(sprintf("✓ Directories ready (OUT=%s)\n", PATH_OUT))

# ============================================================
# ---- 日志系统初始化（CH 对象）----
# ============================================================
LOG_FILE <- file.path(PATH_OUT, sprintf("pipeline_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))

# 初始化 CH（如果之前没定义）
if (!exists("CH")) CH <- new.env(parent = emptyenv())

CH$log_file <- LOG_FILE

CH$log <- function(level, message, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- sprintf("[%s] %-7s | %s", timestamp, level, sprintf(message, ...))
  
  if (VERBOSE) {
    if (level == "ERROR") message(msg)
    else if (level == "WARN") warning(msg, call. = FALSE)
    else cat(msg, "\n")
  }
  
  cat(msg, "\n", file = CH$log_file, append = TRUE)
  
  # ★ 新增：ERROR 时自动 dump traceback 到日志
  if (level == "ERROR") {
    cat("\n--- TRACEBACK ---\n", file = CH$log_file, append = TRUE)
    traceback(max.lines = 10)
    # traceback 只打印到控制台，想写文件需要 sink
  }
  
  if (level == "ERROR" && STOP_ON_FAIL) {
    stop("Pipeline halted.")
  }
}

CH$report <- function(step_name, status = "SUCCESS") {
  CH$log("INFO", sprintf("Step '%s' %s", step_name, status))
}

# ---- 全局错误/警告捕获 ----
CH$log("INFO", "Pipeline initialized. Log file: %s", LOG_FILE)
CH$log("INFO", "R version: %s, Platform: %s", R.version.string, R.version$platform)