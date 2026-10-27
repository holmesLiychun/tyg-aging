# ============================================================
# 00_run_all.R — 主控脚本
# ============================================================
if (!exists("CH")) source("setup.R")
install.packages("here")
steps <- list(
  "1_read"           = "Step1_read.R",
  "2_audit"          = "Step2_audit.R",
  "2b_table1"          = "Step2b_table1.R",
  "3_survival"       = "Step3_survival.R",
  "4_BA"             = "Step4_BA_mod.R",
  "4bto10"      = "Step4b_to_10_KDMcore6.R"   # <— 加这一行
)
# START_AT 可以设为 "10_export" 单独跑导出，不用重跑全部
START_AT <- "1_read"



started <- FALSE
for (step_name in names(steps)) {
  if (!started) {
    if (step_name == START_AT) started <- TRUE else next
  }
  CH$log("INFO", sprintf("========== Running %s ==========", step_name))
  source(steps[[step_name]])
  CH$report(step_name)
}
