# ============================================================
# R/setup.R — 包管理 + CH 检查点系统
# ============================================================
CH <- list(
  checks = list(),
  reset = function() { CH$checks <<- list(); invisible(NULL) },
  log = function(level, msg, step = NA) {
    CH$checks <<- append(CH$checks, list(list(
      level = level, msg = msg, step = step,
      time = format(Sys.time(), "%H:%M:%S")
    )))
    if (level == "FAIL") cat(sprintf("  ✗ FAIL: %s\n", msg))
    else if (level == "WARN") cat(sprintf("  ⚠ WARN: %s\n", msg))
    else if (exists("VERBOSE") && VERBOSE) cat(sprintf("  ✓ %s\n", msg))
  },
  report = function(step_name = NA) {
    n_fail <- sum(sapply(CH$checks, function(x) x$level == "FAIL"))
    n_warn <- sum(sapply(CH$checks, function(x) x$level == "WARN"))
    n_ok   <- sum(sapply(CH$checks, function(x) x$level == "OK"))
    cat(sprintf("\n--- CH Report [%s]: %d OK, %d WARN, %d FAIL ---\n",
                ifelse(is.na(step_name), "overall", step_name),
                n_ok, n_warn, n_fail))
    if (n_fail > 0) {
      cat("FAIL details:\n")
      for (c in CH$checks) if (c$level == "FAIL") cat(sprintf("  %s\n", c$msg))
      if (exists("STOP_ON_FAIL") && STOP_ON_FAIL) stop("STOP_ON_FAIL=TRUE, aborting.")
    }
    invisible(CH$checks)
  }
)

# 包管理
pkg <- c("haven", "dplyr", "tidyr", "lme4", "lmerTest", "survival",
         "gee", "rms", "ggplot2", "splines", "boot", "ipw", "lubridate")
for (p in pkg) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
    library(p, character.only = TRUE)
  }
}

# 目录
if (exists("PATH_OUT")) {
  for (d in c(PATH_OUT, paste0(PATH_OUT, "/checkpoints"),
              paste0(PATH_OUT, "/tables"), paste0(PATH_OUT, "/figures"))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
}

# 检查点 I/O
save_ck <- function(obj_name, step) {
  if (!exists("PATH_RDATA")) return()
  dir.create(PATH_RDATA, recursive = TRUE, showWarnings = FALSE)
  fname <- sprintf("%s/ck_%s_%s.rds", PATH_RDATA, step, obj_name)
  saveRDS(get(obj_name), fname)
  CH$log("OK", sprintf("Saved checkpoint: %s", basename(fname)), step)
}
load_ck <- function(obj_name, step) {
  if (!exists("PATH_RDATA")) return(FALSE)
  dir.create(PATH_RDATA, recursive = TRUE, showWarnings = FALSE)
  fname <- sprintf("%s/ck_%s_%s.rds", PATH_RDATA, step, obj_name)
  if (file.exists(fname)) {
    assign(obj_name, readRDS(fname), envir = .GlobalEnv)
    CH$log("OK", sprintf("Loaded checkpoint: %s", basename(fname)), step)
    return(TRUE)
  }
  return(FALSE)
}

cat("✓ CH system ready\n")

# ---- 覆盖 CH$log：彻底修复 % 字符导致的 sprintf 警告 ----
CH$log <- function(level, message, ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- paste0("[", ts, "] ", level, " | ", as.character(message), "\n")
  cat(log_line)
}