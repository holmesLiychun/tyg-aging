# ============================================================================
# 模块 02：通用工具函数
# ----------------------------------------------------------------------------
# 内容：日志、计时、调试断言、变量诊断、检查点存取、自定义运算符。
# 依赖：01_config.R（LOG_FILE / OUT_DIR / DEBUG / VERBOSE / CHECKPOINT）
# ============================================================================

# ---- 2.1 运算符说明 ----
# %||% / %nin% 已在 00_main.R §0.0 定义，此处不再重复定义，避免覆盖

# ---- 2.2 日志系统 -----------------------------------------------------------
# 日志同时输出到控制台与日志文件，避免 Rscript 模式下丢失过程信息
log_conn <- NULL

log_init <- function() {
  log_conn <<- file(LOG_FILE, open = "wt")
  on.exit(NULL)
  sink(log_conn, append = TRUE, split = TRUE)
  cat(sprintf("===== NHANES 分析日志 开始 %s =====\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

log_write <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

log_close <- function() {
  if (!is.null(log_conn) && isOpen(log_conn)) {
    cat(sprintf("===== 日志结束 %s =====\n", format(Sys.time(), "%H:%M:%S")))
    sink()
    close(log_conn)
    log_conn <<- NULL
  }
}

# 注册退出钩子，确保异常退出时日志文件也能正常关闭
reg.finalizer(environment(), function(e) log_close(), onexit = TRUE)

# ---- 2.3 计时器 -------------------------------------------------------------
TIMERS <- list()

timer_start <- function(name) {
  TIMERS[[name]] <<- Sys.time()
  if (VERBOSE) log_write(sprintf("▶ 开始: %s", name))
}

timer_end <- function(name) {
  if (!is.null(TIMERS[[name]])) {
    dt <- difftime(Sys.time(), TIMERS[[name]], units = "secs")
    log_write(sprintf("✓ 结束: %s (耗时 %.1f 秒)", name, as.numeric(dt)))
    TIMERS[[name]] <<- NULL
  }
}

# ---- 2.4 调试断言（只在 DEBUG=TRUE 时生效，否则是空操作）--------------------
# dbg_assert：条件不成立时打印警告（不中断，便于先跑通再看问题）
dbg_assert <- function(condition, msg) {
  if (!DEBUG) return(invisible(NULL))
  if (!isTRUE(condition)) {
    log_write(sprintf("  ⚠ 断言失败: %s", msg))
  }
  invisible(NULL)
}

# dbg_stop：条件不成立时直接报错中断（用于不可容忍的数据问题）
dbg_stop <- function(condition, msg) {
  if (!isTRUE(condition)) {
    stop(msg, call. = FALSE)
  }
}

# ---- 2.5 变量诊断（每个模块结束自动调用）------------------------------------
# 打印 n、缺失率、均值/SD，快速发现数据异常
diag_vars <- function(data, vars, title = "变量诊断") {
  if (!VERBOSE) return(invisible(NULL))
  if (!is.data.frame(data)) return(invisible(NULL))
  vars <- intersect(vars, names(data))
  if (length(vars) == 0) {
    log_write(sprintf("  [%s] 无可诊断变量", title))
    return(invisible(NULL))
  }

  log_write(sprintf("  [%s] n=%d", title, nrow(data)))
  for (v in vars) {
    x <- data[[v]]
    n_na <- sum(is.na(x))
    miss_pct <- round(100 * n_na / max(nrow(data), 1), 1)
    if (is.numeric(x) && n_na < nrow(data)) {
      rng <- sprintf("范围[%.1f, %.1f]", min(x, na.rm = TRUE), max(x, na.rm = TRUE))
      stats <- sprintf("mean=%.2f sd=%.2f %s", mean(x, na.rm = TRUE),
                       sd(x, na.rm = TRUE), rng)
    } else {
      stats <- sprintf("非数值型 / 全缺失 (levels=%d)",
                       if (is.factor(x)) nlevels(x) else 0)
    }
    log_write(sprintf("    %-22s NA=%d (%.1f%%) | %s", v, n_na, miss_pct, stats))
  }
  invisible(NULL)
}

# ---- 2.6 检查点存取（CHECKPOINT=TRUE 时每个模块结束保存一次 RDS）------------
# 命名规则：checkpoint_<模块编号>_<描述>.rds
cp_save <- function(obj, module, desc) {
  if (!CHECKPOINT) return(invisible(NULL))
  fname <- sprintf("checkpoint_%s_%s.rds", module, desc)
  fpath <- file.path(OUT_DIR, "cache", fname)
  tryCatch({
    saveRDS(obj, fpath)
    log_write(sprintf("  ☑ 检查点已保存: %s (%d 行 x %d 列)", fname, nrow(obj), ncol(obj)))
  }, error = function(e) {
    log_write(sprintf("  ⚠ 检查点保存失败: %s -> %s", fname, e$message))
  })
  invisible(NULL)
}

cp_load <- function(module, desc) {
  fname <- sprintf("checkpoint_%s_%s.rds", module, desc)
  fpath <- file.path(OUT_DIR, "cache", fname)
  if (file.exists(fpath)) {
    return(readRDS(fpath))
  }
  NULL
}

# ---- 2.7 输出文件清单（带行数/列数/变量名，便于核对产物）--------------------
log_file_list <- function() {
  files <- list.files(OUT_DIR, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir == FALSE]
  if (length(files) == 0) {
    log_write("  （输出目录为空）")
    return(invisible(NULL))
  }
  for (f in files) {
    rel <- gsub(OUT_DIR, "", f, fixed = TRUE)
    size_kb <- round(file.info(f)$size / 1024, 1)
    cat(sprintf("  %-55s %8.1f KB\n", rel, size_kb))
  }
  cat(sprintf("  共 %d 个文件\n", length(files)))
  invisible(NULL)
}

# ---- 2.8 包加载器（首次自动安装，避免原稿手动 install.packages 的重复）-----
load_or_install <- function(pkg) {
  for (p in pkg) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
      log_write(sprintf("  安装缺失包: %s", p))
      install.packages(p, repos = "https://cloud.r-project.org")
      library(p, character.only = TRUE)
    }
  }
}

# 首次加载时即安装并加载本项目所有依赖包
load_or_install(c("haven", "dplyr", "tidyr", "survey", "rms", "ggplot2",
                  "splines", "devtools", "mice", "tableone", "table1",
                  "flextable", "officer", "forestploter", "mediation",
                  "grid", "readr"))

# BioAge 从 GitHub 安装
if (!require("BioAge", quietly = TRUE)) {
  log_write("  安装缺失包: BioAge (GitHub)")
  devtools::install_github("dayoonkwon/BioAge")
  library(BioAge)
}

# 全局选项：lonely PSU 用 adjust 避免报错、抑制 dplyr 的 summarise 提示
options(survey.lonely.psu = "adjust",
        dplyr.summarise.inform = FALSE,
        scipen = 10)

log_write("✓ 模块 02 工具函数加载完成")
