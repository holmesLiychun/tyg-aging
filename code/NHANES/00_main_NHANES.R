# ============================================================================
# NHANES 2005-2018 (7 周期) 胰岛素抵抗(TyG) → 生物年龄加速 分析流程 v9
# ----------------------------------------------------------------------------
# 主题：主调度脚本（入口文件）
#
# v9 重构目标
#   1. 模块化拆分 —— 每个文件一个主题，消除原稿中函数重复定义、
#      模块顺序错乱（B 模块中间插入 B.5、B+ 散落）、library 散落各处的问题
#   2. 内置 debug 能力 —— 统一开关 DEBUG / VERBOSE / CHECKPOINT / STOP_ON_ERROR，
#      每模块有断点续跑（skip if cached）与诊断函数
#   3. 健壮的错误处理 —— 自定 %||% / %nin% 运算符、tryCatch 包裹风险步骤、
#      get_design 不再直接 stop，改为返回诊断信息
#   4. 全量打点 —— 每步记录耗时、未加权 n、加权 N、变量缺失率，
#      输出一份带时间戳的运行日志与文件清单（含行数/列数/变量名）
#   5. 变量命名规范 —— 见 01_config.R 顶部的命名约定表
#
# 运行方式
#   Rscript 00_main.R            # 控制台运行，跑全部模块
#   source("00_main.R")          # RStudio 交互运行（完成后对象留在环境）
#
# 控制选项见 01_config.R 中的 DEBUG / VERBOSE / CHECKPOINT / SKIP_IF_CACHED
# ============================================================================
install.packages("here")
# ---- 0. 清理与基础设置 -------------------------------------------------------
# rm(list = ls(all = TRUE))      # 交互式运行若想清空环境可取消注释

# ---- 0.0 运算符预定义（必须在最前，PROJECT_ROOT 那行就用到了 %||%）----------
# %||%：左值为 NULL/NA 时取右值（Elvis 运算符）
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b
# %nin%：不在集合中（与 %in% 相反）
`%nin%` <- function(x, table) !match(x, table, nomatch = 0L)

# 记录整体开始时间（全程计时用）
t0_global <- Sys.time()

# 定位本项目根目录：本文件所在目录
PROJECT_ROOT <- here::here()

# 加载配置：路径、周期、算法参数、调试开关、命名约定
source(file.path(PROJECT_ROOT, "code/NHANES", "01_config.R"))

# 初始化日志文件与计时器
source(file.path(PROJECT_ROOT, "code/NHANES", "02_utils.R"))
log_init()
timer_start("全分析流程")

# 控制台欢迎信息
cat(sprintf("\n========== NHANES 分析流程 v9 ==========\n"))
cat(sprintf("项目根目录 : %s\n", PROJECT_ROOT))
cat(sprintf("输出目录   : %s\n", OUT_DIR))
cat(sprintf("调试开关   : DEBUG=%s | VERBOSE=%s | CHECKPOINT=%s | SKIP_IF_CACHED=%s\n",
            DEBUG, VERBOSE, CHECKPOINT, SKIP_IF_CACHED))
cat(sprintf("运行开始   : %s\n", format(t0_global, "%Y-%m-%d %H:%M:%S")))
cat(sprintf("周期数     : %d (%s)\n\n", length(CYCLES), paste(CYCLES, collapse = ", ")))

# ---- 1. 数据导入与衍生变量 ---------------------------------------------------
# 模块 10：读取 7 周期 XPT，标准化为统一变量名，计算 TyG / HOMA-IR / 糖尿病
#           状态等衍生变量。结果缓存为 nhanes_raw_derived.rds，断点续跑。
source(file.path(PROJECT_ROOT, "code/NHANES", "10_import.R"))
# 导出对象：nhanes（全周期合并后的主数据框）

# ---- 2. 生物年龄算法计算 ----------------------------------------------------
# 模块 20：KDM-11（主分析，7 周期，不含 CRP）、Modified PhenoAge（7 周期，
#           不含 CRP）、PhenoAge-9 原版（5 周期，含 CRP）、KDM-9（5 周期，
#           含 CRP）。按性别算加速残差。
source(file.path(PROJECT_ROOT,"code/NHANES", "20_bioage.R"))
source(file.path(PROJECT_ROOT, "code/NHANES","25_phenoage.R")) #6种PhenoAge
# 导出对象：nhanes（新增 kdm_ba/kdm_accel、pheno_ba/pheno_accel、
#           pheno9_ba/pheno9_accel、kdm9_ba/kdm9_accel）

# ---- 3. 分析人群筛选与流程图 ------------------------------------------------
# 模块 30：按 成人≥20 → MEC合格 → 空腹合格 → TyG非缺失 → marker完整 的
#           步骤统计各阶段未加权/加权样本量，输出 STROBE 式流程图数据。
source(file.path(PROJECT_ROOT,"code/NHANES", "30_flowchart.R"))
# 导出对象：ana_kdm、ana_pheno、ana_pheno9、ana_kdm9、flow_data

# ---- 4. 描述性统计与 Table 1 -----------------------------------------------
# 模块 40：加权分布、TyG 三分位分布、各算法 BA 分布（按性别）、CRP 分布、
#           Table 1 基线特征（加权，含 SMD）。
source(file.path(PROJECT_ROOT,"code/NHANES", "40_descriptive.R"))

# ---- 5. 主回归分析 ----------------------------------------------------------
# 模块 50：TyG 三分位对 4 个生物年龄加速残差的加权线性回归（Model 1-5），
#           提取 β/SE/P 值。
source(file.path(PROJECT_ROOT, "code/NHANES","50_regression.R"))
# 导出对象：main_effects（4 算法 × 5 模型的结果长表）

# ---- 6. 敏感性分析 ----------------------------------------------------------
# 模块 60：截尾（±3SD）、种族分层与交互、4 套算法 TyG 关联、解耦矩阵（4
#           暴露 × 4 结局）、中介分析（TyG → 炎症 → 生物年龄加速）。
source(file.path(PROJECT_ROOT,"code/NHANES", "60_sensitivity.R"))
source(file.path(PROJECT_ROOT, "code/NHANES","65_interaction.R")) #交互项分析
source(file.path(PROJECT_ROOT, "code/NHANES","65b_stratified_analysis.R")) #分层分析
source(file.path(PROJECT_ROOT, "code/NHANES","65c_stratified_figures.R")) #分层分析图
source(file.path(PROJECT_ROOT, "code/NHANES","66_nri_idi.R")) #nri、idi分析

# ---- 7. 非线性与可视化 ------------------------------------------------------
# 模块 70：RCS 受限立方样条（主分析 7 周期 + 敏感性 5 周期）、森林图、
#           4×4 解耦矩阵热图、β 轨迹图。
source(file.path(PROJECT_ROOT,"code/NHANES", "70_visualization.R"))

# ---- 8. 收尾：保存对象 + 输出文件清单 + 结束日志 -----------------------------
timer_end("全分析流程")
t1_global <- Sys.time()
elapsed <- difftime(t1_global, t0_global, units = "mins")

log_write(sprintf("\n========== 全部分析完成 =========="))
log_write(sprintf("主分析  : KDM-11 (不含 CRP) → %d 行, 7 周期", nrow(ana_kdm)))
if (exists("ana_pheno9")) {
  log_write(sprintf("敏感性  : PhenoAge-9 原版 (含 CRP) → %d 行, 5 周期", nrow(ana_pheno9)))
}
log_write(sprintf("输出目录 : %s", OUT_DIR))
log_write(sprintf("运行耗时 : %.1f 分钟", as.numeric(elapsed)))
log_write(sprintf("结束时间 : %s", format(t1_global, "%Y-%m-%d %H:%M:%S")))

# 控制台同步输出一份文件清单（带行数/列数/关键变量名）
cat("\n========== 输出文件清单 ==========\n")
log_file_list()

cat(sprintf("\n=== 流程结束，日志写入 %s ===\n", LOG_FILE))
