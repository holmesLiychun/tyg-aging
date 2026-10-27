# ============================================================================
# 模块 01：配置（路径 / 周期 / 算法参数 / 调试开关 / 命名约定）
# ----------------------------------------------------------------------------
# 唯一存放“可变参数”的地方。改路径、改周期、改算法都只动本文件。
# ============================================================================

# ---- 1.1 路径设置 -----------------------------------------------------------
# 原始 XPT 文件根目录（按周期分子文件夹）
NHANES_ROOT <- "~:/NHANES/"

# 结果输出目录（含子目录：tables/、figures/、cache/、logs/）
OUT_DIR <- here::here()


# 确保输出目录与子目录存在
for (d in c(OUT_DIR,
            file.path(OUT_DIR, "tables"),
            file.path(OUT_DIR, "figures"),
            file.path(OUT_DIR, "cache"),
            file.path(OUT_DIR, "logs"))) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# 缓存 RDS 路径（模块 10 的断点续跑文件）
CACHE_RAW_DERIVED <- file.path(OUT_DIR, "data/example", "nhanes_raw_derived.rds")

# 日志文件路径（带时间戳，避免覆盖）
LOG_FILE <- file.path(
  OUT_DIR, "logs",
  sprintf("run_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

# ---- 1.2 周期与 CRP 可用性 --------------------------------------------------
# 全部 7 个周期（KDM-11 / Modified PhenoAge 主分析用）
CYCLES <- c("2005-2006", "2007-2008", "2009-2010",
            "2011-2012", "2013-2014",
            "2015-2016", "2017-2018")

# 有 CRP 的 5 个周期（PhenoAge-9 原版、KDM-9、中介分析用）
# 2011-2012 与 2013-2014 无 CRP，PhenoAge 原版无法计算
CRP_CYCLES <- c("2005-2006", "2007-2008", "2009-2010", "2015-2016", "2017-2018")

# ---- 1.3 算法 marker 定义 ---------------------------------------------------
# KDM-11：主分析算法（不含 CRP，可覆盖全部 7 周期）
KDM11_MARKERS <- c("albumin", "alp", "lncreat", "hba1c", "sbp", "bun", "uap",
                   "lymph", "mcv", "wbc", "totchol")

# Modified PhenoAge：剔除 CRP 后的 9-marker 版本（7 周期）
PHENO_NOCRP_MARKERS <- c("albumin", "lncreat", "glucose", "hba1c", "sbp", "bmi",
                        "totchol", "trig", "age")

# PhenoAge-9 原版：含 CRP（仅 5 周期）
PHENO9_MARKERS <- c("albumin", "lncreat", "glucose", "hba1c", "sbp", "bmi",
                    "totchol", "trig", "crp")

# KDM-9：含 CRP + 含 glucose（13 marker，仅 5 周期）
KDM9_MARKERS <- c("albumin", "alp", "lncreat", "hba1c", "sbp", "bun", "uap",
                  "lymph", "mcv", "wbc", "totchol", "glucose", "lncrp")

# KDM-noG：含 CRP、不含 glucose（即 KDM 完整版剔除 glucose，12 marker）
KDMNOG_MARKERS <- c("albumin", "alp", "lncreat", "hba1c", "sbp", "bun", "uap",
                    "lymph", "mcv", "wbc", "totchol", "lncrp")

# Pheno-noG-6：含 CRP、不含 glucose（6 marker，仅 5 周期）
PHENO_NOG6_MARKERS <- c("albumin", "lncreat", "hba1c", "sbp", "totchol", "lncrp")

# ---- 1.4 调试与运行开关 ------------------------------------------------------
# DEBUG      : TRUE 时打印完整诊断表、变量分布、相关性矩阵等细节
#              会显著增加运行时间与日志体积
# VERBOSE    : TRUE 时打印每步的样本量/耗时/变量缺失率
# CHECKPOINT : TRUE 时每个模块结束保存一次 RDS 检查点，便于断点续跑与回溯
# SKIP_IF_CACHED : TRUE 时若检测到缓存文件已存在则跳过对应模块（加速重跑）
DEBUG          <- TRUE
VERBOSE        <- TRUE
CHECKPOINT     <- TRUE
SKIP_IF_CACHED <- TRUE  #

# ---- 1.5 变量命名约定（约定优于配置）----------------------------------------
# 保持统一后缀，grep 就能定位一类变量：
#   *_raw   原始读取未处理
#   *_ba    生物年龄估计值（单位：岁）
#   *_accel 生物年龄加速残差（BA 对 age 回归的残差，Mean≈0）
#   *_z     标准化 Z-score（SD=1，用于跨算法 β 比较）
#   *_tert  TyG 三分位因子（T1 low / T2 mid / T3 high，T1 为参考）
#   *_winsor 截尾后变量（±3SD）
#   ln*      对数转换（ln(x+1)，避免 0 取对数）
#
# 算法相关变量一览：
#   kdm_ba / kdm_accel           KDM-11（主分析，7 周期，不含 CRP）
#   pheno_ba / pheno_accel       Modified PhenoAge（7 周期，不含 CRP）
#   pheno9_ba / pheno9_accel     PhenoAge-9 原版（5 周期，含 CRP，敏感性）
#   kdm9_ba / kdm9_accel         KDM-9（5 周期，含 CRP，敏感性）
# ============================================================================
# Toggle: "example" for testing without NHANES access, "real" for actual data

