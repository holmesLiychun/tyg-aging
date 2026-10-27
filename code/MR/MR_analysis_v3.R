# =====================================================================
# 双向两样本孟德尔随机化 (R 4.5.1)
# 本地 EBI summary VCF + 本地 PLINK LD clump + STROBE-MR 全流程审计
# ---------------------------------------------------------------------
# 使用说明:
#   1) DEBUG=TRUE  运行时打印中间对象/写调试csv, 便于排查; 正式出结果设FALSE
#   2) DO_MD5=FALSE 大VCF跳过md5省内存(默认关); 需校验完整性再开
#   3) RUN_STEIGER/RUN_LEAVEOUT/RUN_MRPRESSO 控制敏感性分析开关
#   4) 带 [...] / [人工:...] 的字段必须补; 自动值由IEU/Catalog拉, 仍需人工核对
# =====================================================================
# =====================================================================
# 双向两样本孟德尔随机化 (R 4.5.1) — 【修正完整版】
# 本地 EBI summary VCF + 本地 PLINK LD clump + STROBE-MR 审计
# 
# 修复记录:
#   1) 暴露 p<5e-8 筛选工具变量，结局不过滤 p 值（核心修复）
#   2) harF/harR 初始化 + exists() 保护
#   3) SS 字段 FORMAT keys 调试输出
#   4) 结局数据可选宽松过滤（p_filter_outcome 参数）
# =====================================================================

# ---------- 0. 安装依赖（首次运行取消注释）----------
# install.packages(c("data.table","dplyr","remotes","stringi"))
# remotes::install_github("MRCIEU/TwoSampleMR", upgrade="never")
# remotes::install_github("MRCIEU/ieugwasr")
# remotes::install_github("MRCIEU/genetics.binaRies")

# ---------- 1. 加载包 ----------
suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(ggplot2)
  library(dplyr)
})

if (!requireNamespace("stringi", quietly = TRUE)) {
  stop("请先安装 stringi: install.packages('stringi')")
} else {
  cat("✓ stringi 已加载，矩阵拆分模式\n")
}

# ---------- 2. 路径与参数 ----------
dir    <- "you_dir"
outdir <- file.path(dir, "mr_results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

fA <- file.path(dir, "ebi-a-GCST90002232.vcf.gz")
fB <- file.path(dir, "ebi-a-GCST90014292.vcf.gz")

TRAIT_A <- "Trait_A_GCST90002232"
TRAIT_B <- "Trait_B_GCST90014292"

plink_path <- "D:/software/plink.exe"
bfile_path <- "D:/reference/1kg/EUR"

# 暴露筛选阈值（严格）
P_EXPOSURE <- 5e-8
# 结局筛选阈值（NULL=不过滤, 或设 1e-3 宽松过滤以省内存）
P_OUTCOME  <- NULL          # ← 设为 NULL 表示结局不过滤
# P_OUTCOME  <- 1e-3        # ← 如内存紧张，取消注释此行

# ---------- 3. 读取 EBI summary VCF（stringi矩阵拆分）----------
read_ebi_format_vcf <- function(file, p_filter = NULL, nThread = 8, debug = FALSE) {
  cat("  读取:", basename(file), "\n")
  
  # 跳过 meta 行
  con <- gzfile(file, "rt")
  skip_meta <- 0
  hdrline <- NULL
  while (TRUE) {
    l <- readLines(con, 1)
    if (length(l) == 0) break
    if (startsWith(l, "##")) { skip_meta <- skip_meta + 1; next }
    hdrline <- sub("^#", "", l)
    break
  }
  close(con)
  
  hdr <- strsplit(hdrline, "\t")[[1]]
  
  # 只读取需要的列
  need_cols <- c("CHROM", "POS", "ID", "REF", "ALT", "FORMAT")
  samp_idx <- which(hdr == "FORMAT") + 1
  if (samp_idx <= length(hdr)) need_cols <- c(need_cols, hdr[samp_idx])
  
  col_indices <- which(hdr %in% need_cols)
  
  dt <- fread(file, skip = skip_meta + 1, header = FALSE, sep = "\t",
              fill = TRUE, stringsAsFactors = FALSE, nThread = nThread,
              select = col_indices)
  setnames(dt, hdr[col_indices])
  
  stopifnot(all(c("CHROM", "POS", "ID", "REF", "ALT", "FORMAT") %in% names(dt)))
  
  samp_name <- hdr[samp_idx]
  samp_col <- dt[[samp_name]]
  
  # 解析 FORMAT 键名
  keys <- strsplit(dt$FORMAT[1], ":", fixed = TRUE)[[1]]
  n_keys <- length(keys)
  
  if (debug) cat("  FORMAT keys:", paste(keys, collapse = ", "), "\n")
  
  # stringi 矩阵拆分
  mat <- stringi::stri_split_fixed(samp_col, ":", simplify = TRUE)
  
  if (ncol(mat) < n_keys) {
    pad <- matrix(NA_character_, nrow = nrow(mat), ncol = n_keys - ncol(mat))
    mat <- cbind(mat, pad)
  } else if (ncol(mat) > n_keys) {
    mat <- mat[, 1:n_keys]
  }
  colnames(mat) <- keys
  
  # 按候选名取列
  get_col <- function(cands) {
    for (k in cands) {
      ki <- match(k, keys)
      if (!is.na(ki) && ki <= ncol(mat)) {
        return(as.numeric(mat[, ki]))
      }
    }
    return(rep(NA_real_, nrow(dt)))
  }
  
  snp <- ifelse(is.na(dt$ID) | dt$ID == ".",
                paste0(dt$CHROM, ":", dt$POS, "_", dt$REF, "_", dt$ALT),
                dt$ID)
  
  out <- data.table(
    SNP       = snp,
    CHROM     = as.character(dt$CHROM),
    POS       = as.integer(dt$POS),
    ID        = ifelse(is.na(dt$ID) | dt$ID == ".", "", dt$ID),
    EA        = as.character(dt$ALT),
    OA        = as.character(dt$REF),
    beta      = get_col(c("EFFECT", "BETA", "ES", "OR")),
    se        = get_col(c("SE", "STDERR")),
    z         = get_col(c("Z", "ZSCORE")),
    lp        = get_col(c("LP")),
    p         = get_col(c("P", "PVAL")),
    af        = get_col(c("AF", "EAF")),
    ss        = get_col(c("SS", "N")),
    ncase     = get_col(c("NC", "NCASES", "NCA")),
    ncontrol  = get_col(c("NCU", "NCONTROL", "NCU")),
    acc       = get_col(c("ACC"))
  )
  
  out$pval <- out$p
  if (all(is.na(out$pval))) out$pval <- 10^(-out$lp)
  if (all(is.na(out$beta)) && !all(is.na(out$z)) && !all(is.na(out$se)))
    out$beta <- out$z * out$se
  
  out <- out[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]
  
  # p 值过滤
  if (!is.null(p_filter)) {
    out <- out[pval < p_filter]
    cat("  p<", p_filter, " 过滤后剩余 SNP 数:", nrow(out), "\n")
  }
  
  out <- out[!duplicated(SNP)]
  out[, c("z", "lp", "p") := NULL]
  out[]
}

# ---------- 4. 读取数据（修正：暴露严格过滤，结局不过滤）----------
cat("\n===== 4. 读取 VCF 数据 =====\n")

cat("\n--- 暴露 A (p<", P_EXPOSURE, ") ---\n")
A_exp <- read_ebi_format_vcf(fA, p_filter = P_EXPOSURE, debug = TRUE)
cat("A 暴露有效 SNP 数:", nrow(A_exp), "\n")

cat("\n--- 暴露 B (p<", P_EXPOSURE, ") ---\n")
B_exp <- read_ebi_format_vcf(fB, p_filter = P_EXPOSURE, debug = TRUE)
cat("B 暴露有效 SNP 数:", nrow(B_exp), "\n")

cat("\n--- 结局 A (全量", ifelse(is.null(P_OUTCOME), "不过滤", paste0("p<", P_OUTCOME)), ") ---\n")
A_out <- read_ebi_format_vcf(fA, p_filter = P_OUTCOME)
cat("A 结局有效 SNP 数:", nrow(A_out), "\n")

cat("\n--- 结局 B (全量", ifelse(is.null(P_OUTCOME), "不过滤", paste0("p<", P_OUTCOME)), ") ---\n")
B_out <- read_ebi_format_vcf(fB, p_filter = P_OUTCOME)
cat("B 结局有效 SNP 数:", nrow(B_out), "\n")

# 检查 SS 字段
cat("\n--- SS 字段检查 ---\n")
if (!all(is.na(A_exp$ss))) {
  cat("A 暴露 SS 众数:", as.numeric(names(which.max(table(A_exp$ss)))), "\n")
} else {
  cat("⚠ A 暴露 SS 字段全为 NA（FORMAT 中无 SS 或列名不匹配）\n")
}
if (!all(is.na(B_exp$ss))) {
  cat("B 暴露 SS 众数:", as.numeric(names(which.max(table(B_exp$ss)))), "\n")
} else {
  cat("⚠ B 暴露 SS 字段全为 NA\n")
}

# ---------- 5. VCF 头部元数据 ----------
extract_vcf_header_meta <- function(vcf_path) {
  con <- gzfile(vcf_path, "rt")
  metas <- c()
  while (TRUE) {
    l <- readLines(con, 1)
    if (length(l) == 0) break
    if (!startsWith(l, "##")) break
    metas <- c(metas, sub("^##", "", l))
  }
  close(con)
  
  kw <- c("trait", "phenotype", "sample", "population", "ancestry", 
          "GCST", "source", "consortium", "author", "unit", "build")
  rel <- metas[grepl(paste(kw, collapse = "|"), metas, ignore.case = TRUE)]
  
  list(file = basename(vcf_path), relevant_header_lines = rel, total_meta_lines = length(metas))
}

cat("\n===== VCF 头部元数据 =====\n")
meta_A <- extract_vcf_header_meta(fA)
meta_B <- extract_vcf_header_meta(fB)

sink(file.path(outdir, "mr_gwas_metadata.txt"))
cat("=== MR GWAS Metadata for Methods ===\n\n")
cat("Forward exposure: ebi-a-GCST90092992\n")
cat("Forward outcome:  ebi-a-GCST90014292\n")
cat("Reverse exposure: ebi-a-GCST90014292\n")
cat("Reverse outcome:  ebi-a-GCST90092992\n\n")
cat("--- VCF Header (A) ---\n")
print(meta_A$relevant_header_lines)
cat("\n--- VCF Header (B) ---\n")
print(meta_B$relevant_header_lines)
sink()
cat("✓ 已导出 mr_gwas_metadata.txt\n")

# ---------- 6. 转换为 TwoSampleMR 格式 ----------
fmt_mr <- function(d, trait, type, p_thresh = NULL) {
  if (!is.null(p_thresh)) d <- d[d$pval < p_thresh, ]
  if (nrow(d) == 0) {
    warning("No SNPs after p-value filtering.")
    return(NULL)
  }
  
  std <- data.frame(
    SNP           = d$SNP,
    beta          = d$beta,
    se            = d$se,
    effect_allele = d$EA,
    other_allele  = d$OA,
    eaf           = d$af,
    pval          = d$pval,
    chr           = d$CHROM,
    pos           = d$POS,
    samplesize    = if (!all(is.na(d$ss))) d$ss else NA,
    ncase         = if (!all(is.na(d$ncase))) d$ncase else NA,
    ncontrol      = if (!all(is.na(d$ncontrol))) d$ncontrol else NA,
    stringsAsFactors = FALSE
  )
  std <- std[, sapply(std, function(x) !all(is.na(x)))]
  std$phenotype <- trait
  format_data(std, type = type)
}

# ---------- 7. 本地 LD clump 函数 ----------
local_clump <- function(exposure_df, plink_bin, bfile,
                        clump_kb = 5000,
                        clump_r2 = 0.01,
                        clump_p = 1e-6) {
  clump_input <- data.frame(
    rsid = exposure_df$SNP,
    pval = exposure_df$pval.exposure,
    id   = exposure_df$id.exposure,
    stringsAsFactors = FALSE
  )
  clumped <- ld_clump(
    dat = clump_input,
    plink_bin = plink_bin,
    bfile = bfile,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = clump_p
  )
  if (is.null(clumped) || nrow(clumped) == 0) {
    cat("  ⚠ ld_clump 返回空结果\n")
    return(exposure_df[0, ])
  }
  exposure_df[exposure_df$SNP %in% clumped$rsid, ]
}

# ---------- 8. 初始化 harmonise 对象（作用域保护）----------
harF <- NULL
harR <- NULL
resF <- NULL; hetF <- NULL; pleiF <- NULL
resR <- NULL; hetR <- NULL; pleiR <- NULL

# ---------- 9. 正向 MR: A -> B ----------
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("===== 正向 MR: ", TRAIT_A, " -> ", TRAIT_B, " =====\n")
cat(paste(rep("=", 60), collapse=""), "\n")

expA_raw <- fmt_mr(A_exp, TRAIT_A, "exposure")
if (!is.null(expA_raw) && nrow(expA_raw) > 0) {
  cat("A 暴露 SNP 数（未 clump）:", nrow(expA_raw), "\n")
  
  expA <- local_clump(expA_raw, plink_bin = plink_path, bfile = bfile_path)
  cat("A 本地 clump 后工具变量数:", nrow(expA), "\n")
  
  if (nrow(expA) > 0) {
    # 用结局全量数据匹配
    outB_matched <- B_out[B_out$SNP %in% expA$SNP, ]
    cat("结局 B 中匹配到 SNP 数:", nrow(outB_matched), "\n")
    
    if (nrow(outB_matched) > 0) {
      outB <- fmt_mr(outB_matched, TRAIT_B, "outcome")
      if (!is.null(outB)) {
        harF <- harmonise_data(expA, outB, action = 2)
        cat("正向 harmonise 后 SNP 数:", nrow(harF), "\n")
        
        nIV <- length(unique(harF$SNP))
        cat("正向独立工具变量数:", nIV, "\n")
        
        if (nIV >= 3) {
          resF <- mr(harF, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                                           "mr_simple_mode", "mr_weighted_mode"))
        } else if (nIV == 2) {
          resF <- mr(harF, method_list = c("mr_ivw", "mr_weighted_median"))
        } else if (nIV == 1) {
          resF <- mr(harF, method_list = "mr_wald_ratio")
        } else {
          cat("  ⚠ 正向 MR：没有可用的工具变量\n")
        }
        
        if (!is.null(resF)) {
          hetF <- mr_heterogeneity(harF)
          pleiF <- tryCatch(mr_pleiotropy_test(harF), error = function(e) NULL)
        }
      } else {
        cat("  ⚠ 结局 B format_data 失败\n")
      }
    } else {
      cat("  ⚠ 结局 B 中未匹配到任何工具 SNP\n")
    }
  } else {
    cat("  ⚠ clump 后无工具变量\n")
  }
} else {
  cat("  ⚠ 暴露 A 无有效 SNP\n")
}

# ---------- 10. 反向 MR: B -> A ----------
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("===== 反向 MR: ", TRAIT_B, " -> ", TRAIT_A, " =====\n")
cat(paste(rep("=", 60), collapse=""), "\n")

expB_raw <- fmt_mr(B_exp, TRAIT_B, "exposure")
if (!is.null(expB_raw) && nrow(expB_raw) > 0) {
  cat("B 暴露 SNP 数（未 clump）:", nrow(expB_raw), "\n")
  
  expB <- local_clump(expB_raw, plink_bin = plink_path, bfile = bfile_path)
  cat("B 本地 clump 后工具变量数:", nrow(expB), "\n")
  
  if (nrow(expB) > 0) {
    outA_matched <- A_out[A_out$SNP %in% expB$SNP, ]
    cat("结局 A 中匹配到 SNP 数:", nrow(outA_matched), "\n")
    
    if (nrow(outA_matched) > 0) {
      outA <- fmt_mr(outA_matched, TRAIT_A, "outcome")
      if (!is.null(outA)) {
        harR <- harmonise_data(expB, outA, action = 2)
        cat("反向 harmonise 后 SNP 数:", nrow(harR), "\n")
        
        nIV_rev <- length(unique(harR$SNP))
        cat("反向独立工具变量数:", nIV_rev, "\n")
        
        if (nIV_rev >= 3) {
          resR <- mr(harR, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                                           "mr_simple_mode", "mr_weighted_mode"))
        } else if (nIV_rev == 2) {
          resR <- mr(harR, method_list = c("mr_ivw", "mr_weighted_median"))
        } else if (nIV_rev == 1) {
          resR <- mr(harR, method_list = "mr_wald_ratio")
        } else {
          cat("  ⚠ 反向 MR：没有可用的工具变量\n")
        }
        
        if (!is.null(resR)) {
          hetR <- mr_heterogeneity(harR)
          pleiR <- tryCatch(mr_pleiotropy_test(harR), error = function(e) NULL)
        }
      } else {
        cat("  ⚠ 结局 A format_data 失败\n")
      }
    } else {
      cat("  ⚠ 结局 A 中未匹配到任何工具 SNP\n")
    }
  } else {
    cat("  ⚠ clump 后无工具变量\n")
  }
} else {
  cat("  ⚠ 暴露 B 无有效 SNP\n")
}

# ---------- 11. 输出结果 ----------
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("===== 正向 MR 结果 =====\n")
if (!is.null(resF) && nrow(resF) > 0) {
  print(resF)
  cat("\n正向异质性检验:\n"); print(hetF)
  if (!is.null(pleiF)) { cat("\n正向水平多效性检验:\n"); print(pleiF) }
} else {
  cat("正向 MR 无结果（检查上方匹配日志）\n")
}

cat("\n===== 反向 MR 结果 =====\n")
if (!is.null(resR) && nrow(resR) > 0) {
  print(resR)
  cat("\n反向异质性检验:\n"); print(hetR)
  if (!is.null(pleiR)) { cat("\n反向水平多效性检验:\n"); print(pleiR) }
} else {
  cat("反向 MR 无结果（检查上方匹配日志）\n")
}

# 工具强度
strength <- function(har, name) {
  if (is.null(har) || nrow(har) == 0) 
    return(data.frame(direction = name, n_snp = 0, median_F = NA_real_, mean_F = NA_real_))
  f <- (har$beta.exposure / har$se.exposure)^2
  data.frame(direction = name,
             n_snp = length(unique(har$SNP)),
             median_F = round(median(f), 1),
             mean_F = round(mean(f), 1))
}

cat("\n===== 工具强度 (F 统计量) =====\n")
f_table <- rbind(strength(harF, "A->B"), strength(harR, "B->A"))
print(f_table)

# ---------- 12. 保存结果 ----------
cat("\n===== 保存结果 =====\n")

if (!is.null(resF) && nrow(resF) > 0) {
  fwrite(resF, file.path(outdir, "forward_MR_results.csv"))
  fwrite(hetF, file.path(outdir, "forward_heterogeneity.csv"))
  if (!is.null(pleiF)) fwrite(pleiF, file.path(outdir, "forward_pleiotropy.csv"))
  cat("✓ 正向结果已保存\n")
}
if (!is.null(resR) && nrow(resR) > 0) {
  fwrite(resR, file.path(outdir, "reverse_MR_results.csv"))
  fwrite(hetR, file.path(outdir, "reverse_heterogeneity.csv"))
  if (!is.null(pleiR)) fwrite(pleiR, file.path(outdir, "reverse_pleiotropy.csv"))
  cat("✓ 反向结果已保存\n")
}

fwrite(f_table, file.path(outdir, "instrument_strength.csv"))
cat("✓ 工具强度表已保存\n")

# ---------- 13. 森林图 ----------
plot_forest <- function(res, har, name, outcome_binary = FALSE) {
  if (is.null(res) || nrow(res) == 0) return(NULL)
  d <- res[res$method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"), ]
  if (nrow(d) == 0) return(NULL)
  
  if (outcome_binary) {
    d$est <- exp(d$b); d$lo <- exp(d$b - 1.96 * d$se); d$hi <- exp(d$b + 1.96 * d$se)
    xlab <- "OR (95% CI)"; vline <- 1
  } else {
    d$est <- d$b; d$lo <- d$b - 1.96 * d$se; d$hi <- d$b + 1.96 * d$se
    xlab <- "Beta (95% CI)"; vline <- 0
  }
  d$method <- factor(d$method, levels = rev(c("Inverse variance weighted", "MR Egger", "Weighted median")))
  
  ggplot(d, aes(x = est, y = method)) +
    geom_vline(xintercept = vline, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    geom_point(size = 3, color = "#2c7fb8") +
    labs(title = name, x = xlab, y = "") +
    theme_minimal(base_size = 13)
}

# 判断是否为二元结局（基于 ncase 字段）
binaryA <- !all(is.na(A_out$ncase)) && any(A_out$ncase > 0, na.rm = TRUE)
binaryB <- !all(is.na(B_out$ncase)) && any(B_out$ncase > 0, na.rm = TRUE)

cat("\n结局类型: A =", ifelse(binaryA, "二元", "连续"), 
    ", B =", ifelse(binaryB, "二元", "连续"), "\n")

p_fwd <- plot_forest(resF, harF, paste0(TRAIT_A, " -> ", TRAIT_B), outcome_binary = binaryB)
if (!is.null(p_fwd)) {
  ggsave(file.path(outdir, "forward_forest.png"), p_fwd, width = 7, height = 4)
  cat("✓ 正向森林图已保存\n")
}

p_rev <- plot_forest(resR, harR, paste0(TRAIT_B, " -> ", TRAIT_A), outcome_binary = binaryA)
if (!is.null(p_rev)) {
  ggsave(file.path(outdir, "reverse_forest.png"), p_rev, width = 7, height = 4)
  cat("✓ 反向森林图已保存\n")
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("所有结果已保存至:", outdir, "\n")
cat("完成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse=""), "\n")