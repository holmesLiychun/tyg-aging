#!/usr/bin/env Rscript
# ============================================================
# MVMR: Fasting Glucose (FG) + Triglycerides (TG) -> PhenoAge
# IEU OpenGWAS VCF + local plink + 1KG EUR
# 修复: VCF FORMAT解析(非INFO), 正确S4调用, 等位基因对齐
# ============================================================

## ========== 0. 路径（已纠正） ==========
PLINK   <- "you_dir/plink.exe"
KG_EUR  <- "you_dir/EUR"     # 需有 EUR.bed/.bim/.fam
OUT_DIR <- "you_dir"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# FG = Fasting glucose; TG = Total triglycerides
GWAS_FG       <- "D:/gwas/vcf/ebi-a-GCST90002232.vcf.gz"   # ← FG
GWAS_TG       <- "D:/gwas/vcf/ebi-a-GCST90092992.vcf.gz"   # ← TG
GWAS_PHENOAGE <- "D:/gwas/vcf/ebi-a-GCST90014292.vcf.gz"   # Outcome

CLUMP_P  <- 5e-8
R2_CLUMP <- 0.001
KB_CLUMP <- 10000

cat("==================================================================\n")
cat("  FG   = GCST90002232 (Fasting glucose)\n")
cat("  TG   = GCST90092992 (Total triglycerides)\n")
cat("  Out  = GCST90014292 (PhenoAge)\n")
cat("==================================================================\n\n")

## ========== 1. 包 ==========
pkgs <- c("data.table", "MendelianRandomization")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(p, character.only = TRUE)
}

## ========== 2. VCF 读取（手动拆 FORMAT 字段，IEU标准） ==========
find_chrom_line <- function(vcf) {
  con <- gzfile(vcf, "r")
  i <- 0; hl <- -1
  while (TRUE) {
    l <- readLines(con, 1)
    if (!length(l)) break
    i <- i + 1
    if (startsWith(l, "#CHROM")) { hl <- i; break }
  }
  close(con)
  hl
}

read_ieu_vcf <- function(vcf, trait) {
  cat(sprintf("[%s] Reading %s ...\n", Sys.time(), basename(vcf)))
  
  hl <- find_chrom_line(vcf)
  if (hl < 0) stop(sprintf("[%s] No #CHROM line found", trait))
  
  d <- fread(vcf, sep = "\t", header = TRUE, skip = hl - 1, fill = TRUE)
  setnames(d, tolower(names(d)))
  if ("#chrom" %in% names(d)) setnames(d, "#chrom", "chr")
  if ("chrom" %in% names(d)) setnames(d, "chrom", "chr")
  
  # 找样本列（FORMAT之后的第一列）
  fixed_cols <- c("chr","pos","id","ref","alt","qual","filter","info","format")
  samp_cols <- setdiff(names(d), fixed_cols)
  if (length(samp_cols) == 0) stop(sprintf("[%s] No sample column found", trait))
  sc <- samp_cols[1]  # 取第一个样本列
  
  cat(sprintf("  Sample column: %s\n", sc))
  
  # 解析 FORMAT 字段，提取每个SNP的key顺序
  fmt_keys <- strsplit(as.character(d$format[1]), ":")[[1]]
  cat(sprintf("  FORMAT keys: %s\n", paste(fmt_keys, collapse=", ")))
  
  # 提取函数
  extract_field <- function(key) {
    ki <- match(key, fmt_keys)
    if (is.na(ki)) return(rep(NA_character_, nrow(d)))
    sapply(seq_len(nrow(d)), function(r) {
      vals <- strsplit(as.character(d[[sc]][r]), ":")[[1]]
      if (length(vals) >= ki) vals[ki] else NA_character_
    })
  }
  
  ES <- as.numeric(extract_field("ES"))
  SE <- as.numeric(extract_field("SE"))
  LP <- as.numeric(extract_field("LP"))
  AF <- as.numeric(extract_field("AF"))
  SS <- as.numeric(extract_field("SS"))
  EZ <- as.numeric(extract_field("EZ"))
  
  # 容错：ES缺失但EZ/SE有 → beta = EZ * SE
  if (all(is.na(ES)) && !all(is.na(EZ)) && !all(is.na(SE))) {
    ES <- EZ * SE
    cat("  Note: ES derived from EZ*SE\n")
  }
  
  # p值
  p <- ifelse(is.na(LP), NA_real_, 10^(-LP))
  if (all(is.na(p)) && !all(is.na(EZ))) {
    p <- 2 * pnorm(abs(EZ), lower.tail = FALSE)
  }
  
  # 检查是否解析成功
  if (all(is.na(ES))) {
    cat("  ⚠️ First 3 FORMAT entries:\n")
    for (k in 1:min(3, nrow(d))) {
      cat(sprintf("    row%d: format=%s | %s=%s\n", k, d$format[k], sc, substr(as.character(d[[sc]][k]), 1, 60)))
    }
    stop(sprintf("[%s] ES not found. Check FORMAT keys above.", trait))
  }
  
  res <- data.table(
    snp   = d$id,
    chr   = suppressWarnings(as.integer(gsub("chr", "", d$chr))),
    bp    = suppressWarnings(as.integer(d$pos)),
    a1    = d$alt,    # ALT = effect allele
    a2    = d$ref,    # REF = other allele
    beta  = ES,
    se    = SE,
    p     = p,
    af    = AF,
    n     = SS
    
  )
  
  # QC
  res <- res[!is.na(beta) & !is.na(se) & !is.na(p) & p > 0 & p < 1]
  res <- res[!is.na(chr) & chr %in% 1:22]
  
  # 如果 snp 是 "." 或 NA，用 chr:pos:ref:alt 构建
  res[snp == "." | is.na(snp), snp := sprintf("%d:%d:%s:%s", chr, bp, a2, a1)]
  
  cat(sprintf("  %s: %d variants (ES=%d, AF=%d, N=%d, p<5e-8=%d)\n",
              trait, nrow(res), sum(!is.na(res$beta)), sum(!is.na(res$af)),
              sum(!is.na(res$n)), sum(res$p < 5e-8, na.rm = TRUE)))
  
  # ===== 内置诊断（不往外传，不拖慢速度）=====
  cat(sprintf("  [诊断] FORMAT 前5行:\n"))
  fmt_sample <- d$format[1:min(5, nrow(d))]
  for (k in seq_along(fmt_sample)) {
    cat(sprintf("    [%d] %s\n", k, fmt_sample[k]))
  }
  
  # 检查 SS 是否在任何行的 FORMAT 里
  if (any(grepl("SS", d$format, fixed = TRUE))) {
    cat("  [诊断] ✅ 部分行 FORMAT 含 'SS'\n")
  } else {
    cat("  [诊断] ❌ 所有行 FORMAT 均不含 'SS' → n 列将全为 NA\n")
  }
  
  # 统计 FORMAT 字符串的唯一值
  fmt_counts <- head(sort(table(d$format), decreasing = TRUE), 3)
  cat(sprintf("  [诊断] FORMAT 唯一值 top-3:\n"))
  for (i in seq_along(fmt_counts)) {
    cat(sprintf("    '%s' (%.1f%%)\n", names(fmt_counts)[i], 
                100 * fmt_counts[i] / nrow(d)))
  }
  return(res)
}

## ========== 3. 读取三个文件 ==========
dat_fg  <- read_ieu_vcf(GWAS_FG, "FG")
dat_tg  <- read_ieu_vcf(GWAS_TG, "TG")
dat_out <- read_ieu_vcf(GWAS_PHENOAGE, "PhenoAge")
## ========== 3.5 诊断：SS/FORMAT 解析检查 ==========
diagnose_vcf_format <- function(dat, trait) {
  cat(sprintf("\n===== [%s] FORMAT 诊断 =====\n", trait))
  
  # 1. n 列（SS 提取结果）
  cat(sprintf("  n 列非 NA 数: %d / %d\n", sum(!is.na(dat$n)), nrow(dat)))
  if (sum(!is.na(dat$n)) > 0) {
    cat(sprintf("  n 值范围: %.0f - %.0f\n", min(dat$n, na.rm=TRUE), max(dat$n, na.rm=TRUE)))
    cat(sprintf("  n 值众数: %.0f\n", as.numeric(names(which.max(table(dat$n, exclude=NULL))))))
  } else {
    cat("  ⚠️ n 列全为 NA — SS 未从 FORMAT 中解析到\n")
  }
  
  # 2. 如果你在 read_ieu_vcf 里保存了原始 format 字符串，检查它
  # 假设 dat 里有 format_raw 列（需要在 read_ieu_vcf 里加一行保存）
  if ("format_raw" %in% names(dat)) {
    fmt_counts <- head(sort(table(dat$format_raw), decreasing=TRUE), 5)
    cat(sprintf("  FORMAT 字符串 top-5:\n"))
    for (i in seq_along(fmt_counts)) {
      cat(sprintf("    [%d] '%s' (出现 %d 次)\n", i, names(fmt_counts)[i], fmt_counts[i]))
    }
    # 检查是否含 SS
    has_ss <- any(grepl("SS", names(fmt_counts)))
    cat(sprintf("  含 'SS' key: %s\n", has_ss))
  } else {
    cat("  ⚠️ 未保存 format_raw 列，无法检查 FORMAT 字符串分布\n")
  }
  
  # 3. p 值分布（确认数据质量）
  cat(sprintf("  p 值分布: min=%.2e, median=%.2e, max=%.2e\n", 
              min(dat$p), median(dat$p), max(dat$p)))
  cat(sprintf("  p<5e-8: %d SNPs\n", sum(dat$p < 5e-8, na.rm=TRUE)))
  cat(sprintf("  p<1e-5: %d SNPs\n", sum(dat$p < 1e-5, na.rm=TRUE)))
  
  # 4. beta/se 完整性
  cat(sprintf("  beta 非 NA: %d, se 非 NA: %d\n", sum(!is.na(dat$beta)), sum(!is.na(dat$se))))
  
  invisible(NULL)
}

# 逐个诊断
diagnose_vcf_format(dat_fg, "FG")
diagnose_vcf_format(dat_tg, "TG")
diagnose_vcf_format(dat_out, "PhenoAge")
## ========== 4. plink clump ==========
clump_plink <- function(dat, trait, plink, bfile, p_thresh = 5e-8, r2 = 0.001, kb = 10000) {
  cat(sprintf("\n[%s] Clumping %s (p<%g)...\n", Sys.time(), trait, p_thresh))
  
  if (!file.exists(plink)) stop(sprintf("plink not found: %s", plink))
  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(bfile, ext))) stop(sprintf("Missing %s%s", bfile, ext))
  }
  
  sig <- dat[p < p_thresh]
  if (nrow(sig) == 0) {
    cat(sprintf("  ⚠️ No SNPs with p<%g. p distribution:\n", p_thresh))
    cat(sprintf("    min p=%.2e, median p=%.2e, max p=%.2e\n", 
                min(dat$p, na.rm=TRUE), median(dat$p, na.rm=TRUE), max(dat$p, na.rm=TRUE)))
    stop(sprintf("[%s] Cannot clump without significant SNPs.", trait))
  }
  
  inp <- tempfile(fileext = ".clump.txt")
  fwrite(sig[, .(SNP = snp, P = p, CHR = chr, BP = bp)], inp, sep = "\t")
  
  out_prefix <- tempfile(pattern = "clump_")
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-p1 %g --clump-p2 1 --clump-r2 %g --clump-kb %d --clump-snp-field SNP --clump-field P --out %s --allow-no-sex --silent',
    plink, bfile, inp, p_thresh, r2, kb, out_prefix
  )
  cat(sprintf("  CMD: %s\n", cmd))
  system(cmd)
  
  clf <- paste0(out_prefix, ".clumped")
  if (!file.exists(clf)) {
    stop(sprintf("[%s] Clump failed - no .clumped output. Check plink/1KG.", trait))
  }
  
  cl <- fread(clf, fill = TRUE)
  if ("SNP" %in% names(cl)) {
    ivs <- unique(na.omit(cl$SNP))
  } else if ("ID" %in% names(cl)) {
    ivs <- unique(na.omit(cl$ID))
  } else {
    stop(sprintf("[%s] Unexpected clump output format", trait))
  }
  
  cat(sprintf("  %s: %d IVs retained (from %d significant SNPs)\n", trait, length(ivs), nrow(sig)))
  return(ivs)
}

iv_fg <- clump_plink(dat_fg, "FG", PLINK, KG_EUR)
iv_tg <- clump_plink(dat_tg, "TG", PLINK, KG_EUR)

iv_union <- unique(c(iv_fg, iv_tg))
cat(sprintf("\nUnion IVs: %d SNPs\n", length(iv_union)))

## ========== 5. 构建 MVMR 数据 ==========
cat("\nBuilding MVMR dataset...\n")
m <- dat_out[snp %in% iv_union, .(snp, chr, bp, a1_out = a1, a2_out = a2, 
                                  beta_out = beta, se_out = se, p_out = p, af_out = af, n_out = n)]
m <- merge(m, dat_fg[snp %in% iv_union, .(snp, a1_fg = a1, a2_fg = a2, 
                                          beta_fg = beta, se_fg = se, p_fg = p, af_fg = af, n_fg = n)], 
           by = "snp", all.x = TRUE)
m <- merge(m, dat_tg[snp %in% iv_union, .(snp, a1_tg = a1, a2_tg = a2, 
                                          beta_tg = beta, se_tg = se, p_tg = p, af_tg = af, n_tg = n)], 
           by = "snp", all.x = TRUE)

m <- m[!is.na(beta_fg) & !is.na(beta_tg) & !is.na(beta_out)]
cat(sprintf("SNPs with complete data: %d\n", nrow(m)))

## ========== 6. 等位基因对齐（以结局 A1 为基准） ==========
cat("Harmonising alleles...\n")
m$bfg  <- m$beta_fg
m$btg  <- m$beta_tg
m$bout <- m$beta_out

flip_fg <- m$a1_fg == m$a2_out & m$a2_fg == m$a1_out
flip_tg <- m$a1_tg == m$a2_out & m$a2_tg == m$a1_out

m$bfg[flip_fg]  <- -m$beta_fg[flip_fg]
m$btg[flip_tg]  <- -m$beta_tg[flip_tg]

align_ok <- (m$a1_fg == m$a1_out | flip_fg) & (m$a1_tg == m$a1_out | flip_tg)

# 回文(A/T 或 C/G) + 无AF → 丢
pal <- (m$a1_out %in% c("A","T") & m$a2_out %in% c("A","T")) | 
  (m$a1_out %in% c("C","G") & m$a2_out %in% c("C","G"))
drop_pal <- pal & (is.na(m$af_out) | is.na(m$af_fg) | is.na(m$af_tg))
keep <- align_ok & !drop_pal

cat(sprintf("  Retained %d SNPs (dropped %d ambiguous/palindromic)\n", sum(keep), sum(!keep)))
m <- m[keep]

## ========== 7. Steiger 方向性过滤 ==========
cat("Steiger filter...\n")
if (!all(is.na(m$af_fg)) && !all(is.na(m$n_fg))) {
  r2_fn <- function(af, b) 2 * af * (1 - af) * b^2
  F_fn  <- function(af, b, n) { r <- r2_fn(af, b); r * (n - 2) / (1 - r) }
  Ff <- F_fn(m$af_fg, m$bfg, m$n_fg)
  Ft <- F_fn(m$af_tg, m$btg, m$n_tg)
  Fo <- F_fn(m$af_out, m$bout, m$n_out)
  steiger_ok <- Ff > Fo & Ft > Fo
} else {
  cf <- (m$bfg / m$se_fg)^2
  ct <- (m$btg / m$se_tg)^2
  co <- (m$bout / m$se_out)^2
  steiger_ok <- cf > co & ct > co
}
cat(sprintf("  Steiger passed: %d / %d\n", sum(steiger_ok), nrow(m)))
m <- m[steiger_ok]

if (nrow(m) < 2) stop("Need ≥2 SNPs after Steiger for MVMR (2 exposures).")

## ========== 8. MVMR（修复参数名） ==========
cat(sprintf("\n[%s] Running MVMR with %d SNPs...\n", Sys.time(), nrow(m)))

mvobj <- mr_mvinput(
  bx      = cbind(m$bfg, m$btg),
  bxse    = cbind(m$se_fg, m$se_tg),
  by      = m$bout,
  byse    = m$se_out,
  exposure = c("FG", "TG"),
  outcome  = "PhenoAge"
)

cat("\n=== MVMR-IVW (Main) ===\n")
ivw <- mr_mvivw(mvobj, model = "default", correl = FALSE, distribution = "normal", alpha = 0.05)
print(ivw)

cat("\n=== MVMR-Egger ===\n")
egger <- tryCatch(mr_mvegger(mvobj, orientate = 1), 
                  error = function(e) { cat("  Egger error:", e$message, "\n"); NULL })
if (!is.null(egger)) print(egger)

cat("\n=== MVMR-Median ===\n")
med <- tryCatch(mr_mvmedian(mvobj), 
                error = function(e) { cat("  Median error:", e$message, "\n"); NULL })
if (!is.null(med)) print(med)

## ========== 8.5 手动回填 GWAS 样本量（解决 NaN） ==========
cat("\n=== Sample Size ===\n")
N_FG  <- 200622   # GCST90002232  Chen J 2021 MAGIC+UKB
N_TG  <- 115082   # GCST90092992  Richardson TG 2022 UKB
N_OUT <- 34463    # GCST90014292  McCartney DL 2021 UKB

m$n_fg  <- N_FG
m$n_tg  <- N_TG
m$n_out <- N_OUT

cat(sprintf("  FG=%d, TG=%d, PhenoAge=%d\n", N_FG, N_TG, N_OUT))
cat(sprintf("  Effective N (min): %d\n", min(N_FG, N_TG, N_OUT)))

## ========== 9. Diagnostics + strength_mvmr ==========
cat("\n=== Diagnostics ===\n")

# VIF
if (nrow(m) >= 3) {
  vif_val <- 1 / (1 - summary(lm(bfg ~ btg, data = m))$r.squared)
  cat(sprintf("  VIF (FG ~ TG): %.2f\n", vif_val))
} else {
  vif_val <- NA
  cat("  VIF: N/A (need ≥3 SNPs)\n")
}

# ---- strength_mvmr (MVMR包条件F) ----
cat("\n--- MVMR Conditional F (strength_mvmr) ---\n")
library(MVMR)

mvdat <- format_mvmr(
  BXGs   = cbind(m$bfg, m$btg),
  BYG    = m$bout,
  seBXGs = cbind(m$se_fg, m$se_tg),
  seBYG  = m$se_out,
  RSID   = m$snp
)

class(mvdat)  # 确认是 "mvmr_format"
cfs <- strength_mvmr(mvdat)
print(cfs)

# ---- 条件 F 阈值判断 ----
cfs_vals <- as.numeric(cfs)
if (any(cfs_vals < 10, na.rm = TRUE)) {
  cat("  ⚠️  WARNING: Conditional F < 10 for some exposure → weak instrument bias possible\n")
} else {
  cat("  ✅ All conditional F > 10 → instrument strength OK\n")
}

## ========== 10. Power Analysis (修正版) ==========
cat("\n=== MVMR Power Analysis (Corrected) ===\n")

pwr_from_se <- function(beta_alt, cond_se, alpha = 0.05) {
  ncp   <- abs(beta_alt) / cond_se
  zc    <- qnorm(1 - alpha / 2)
  # 双侧检验正确公式
  power <- pnorm(ncp - zc) + pnorm(-zc - ncp)
  # 80% power 需要的样本倍数
  z80   <- qnorm(0.80)
  se_need <- abs(beta_alt) / (z80 + zc)
  n_mult  <- (cond_se / se_need)^2
  return(list(ncp = ncp, power = power, se_need = se_need, n_mult = n_mult))
}

# 不同效应量下的 power 扫描
beta_grid <- c(0.01, 0.02, 0.05, 0.10, 0.20)
cat(sprintf("\n  Power under different assumed effects:\n"))
cat(sprintf("  %-10s %10s %10s %10s\n", "β_alt", "FG_Power", "TG_Power", "CondSE_FG"))
cat(sprintf("  %-10s %10s %10s %10s\n", "----------", "--------", "--------", "--------"))
for (b in beta_grid) {
  pf_b <- pwr_from_se(b, ivw$StdError[1])
  pt_b <- pwr_from_se(b, ivw$StdError[2])
  cat(sprintf("  β=%-8.3f %9.1f%% %9.1f%% %10.4f\n", 
              b, 100*pf_b$power, 100*pt_b$power, ivw$StdError[1]))
}

# 主效应量
beta_alt_fg <- 0.05
beta_alt_tg <- 0.05
pf <- pwr_from_se(beta_alt_fg, ivw$StdError[1])
pt <- pwr_from_se(beta_alt_tg, ivw$StdError[2])

cat(sprintf("\n  ┌─────────────────────────────────────────────────────────┐\n"))
cat(sprintf("  │  Exposure │ CondSE  │ NCP   │ Power   │ 80%% N-mult │\n"))
cat(sprintf("  ├─────────────────────────────────────────────────────────┤\n"))
cat(sprintf("  │  FG       │ %6.4f │ %5.2f │ %5.1f%%   │ %.0fx       │\n",
            ivw$StdError[1], pf$ncp, 100*pf$power, pf$n_mult))
cat(sprintf("  │  TG       │ %6.4f │ %5.2f │ %5.1f%%   │ %.0fx       │\n",
            ivw$StdError[2], pt$ncp, 100*pt$power, pt$n_mult))
cat(sprintf("  └─────────────────────────────────────────────────────────┘\n"))

cat(sprintf("\n  Interpretation:\n"))
cat(sprintf("  ⚠️  FG Power = %.1f%% → SEVERELY UNDERPOWERED\n", 100*pf$power))
cat(sprintf("  ⚠️  TG Power = %.1f%% → SEVERELY UNDERPOWERED\n", 100*pt$power))
cat(sprintf("  ⚠️  Negative results (p=%.3f, p=%.3f) do NOT mean no causality\n", 
            ivw$Pvalue[1], ivw$Pvalue[2]))
cat(sprintf("  → Need ~%.0fx (FG) / ~%.0fx (TG) more samples for 80%% power at β=0.05\n",
            pf$n_mult, pt$n_mult))
cat(sprintf("  → Or: need PhenoAge GWAS N ≥ %d for FG, ≥ %d for TG\n",
            round(34463 * pf$n_mult), round(34463 * pt$n_mult)))
## ========== 11. 留一法 ==========
cat("\n=== Leave-One-Out ===\n")
loo_list <- list()
for (i in seq_len(nrow(m))) {
  idx <- setdiff(seq_len(nrow(m)), i)
  obj_loo <- mr_mvinput(
    bx       = cbind(m$bfg, m$btg)[idx, , drop = FALSE],
    bxse     = cbind(m$se_fg, m$se_tg)[idx, , drop = FALSE],
    by       = m$bout[idx],
    byse     = m$se_out[idx],
    exposure = c("FG", "TG"),
    outcome  = "PhenoAge"
  )
  r <- tryCatch(mr_mvivw(obj_loo), error = function(e) NULL)
  if (!is.null(r)) {
    loo_list[[i]] <- data.frame(
      leftout = m$snp[i],
      beta_FG = r$Estimate[1],
      beta_TG = r$Estimate[2],
      p_FG    = r$Pvalue[1],
      p_TG    = r$Pvalue[2]
    )
  }
}
if (length(loo_list) > 0) {
  loo_df <- do.call(rbind, loo_list)
  fwrite(loo_df, file.path(OUT_DIR, "mvmr_loo.csv"), sep = "\t")
  cat(sprintf("  LOO saved (%d iterations)\n", nrow(loo_df)))
}

## ========== 12. 导出所有结果 ==========
cat("\n=== Exporting Results ===\n")

# IVW 主结果
main_res <- data.frame(
  Exposure = c("FG", "TG"),
  Beta     = ivw$Estimate,
  SE       = ivw$StdError,
  CI_lower = ivw$CILower,
  CI_upper = ivw$CIUpper,
  Pvalue   = ivw$Pvalue
)
fwrite(main_res, file.path(OUT_DIR, "mvmr_ivw_main.csv"), sep = "\t")

# Power 结果
power_tab <- data.frame(
  Exposure      = c("FG", "TG"),
  N_exposure    = c(N_FG, N_TG),
  N_outcome     = c(N_OUT, N_OUT),
  IVW_Beta      = ivw$Estimate,
  CondSE        = ivw$StdError,
  Assumed_Beta  = c(beta_alt_fg, beta_alt_tg),
  NCP           = c(pf$ncp, pt$ncp),
  Power         = c(pf$power, pt$power),
  N_for_80pc    = c(pf$n_mult, pt$n_mult)
)
fwrite(power_tab, file.path(OUT_DIR, "mvmr_power.csv"), sep = "\t")

# 条件 F
cfs_df <- data.frame(
  Exposure = c("FG", "TG"),
  Conditional_F = cfs_vals[1:2]
)
fwrite(cfs_df, file.path(OUT_DIR, "mvmr_conditionalF.csv"), sep = "\t")

# Audit
audit <- data.frame(
  Step  = c("FG_IVs", "TG_IVs", "Union", "Complete", "Aligned_Steiger", "Final_MVMR", "VIF", "Min_N", "FG_Power", "TG_Power", "FG_CondF", "TG_CondF"),
  Value = c(length(iv_fg), length(iv_tg), length(iv_union), 
            nrow(m) + sum(!keep), nrow(m), nrow(m), 
            ifelse(is.na(vif_val), "NA", round(vif_val, 2)),
            min(N_FG, N_TG, N_OUT),
            round(pf$power, 3), round(pt$power, 3),
            round(cfs_vals[1], 1), round(cfs_vals[2], 1))
)
fwrite(audit, file.path(OUT_DIR, "mvmr_audit.csv"), sep = "\t")

cat(sprintf("\n✅ Done! All results in %s\n", OUT_DIR))
cat("  - mvmr_ivw_main.csv        : IVW estimates\n")
cat("  - mvmr_power.csv           : Power analysis\n")
cat("  - mvmr_conditionalF.csv   : Conditional F-statistics\n")
cat("  - mvmr_loo.csv             : Leave-one-out sensitivity\n")
cat("  - mvmr_audit.csv           : QC audit trail\n")