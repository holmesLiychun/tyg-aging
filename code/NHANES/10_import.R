# ============================================================================
# 模块 10：数据导入与变量衍生（match 直接赋值版 —— 无 join 无 merge）
# ============================================================================

timer_start("模块10 数据导入与变量衍生")

if (SKIP_IF_CACHED && file.exists(CACHE_RAW_DERIVED)) {
  nhanes <- readRDS(CACHE_RAW_DERIVED)
  log_write(sprintf("✓ 从缓存加载主数据框: %d 行 x %d 列", nrow(nhanes), ncol(nhanes)))
  diag_vars(nhanes, c("age", "tyg", "homa_ir", "dm"), "缓存数据关键变量")
  timer_end("模块10 数据导入与变量衍生")
  N_CYC <- length(CYCLES)
  log_write("===== 模块 10 结束（已加载缓存）=====\n")
} else {
  
  # ---- 文件路径（不变）----
  FILES <- list(
    "2005-2006" = list(demo = "Demographics/DEMO_D.XPT", bmx = "Examination/BMX_D.XPT",
                       bpx = "Examination/BPX_D.XPT", bio = "Laboratory/BIOPRO_D.XPT",
                       glu = "Laboratory/GLU_D.XPT", cbc = "Laboratory/CBC_D.XPT",
                       ghb = "Laboratory/GHB_D.XPT", crp = "Laboratory/CRP_D.XPT",
                       hdl = "Laboratory/HDL_D.XPT", albcr = "Laboratory/ALB_CR_D.XPT",
                       smq = "Questionnaire/SMQ_D.XPT", diq = "Questionnaire/DIQ_D.XPT",
                       alq = "Questionnaire/ALQ_D.XPT", paq = "Questionnaire/PAQ_D.XPT"),
    "2007-2008" = list(demo = "Demographics/DEMO_E.XPT", bmx = "Examination/BMX_E.XPT",
                       bpx = "Examination/BPX_E.XPT", bio = "Laboratory/BIOPRO_E.XPT",
                       glu = "Laboratory/GLU_E.XPT", cbc = "Laboratory/CBC_E.XPT",
                       ghb = "Laboratory/GHB_E.XPT", crp = "Laboratory/CRP_E.XPT",
                       hdl = "Laboratory/HDL_E.XPT", albcr = "Laboratory/ALB_CR_E.XPT",
                       smq = "Questionnaire/SMQ_E.XPT", diq = "Questionnaire/DIQ_E.XPT",
                       alq = "Questionnaire/ALQ_E.XPT", paq = "Questionnaire/PAQ_E.XPT"),
    "2009-2010" = list(demo = "Demographics/DEMO_F.XPT", bmx = "Examination/BMX_F.XPT",
                       bpx = "Examination/BPX_F.XPT", bio = "Laboratory/BIOPRO_F.XPT",
                       glu = "Laboratory/GLU_F.XPT", cbc = "Laboratory/CBC_F.XPT",
                       ghb = "Laboratory/GHB_F.XPT", crp = "Laboratory/CRP_F.XPT",
                       hdl = "Laboratory/HDL_F.XPT", albcr = "Laboratory/ALB_CR_F.XPT",
                       smq = "Questionnaire/SMQ_F.XPT", diq = "Questionnaire/DIQ_F.XPT",
                       alq = "Questionnaire/ALQ_F.XPT", paq = "Questionnaire/PAQ_F.XPT"),
    "2011-2012" = list(demo = "Demographics/DEMO_G.XPT", bmx = "Examination/BMX_G.XPT",
                       bpx = "Examination/BPX_G.XPT", bio = "Laboratory/BIOPRO_G.XPT",
                       glu = "Laboratory/GLU_G.XPT", cbc = "Laboratory/CBC_G.XPT",
                       ghb = "Laboratory/GHB_G.XPT", crp = NULL,
                       hdl = "Laboratory/HDL_G.XPT", albcr = "Laboratory/ALB_CR_G.XPT",
                       smq = "Questionnaire/SMQ_G.XPT", diq = "Questionnaire/DIQ_G.XPT",
                       alq = "Questionnaire/ALQ_G.XPT", paq = "Questionnaire/PAQ_G.XPT"),
    "2013-2014" = list(demo = "Demographics/DEMO_H.XPT", bmx = "Examination/BMX_H.XPT",
                       bpx = "Examination/BPX_H.XPT", bio = "Laboratory/BIOPRO_H.XPT",
                       glu = "Laboratory/GLU_H.XPT", cbc = "Laboratory/CBC_H.XPT",
                       ghb = "Laboratory/GHB_H.XPT", crp = NULL,
                       hdl = "Laboratory/HDL_H.XPT", albcr = "Laboratory/ALB_CR_H.XPT",
                       smq = "Questionnaire/SMQ_H.XPT", diq = "Questionnaire/DIQ_H.XPT",
                       alq = "Questionnaire/ALQ_H.XPT", paq = "Questionnaire/PAQ_H.XPT"),
    "2015-2016" = list(demo = "Demographics/DEMO_I.XPT", bmx = "Examination/BMX_I.XPT",
                       bpx = "Examination/BPX_I.XPT", bio = "Laboratory/BIOPRO_I.XPT",
                       glu = "Laboratory/GLU_I.XPT", cbc = "Laboratory/CBC_I.XPT",
                       ghb = "Laboratory/GHB_I.XPT", crp = "Laboratory/HSCRP_I.XPT",
                       hdl = "Laboratory/HDL_I.XPT", albcr = "Laboratory/ALB_CR_I.XPT",
                       smq = "Questionnaire/SMQ_I.XPT", diq = "Questionnaire/DIQ_I.XPT",
                       alq = "Questionnaire/ALQ_I.XPT", paq = "Questionnaire/PAQ_I.XPT"),
    "2017-2018" = list(demo = "Demographics/DEMO_J.XPT", bmx = "Examination/BMX_J.XPT",
                       bpx = "Examination/BPX_J.XPT", bio = "Laboratory/BIOPRO_J.XPT",
                       glu = "Laboratory/GLU_J.XPT", cbc = "Laboratory/CBC_J.XPT",
                       ghb = "Laboratory/GHB_J.XPT", crp = "Laboratory/HSCRP_J.XPT",
                       hdl = "Laboratory/HDL_J.XPT", albcr = "Laboratory/ALB_CR_J.XPT",
                       smq = "Questionnaire/SMQ_J.XPT", diq = "Questionnaire/DIQ_J.XPT",
                       alq = "Questionnaire/ALQ_J.XPT", paq = "Questionnaire/PAQ_J.XPT")
  )
  
  read_xpt_safe <- function(fp) {
    if (!file.exists(fp)) { log_write(sprintf("  WARN 缺文件: %s", fp)); return(NULL) }
    tryCatch({
      df <- as.data.frame(haven::read_xpt(fp))
      df[] <- lapply(df, function(x) if (haven::is.labelled(x)) { if (is.numeric(x)) as.numeric(x) else if (is.character(x)) as.character(x) else x } else x)
      df
    }, error = function(e) { log_write(sprintf("  ERROR: %s", e$message)); NULL })
  }
  
  read_cycle <- function(cy, f) {
    cp <- file.path(NHANES_ROOT, cy)
    demo <- read_xpt_safe(file.path(cp, f$demo))
    if (is.null(demo)) { log_write(sprintf("  ✗ %s demo 缺失", cy)); return(NULL) }
    
    # 基准 data.frame，只含 SEQN + demo 变量
    out <- data.frame(
      SEQN       = demo$SEQN,
      RIDAGEYR   = as.numeric(demo$RIDAGEYR),
      RIAGENDR   = as.numeric(demo$RIAGENDR),
      RIDRETH1   = as.numeric(demo$RIDRETH1),
      DMDEDUC2   = as.numeric(demo$DMDEDUC2),
      INDFMPIR   = as.numeric(demo$INDFMPIR),
      RIDEXPRG   = as.numeric(demo$RIDEXPRG),
      RIDSTATR   = as.numeric(demo$RIDSTATR),
      WTMEC2YR   = as.numeric(demo$WTMEC2YR),
      SDMVSTRA   = as.numeric(demo$SDMVSTRA),
      SDMVPSU    = as.numeric(demo$SDMVPSU),
      cycle      = cy,
      stringsAsFactors = FALSE
    )
    n <- nrow(out)
    
    # 预建所有变量列为 NA
    for (v in c("crp","hdl","uacr","SMQ040","DIQ010","DIQ050","wtsfa",
                "sbp","dbp","bmi","waist","albumin","alp","bun","creatinine",
                "ggt","uric_acid","totchol","trig","glucose","insulin",
                "wbc","lymph","mcv","rdw","plate","hba1c","pad200","pad320",
                "paq560","lvig_d","lvig_m","lmod_d","lmod_m","wvig_d","wvig_m",
                "wmod_d","wmod_m","trans_d","trans_m","paq_any_mv",
                "ALQ101","ALQ110","ALQ120Q","ALQ120U","ALQ130","ALQ111","ALQ121")) {
      out[[v]] <- NA_real_
    }
    
    # ---- 辅助：按 SEQN 匹配赋值 ----
    assign_by_seqn <- function(src_df, var_map) {
      # var_map: 命名向量，names = 目标列名，values = 源列名
      if (is.null(src_df)) return()
      idx <- match(out$SEQN, src_df$SEQN)
      ok <- !is.na(idx)
      for (target in names(var_map)) {
        src_col <- var_map[[target]]
        if (src_col %in% names(src_df)) {
          val <- haven::zap_labels(src_df[[src_col]][idx[ok]])
          out[[target]][ok] <<- as.numeric(val)
        }
      }
    }
    
    # 2) 血压
    bpx <- read_xpt_safe(file.path(cp, f$bpx))
    if (!is.null(bpx)) {
      sbp_cols <- intersect(c("BPXSY1","BPXSY2","BPXSY3","BPXSY4"), names(bpx))
      dbp_cols <- intersect(c("BPXDI1","BPXDI2","BPXDI3","BPXDI4"), names(bpx))
      idx <- match(out$SEQN, bpx$SEQN)
      ok <- !is.na(idx)
      out$sbp[ok] <- as.numeric(rowMeans(bpx[idx[ok], sbp_cols, drop=FALSE], na.rm=TRUE))
      out$dbp[ok] <- as.numeric(rowMeans(bpx[idx[ok], dbp_cols, drop=FALSE], na.rm=TRUE))
    }
    
    # 3) 体格
    bmx <- read_xpt_safe(file.path(cp, f$bmx))
    if (!is.null(bmx)) {
      idx <- match(out$SEQN, bmx$SEQN); ok <- !is.na(idx)
      out$bmi[ok]   <- as.numeric(haven::zap_labels(bmx$BMXBMI[idx[ok]]))
      out$waist[ok] <- as.numeric(haven::zap_labels(bmx$BMXWAIST[idx[ok]]))
    }
    
    # 4) 生化
    bio <- read_xpt_safe(file.path(cp, f$bio))
    if (!is.null(bio)) {
      pick <- function(df, cands) {
        for (v in cands) if (v %in% names(df) && any(!is.na(df[[v]]))) return(df[[v]])
        cands[1]
      }
      idx <- match(out$SEQN, bio$SEQN); ok <- !is.na(idx)
      for (target in c("albumin","alp","bun","creatinine","ggt","uric_acid","totchol","trig")) {
        # 映射
      }
      # albumin
      if ("LBXSAL" %in% names(bio)) out$albumin[ok] <- as.numeric(haven::zap_labels(bio$LBXSAL[idx[ok]]))
      else if ("LBDSALSI" %in% names(bio)) out$albumin[ok] <- as.numeric(haven::zap_labels(bio$LBDSALSI[idx[ok]]))
      if ("LBXSAPSI" %in% names(bio)) out$alp[ok] <- as.numeric(haven::zap_labels(bio$LBXSAPSI[idx[ok]]))
      else if ("LBXSALP" %in% names(bio)) out$alp[ok] <- as.numeric(haven::zap_labels(bio$LBXSALP[idx[ok]]))
      if ("LBXSBU" %in% names(bio)) out$bun[ok] <- as.numeric(haven::zap_labels(bio$LBXSBU[idx[ok]]))
      else if ("LBXBUN" %in% names(bio)) out$bun[ok] <- as.numeric(haven::zap_labels(bio$LBXBUN[idx[ok]]))
      if ("LBXSCR" %in% names(bio)) out$creatinine[ok] <- as.numeric(haven::zap_labels(bio$LBXSCR[idx[ok]]))
      if ("LBXSGTSI" %in% names(bio)) out$ggt[ok] <- as.numeric(haven::zap_labels(bio$LBXSGTSI[idx[ok]]))
      else if ("LBXSGT" %in% names(bio)) out$ggt[ok] <- as.numeric(haven::zap_labels(bio$LBXSGT[idx[ok]]))
      if ("LBXSUA" %in% names(bio)) out$uric_acid[ok] <- as.numeric(haven::zap_labels(bio$LBXSUA[idx[ok]]))
      if ("LBXSCH" %in% names(bio)) out$totchol[ok] <- as.numeric(haven::zap_labels(bio$LBXSCH[idx[ok]]))
      else if ("LBXTC" %in% names(bio)) out$totchol[ok] <- as.numeric(haven::zap_labels(bio$LBXTC[idx[ok]]))
      if ("LBXSTR" %in% names(bio)) out$trig[ok] <- as.numeric(haven::zap_labels(bio$LBXSTR[idx[ok]]))
      else if ("LBXTR" %in% names(bio)) out$trig[ok] <- as.numeric(haven::zap_labels(bio$LBXTR[idx[ok]]))
    }
    
    # 5) 血糖/胰岛素
    glu <- read_xpt_safe(file.path(cp, f$glu))
    if (!is.null(glu)) {
      idx <- match(out$SEQN, glu$SEQN); ok <- !is.na(idx)
      if ("LBXGLU" %in% names(glu)) out$glucose[ok] <- as.numeric(haven::zap_labels(glu$LBXGLU[idx[ok]]))
      if ("LBXIN" %in% names(glu)) out$insulin[ok] <- as.numeric(haven::zap_labels(glu$LBXIN[idx[ok]]))
      else if ("LBXINSI" %in% names(glu)) out$insulin[ok] <- as.numeric(haven::zap_labels(glu$LBXINSI[idx[ok]]))
      if ("WTSAF2YR" %in% names(glu)) out$wtsfa[ok] <- as.numeric(haven::zap_labels(glu$WTSAF2YR[idx[ok]]))
    }
    
    # 6) 血常规
    cbc <- read_xpt_safe(file.path(cp, f$cbc))
    if (!is.null(cbc)) {
      idx <- match(out$SEQN, cbc$SEQN); ok <- !is.na(idx)
      if ("LBXWBCSI" %in% names(cbc)) out$wbc[ok] <- as.numeric(haven::zap_labels(cbc$LBXWBCSI[idx[ok]]))
      else if ("LBXWBC" %in% names(cbc)) out$wbc[ok] <- as.numeric(haven::zap_labels(cbc$LBXWBC[idx[ok]]))
      if ("LBXLYPCT" %in% names(cbc)) out$lymph[ok] <- as.numeric(haven::zap_labels(cbc$LBXLYPCT[idx[ok]]))
      else if ("LBXLYPCTS" %in% names(cbc)) out$lymph[ok] <- as.numeric(haven::zap_labels(cbc$LBXLYPCTS[idx[ok]]))
      if ("LBXMCVSI" %in% names(cbc)) out$mcv[ok] <- as.numeric(haven::zap_labels(cbc$LBXMCVSI[idx[ok]]))
      else if ("LBXMCV" %in% names(cbc)) out$mcv[ok] <- as.numeric(haven::zap_labels(cbc$LBXMCV[idx[ok]]))
      if ("LBXRDW" %in% names(cbc)) out$rdw[ok] <- as.numeric(haven::zap_labels(cbc$LBXRDW[idx[ok]]))
      else if ("LBXRDWTS" %in% names(cbc)) out$rdw[ok] <- as.numeric(haven::zap_labels(cbc$LBXRDWTS[idx[ok]]))
      if ("LBXPLTSI" %in% names(cbc)) out$plate[ok] <- as.numeric(haven::zap_labels(cbc$LBXPLTSI[idx[ok]]))
      else if ("LBXPLT" %in% names(cbc)) out$plate[ok] <- as.numeric(haven::zap_labels(cbc$LBXPLT[idx[ok]]))
    }
    
    # 7) 糖化
    ghb <- read_xpt_safe(file.path(cp, f$ghb))
    if (!is.null(ghb)) {
      idx <- match(out$SEQN, ghb$SEQN); ok <- !is.na(idx)
      if ("LBXGH" %in% names(ghb)) out$hba1c[ok] <- as.numeric(haven::zap_labels(ghb$LBXGH[idx[ok]]))
    }
    
    # 8) CRP —— 改为 match 赋值（与其他段一致）
    if (!is.null(f$crp)) {
      crp_raw <- read_xpt_safe(file.path(cp, f$crp))
      if (!is.null(crp_raw)) {
        idx <- match(out$SEQN, crp_raw$SEQN)
        ok <- !is.na(idx)
        out$crp[ok] <- if ("LBXHSCRP" %in% names(crp_raw)) {
          as.numeric(haven::zap_labels(crp_raw$LBXHSCRP[idx[ok]]))
        } else if ("LBXCRP" %in% names(crp_raw)) {
          as.numeric(haven::zap_labels(crp_raw$LBXCRP[idx[ok]]))
        }
      }
    }
    
    # 9) HDL
    hdl_raw <- read_xpt_safe(file.path(cp, f$hdl))
    if (!is.null(hdl_raw) && "LBDHDD" %in% names(hdl_raw)) {
      idx <- match(out$SEQN, hdl_raw$SEQN); ok <- !is.na(idx)
      out$hdl[ok] <- as.numeric(haven::zap_labels(hdl_raw$LBDHDD[idx[ok]]))
    }
    
    # 10) 尿白蛋白
    albcr <- read_xpt_safe(file.path(cp, f$albcr))
    if (!is.null(albcr) && "URDACT" %in% names(albcr)) {
      idx <- match(out$SEQN, albcr$SEQN); ok <- !is.na(idx)
      out$uacr[ok] <- as.numeric(haven::zap_labels(albcr$URDACT[idx[ok]]))
    }
    
    # 11) 吸烟
    smq <- read_xpt_safe(file.path(cp, f$smq))
    if (!is.null(smq) && "SMQ040" %in% names(smq)) {
      idx <- match(out$SEQN, smq$SEQN); ok <- !is.na(idx)
      out$SMQ040[ok] <- as.numeric(haven::zap_labels(smq$SMQ040[idx[ok]]))
    }
    
    # 12) 糖尿病问卷
    diq <- read_xpt_safe(file.path(cp, f$diq))
    if (!is.null(diq)) {
      idx <- match(out$SEQN, diq$SEQN); ok <- !is.na(idx)
      if ("DIQ010" %in% names(diq)) out$DIQ010[ok] <- as.numeric(haven::zap_labels(diq$DIQ010[idx[ok]]))
      if ("DIQ050" %in% names(diq)) out$DIQ050[ok] <- as.numeric(haven::zap_labels(diq$DIQ050[idx[ok]]))
    }
    
    # 13) 饮酒
    alq <- read_xpt_safe(file.path(cp, f$alq))
    if (!is.null(alq)) {
      idx <- match(out$SEQN, alq$SEQN); ok <- !is.na(idx)
      for (v in c("ALQ101","ALQ110","ALQ120Q","ALQ120U","ALQ130","ALQ111","ALQ121")) {
        if (v %in% names(alq)) out[[v]][ok] <- as.numeric(haven::zap_labels(alq[[v]][idx[ok]]))
      }
    }
    
    # 14) 体力活动
    paq <- read_xpt_safe(file.path(cp, f$paq))
    if (!is.null(paq)) {
      idx <- match(out$SEQN, paq$SEQN); ok <- !is.na(idx)
      if (cy == "2005-2006") {
        if ("PAD200" %in% names(paq)) out$pad200[ok] <- as.numeric(haven::zap_labels(paq$PAD200[idx[ok]]))
        if ("PAD320" %in% names(paq)) out$pad320[ok] <- as.numeric(haven::zap_labels(paq$PAD320[idx[ok]]))
        if ("PAQ560" %in% names(paq)) out$paq560[ok] <- as.numeric(haven::zap_labels(paq$PAQ560[idx[ok]]))
        out$paq_any_mv <- ifelse((!is.na(out$pad200) & out$pad200==1) |
                                   (!is.na(out$pad320) & out$pad320==1) |
                                   (!is.na(out$paq560) & out$paq560>=1), 1, 0)
      } else {
        if ("PAQ655" %in% names(paq)) out$lvig_d[ok] <- as.numeric(haven::zap_labels(paq$PAQ655[idx[ok]]))
        if ("PAD660" %in% names(paq)) out$lvig_m[ok] <- as.numeric(haven::zap_labels(paq$PAD660[idx[ok]]))
        if ("PAQ670" %in% names(paq)) out$lmod_d[ok] <- as.numeric(haven::zap_labels(paq$PAQ670[idx[ok]]))
        if ("PAD675" %in% names(paq)) out$lmod_m[ok] <- as.numeric(haven::zap_labels(paq$PAD675[idx[ok]]))
        if ("PAQ610" %in% names(paq)) out$wvig_d[ok] <- as.numeric(haven::zap_labels(paq$PAQ610[idx[ok]]))
        if ("PAD615" %in% names(paq)) out$wvig_m[ok] <- as.numeric(haven::zap_labels(paq$PAD615[idx[ok]]))
        if ("PAQ625" %in% names(paq)) out$wmod_d[ok] <- as.numeric(haven::zap_labels(paq$PAQ625[idx[ok]]))
        if ("PAD630" %in% names(paq)) out$wmod_m[ok] <- as.numeric(haven::zap_labels(paq$PAD630[idx[ok]]))
        if ("PAQ640" %in% names(paq)) out$trans_d[ok] <- as.numeric(haven::zap_labels(paq$PAQ640[idx[ok]]))
        if ("PAD645" %in% names(paq)) out$trans_m[ok] <- as.numeric(haven::zap_labels(paq$PAD645[idx[ok]]))
      }
    }
    
    log_write(sprintf("  ✓ %s: %d 行 x %d 列", cy, nrow(out), ncol(out)))
    out
  }
  
  # ---- 逐周期读取并合并 ----
  data_list <- list()
  for (cy in CYCLES) {
    data_list[[cy]] <- tryCatch(read_cycle(cy, FILES[[cy]]),
                                error = function(e) { log_write(sprintf("  ✗ %s 读取失败: %s", cy, e$message)); NULL })
  }
  data_list <- Filter(Negate(is.null), data_list)
  dbg_stop(length(data_list) > 0, "没有任何周期读取成功")
  
  nhanes <- dplyr::bind_rows(data_list)
  log_write(sprintf("合并完成: %d 行 x %d 列 (%d 个周期)", nrow(nhanes), ncol(nhanes), length(data_list)))
  
  # ---- 变量衍生（mutate）----
  nhanes <- nhanes %>% mutate(
    age  = as.numeric(RIDAGEYR),
    sex  = factor(ifelse(RIAGENDR == 1, "Male", "Female"), levels = c("Male", "Female")),
    race = factor(case_when(
      RIDRETH1 == 1 ~ "Mexican American", RIDRETH1 == 2 ~ "Other Hispanic",
      RIDRETH1 == 3 ~ "Non-Hispanic White", RIDRETH1 == 4 ~ "Non-Hispanic Black",
      RIDRETH1 == 5 ~ "Other", TRUE ~ NA_character_),
      levels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other")),
    edu_cat = factor(case_when(
      DMDEDUC2 %in% c(1,2) ~ "Less than high school", DMDEDUC2 == 3 ~ "High school",
      DMDEDUC2 == 4 ~ "Some college", DMDEDUC2 == 5 ~ "College or above", TRUE ~ NA_character_),
      levels = c("Less than high school", "High school", "Some college", "College or above")),
    edu = as.numeric(DMDEDUC2), pir = as.numeric(INDFMPIR), preg = RIDEXPRG, exam_stat = RIDSTATR,
    
    # ★ wtsfa 来自 glu 文件，已在 read_cycle 里赋值，这里直接引用
    weight_fasting = as.numeric(ifelse(!is.na(wtsfa) & wtsfa > 0, wtsfa / length(CYCLES), NA_real_)),
    weight_mec    = as.numeric(ifelse(!is.na(WTMEC2YR) & WTMEC2YR > 0, WTMEC2YR / length(CYCLES), NA_real_)),
    weight = weight_fasting,
    strata = interaction(cycle, SDMVSTRA, drop = TRUE),
    psu    = interaction(cycle, SDMVPSU, drop = TRUE),
    fasting = ifelse(!is.na(wtsfa) & wtsfa > 0, 1L, 0L),
    
    lncreat = log(creatinine + 1),
    lncrp   = log(crp + 1),
    insulin_uIU = ifelse(!is.na(insulin), ifelse(insulin > 200, insulin / 6.945, insulin), NA_real_),
    glucose_mmol = glucose / 18.0,
    tyg = ifelse(fasting == 1 & glucose > 0 & trig > 0, log((glucose * trig) / 2), NA_real_),
    homa_ir = ifelse(fasting == 1 & glucose_mmol > 0 & insulin_uIU > 0, glucose_mmol * insulin_uIU / 22.5, NA_real_),
    
    smoke = factor(case_when(SMQ040 == 1 ~ "Current", SMQ040 == 2 ~ "Former", SMQ040 == 3 ~ "Never", TRUE ~ NA_character_),
                   levels = c("Never", "Former", "Current")),
    alc_ever = case_when(!is.na(ALQ110) ~ ifelse(ALQ110==1,"ever","never"), !is.na(ALQ111) ~ ifelse(ALQ111==1,"ever","never"), TRUE ~ NA_character_),
    alc_pastyr = case_when(!is.na(ALQ101) ~ ifelse(ALQ101==1,"yes","no"), !is.na(ALQ121) ~ ifelse(ALQ121 %in% 1:6,"yes","no"), TRUE ~ NA_character_),
    alcohol = factor(case_when(alc_ever=="never" ~ "Never", alc_ever=="ever" & alc_pastyr=="no" ~ "Former", alc_ever=="ever" & alc_pastyr=="yes" ~ "Current", TRUE ~ NA_character_),
                     levels = c("Never","Former","Current")),
    avg_drinks_day = ALQ130,
    
    leis_vig_met = ifelse(!is.na(lvig_d)&!is.na(lvig_m)&lvig_d>0&lvig_m>0, lvig_d*lvig_m*8, 0),
    leis_mod_met = ifelse(!is.na(lmod_d)&!is.na(lmod_m)&lmod_d>0&lmod_m>0, lmod_d*lmod_m*4, 0),
    work_vig_met = ifelse(!is.na(wvig_d)&!is.na(wvig_m)&wvig_d>0&wvig_m>0, wvig_d*wvig_m*8, 0),
    work_mod_met = ifelse(!is.na(wmod_d)&!is.na(wmod_m)&wmod_d>0&wmod_m>0, wmod_d*wmod_m*4, 0),
    trans_met    = ifelse(!is.na(trans_d)&!is.na(trans_m)&trans_d>0&trans_m>0, trans_d*trans_m*4, 0),
    paq_leisure_mets = leis_vig_met + leis_mod_met,
    paq_total_mets   = paq_leisure_mets + work_vig_met + work_mod_met + trans_met,
    paq_active_unified = factor(case_when(
      !is.na(paq_leisure_mets) & paq_leisure_mets >= 600 ~ "Active",
      !is.na(paq_leisure_mets) & paq_leisure_mets < 600 ~ "Insufficient",
      cycle=="2005-2006" & paq_any_mv==1 ~ "Active", cycle=="2005-2006" & paq_any_mv==0 ~ "Insufficient",
      TRUE ~ NA_character_), levels = c("Insufficient","Active")),
    
    dm_diag = ifelse(!is.na(DIQ010) & DIQ010==1, 1L, 0L),
    dm_med  = ifelse(!is.na(DIQ050) & DIQ050==1, 1L, 0L),
    dm_fpg  = ifelse(!is.na(glucose) & glucose >= 126, 1L, 0L),
    dm_hba1c= ifelse(!is.na(hba1c) & hba1c >= 6.5, 1L, 0L),
    dm = case_when(dm_diag==1|dm_med==1|dm_fpg==1|dm_hba1c==1 ~ 1L,
                   (!is.na(dm_diag)|!is.na(dm_med)|!is.na(dm_fpg)|!is.na(dm_hba1c)) ~ 0L, TRUE ~ NA_integer_),
    dm_status = factor(dm, levels=c(0,1), labels=c("Non-diabetes","Diabetes")),
    
    lncreat_umol = log(creatinine * 88.42 + 1),
    lncrp_mgdl   = log(crp / 10 + 1),
    albumin_gL   = albumin * 10,
    creat_umol   = creatinine * 88.42,
    crp_mgdl     = crp / 10,
    sex_male     = ifelse(sex == "Male", 1, 0),
    tyg_tert = factor(ntile(tyg, 3), levels=1:3, labels=c("T1 (low)","T2 (mid)","T3 (high)")),
    tyg_z = as.numeric(scale(tyg))
  )
  
  N_CYC <- length(CYCLES)
  
  # ---- 诊断 ----
  diag_vars(nhanes, c("age","tyg","homa_ir","crp","wtsfa"), "关键变量")
  log_write("\n--- 按周期诊断 ---")
  print(nhanes %>% group_by(cycle) %>% summarise(
    n=n(), fasting_n=sum(fasting==1,na.rm=TRUE),
    crp_n=sum(!is.na(crp)), tyg_n=sum(!is.na(tyg))))
  
  # ---- 缓存 ----
  saveRDS(nhanes, CACHE_RAW_DERIVED)
  log_write(sprintf("✓ 已保存缓存: %s", basename(CACHE_RAW_DERIVED)))
  cp_save(nhanes, "10", "main")
  timer_end("模块10 数据导入与变量衍生")
  log_write("===== 模块 10 结束 =====\n")
}