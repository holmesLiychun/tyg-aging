# ============================================================================
# 模块 70：非线性与可视化（全算法版 + PNG/PDF 双导出）
# ★ 标签更新：KDM-9→KDM+CRP, KDM-11→KDM+TG, PhenoAge-9→PhenoAge-original
# ============================================================================

timer_start("模块70 非线性与可视化")

library(ggplot2)
library(patchwork)

# ============================================================================
# ★ 标签映射（统一替换，方便后续维护）
# ============================================================================
algo_label_map <- c(
  "KDM-9"       = "KDM+CRP",
  "KDM-11"      = "KDM+TG",
  "PhenoAge-9"  = "PhenoAge-original"
)

# ============================================================================
# ★ 双格式保存封装函数
# ============================================================================
save_fig <- function(p, name, width = 8, height = 6, dpi = 300,
                     dir = file.path(OUT_DIR, "figures"),
                     tiff = FALSE) {
  if (is.null(p)) return(invisible(NULL))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  
  ggsave(file.path(dir, paste0(name, ".png")), plot = p,
         width = width, height = height, units = "in", dpi = dpi)
  ggsave(file.path(dir, paste0(name, ".pdf")), plot = p,
         width = width, height = height, units = "in")
  
  if (tiff) {
    tiff(file.path(dir, paste0(name, ".tiff")),
         width = width, height = height, units = "in",
         res = 600, compression = "lzw")
    print(p)
    dev.off()
  }
  log_write(sprintf("  ✓ 已保存 %s (.png/.pdf)", name))
}

# ---- 70.1 RCS 绘图函数（通用）-----------------------------------------------
draw_rcs <- function(design_obj, data_obj, endpoint_var, title_str, out_name,
                     extra_covs = NULL) {
  log_write(sprintf(">>> 拟合 %s RCS...", title_str))
  
  cov_base <- c("age", "sex", "race", "bmi", "smoke", "alcohol",
                "paq_active_unified", "edu_cat")
  covs <- c(cov_base, extra_covs)
  covs <- covs[covs %in% names(data_obj)]
  
  form_rcs <- as.formula(paste0(endpoint_var, " ~ splines::ns(tyg, df=3) + ",
                                paste(covs, collapse = " + ")))
  fit_rcs <- tryCatch(
    survey::svyglm(form_rcs, design = design_obj, family = gaussian()),
    error = function(e) {
      log_write(sprintf("  ERROR RCS 拟合失败: %s", e$message))
      NULL
    })
  if (is.null(fit_rcs)) return(NULL)
  
  log_write("  ✓ 模型拟合成功")
  
  p_overall <- survey::regTermTest(fit_rcs, ~splines::ns(tyg, df = 3))
  p_nonlin  <- survey::regTermTest(fit_rcs, ~splines::ns(tyg, df = 3) - tyg)
  
  fmt_p <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("<0.001")
    if (p < 0.01) return(sprintf("%.3f", p))
    sprintf("%.2f", p)
  }
  log_write(sprintf("  Overall P=%s, Non-linear P=%s",
                    fmt_p(p_overall$p), fmt_p(p_nonlin$p)))
  
  cont_grid <- c("age", "bmi")
  cont_grid <- cont_grid[cont_grid %in% names(data_obj)]
  cat_grid  <- covs[!covs %in% c(cont_grid, "tyg")]
  
  rng <- quantile(data_obj$tyg, c(0.05, 0.95), na.rm = TRUE)
  pred_data <- data.frame(tyg = seq(rng[1], rng[2], length.out = 100))
  for (v in cont_grid) pred_data[[v]] <- median(data_obj[[v]], na.rm = TRUE)
  for (v in cat_grid) {
    if (is.factor(data_obj[[v]])) {
      tab <- table(data_obj[[v]])
      pred_data[[v]] <- factor(names(tab)[which.max(tab)],
                               levels = levels(data_obj[[v]]))
    }
  }
  
  mm <- model.matrix(delete.response(terms(fit_rcs)), pred_data)
  pred_data$fit   <- as.vector(mm %*% coef(fit_rcs))
  pred_data$se    <- sqrt(diag(mm %*% vcov(fit_rcs) %*% t(mm)))
  pred_data$lower <- pred_data$fit - 1.96 * pred_data$se
  pred_data$upper <- pred_data$fit + 1.96 * pred_data$se
  
  log_write(sprintf("  fit range: [%.3f, %.3f]",
                    min(pred_data$fit), max(pred_data$fit)))
  
  p <- ggplot(pred_data, aes(x = tyg, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.25) +
    geom_line(color = "steelblue", size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "black", size = 0.5) +
    labs(
      title = title_str,
      x = "TyG Index", y = "Biological Age Acceleration (years)",
      subtitle = sprintf("Overall P=%s | Non-linear P=%s",
                         fmt_p(p_overall$p), fmt_p(p_nonlin$p))) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkred"))
  
  save_fig(p, paste0("RCS_", out_name), width = 8, height = 6)
  
  list(fit = fit_rcs, p = p, p_overall = p_overall$p, p_nonlin = p_nonlin$p)
}

# ============================================================================
# 70.2 主分析 RCS（7 算法，KDM_core6 排第一）
# ============================================================================
log_write(">>> 70.2 主分析 RCS（7 算法）")

rcs_results <- list()
rcs_plots <- list()

rcs_config <- list(
  "KDM_core6"          = list(ds = "ds_kdmcore6", data = "ana_kdmcore6", ep = "kdm_core6_accel",
                              label = "TyG vs KDM_core6 BA Acceleration (7 cycles)", out = "KDMcore6"),
  "KDM+TG"             = list(ds = "ds_kdm",    data = "ana_kdm",    ep = "kdm_accel",
                              label = "TyG vs KDM+TG BA Acceleration", out = "KDM11"),
  "KDM+CRP"            = list(ds = "ds_kdm9",   data = "ana_kdm9",   ep = "kdm9_accel",
                              label = "TyG vs KDM+CRP BA Acceleration", out = "KDM9"),
  "PhenoAge-original"  = list(ds = "ds_pheno9", data = "ana_pheno9", ep = "pheno9_accel",
                              label = "TyG vs PhenoAge-original Acceleration", out = "Pheno9"),
  "Modified PhenoAge"  = list(ds = "ds_pheno",  data = "ana_pheno",  ep = "pheno_accel",
                              label = "TyG vs Modified PhenoAge Acceleration", out = "PhenoAge"),
  "Pheno-noG"          = list(ds = "ds_pheno_nog", data = "ana_pheno_nog", ep = "pheno_nog_accel",
                              label = "TyG vs Pheno-noG Acceleration", out = "PhenoNoG"),
  "Pheno-noCRP"        = list(ds = "ds_phenocrp",  data = "ana_phenocrp",  ep = "phenocrp_accel",
                              label = "TyG vs Pheno-noCRP Acceleration", out = "PhenoNoCRP")
)

for (nm in names(rcs_config)) {
  cfg <- rcs_config[[nm]]
  if (!exists(cfg$ds) || !exists(cfg$data)) {
    log_write(sprintf("  ⚠ %s: %s 或 %s 不存在，跳过", nm, cfg$ds, cfg$data))
    next
  }
  ds_obj <- get(cfg$ds)
  data_obj <- get(cfg$data)
  if (!cfg$ep %in% names(data_obj)) {
    log_write(sprintf("  ⚠ %s: %s 中无 %s 列，跳过", nm, cfg$data, cfg$ep))
    next
  }
  rcs_results[[nm]] <- draw_rcs(ds_obj, data_obj, cfg$ep, cfg$label, cfg$out)
  if (!is.null(rcs_results[[nm]]$p)) {
    rcs_plots[[nm]] <- rcs_results[[nm]]$p
  }
}

log_write("========== RCS P 值汇总 ==========")
fmt <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
for (nm in names(rcs_results)) {
  if (!is.null(rcs_results[[nm]])) {
    log_write(sprintf("%-22s: Overall P=%s, Non-linear P=%s",
                      nm, fmt(rcs_results[[nm]]$p_overall),
                      fmt(rcs_results[[nm]]$p_nonlin)))
  }
}

# ============================================================================
# 70.2c ★ 合并 7 张 RCS 图为 PNG + PDF（4×2 布局）
# ============================================================================
if (length(rcs_plots) >= 2) {
  plot_order <- c("KDM_core6", "KDM+TG", "KDM+CRP", "PhenoAge-original",
                  "Modified PhenoAge", "Pheno-noG", "Pheno-noCRP")
  plots_to_combine <- rcs_plots[plot_order[plot_order %in% names(rcs_plots)]]
  
  if (length(plots_to_combine) >= 6) {
    log_write(sprintf(">>> 合并 %d 张 RCS 图 (PNG+PDF)...", length(plots_to_combine)))
    
    combined_plot <- wrap_plots(plots_to_combine, ncol = 4, nrow = 2,
                                guides = "collect") &
      theme(legend.position = "bottom",
            plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 8))
    
    save_fig(combined_plot, "Figure_RCS_7in1", width = 14, height = 7)
  } else {
    log_write(sprintf("  ⚠ 仅收集到 %d 张图，跳过合并", length(plots_to_combine)))
  }
}

# ============================================================================
# 70.2b 拐点计算
# ============================================================================
log_write(">>> 70.2b 拐点计算（数值导数法）")

inflection_results <- list()

for (nm in names(rcs_results)) {
  if (is.null(rcs_results[[nm]]$fit)) {
    log_write(sprintf("  %s: 无模型对象，跳过", nm))
    next
  }
  
  fit <- rcs_results[[nm]]$fit
  cfg <- rcs_config[[nm]]
  data_obj <- get(cfg$data)
  
  cov_base <- c("age", "sex", "race", "bmi", "smoke", "alcohol",
                "paq_active_unified", "edu_cat")
  covs <- cov_base[cov_base %in% names(data_obj)]
  cont_grid <- c("age", "bmi")
  cont_grid <- cont_grid[cont_grid %in% names(data_obj)]
  cat_grid  <- covs[!covs %in% c(cont_grid, "tyg")]
  
  rng <- quantile(data_obj$tyg, c(0.05, 0.95), na.rm = TRUE)
  n_grid <- 2000
  pred_grid <- data.frame(tyg = seq(rng[1], rng[2], length.out = n_grid))
  for (v in cont_grid) pred_grid[[v]] <- median(data_obj[[v]], na.rm = TRUE)
  for (v in cat_grid) {
    if (is.factor(data_obj[[v]])) {
      tab <- table(data_obj[[v]])
      pred_grid[[v]] <- factor(names(tab)[which.max(tab)],
                               levels = levels(data_obj[[v]]))
    }
  }
  
  mm <- model.matrix(delete.response(terms(fit)), pred_grid)
  pred_grid$fit <- as.vector(mm %*% coef(fit))
  
  h <- pred_grid$tyg[2] - pred_grid$tyg[1]
  deriv <- c(NA, (pred_grid$fit[3:n_grid] - pred_grid$fit[1:(n_grid-2)]) / (2*h), NA)
  pred_grid$deriv <- deriv
  
  valid <- !is.na(deriv)
  idx_infl <- which.min(abs(deriv[valid]))
  infl_tyg <- pred_grid$tyg[valid][idx_infl]
  infl_deriv <- deriv[valid][idx_infl]
  infl_fit <- pred_grid$fit[valid][idx_infl]
  
  left_mask  <- pred_grid$tyg <= infl_tyg & pred_grid$tyg >= rng[1]
  right_mask <- pred_grid$tyg >= infl_tyg & pred_grid$tyg <= rng[2]
  
  if (sum(left_mask, na.rm = TRUE) > 5) {
    left_lm <- lm(fit ~ tyg, data = pred_grid[left_mask, ])
    left_slope <- coef(left_lm)[2]; left_ci <- confint(left_lm)[2, ]
  } else { left_slope <- NA; left_ci <- c(NA, NA) }
  
  if (sum(right_mask, na.rm = TRUE) > 5) {
    right_lm <- lm(fit ~ tyg, data = pred_grid[right_mask, ])
    right_slope <- coef(right_lm)[2]; right_ci <- confint(right_lm)[2, ]
  } else { right_slope <- NA; right_ci <- c(NA, NA) }
  
  pct_rank <- mean(data_obj$tyg <= infl_tyg, na.rm = TRUE)
  boundary_flag <- ifelse(pct_rank < 0.10 || pct_rank > 0.90, " ⚠ 靠近边界", "")
  
  log_write(sprintf(
    "  %s: 拐点 TyG=%.3f (%.1f%% 分位)%s | 左斜率=%.3f | 右斜率=%.3f",
    nm, infl_tyg, pct_rank * 100, boundary_flag,
    ifelse(is.na(left_slope), NA, round(left_slope, 3)),
    ifelse(is.na(right_slope), NA, round(right_slope, 3))))
  
  inflection_results[[nm]] <- data.frame(
    Algorithm = nm,
    Inflection_TyG = round(infl_tyg, 3),
    Percentile = round(pct_rank * 100, 1),
    Left_Slope = round(left_slope, 4),
    Left_CI_lower = round(left_ci[1], 4), Left_CI_upper = round(left_ci[2], 4),
    Right_Slope = round(right_slope, 4),
    Right_CI_lower = round(right_ci[1], 4), Right_CI_upper = round(right_ci[2], 4),
    Fit_at_infl = round(infl_fit, 4), stringsAsFactors = FALSE
  )
}

if (length(inflection_results) > 0) {
  infl_df <- do.call(rbind, inflection_results)
  write.csv(infl_df, file.path(OUT_DIR, "table_rcs_inflection_points.csv"), row.names = FALSE)
  log_write("✓ 拐点汇总表已保存")
}

# ============================================================================
# 70.3 敏感性 RCS 对比图（KDM + PhenoAge）
# ============================================================================
log_write(">>> 70.3 敏感性 RCS 对比图")

d5 <- NULL
ds_5full <- NULL

if (exists("ana_kdm9") && is.data.frame(ana_kdm9) && nrow(ana_kdm9) > 0) {
  d5 <- ana_kdm9
  log_write(sprintf("  基底: ana_kdm9, n=%d", nrow(d5)))
  
  if (exists("ana_kdmcore6") && "kdm_core6_accel" %in% names(ana_kdmcore6)) {
    d5 <- merge(d5, ana_kdmcore6[, c("SEQN", "kdm_core6_accel")], by = "SEQN", all.x = TRUE)
    log_write("  ✓ kdm_core6_accel 已合并")
  }
  
  inflam_map <- list()
  if ("wbc" %in% names(d5)) inflam_map$wbc <- "wbc"
  else if ("LBXWBCSI" %in% names(d5)) { d5$wbc <- d5$LBXWBCSI; inflam_map$wbc <- "wbc" }
  if ("lymph" %in% names(d5)) inflam_map$lymph <- "lymph"
  else if ("LBXLYPCT" %in% names(d5)) { d5$lymph <- d5$LBXLYPCT; inflam_map$lymph <- "lymph" }
  if ("lncrp_mgdl" %in% names(d5)) inflam_map$lncrp <- "lncrp_mgdl"
  else if ("lncrp" %in% names(d5)) inflam_map$lncrp <- "lncrp"
  else if ("LBXCRP" %in% names(d5)) { d5$lncrp_mgdl <- d5$LBXCRP; inflam_map$lncrp <- "lncrp_mgdl" }
  
  covars_needed <- c("pir", "sbp", "dbp", "dm_status", "psu", "strata", "weight_fasting")
  for (v in covars_needed) {
    if (!v %in% names(d5) && exists("nhanes") && v %in% names(nhanes)) {
      d5 <- merge(d5, nhanes[, c("SEQN", v)], by = "SEQN", all.x = TRUE)
    }
  }
  
  required_all <- c("tyg", "kdm_core6_accel", "kdm_accel", "kdm9_accel",
                    "pheno9_accel", "pheno_accel", "pheno_nog_accel", "phenocrp_accel",
                    "age", "sex", "race", "edu_cat", "pir", "bmi", "sbp", "dbp",
                    "smoke", "alcohol", "paq_active_unified", "dm_status",
                    "wbc", "lymph", "lncrp_mgdl", "weight_fasting")
  required_present <- intersect(required_all, names(d5))
  d5 <- d5[complete.cases(d5[, required_present]), ]
  log_write(sprintf("  完整案例: n=%d", nrow(d5)))
  
  if (all(c("psu", "strata") %in% names(d5))) {
    ds_5full <- survey::svydesign(id = ~psu, strata = ~strata,
                                  weights = ~weight_fasting, nest = TRUE, data = d5)
  } else {
    ds_5full <- survey::svydesign(id = ~1, weights = ~weight_fasting, nest = TRUE, data = d5)
  }
} else {
  log_write("  ⚠ ana_kdm9 不存在，跳过 70.3")
}

p_kdm <- NULL
p_pheno <- NULL

if (!is.null(d5) && !is.null(ds_5full) && nrow(d5) > 100) {
  fmt_p <- function(p) { if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p) }
  
  cov_base_formula <- "age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"
  cov_inflam_formula <- paste(c("wbc", "lymph", "lncrp_mgdl"), collapse = " + ")
  
  fit_kdm11_s <- tryCatch(survey::svyglm(as.formula(paste("kdm_accel ~ splines::ns(tyg, df=3) +", cov_base_formula)),
                                         design = ds_5full, family = gaussian()), error = function(e) NULL)
  fit_kdm9_s  <- tryCatch(survey::svyglm(as.formula(paste("kdm9_accel ~ splines::ns(tyg, df=3) +", cov_base_formula, "+", cov_inflam_formula)),
                                         design = ds_5full, family = gaussian()), error = function(e) NULL)
  fit_p9_s    <- tryCatch(survey::svyglm(as.formula(paste("pheno9_accel ~ splines::ns(tyg, df=3) +", cov_base_formula, "+", cov_inflam_formula)),
                                         design = ds_5full, family = gaussian()), error = function(e) NULL)
  fit_pp_s    <- tryCatch(survey::svyglm(as.formula(paste("pheno_accel ~ splines::ns(tyg, df=3) +", cov_base_formula)),
                                         design = ds_5full, family = gaussian()), error = function(e) NULL)
  fit_png_s   <- tryCatch(survey::svyglm(as.formula(paste("pheno_nog_accel ~ splines::ns(tyg, df=3) +", cov_base_formula)),
                                         design = ds_5full, family = gaussian()), error = function(e) NULL)
  
  tq <- quantile(d5$tyg, c(0.05, 0.95), na.rm = TRUE)
  tyg_seq <- seq(tq[1], tq[2], length.out = 200)
  pred_grid <- data.frame(tyg = tyg_seq)
  for (v in c("age", "bmi", "sbp", "dbp", "pir", "wbc", "lymph", "lncrp_mgdl"))
    pred_grid[[v]] <- if (v %in% names(d5)) median(d5[[v]], na.rm = TRUE) else 0
  for (v in c("sex", "race", "edu_cat", "smoke", "alcohol", "paq_active_unified", "dm_status")) {
    if (v %in% names(d5)) {
      if (is.factor(d5[[v]])) {
        tab <- table(d5[[v]]); pred_grid[[v]] <- factor(names(tab)[which.max(tab)], levels = levels(d5[[v]]))
      } else pred_grid[[v]] <- median(d5[[v]], na.rm = TRUE)
    }
  }
  
  # ---- KDM 对比图（括号描述保留原样）----
  if (!is.null(fit_kdm11_s) && !is.null(fit_kdm9_s)) {
    X1 <- model.matrix(delete.response(fit_kdm11_s$terms), data = pred_grid)
    X2 <- model.matrix(delete.response(fit_kdm9_s$terms), data = pred_grid)
    se1 <- sqrt(diag(X1 %*% vcov(fit_kdm11_s) %*% t(X1)))
    se2 <- sqrt(diag(X2 %*% vcov(fit_kdm9_s) %*% t(X2)))
    
    rcs_kdm_df <- rbind(
      data.frame(tyg = tyg_seq, fit = as.vector(X1 %*% coef(fit_kdm11_s)),
                 lower = as.vector(X1 %*% coef(fit_kdm11_s)) - 1.96*se1,
                 upper = as.vector(X1 %*% coef(fit_kdm11_s)) + 1.96*se1,
                 Model = factor("KDM+TG (no CRP)", levels = c("KDM+TG (no CRP)", "KDM+CRP (with CRP)"))),
      data.frame(tyg = tyg_seq, fit = as.vector(X2 %*% coef(fit_kdm9_s)),
                 lower = as.vector(X2 %*% coef(fit_kdm9_s)) - 1.96*se2,
                 upper = as.vector(X2 %*% coef(fit_kdm9_s)) + 1.96*se2,
                 Model = factor("KDM+CRP (with CRP)", levels = c("KDM+TG (no CRP)", "KDM+CRP (with CRP)")))
    )
    
    p_kdm <- ggplot(rcs_kdm_df, aes(x = tyg, y = fit, color = Model, fill = Model)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
      geom_line(size = 1.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      scale_color_manual(values = c("#0072B2", "#D55E00")) +
      scale_fill_manual(values = c("#0072B2", "#D55E00")) +
      labs(title = "TyG Index and KDM Biological Age Acceleration",
           subtitle = sprintf("RCS (df=3); N=%d; Model 5 adjusted", nrow(d5)),
           x = "TyG Index", y = "Biological Age Acceleration (years)") +
      theme_classic(base_size = 12) +
      theme(legend.position = c(0.80, 0.15),
            legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3))
    
    save_fig(p_kdm, "RCS_kdm11_vs_kdm9", width = 8, height = 5.3)
  }
  
  # ---- PhenoAge 对比图（括号描述保留原样）----
  pheno_fits <- list("PhenoAge-original (with CRP)" = fit_p9_s,
                     "Modified PhenoAge (no CRP)" = fit_pp_s,
                     "Pheno-noG (no glucose/trig)" = fit_png_s)
  valid_pheno <- pheno_fits[!sapply(pheno_fits, is.null)]
  
  if (length(valid_pheno) >= 2) {
    rcs_pheno_list <- list()
    for (i in seq_along(valid_pheno)) {
      fit_i <- valid_pheno[[i]]
      X_i <- model.matrix(delete.response(fit_i$terms), data = pred_grid)
      se_i <- sqrt(diag(X_i %*% vcov(fit_i) %*% t(X_i)))
      rcs_pheno_list[[names(valid_pheno)[i]]] <- data.frame(
        tyg = tyg_seq, fit = as.vector(X_i %*% coef(fit_i)),
        lower = as.vector(X_i %*% coef(fit_i)) - 1.96*se_i,
        upper = as.vector(X_i %*% coef(fit_i)) + 1.96*se_i,
        Model = factor(names(valid_pheno)[i], levels = names(valid_pheno)))
    }
    rcs_pheno_df <- do.call(rbind, rcs_pheno_list)
    
    p_pheno <- ggplot(rcs_pheno_df, aes(x = tyg, y = fit, color = Model, fill = Model)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
      geom_line(size = 1.2) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
      scale_color_brewer(palette = "Set1") +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "TyG Index and Phenotypic Age Acceleration",
           subtitle = sprintf("RCS (df=3); N=%d; Model 5 adjusted", nrow(d5)),
           x = "TyG Index", y = "Phenotypic Age Acceleration (years)") +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom")
    
    save_fig(p_pheno, "RCS_phenoage_comparison", width = 9, height = 5.5)
  }
} else {
  log_write("  ⚠ d5 或 ds_5full 构造失败，跳过 70.3")
}

# ============================================================================
# 70.3b ★ KDM + PhenoAge 并排合并 PNG + PDF
# ============================================================================
if (!is.null(p_kdm) && !is.null(p_pheno)) {
  log_write(">>> 合并 KDM + PhenoAge 对比图 (PNG+PDF)...")
  combined_two <- wrap_plots(p_kdm, p_pheno, ncol = 2, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9))
  
  save_fig(combined_two, "Figure_RCS_KDM_Phenoage_combined", width = 14, height = 5.5)
}

# ============================================================================
# 70.4 森林图（PNG + PDF 双导出，标签映射替换）
# ============================================================================
log_write(">>> 70.4 森林图")

forest_data <- NULL

if (exists("tyg_assoc") && is.data.frame(tyg_assoc) && nrow(tyg_assoc) > 0) {
  req <- c("Algorithm", "Contrast", "Beta", "SE", "P")
  if (all(req %in% names(tyg_assoc))) {
    forest_data <- tyg_assoc %>%
      filter(Contrast %in% c("T2 vs T1", "T3 vs T1")) %>%
      mutate(
        Contrast_clean = Contrast,
        Lower = Beta - 1.96 * SE, Upper = Beta + 1.96 * SE,
        Sig = case_when(
          P < 0.001 ~ "***", P < 0.01 ~ "**", P < 0.05 ~ "*", TRUE ~ ""),
        Effect_CI = sprintf("%.2f (%.2f, %.2f) %s", Beta, Lower, Upper, Sig),
        Endpoint = recode(Algorithm, !!!algo_label_map, .default = Algorithm)
      ) %>%
      dplyr::select(Endpoint, Contrast_clean, Beta, Lower, Upper, Effect_CI)
  } else log_write("  ⚠ tyg_assoc 缺少必要列")
} else log_write("  ⚠ tyg_assoc 不存在或为空")

if (!is.null(forest_data) && nrow(forest_data) > 0) {
  forest_data <- forest_data %>%
    mutate(
      Contrast_clean = factor(Contrast_clean, levels = c("T2 vs T1", "T3 vs T1")),
      Endpoint = factor(Endpoint, levels = unique(Endpoint))) %>%
    arrange(Endpoint, Contrast_clean)
  
  dt <- data.frame(
    Outcome    = as.character(forest_data$Endpoint),
    Comparison = as.character(forest_data$Contrast_clean),
    ` `        = rep(" ", nrow(forest_data)),
    `Beta (95% CI)` = forest_data$Effect_CI,
    stringsAsFactors = FALSE, check.names = FALSE)
  
  xmin <- floor(min(forest_data$Lower, na.rm = TRUE))
  xmax <- ceiling(max(forest_data$Upper, na.rm = TRUE))
  tk   <- seq(xmin, xmax, by = 2)
  if (length(tk) == 0) tk <- c(xmin, xmax)
  
  fp <- forestploter::forest(dt, est = forest_data$Beta, lower = forest_data$Lower,
                             upper = forest_data$Upper, ci_column = 3,
                             ref_line = 0, col_widths = c(2.5, 2, 6, 3))
  fp <- forestploter::set_xaxis(fp, xlim = c(xmin, xmax), ticks_at = tk)
  fp <- forestploter::set_labs(fp,
                               title = "TyG Tertile and Biological Age Acceleration\nT2 vs T1 and T3 vs T1, Model 5 adjusted",
                               xlab = "Beta (years)")
  fp <- forestploter::set_style(fp,
                                core = list(padding = unit(c(5, 3), "mm"), gp = gpar(fontsize = 8)),
                                colhead = list(padding = unit(c(25, 3), "mm"), gp = gpar(fontsize = 8, fontface = "bold")))
  
  png(file.path(OUT_DIR, "figures", "forest_plot_tyg_ba.png"),
      width = 2000, height = 1000, res = 150)
  grid.newpage(); grid.draw(fp); dev.off()
  
  pdf(file.path(OUT_DIR, "figures", "forest_plot_tyg_ba.pdf"),
      width = 13, height = 7)
  grid.newpage(); grid.draw(fp); dev.off()
  
  log_write("✓ 森林图已保存: .png / .pdf")
} else log_write("  无可绘制的森林图数据")

# ============================================================================
# 70.5 解耦矩阵热图 + β 轨迹图（★ Outcome 标签同步替换）
# ============================================================================
log_write(">>> 70.5 解耦矩阵热图与 β 轨迹图")

if (exists("mat_final") && is.data.frame(mat_final) && nrow(mat_final) > 0) {
  
  # ★ 对 mat_final$Outcome 做标签映射
  mat_labeled <- mat_final %>%
    mutate(Outcome = recode(Outcome, !!!algo_label_map, .default = Outcome))
  
  heat_data <- mat_labeled %>%
    filter(Model == "Model 5") %>%
    dplyr::select(Exposure, Outcome, Beta, P_value) %>%
    mutate(
      Exposure_short = gsub(" \\(per SD\\)", "", Exposure),
      sig = case_when(
        P_value < 0.001 ~ "***", P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*", TRUE ~ ""),
      label = sprintf("%.1f%s", Beta * 100, sig))
  
  p_heat <- ggplot(heat_data, aes(x = Exposure_short, y = Outcome, fill = Beta)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), color = "black", size = 3.5) +
    scale_fill_gradient2(low = "#2c7fb8", mid = "#f7f7f7", high = "#d95f0e",
                         midpoint = 0, name = expression(beta~"(per SD)")) +
    labs(title = "Decoupling Matrix: TyG-related Exposures x Biological Age Acceleration",
         subtitle = "Model 5 (fully adjusted); beta = years per 1 SD",
         x = "Exposure", y = "Outcome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 13, face = "bold"))
  
  save_fig(p_heat, "decoupling_heatmap", width = 9, height = 6)
  
  # ★ 轨迹图同样用映射后的 Outcome
  traj_data <- mat_labeled %>%
    filter(ExposureVar %in% c("tyg_z", "tyg", "TyG", "tyg_sd")) %>%
    mutate(Model_num = as.numeric(gsub("Model ", "", Model)))
  
  if (nrow(traj_data) > 0) {
    p_traj <- ggplot(traj_data, aes(x = Model_num, y = Beta, color = Outcome, group = Outcome)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1, alpha = 0.6) +
      scale_x_continuous(breaks = sort(unique(traj_data$Model_num)),
                         labels = paste0("M", sort(unique(traj_data$Model_num)))) +
      labs(title = "Effect Trajectory: TyG -> Biological Age Acceleration",
           subtitle = "Beta across Model 1 (crude) to Model 5 (fully adjusted)",
           x = "Adjustment Model", y = expression(beta~"(per SD)")) +
      theme_minimal() + theme(legend.title = element_blank())
    
    save_fig(p_traj, "tyg_effect_trajectory", width = 8, height = 5)
  }
} else log_write("  解耦矩阵数据不存在，跳过热图与轨迹图")

# ============================================================================
# 70.6 线性模型汇总（7 算法，标签更新）
# ============================================================================
log_write(">>> 70.6 线性模型汇总（全算法）")

linear_summary_rows <- list()

if (exists("ds_5full") && !is.null(ds_5full) && exists("d5") && !is.null(d5)) {
  algo_linear <- list(
    "KDM_core6"          = "kdm_core6_accel",
    "KDM+TG"             = "kdm_accel",
    "KDM+CRP"            = "kdm9_accel",
    "PhenoAge-original"  = "pheno9_accel",
    "Modified PhenoAge"  = "pheno_accel",
    "Pheno-noG"          = "pheno_nog_accel",
    "Pheno-noCRP"        = "phenocrp_accel"
  )
  
  for (nm in names(algo_linear)) {
    ep <- algo_linear[[nm]]
    if (!ep %in% names(d5)) { log_write(sprintf("  %s: %s 不在 d5，跳过", nm, ep)); next }
    
    form_lin <- as.formula(paste(
      ep, "~ tyg + age + sex + race + edu_cat + pir + bmi + sbp + dbp + smoke + alcohol + paq_active_unified + dm_status"))
    fit_lin <- tryCatch(survey::svyglm(form_lin, design = ds_5full, family = gaussian()),
                        error = function(e) NULL)
    if (is.null(fit_lin)) next
    
    sm <- summary(fit_lin)$coefficients
    if (!"tyg" %in% rownames(sm)) next
    tr <- sm["tyg", ]; pv <- tr["Pr(>|t|)"]
    
    log_write(sprintf("  %s: TyG Beta=%.3f, P=%s", nm, tr["Estimate"],
                      if (pv < 0.001) "<0.001" else sprintf("%.3f", pv)))
    linear_summary_rows[[nm]] <- data.frame(
      Algorithm = nm, Endpoint = ep,
      TyG_Beta = round(tr["Estimate"], 4), TyG_SE = round(tr["Std. Error"], 4),
      TyG_P = pv, stringsAsFactors = FALSE)
  }
}

if (length(linear_summary_rows) > 0) {
  linear_summary_df <- do.call(rbind, linear_summary_rows)
  linear_summary_df$TyG_P_formatted <- ifelse(linear_summary_df$TyG_P < 0.001, "<0.001",
                                              round(linear_summary_df$TyG_P, 3))
  write.csv(linear_summary_df, file.path(OUT_DIR, "table_linear_tyg_summary.csv"), row.names = FALSE)
  log_write("✓ 线性汇总表已保存")
}

# 70.7
if (exists("nhanes")) {
  saveRDS(nhanes, file.path(OUT_DIR, "nhanes_bioage_final_v10.rds"))
  log_write("✓ 最终结果 RDS 已保存")
}

timer_end("模块70 非线性与可视化")
log_write("===== 模块 70 结束 =====\n")