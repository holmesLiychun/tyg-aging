# ============================================================
# R/ba_algorithms.R — KDM(校准) + SPA(残差) 四版
# ============================================================

# ---------- KDM 校准法（分性别，按 KD 公式）----------
calc_kdm <- function(d, markers, sex_specific = TRUE, label = "") {
  out <- rep(NA_real_, nrow(d))
  sexes <- if (sex_specific) sort(unique(d$gender)) else "pool"
  
  for (g in sexes) {
    if (g != "pool") idx <- which(d$gender == g) else idx <- 1:nrow(d)
    sub <- d[idx, , drop = FALSE]
    req <- c("age", markers)
    ok <- complete.cases(sub[, req])
    if (sum(ok) < 50) { CH$log("WARN", sprintf("%s KDM %s n=%d skip", label, g, sum(ok)), "4_BA"); next }
    
    sg <- sub[ok, , drop = FALSE]
    agev <- sg$age
    # 每个 marker: age ~ x  -> k截距,q斜率,s=rmse
    k_i <- q_i <- s_i <- numeric(length(markers)); names(k_i) <- names(q_i) <- names(s_i) <- markers
    for (m in markers) {
      fit <- lm(agev ~ sg[[m]])
      k_i[m] <- coef(fit)[1]; q_i[m] <- coef(fit)[2]
      s_i[m] <- sigma(fit)
      if (s_i[m] == 0 || is.na(s_i[m]) || is.na(q_i[m])) s_i[m] <- NA
    }
    use <- markers[!is.na(s_i) & s_i > 0 & !is.na(q_i)]
    if (length(use) < 1) next
    X <- as.matrix(sg[, use, drop = FALSE])
    # 全模型 rmse -> s_BA
    full <- lm(agev ~ X); s_ba2 <- var(residuals(full), na.rm = TRUE)
    if (!is.finite(s_ba2) || s_ba2 <= 0) s_ba2 <- var(agev, na.rm = TRUE)
    
    num <- den <- 0
    for (m in use) {
      wi <- (q_i[m]^2)/(s_i[m]^2)
      age_hat <- k_i[m] + q_i[m]*X[, m]
      num <- num + wi*age_hat
      den <- den + (k_i[m]/s_i[m]^2)   # KD 标准分母项
    }
    # 标准 KD（按 age~x 参数）：BA = (Σ k/s^2 + CA/s_BA^2)/(Σ (q/s^2? 见下注) + 1/s_BA^2)
    # 实操用加权年龄估计版本，稳定且不溢出：
    num2 <- 0; den2 <- 0
    for (m in use) {
      wi <- (q_i[m]^2)/(s_i[m]^2)
      num2 <- num2 + wi*(k_i[m] + q_i[m]*X[, m])
      den2 <- den2 + wi
    }
    ba <- (num2 + agev/s_ba2)/(den2 + 1/s_ba2)
    out[idx[ok]] <- ba
  }
  CH$log("OK", sprintf("%s KDM done: %d valid", label, sum(!is.na(out))), "4_BA")
  out
}

# ---------- SPA 残差法（age - predicted_age）----------
calc_spa <- function(d, markers, sex_specific = FALSE, label = "") {
  out <- rep(NA_real_, nrow(d))
  train_r2 <- NA_real_
  sexes <- if (sex_specific) sort(unique(d$gender)) else "pool"
  for (g in sexes) {
    if (g != "pool") idx <- which(d$gender == g) else idx <- 1:nrow(d)
    sub <- d[idx, , drop = FALSE]
    req <- c("age", markers)
    ok <- complete.cases(sub[, req])
    if (sum(ok) < 50) next
    fit <- lm(as.formula(paste("age ~", paste(markers, collapse=" + "))), data = sub[ok, ])
    out[idx[ok]] <- sub$age[ok] - predict(fit, newdata = sub[ok, markers])
    if (is.na(train_r2)) train_r2 <- summary(fit)$r.squared
  }
  CH$log("OK", sprintf("%s SPA done: %d valid, train-R2=%.3f", label, sum(!is.na(out)), train_r2), "4_BA")
  out
}

# ---------- 统一残差（避免行号错位）----------
add_resid <- function(dt, rawname, ver) {
  df0 <- data.frame(y = dt[[rawname]], age = dt$age)
  df0 <- df0[!is.na(df0$y) & !is.na(df0$age), ]
  if (nrow(df0) < 30) { dt[[paste0("BA_resid_", ver)]] <- NA_real_; return(dt) }
  r <- residuals(lm(y ~ age, data = df0))
  out <- rep(NA_real_, nrow(dt))
  out[!is.na(dt[[rawname]]) & !is.na(dt$age)] <- r
  dt[[paste0("BA_resid_", ver)]] <- out
  
  # 分性别
  rs <- rep(NA_real_, nrow(dt))
  for (g in sort(unique(dt$gender))) {
    i <- which(dt$gender == g & !is.na(dt[[rawname]]))
    if (length(i) < 30) next
    rf <- residuals(lm(dt[[rawname]][i] ~ dt$age[i]))
    rs[i] <- rf
  }
  dt[[paste0("BA_resid_sex_", ver)]] <- rs
  dt
}

# ---------- KDM_core6 固定系数投影法（跨队列统一）----------
calc_kdm_core6_fixed <- function(d, coef_list, label = "KDM_core6") {
  # coef_list: 从 NHANES 导出的 RDS，含 Male/Female 的 intercept + coefficients
  # 期望 d 含列: lncreat, bl_bun, bl_ua, bl_wbc, bl_mcv, systo, gender, age
  
  out <- rep(NA_real_, nrow(d))
  
  for (g in sort(unique(d$gender))) {
    idx <- which(d$gender == g)
    if (length(idx) == 0) next
    
    # 匹配性别编码
    sex_key <- ifelse(g == 1, "Male", "Female")  # 按你 CHARLS gender 编码调整
    if (!sex_key %in% names(coef_list)) {
      # 尝试反过来
      sex_key <- ifelse(g == 0, "Male", "Female")
    }
    if (!sex_key %in% names(coef_list)) {
      CH$log("WARN", sprintf("%s: sex=%d not found in coef_list", label, g), "4_BA")
      next
    }
    
    cf <- coef_list[[sex_key]]
    sub <- d[idx, , drop = FALSE]
    
    # 检查所有 marker 存在
    markers <- names(cf$coefficients)
    if (!all(markers %in% names(sub))) {
      missing_m <- setdiff(markers, names(sub))
      CH$log("WARN", sprintf("%s %s: missing markers %s", label, sex_key, 
                             paste(missing_m, collapse=",")), "4_BA")
      next
    }
    
    # 完整病例
    ok <- complete.cases(sub[, markers])
    if (sum(ok) < 10) next
    
    X <- as.matrix(sub[ok, markers, drop = FALSE])
    
    # 固定投影: BA = intercept + X %*% β
    ba_pred <- cf$intercept + as.vector(X %*% cf$coefficients)
    
    # 加速残差: age - predicted_age
    out[idx[ok]] <- sub$age[ok] - ba_pred
  }
  
  n_valid <- sum(!is.na(out))
  CH$log("OK", sprintf("%s fixed-projection: %d valid (NHANES coefficients)", label, n_valid), "4_BA")
  out
}

cat("✓ BA algorithms loaded (KDM + SPA, 4 versions)\n")