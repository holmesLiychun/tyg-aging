# ============================================================================
# 模块 65c：分层结果森林图（修正路径版）
# ============================================================================

# ---- 定义路径（如果还没定义）----
if (!exists("PROJECT_ROOT")) {
  stop("请先定义 PROJECT_ROOT，例如: PROJECT_ROOT <- '你的项目路径'")
}

# 用你的实际路径替换下面这行，或者确认 10_setup.R 里叫什么
fig_dir <- ifelse(exists("FIG_DIR"), FIG_DIR, 
                  file.path(OUT_DIR, "figures"))
out_dir <- ifelse(exists("OUT_DIR"), OUT_DIR, 
                  file.path(OUT_DIR, "output"))

# 确保目录存在
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 读取分层结果 ----
csv_path <- file.path(out_dir, "table_S4_stratified_analysis.csv")
if (!file.exists(csv_path)) {
  # 如果65b保存到了别的位置，手动指定
  csv_path <- file.path(PROJECT_ROOT, "table_S4_stratified_analysis.csv")
}

strat_df <- read.csv(csv_path)

# 只保留 PhenoAge-9
strat_p9 <- strat_df[strat_df$Algorithm == "PhenoAge-9", ]



# 计算 95% CI
strat_p9$CI_lower <- strat_p9$Beta - 1.96 * strat_p9$SE
strat_p9$CI_upper <- strat_p9$Beta + 1.96 * strat_p9$SE

# 设置绘图顺序
strat_p9$Label <- paste(strat_p9$Stratification, strat_p9$Level, sep = " | ")
strat_p9$Label <- factor(strat_p9$Label, 
                         levels = rev(c("Diabetes status | Diabetes",
                                        "Diabetes status | Non-diabetes", 
                                        "BMI category | ≥30",
                                        "BMI category | 25-29.9",
                                        "BMI category | <25 (ref)")))

# 2. 极致拉宽画布，左侧留足标签空间，右侧留足P值空间
# 1. 极限拉宽画布（物理尺寸），给左右文本留足绝对空间
# 极限拉宽画布
png(file.path(OUT_DIR, "figures/forest_stratified_PhenoAge9_FINAL_FIX.png"), 
    width = 2400, height = 700, res = 150, pointsize = 12)

# 边距：右侧给足空间
par(mar = c(4, 24, 4, 20), mgp = c(3, 0.5, 0))
n <- nrow(strat_p9)

plot(0, 0, type = "n", 
     xlim = c(-10, 10), ylim = c(0.5, n + 1.5), # Y轴上限稍微调大，给顶部留空间
     xlab = "Beta (years per 1-SD TyG increase)", 
     ylab = "", yaxt = "n", 
     main = "Stratified Analysis: TyG and PhenoAge-9 Acceleration",
     cex.main = 1.5, cex.lab = 1.2)

abline(v = 0, lty = 2, col = "gray60", lwd = 1.5)
cols <- ifelse(strat_p9$Beta > 0, "#E74C3C", "#3498DB")

for (i in 1:n) {
  y_pos <- n - i + 1
  
  # 画CI线段和点
  segments(strat_p9$CI_lower[i], y_pos, 
           strat_p9$CI_upper[i], y_pos, 
           col = cols[i], lwd = 2.5)
  points(strat_p9$Beta[i], y_pos, pch = 18, cex = 2.0, col = cols[i])
  
  # 左侧标签
  text(-9.5, y_pos, strat_p9$Label[i], adj = c(0, 0.5), cex = 0.85)
  
  # Beta值
  beta_text <- sprintf("%.2f", strat_p9$Beta[i])
  text(strat_p9$CI_upper[i] + 0.5, y_pos, beta_text, 
       adj = c(0, 0.5), cex = 0.75, col = cols[i])
  
  # 生成P值文本
  p_text <- ifelse(strat_p9$P_value[i] < 0.001, "P<0.001", 
                   sprintf("P=%.4f", strat_p9$P_value[i]))
  p_col <- ifelse(strat_p9$P_value[i] < 0.05, "#C0392B", "gray40")
  
  # 【核心修复】：使用固定行距偏移，每一行P值都在其对应的y_pos上，彻底告别重叠
  text(9.5, y_pos + 0.15, p_text, adj = c(1, 0.5), cex = 0.7, col = p_col)
}

dev.off()
log_write("✓ P值重叠彻底修复版已保存")