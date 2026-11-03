# ============================================================================
# 70.4b ★ DMRC 解耦前后对比森林图（投稿级优化·最终版 v3）
# ============================================================================
# 🆕 修改日志：
#   - Outcome 标签重映射（case_when 写法，修复 recode .default 长度报错）
#   - 标题改为 CKM/DMRC 框架友好表述
#   - caption 补充 DMRC 全称
#   - y 轴使用 Outcome_label 因子
# ============================================================================

timer_start("70.4b DMRC森林图")
log_write(">>> 70.4b DMRC 解耦前后对比森林图 v3")

if (exists("dmrc_combined") && nrow(dmrc_combined) > 0) {
  
  # ---- 0a. 统一 FDR 列名 ----
  dmrc_plot <- dmrc_combined %>%
    mutate(
      P_fdr_unified = case_when(
        Exposure == "Before decoupling (Full TyG)" ~ P_before_fdr.x,
        Exposure == "After decoupling (TG+glucose removed)" ~ P_after_fdr.x,
        TRUE ~ P_value
      ),
      sig_star = case_when(
        P_fdr_unified < 0.001 ~ "***",
        P_fdr_unified < 0.01  ~ "**",
        P_fdr_unified < 0.05  ~ "*",
        TRUE                 ~ ""
      )
    )
  
  # ---- 0a+. Outcome 标签重映射（投稿版 · case_when 写法）----
  # ❓ 作者请确认：是否还有其他 Outcome 需要映射？
  dmrc_plot <- dmrc_plot %>%
    mutate(
      Outcome_label = case_when(
        Outcome == "KDM-9"       ~ "KDM+CRP",
        Outcome == "KDM-11"      ~ "KDM+TG",
        Outcome == "PhenoAge-9"  ~ "PhenoAge-original",
        # 🆕 以下为推测的其他常见 Outcome，请作者按需取消注释或补充
        # Outcome == "Pheno-noG"   ~ "PhenoAge-noGlucose",
        # Outcome == "Pheno-noCRP" ~ "PhenoAge-noCRP",
        # Outcome == "KDM-BA"      ~ "KDM-BiologicalAge",
        TRUE                      ~ as.character(Outcome)
      )
    )
  
  # ---- 0b. 排序因子（按 Before β 绝对值降序，使用 Outcome_label）----
  outcome_order <- dmrc_plot %>%
    filter(Exposure == "Before decoupling (Full TyG)") %>%
    arrange(desc(abs(Beta))) %>%
    pull(Outcome_label)
  
  dmrc_plot$Outcome_label <- factor(dmrc_plot$Outcome_label, levels = rev(outcome_order))
  
  # ---- 1. 主图 ----
  p_dmrc <- ggplot(dmrc_plot, 
                   aes(x = Beta, y = Outcome_label, color = Exposure)) +
    
    geom_segment(aes(x = CI_low, xend = CI_high, yend = Outcome_label),
                 linewidth = 1.1, 
                 position = position_dodge(width = 0.35),
                 alpha = 0.85) +
    
    geom_point(size = 3.2, 
               position = position_dodge(width = 0.35)) +
    
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "grey50", linewidth = 0.5) +
    
    # 右侧文本列（仅 After 侧，x=1.25 防拥挤）
    geom_text(
      data = dmrc_plot %>% filter(Exposure == "After decoupling (TG+glucose removed)"),
      aes(x = 1.25, 
          label = sprintf("%s %.2f [%.2f, %.2f]", sig_star, Beta, CI_low, CI_high),
          color = NULL),
      hjust = 0, size = 3.1, fontface = "bold", color = "black"
    ) +
    
    scale_color_manual(
      values = c("Before decoupling (Full TyG)" = "#D55E00",
                 "After decoupling (TG+glucose removed)" = "#0072B2"),
      name = NULL,
      labels = c("Before decoupling\n(Full TyG)", 
                 "After decoupling\n(TG+glucose removed)")
    ) +
    
    scale_x_continuous(limits = c(-2.6, 2.2), 
                       breaks = seq(-2, 1.5, by = 0.5)) +
    
    labs(
      title = "Shared Metabolic Components Inflate TyG–Biological Age Associations: A DMRC Sensitivity Analysis",
      subtitle = "Multivariable-adjusted (Model 5); NHANES discovery (n=15,433) and CHARLS replication; FDR-corrected",
      x = expression(beta ~ "(years per 1-SD increase in TyG)"),
      y = NULL,
      caption = paste0(
        "DMRC, Decoupled Metabolic Residualization and Covariate adjustment; ",
        "Residualization: TyG ~ triglycerides + glucose + age, sex-stratified.\n",
        "*FDR < 0.05, **FDR < 0.01, ***FDR < 0.001. ",
        "Pheno-noG: PhenoAge excluding glucose; Pheno-noCRP: PhenoAge excluding CRP."
      )
    ) +
    
    theme_classic(base_size = 13) +
    theme(
      legend.position = c(0.12, 0.12),
      legend.background = element_rect(fill = "white", 
                                       color = "grey70", linewidth = 0.3),
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = 9),
      
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, 
                                   color = "grey30"),
      plot.caption = element_text(hjust = 0, size = 8, 
                                  color = "grey50", face = "italic"),
      
      axis.title.y = element_blank(),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(face = "bold", size = 11, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      
      axis.line = element_line(color = "grey40", linewidth = 0.3)
    )
  
  save_fig(p_dmrc, "Forest_DMRC_Before_After_Final", width = 11, height = 6.5)
  log_write("✓ DMRC 森林图最终版 已保存")
  
} else {
  log_write("  ⚠ dmrc_combined 不存在，跳过")
}

timer_end("70.4b DMRC森林图")
log_write("===== 模块 70.4b 结束 =====\n")
