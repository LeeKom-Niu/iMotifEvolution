# 07_score_all_versions.R
# 目的：将 iMotif 的稳定性分数与分组信息关联，比较各组（C1、C2、C3、Other）的分数分布。
# 图形类型：箱线图、小提琴图+箱线图、箱线图+抖动散点图（三种全部生成）
# 显著性标记：采用原始脚本样式，三个比较组共用同一 y 轴高度（可能会重叠）
# 符合OUP插图指南 - 字体≥7pt，线条粗细0.25-1pt，色盲友好配色

library(dplyr)
library(ggplot2)
library(ggsignif)

# ========== 参数设置 ==========
input_bed <- "imotif_clean.bed"
input_group <- "imotif_group.tsv"
output_boxplot_pdf <- "stability_boxplot_sig.pdf"
output_violin_pdf <- "stability_violin_sig.pdf"
output_boxplot_jitter_pdf <- "stability_boxplot_jitter.pdf"
output_boxplot_tiff <- "stability_boxplot_sig.tiff"
output_violin_tiff <- "stability_violin_sig.tiff"
output_boxplot_jitter_tiff <- "stability_boxplot_jitter.tiff"
output_stats <- "stability_group_summary.tsv"
# 色盲友好配色方案
colors <- c("C1" = "#D62728", 
            "C2" = "#1F77B4", 
            "C3" = "#2CA02C", 
            "Other" = "#B0BDB0")
# ==============================

# 定义符合OUP指南的主题函数
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 轴线粗细
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # 刻度标签大小
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      # 图例
      legend.position = "none",
      # 网格线（使用细虚线）
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "#CCCCCC"),
      panel.grid.minor = element_blank(),
      # 面板边框
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # 图形边距
      plot.margin = margin(25, 25, 25, 25),
      # 确保所有文本元素使用相同字体
      text = element_text(family = base_family)
    )
}

# 1. 读取稳定性分数文件
stability <- read.table(input_bed, header = FALSE,
                        col.names = c("chr", "start", "end", "name", "score", "strand"),
                        stringsAsFactors = FALSE)
stability <- stability[, c("name", "score")]

# 2. 读取分组信息
group_df <- read.table(input_group, header = TRUE, stringsAsFactors = FALSE)

# 3. 合并数据
combined <- left_join(group_df, stability, by = c("iMotif" = "name"))
if (any(is.na(combined$score))) {
  warning("部分 iMotif 在稳定性文件中未找到，这些行将被过滤。")
  combined <- combined %>% filter(!is.na(score))
}

# 4. 设置分组顺序
combined$group <- factor(combined$group, levels = c("C1", "C2", "C3", "Other"))

# 5. 计算y轴最大值用于显著性标记位置（原始脚本样式：所有比较组共用同一高度）
y_max <- max(combined$score, na.rm = TRUE)
y_pos <- y_max * 1.1

# 6. 定义比较组
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))

# 7. 绘制箱线图 + 显著性标记（原始样式）
p_box <- ggplot(combined, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Stability score") +
  theme_oup() +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = rep(y_pos, length(comparisons)),
              tip_length = 0.01,
              textsize = 3.5,
              linewidth = 0.3)

# 8. 绘制小提琴图 + 箱线图 + 显著性标记（原始样式）
p_violin <- ggplot(combined, aes(x = group, y = score, fill = group)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA, 
               fill = "white", color = "black", alpha = 0.7,
               linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Stability score") +
  theme_oup() +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = rep(y_pos, length(comparisons)),
              tip_length = 0.01,
              textsize = 3.5,
              linewidth = 0.3)

# 9. 绘制箱线图 + 抖动散点图 + 显著性标记（原始样式，散点透明度提高避免遮盖）
p_box_jitter <- ggplot(combined, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4, alpha = 0.8) +
  geom_jitter(aes(color = group), width = 0.2, size = 0.6, alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  labs(x = NULL, y = "Stability score") +
  theme_oup() +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = rep(y_pos, length(comparisons)),
              tip_length = 0.01,
              textsize = 3.5,
              linewidth = 0.3)

# 10. 保存为AI可编辑的PDF
ggsave(output_boxplot_pdf, p_box, width = 5, height = 4, 
       device = cairo_pdf, dpi = 300)
ggsave(output_violin_pdf, p_violin, width = 5, height = 4, 
       device = cairo_pdf, dpi = 300)
ggsave(output_boxplot_jitter_pdf, p_box_jitter, width = 5, height = 4, 
       device = cairo_pdf, dpi = 300)

# 保存为TIFF用于印刷
ggsave(output_boxplot_tiff, p_box, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")
ggsave(output_violin_tiff, p_violin, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")
ggsave(output_boxplot_jitter_tiff, p_box_jitter, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")

# 11. 输出统计摘要
summary_stats <- combined %>%
  group_by(group) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    sd_score = sd(score, na.rm = TRUE),
    n = n()
  )
write.table(summary_stats, output_stats, sep = "\t", quote = FALSE, row.names = FALSE)

cat("完成！图片已保存为:\n")
cat("  PDF: ", output_boxplot_pdf, ", ", output_violin_pdf, ", ", output_boxplot_jitter_pdf, "\n")
cat("  TIFF: ", output_boxplot_tiff, ", ", output_violin_tiff, ", ", output_boxplot_jitter_tiff, "\n")
cat("统计摘要保存为", output_stats, "\n")
