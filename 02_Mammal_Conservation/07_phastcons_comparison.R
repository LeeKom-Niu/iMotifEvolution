# =====================================================
# 09_phastcons_comparison.R
# 功能：比较C1、C2、C3、Other组的PhastCons保守性分数
# 修改：符合OUP插图指南 - 字体≥7pt，线条粗细0.25-1pt，色盲友好配色
# =====================================================

library(ggplot2)
library(dplyr)
library(ggsignif)

# ========== 参数设置 ==========
imotif_bed_file <- "imotif_clean.bed4"
group_file <- "imotif_group.tsv"
bigwig_file <- "hg38.phastCons470way.bw"
output_scores <- "phastcons_scores.txt"
output_boxplot_pdf <- "phastcons_boxplot_sig.pdf"
output_violin_pdf <- "phastcons_violin_sig.pdf"
output_boxplot_tiff <- "phastcons_boxplot_sig.tiff"
output_violin_tiff <- "phastcons_violin_sig.tiff"
output_stats <- "phastcons_summary.tsv"
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

# 1. 检查bigWig文件是否存在
if (!file.exists(bigwig_file)) {
  stop("错误：找不到PhastCons bigWig文件，请先下载：",
       "wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons470way/hg38.phastCons470way.bw")
}

# 2. 使用bigWigAverageOverBed提取分数
cat("正在使用bigWigAverageOverBed提取PhastCons分数...\n")
cmd <- paste("bigWigAverageOverBed", bigwig_file, imotif_bed_file, output_scores)
system(cmd)

# 3. 读取分数文件
scores <- read.table(output_scores, header = FALSE,
                     col.names = c("name", "size", "covered", "sum", "mean0", "mean"),
                     stringsAsFactors = FALSE)

# 4. 读取分组信息并合并
group_df <- read.table(group_file, header = TRUE, stringsAsFactors = FALSE)
combined <- left_join(group_df, scores, by = c("iMotif" = "name"))
if (any(is.na(combined$mean))) {
  warning("部分 iMotif 在分数文件中未找到，它们将被过滤。")
  combined <- combined %>% filter(!is.na(mean))
}

# 5. 设置分组顺序
combined$group <- factor(combined$group, levels = c("C1", "C2", "C3", "Other"))

# 6. 计算y轴最大值用于显著性标记位置
y_max <- max(combined$mean, na.rm = TRUE)
y_pos <- y_max * 1.1

# 7. 定义比较组
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))

# 8. 绘制箱线图 + 显著性标记
p_box <- ggplot(combined, aes(x = group, y = mean, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +   # 箱线边框粗细0.4pt
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "PhastCons score") +
  theme_oup() +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = rep(y_pos, length(comparisons)),
              tip_length = 0.01,
              textsize = 3.5,  # 约9.3pt
              linewidth = 0.3)

# 9. 绘制小提琴图 + 箱线图 + 显著性标记
p_violin <- ggplot(combined, aes(x = group, y = mean, fill = group)) +
  geom_violin(alpha = 0.5, scale = "width", linewidth = 0.4) +
  geom_boxplot(width = 0.2, outlier.shape = NA, 
               fill = "white", color = "black", alpha = 0.7,
               linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "PhastCons score") +
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

# 保存为TIFF用于印刷
ggsave(output_boxplot_tiff, p_box, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")
ggsave(output_violin_tiff, p_violin, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")

# 11. 输出统计摘要
stats <- combined %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean = mean(mean, na.rm = TRUE),
    sd = sd(mean, na.rm = TRUE),
    median = median(mean, na.rm = TRUE),
    Q1 = quantile(mean, 0.25, na.rm = TRUE),
    Q3 = quantile(mean, 0.75, na.rm = TRUE)
  )
write.table(stats, output_stats, sep = "\t", quote = FALSE, row.names = FALSE)

cat("PhastCons 分析完成！\n")
cat("图片已保存为:\n")
cat("  PDF: ", output_boxplot_pdf, ", ", output_violin_pdf, "\n")
cat("  TIFF: ", output_boxplot_tiff, ", ", output_violin_tiff, "\n")
cat("统计摘要已保存为", output_stats, "\n")
