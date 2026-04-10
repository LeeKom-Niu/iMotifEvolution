# 05_cal_coverage_boxplot_sig.R
# 功能：绘制覆盖度箱线图，并添加C1、C2、C3组间显著性标记（***）
# 修改：符合OUP插图指南 - 字体≥7pt，线条粗细0.25-1pt，色盲友好配色

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# ========== 参数设置 ==========
coverage_matrix_file <- "imotif_coverage_matrix.txt"
group_file <- "imotif_group.tsv"
species_list_file <- "sp240_nonhuman.txt"
order_file <- "Zoonomia_sp.txt"
output_plot_pdf <- "coverage_by_category_boxplot_sig.pdf"  # PDF用于AI编辑
output_plot_tiff <- "coverage_by_category_boxplot_sig.tiff" # TIFF用于印刷
output_stats <- "coverage_summary_stats.tsv"
output_tests <- "coverage_pairwise_tests.tsv"
# 色盲友好配色方案 (ColorBrewer Set2 调整版)
fill_colors <- c("C1" = "#D62728",  # 红色
                 "C2" = "#1F77B4",  # 蓝色
                 "C3" = "#2CA02C",  # 绿色
                 "Other" = "#B0BDB0") # 灰色
# =============================

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
      # 分面标签
      strip.background = element_rect(fill = "#F0F0F0", linewidth = 0.5),
      strip.text = element_text(size = rel(1.2), face = "plain"),
      # 网格线（使用细虚线）
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "#CCCCCC"),
      panel.grid.minor = element_blank(),
      # 面板边框
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # 图形边距
      plot.margin = margin(5, 5, 5, 5),
      # 确保所有文本元素使用相同字体
      text = element_text(family = base_family)
    )
}

# 1. 读取覆盖度矩阵
mat <- as.matrix(read.table(coverage_matrix_file, header = TRUE, row.names = 1))

# 2. 读取分组信息
group_df <- read.table(group_file, header = TRUE, stringsAsFactors = FALSE)
group <- group_df$group
names(group) <- group_df$iMotif

# 3. 读取物种列表
species_list <- readLines(species_list_file)

# 4. 读取物种分类信息
order_info <- read.table(order_file, header = TRUE, stringsAsFactors = FALSE)
colnames(order_info) <- tolower(colnames(order_info))
primate_sp <- order_info$species[order_info$order == "Primates"]
primate_sp <- intersect(primate_sp, species_list)
nonprimate_sp <- setdiff(species_list, primate_sp)

cat("灵长类物种数:", length(primate_sp), "\n")
cat("非灵长类物种数:", length(nonprimate_sp), "\n")

# 5. 计算平均覆盖度
avg_all <- rowMeans(mat, na.rm = TRUE)
avg_pri <- rowMeans(mat[, primate_sp, drop = FALSE], na.rm = TRUE)
avg_non <- rowMeans(mat[, nonprimate_sp, drop = FALSE], na.rm = TRUE)

avg_df <- data.frame(iMotif = rownames(mat), 
                     Mammalian = avg_all, 
                     Primate = avg_pri, 
                     Nonprimate = avg_non,
                     group = group)

# 6. 长格式
avg_long <- pivot_longer(avg_df, cols = c(Mammalian, Primate, Nonprimate), 
                         names_to = "Category", values_to = "Coverage")
avg_long$group <- factor(avg_long$group, levels = c("C1", "C2", "C3", "Other"))
avg_long$Category <- factor(avg_long$Category, levels = c("Mammalian", "Primate", "Nonprimate"))

# 7. 计算每个facet的最大y值，用于放置显著性标记
y_max <- avg_long %>%
  group_by(Category) %>%
  summarise(max_val = max(Coverage, na.rm = TRUE)) %>%
  mutate(y_pos = max_val * 1.1)

# 8. 定义比较组
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))

# 9. 绘制箱线图，应用OUP主题
p <- ggplot(avg_long, aes(x = group, y = Coverage, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +  # 箱线边框粗细0.4pt
  facet_wrap(~Category, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = fill_colors) +
  labs(x = NULL, y = "Average coverage") +
  theme_oup()

# 为每个facet手动添加显著性标记
for (cat in levels(avg_long$Category)) {
  y_pos_val <- y_max$y_pos[y_max$Category == cat]
  p <- p + geom_signif(
    data = filter(avg_long, Category == cat),
    aes(x = group, y = Coverage),
    comparisons = comparisons,
    map_signif_level = TRUE,
    y_position = rep(y_pos_val, length(comparisons)),
    tip_length = 0.01,
    textsize = 3.5,  # 约9.3pt，符合≥7pt要求
    linewidth = 0.3,  # 显著性标记线粗细
    inherit.aes = FALSE
  )
}

# 10. 保存为AI可编辑的PDF
ggsave(output_plot_pdf, p, width = 10, height = 5, 
       device = cairo_pdf, dpi = 300)

# 保存为TIFF用于印刷
ggsave(output_plot_tiff, p, width = 10, height = 5, 
       device = "tiff", dpi = 600, compression = "lzw")

cat("绘图完成！输出文件:", output_plot_pdf, "和", output_plot_tiff, "\n\n")

# ========== 后续统计摘要和检验保持不变 ==========
cat("========== 统计摘要 ==========\n")
stats_summary <- avg_long %>%
  group_by(Category, group) %>%
  summarise(
    n = n(),
    mean = mean(Coverage, na.rm = TRUE),
    sd = sd(Coverage, na.rm = TRUE),
    median = median(Coverage, na.rm = TRUE),
    Q1 = quantile(Coverage, 0.25, na.rm = TRUE),
    Q3 = quantile(Coverage, 0.75, na.rm = TRUE),
    .groups = "drop"
  )
print(stats_summary)
write.table(stats_summary, output_stats, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n统计摘要已保存至:", output_stats, "\n")

cat("\n========== 组间 Wilcoxon 检验 p 值 ==========\n")
categories <- levels(avg_long$Category)
group_pairs <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))
test_results <- data.frame()

for (cat in categories) {
  cat("\nCategory:", cat, "\n")
  subdata <- filter(avg_long, Category == cat)
  for (pair in group_pairs) {
    g1 <- pair[1]; g2 <- pair[2]
    x <- subdata$Coverage[subdata$group == g1]
    y <- subdata$Coverage[subdata$group == g2]
    wt <- wilcox.test(x, y, na.rm = TRUE)
    p_val <- wt$p.value
    cat(sprintf("%s vs %s: p = %.3e\n", g1, g2, p_val))
    test_results <- rbind(test_results, data.frame(
      Category = cat, group1 = g1, group2 = g2, p_value = p_val
    ))
  }
}
write.table(test_results, output_tests, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n检验结果已保存至:", output_tests, "\n")
