# =====================================================
# imotif_category_proportion.R
# 功能：仅绘制 i‑motif 区域的进化类别堆叠比例图
# 修改：移除随机区域（Shuffled regions），保留 i‑motif 数据
# 输出：PDF 和 TIFF 格式，符合 OUP 插图指南
# =====================================================

library(ggplot2)
library(dplyr)

# 定义符合 OUP 指南的主题函数
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 轴线粗细
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # 刻度标签大小
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(size = rel(1.2)),
      axis.text.y = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      # 图例
      legend.position = "right",
      legend.title = element_text(size = rel(1), face = "plain"),
      legend.text = element_text(size = rel(0.9)),
      legend.key.size = unit(0.5, "cm"),
      # 网格线
      panel.grid = element_blank(),
      # 面板边框
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # 图形边距
      plot.margin = margin(25, 30, 25, 30),
      # 确保所有文本元素使用相同字体
      text = element_text(family = base_family)
    )
}

# 1. 仅读取 i‑motif 分组数据
imotif_group <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)

# 2. 统计各组数量，计算百分比
counts <- imotif_group %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = n / sum(n) * 100)

# 3. 设置因子顺序
counts$group <- factor(counts$group, levels = c("C1", "C2", "C3", "Other"))

# 4. 添加固定 x 轴标签（单组）
counts$source <- "iMotif"

# 5. 绘制堆叠条形图（单柱）
p <- ggplot(counts, aes(x = source, y = percentage, fill = group)) +
  geom_bar(stat = "identity", width = 0.4,
           color = "black", linewidth = 0.3) +
  # 标签居中显示
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 3,  # 约 8pt
            family = "Arial") +
  scale_fill_manual(values = c("C1" = "#D62728",
                               "C2" = "#1F77B4",
                               "C3" = "#2CA02C",
                               "Other" = "#B0BDB0"),
                    name = "Category") +
  labs(x = NULL, y = "Percentage (%)") +
  theme_oup() +
  theme(
    axis.text.x = element_text(size = 10),  # 约 10pt
    axis.text.y = element_text(size = 8),   # 约 8pt
    legend.position = "right"
  )

# 6. 保存为 AI 可编辑的 PDF
ggsave("imotif_category_proportion.pdf", p,
       width = 4, height = 4,   # 调整为正方形，更适合单组
       device = cairo_pdf, dpi = 300)

# 7. 保存为 TIFF 用于印刷
ggsave("imotif_category_proportion.tiff", p,
       width = 4, height = 4,
       device = "tiff", dpi = 600, compression = "lzw")

# 8. 打印并保存统计表
print(counts)
write.table(counts, "imotif_category_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("完成！图片已保存为:\n")
cat("  PDF: imotif_category_proportion.pdf\n")
cat("  TIFF: imotif_category_proportion.tiff\n")
cat("统计表保存为 imotif_category_counts.tsv\n")
