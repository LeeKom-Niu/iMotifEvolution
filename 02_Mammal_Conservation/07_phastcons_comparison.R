
library(ggplot2)
library(dplyr)
library(ggsignif)
imotif_bed_file <- "imotif_clean.bed4"
group_file <- "imotif_group.tsv"
bigwig_file <- "hg38.phastCons470way.bw"
output_scores <- "phastcons_scores.txt"
output_boxplot_pdf <- "phastcons_boxplot_sig.pdf"
output_violin_pdf <- "phastcons_violin_sig.pdf"
output_boxplot_tiff <- "phastcons_boxplot_sig.tiff"
output_violin_tiff <- "phastcons_violin_sig.tiff"
output_stats <- "phastcons_summary.tsv"
colors <- c("C1" = "
            "C2" = "
            "C3" = "
            "Other" = "
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      plot.margin = margin(25, 25, 25, 25),
      text = element_text(family = base_family)
    )
}
if (!file.exists(bigwig_file)) {
  stop("错误：找不到PhastCons bigWig文件，请先下载：",
       "wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons470way/hg38.phastCons470way.bw")
}
cat("正在使用bigWigAverageOverBed提取PhastCons分数...\n")
cmd <- paste("bigWigAverageOverBed", bigwig_file, imotif_bed_file, output_scores)
system(cmd)
scores <- read.table(output_scores, header = FALSE,
                     col.names = c("name", "size", "covered", "sum", "mean0", "mean"),
                     stringsAsFactors = FALSE)
group_df <- read.table(group_file, header = TRUE, stringsAsFactors = FALSE)
combined <- left_join(group_df, scores, by = c("iMotif" = "name"))
if (any(is.na(combined$mean))) {
  warning("部分 iMotif 在分数文件中未找到，它们将被过滤。")
  combined <- combined %>% filter(!is.na(mean))
}
combined$group <- factor(combined$group, levels = c("C1", "C2", "C3", "Other"))
y_max <- max(combined$mean, na.rm = TRUE)
y_pos <- y_max * 1.1
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))
p_box <- ggplot(combined, aes(x = group, y = mean, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "PhastCons score") +
  theme_oup() +
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = rep(y_pos, length(comparisons)),
              tip_length = 0.01,
              textsize = 3.5,
              linewidth = 0.3)
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
ggsave(output_boxplot_pdf, p_box, width = 5, height = 4, 
       device = cairo_pdf, dpi = 300)
ggsave(output_violin_pdf, p_violin, width = 5, height = 4, 
       device = cairo_pdf, dpi = 300)
ggsave(output_boxplot_tiff, p_box, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")
ggsave(output_violin_tiff, p_violin, width = 5, height = 4, 
       device = "tiff", dpi = 600, compression = "lzw")
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
