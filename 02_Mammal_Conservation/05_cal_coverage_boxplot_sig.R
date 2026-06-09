
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
coverage_matrix_file <- "imotif_coverage_matrix.txt"
group_file <- "imotif_group.tsv"
species_list_file <- "sp240_nonhuman.txt"
order_file <- "Zoonomia_sp.txt"
output_plot_pdf <- "coverage_by_category_boxplot_sig.pdf"
output_plot_tiff <- "coverage_by_category_boxplot_sig.tiff"
output_stats <- "coverage_summary_stats.tsv"
output_tests <- "coverage_pairwise_tests.tsv"
fill_colors <- c("C1" = "
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
      strip.background = element_rect(fill = "
      strip.text = element_text(size = rel(1.2), face = "plain"),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      plot.margin = margin(5, 5, 5, 5),
      text = element_text(family = base_family)
    )
}
mat <- as.matrix(read.table(coverage_matrix_file, header = TRUE, row.names = 1))
group_df <- read.table(group_file, header = TRUE, stringsAsFactors = FALSE)
group <- group_df$group
names(group) <- group_df$iMotif
species_list <- readLines(species_list_file)
order_info <- read.table(order_file, header = TRUE, stringsAsFactors = FALSE)
colnames(order_info) <- tolower(colnames(order_info))
primate_sp <- order_info$species[order_info$order == "Primates"]
primate_sp <- intersect(primate_sp, species_list)
nonprimate_sp <- setdiff(species_list, primate_sp)
cat("灵长类物种数:", length(primate_sp), "\n")
cat("非灵长类物种数:", length(nonprimate_sp), "\n")
avg_all <- rowMeans(mat, na.rm = TRUE)
avg_pri <- rowMeans(mat[, primate_sp, drop = FALSE], na.rm = TRUE)
avg_non <- rowMeans(mat[, nonprimate_sp, drop = FALSE], na.rm = TRUE)
avg_df <- data.frame(iMotif = rownames(mat), 
                     Mammalian = avg_all, 
                     Primate = avg_pri, 
                     Nonprimate = avg_non,
                     group = group)
avg_long <- pivot_longer(avg_df, cols = c(Mammalian, Primate, Nonprimate), 
                         names_to = "Category", values_to = "Coverage")
avg_long$group <- factor(avg_long$group, levels = c("C1", "C2", "C3", "Other"))
avg_long$Category <- factor(avg_long$Category, levels = c("Mammalian", "Primate", "Nonprimate"))
y_max <- avg_long %>%
  group_by(Category) %>%
  summarise(max_val = max(Coverage, na.rm = TRUE)) %>%
  mutate(y_pos = max_val * 1.1)
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))
p <- ggplot(avg_long, aes(x = group, y = Coverage, fill = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.4) +
  facet_wrap(~Category, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = fill_colors) +
  labs(x = NULL, y = "Average coverage") +
  theme_oup()
for (cat in levels(avg_long$Category)) {
  y_pos_val <- y_max$y_pos[y_max$Category == cat]
  p <- p + geom_signif(
    data = filter(avg_long, Category == cat),
    aes(x = group, y = Coverage),
    comparisons = comparisons,
    map_signif_level = TRUE,
    y_position = rep(y_pos_val, length(comparisons)),
    tip_length = 0.01,
    textsize = 3.5,
    linewidth = 0.3,
    inherit.aes = FALSE
  )
}
ggsave(output_plot_pdf, p, width = 10, height = 5, 
       device = cairo_pdf, dpi = 300)
ggsave(output_plot_tiff, p, width = 10, height = 5, 
       device = "tiff", dpi = 600, compression = "lzw")
cat("绘图完成！输出文件:", output_plot_pdf, "和", output_plot_tiff, "\n\n")
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
