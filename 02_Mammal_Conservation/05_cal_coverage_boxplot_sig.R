# 05_cal_coverage_boxplot_sig.R
# Function: Plot coverage boxplot with significance markers between C1, C2, C3 groups (***)
# Modification: OUP compliant - font >=7pt, line width 0.25-1pt, colorblind-friendly palette

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# ========== Parameter Settings ==========
coverage_matrix_file <- "imotif_coverage_matrix.txt"
group_file <- "imotif_group.tsv"
species_list_file <- "sp240_nonhuman.txt"
order_file <- "Zoonomia_sp.txt"
output_plot_pdf <- "coverage_by_category_boxplot_sig.pdf"  # PDF for AI editing
output_plot_tiff <- "coverage_by_category_boxplot_sig.tiff" # TIFF for printing
output_stats <- "coverage_summary_stats.tsv"
output_tests <- "coverage_pairwise_tests.tsv"
# Colorblind-friendly color scheme (ColorBrewer Set2 adjusted)
fill_colors <- c("C1" = "#D62728",  # Red
                 "C2" = "#1F77B4",  # Blue
                 "C3" = "#2CA02C",  # Green
                 "Other" = "#B0BDB0") # Gray
# =============================

# Define OUP-compliant theme function
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis line widths
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # Tick label size
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      # Legend
      legend.position = "none",
      # Facet labels
      strip.background = element_rect(fill = "#F0F0F0", linewidth = 0.5),
      strip.text = element_text(size = rel(1.2), face = "plain"),
      # Grid lines (use thin dotted lines)
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.3, linetype = "dotted", color = "#CCCCCC"),
      panel.grid.minor = element_blank(),
      # Panel border
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # Plot margin
      plot.margin = margin(5, 5, 5, 5),
      # Ensure all text elements use the same font
      text = element_text(family = base_family)
    )
}

# 1. Read coverage matrix
mat <- as.matrix(read.table(coverage_matrix_file, header = TRUE, row.names = 1))

# 2. Read group information
group_df <- read.table(group_file, header = TRUE, stringsAsFactors = FALSE)
group <- group_df$group
names(group) <- group_df$iMotif

# 3. Read species list
species_list <- readLines(species_list_file)

# 4. Read species classification info
order_info <- read.table(order_file, header = TRUE, stringsAsFactors = FALSE)
colnames(order_info) <- tolower(colnames(order_info))
primate_sp <- order_info$species[order_info$order == "Primates"]
primate_sp <- intersect(primate_sp, species_list)
nonprimate_sp <- setdiff(species_list, primate_sp)

cat("Primate species:", length(primate_sp), "\n")
cat("Non-primate species:", length(nonprimate_sp), "\n")

# 5. Calculate average coverage
avg_all <- rowMeans(mat, na.rm = TRUE)
avg_pri <- rowMeans(mat[, primate_sp, drop = FALSE], na.rm = TRUE)
avg_non <- rowMeans(mat[, nonprimate_sp, drop = FALSE], na.rm = TRUE)

avg_df <- data.frame(iMotif = rownames(mat), 
                     Mammalian = avg_all, 
                     Primate = avg_pri, 
                     Nonprimate = avg_non,
                     group = group)

# 6. Long format
avg_long <- pivot_longer(avg_df, cols = c(Mammalian, Primate, Nonprimate), 
                         names_to = "Category", values_to = "Coverage")
avg_long$group <- factor(avg_long$group, levels = c("C1", "C2", "C3", "Other"))
avg_long$Category <- factor(avg_long$Category, levels = c("Mammalian", "Primate", "Nonprimate"))

# 7. Calculate max y per facet for significance marker placement
y_max <- avg_long %>%
  group_by(Category) %>%
  summarise(max_val = max(Coverage, na.rm = TRUE)) %>%
  mutate(y_pos = max_val * 1.1)

# 8. Define comparison groups
comparisons <- list(c("C1", "C2"), c("C1", "C3"), c("C2", "C3"))

# 9. Plot boxplot with OUP theme
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

# 10. Save as AI-editable PDF
ggsave(output_plot_pdf, p, width = 10, height = 5, 
       device = cairo_pdf, dpi = 300)

# Save as TIFF for printing
ggsave(output_plot_tiff, p, width = 10, height = 5, 
       device = "tiff", dpi = 600, compression = "lzw")

cat("Plotting complete! Output files:", output_plot_pdf, "and", output_plot_tiff, "\n\n")

# ========== Subsequent statistical summary and tests kept unchanged ==========
cat("========== Statistical Summary ==========\n")
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
cat("\nStatistical summary saved to:", output_stats, "\n")

cat("\n========== Between-group Wilcoxon test p-values ==========\n")
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
cat("\nTest results saved to:", output_tests, "\n")
