# =====================================================
# imotif_category_proportion.R
# Function: Plot stacked proportion of i-motif region evolutionary categories
# Modification: Remove shuffled regions, keep only i-motif data
# Output: PDF and TIFF format, OUP compliant
# =====================================================

library(ggplot2)
library(dplyr)

# Define OUP-compliant theme function
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis line widths
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # Tick label size
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(size = rel(1.2)),
      axis.text.y = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      # Legend
      legend.position = "right",
      legend.title = element_text(size = rel(1), face = "plain"),
      legend.text = element_text(size = rel(0.9)),
      legend.key.size = unit(0.5, "cm"),
      # Grid lines
      panel.grid = element_blank(),
      # Panel border
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # Plot margin
      plot.margin = margin(25, 30, 25, 30),
      # Ensure all text elements use the same font
      text = element_text(family = base_family)
    )
}

# 1. Read i-motif group data only
imotif_group <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)

# 2. Count per group, calculate percentages
counts <- imotif_group %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = n / sum(n) * 100)

# 3. Set factor order
counts$group <- factor(counts$group, levels = c("C1", "C2", "C3", "Other"))

# 4. Add fixed x-axis label (single group)
counts$source <- "iMotif"

# 5. Plot stacked bar chart (single bar)
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

# 6. Save as AI-editable PDF
ggsave("imotif_category_proportion.pdf", p,
       width = 4, height = 4,   # 调整为正方形，更适合单组
       device = cairo_pdf, dpi = 300)

# 7. Save as TIFF for printing
ggsave("imotif_category_proportion.tiff", p,
       width = 4, height = 4,
       device = "tiff", dpi = 600, compression = "lzw")

# 8. Print and save statistics table
print(counts)
write.table(counts, "imotif_category_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("Complete! Figures saved as:\n")
cat("  PDF: imotif_category_proportion.pdf\n")
cat("  TIFF: imotif_category_proportion.tiff\n")
cat("Statistics table saved as imotif_category_counts.tsv\n")
