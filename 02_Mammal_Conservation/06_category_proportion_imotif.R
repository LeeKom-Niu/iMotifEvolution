
library(ggplot2)
library(dplyr)
theme_oup <- function(base_size = 8, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text = element_text(size = rel(1), color = "black"),
      axis.text.x = element_text(size = rel(1.2)),
      axis.text.y = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.2), face = "plain"),
      legend.position = "right",
      legend.title = element_text(size = rel(1), face = "plain"),
      legend.text = element_text(size = rel(0.9)),
      legend.key.size = unit(0.5, "cm"),
      panel.grid = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      plot.margin = margin(25, 30, 25, 30),
      text = element_text(family = base_family)
    )
}
imotif_group <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)
counts <- imotif_group %>%
  group_by(group) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percentage = n / sum(n) * 100)
counts$group <- factor(counts$group, levels = c("C1", "C2", "C3", "Other"))
counts$source <- "iMotif"
p <- ggplot(counts, aes(x = source, y = percentage, fill = group)) +
  geom_bar(stat = "identity", width = 0.4,
           color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)),
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 3,
            family = "Arial") +
  scale_fill_manual(values = c("C1" = "
                               "C2" = "
                               "C3" = "
                               "Other" = "
                    name = "Category") +
  labs(x = NULL, y = "Percentage (%)") +
  theme_oup() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    legend.position = "right"
  )
ggsave("imotif_category_proportion.pdf", p,
       width = 4, height = 4,
       device = cairo_pdf, dpi = 300)
ggsave("imotif_category_proportion.tiff", p,
       width = 4, height = 4,
       device = "tiff", dpi = 600, compression = "lzw")
print(counts)
write.table(counts, "imotif_category_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("完成！图片已保存为:\n")
cat("  PDF: imotif_category_proportion.pdf\n")
cat("  TIFF: imotif_category_proportion.tiff\n")
cat("统计表保存为 imotif_category_counts.tsv\n")
