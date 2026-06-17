library(ggplot2)
library(patchwork)   # If not installed, run install.packages("patchwork")

# Read data
imotif <- read.table("imotif_af_info.tsv", header = FALSE,
                     col.names = c("chr", "pos", "AF"))
bg <- read.table("background_af_info.tsv", header = FALSE,
                 col.names = c("chr", "pos", "AF"))

# Filter
imotif <- imotif[!is.na(imotif$AF) & imotif$AF >= 0 & imotif$AF <= 1, ]
bg <- bg[!is.na(bg$AF) & bg$AF >= 0 & bg$AF <= 1, ]

# Calculate MAF
imotif$MAF <- pmin(imotif$AF, 1 - imotif$AF)
bg$MAF <- pmin(bg$AF, 1 - bg$AF)

imotif$Group <- "Human-specific iMs"
bg$Group <- "Genomic background"
combined <- rbind(imotif, bg)

# Summary statistics
cat("iM loci count:", nrow(imotif), "\n")
cat("Background loci count:", nrow(bg), "\n")
cat("iM mean MAF:", mean(imotif$MAF), "\n")
cat("Background mean MAF:", mean(bg$MAF), "\n")
wt <- wilcox.test(imotif$MAF, bg$MAF)
cat("Wilcoxon p =", wt$p.value, "\n")

# Save statistical summary
sink("MAF_summary.txt")
cat("=== MAF Statistics ===\n")
cat("iM loci (n =", nrow(imotif), "):\n")
cat("  Mean =", mean(imotif$MAF), "\n")
cat("  Median =", median(imotif$MAF), "\n")
cat("  SD =", sd(imotif$MAF), "\n")
cat("Background loci (n =", nrow(bg), "):\n")
cat("  Mean =", mean(bg$MAF), "\n")
cat("  Median =", median(bg$MAF), "\n")
cat("  SD =", sd(bg$MAF), "\n")
cat("\nWilcoxon rank-sum test p-value =", wt$p.value, "\n")
sink()

# Plot 1: Density curve overlay (recommended main figure)
p1 <- ggplot(combined, aes(x = MAF, color = Group, fill = Group)) +
  geom_density(alpha = 0.2, size = 1.2) +
  scale_color_manual(values = c("Human-specific iMs" = "red",
                                "Genomic background" = "black")) +
  scale_fill_manual(values = c("Human-specific iMs" = "red",
                               "Genomic background" = "grey70")) +
  labs(x = "Minor Allele Frequency (MAF)", y = "Density",
       title = "MAF distribution: human-specific iMs vs. genomic background",
       subtitle = paste0("Wilcoxon P = ", format.pval(wt$p.value, digits = 3),
                         " | Mean MAF: iM = ", round(mean(imotif$MAF), 4),
                         ", Background = ", round(mean(bg$MAF), 4))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Plot 2: Faceted histogram (supplementary details)
p2 <- ggplot(combined, aes(x = MAF, fill = Group)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.02, boundary = 0, alpha = 0.7) +
  facet_grid(Group ~ ., scales = "free_y") +
  scale_fill_manual(values = c("Human-specific iMs" = "red",
                               "Genomic background" = "grey70")) +
  labs(x = "MAF", y = "Density",
       title = "Separate MAF histograms") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 12))

# Combine output
combined_plot <- p1 + p2 + plot_annotation(tag_levels = "A")
ggsave("daf_combined.pdf", combined_plot, width = 14, height = 6)
ggsave("daf_combined.png", combined_plot, width = 14, height = 6, dpi = 300)

# Save density curve separately (for main paper figure)
ggsave("daf_density_overlay.pdf", p1, width = 8, height = 5)
ggsave("daf_density_overlay.png", p1, width = 8, height = 5, dpi = 300)
