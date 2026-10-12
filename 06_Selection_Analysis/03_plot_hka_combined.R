#!/usr/bin/env Rscript
# 03_plot_hka_combined.R
# Purpose: plot HKA test results in promoter regions (three comparisons side by side)
# Input: output/fisher_results_promoter.tsv
# Output: output/Figure_HKA_combined.pdf and .png
# Change: remove titles, add bold uppercase letters A / B / C outside above each subplot

library(data.table)
library(ggplot2)
library(patchwork)  # for combining figures (if not installed, run install.packages("patchwork"))

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
# Root directory for input/output files.
# On the original cluster: Sys.setenv(IMOTIF_BASE = "/datapool/home/2023200496/niulk/my_project")
BASE_DIR <- Sys.getenv("IMOTIF_BASE", unset = ".")

# ----- set working directory -----
workdir <- file.path(BASE_DIR, "19_modified_analysis/09_HKA_PhyloFit")
output_dir <- file.path(workdir, "output")

# ----- read data -----
res_file <- file.path(output_dir, "fisher_results_promoter.tsv")
if (!file.exists(res_file)) {
    stop("File not found: ", res_file)
}
dt <- fread(res_file)

# ensure group order
dt[, Group := factor(Group, levels = c("hominid", "homininae", "hominini"))]

# ----- define colors -----
group_colors <- c("hominid" = "#E6A8B6",
                  "homininae" = "#7A4F2C",
                  "hominini" = "#66C2A5")

# ----- function to create a single subplot -----
# the panel_label argument adds A / B / C outside above each subplot
plot_comparison <- function(data, comparison_name, y_label, panel_label) {
    d <- data[Comparison == comparison_name]
    p <- ggplot(d, aes(x = Group, y = OR, color = Group)) +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                      width = 0.2, linewidth = 0.8, color = "black") +
        geom_point(size = 4, shape = 16) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
        scale_color_manual(values = group_colors) +
        labs(y = y_label, x = NULL) +
        theme_bw() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 12, color = "black"),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "none",
            # remove subplot title (uncomment if you want to add it back later)
            # plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
            # increase top margin to leave space for the outer A / B / C
            plot.margin = margin(25, 8, 5, 5)
        ) +
        scale_y_continuous(limits = c(0, max(d$CI_upper) * 1.2),
                           expand = expansion(mult = c(0.05, 0.05))) +
        # add bold letter outside above the panel
        geom_text(
            aes(x = -Inf, y = Inf, label = panel_label),
            hjust = -0.3,      # left-aligned
            vjust = -0.6,      # negative value = push up outside the panel
            size = 6,          # about 17pt
            fontface = "bold",
            color = "black",
            inherit.aes = FALSE
        ) +
        # allow text to be drawn outside the panel (otherwise it gets clipped)
        coord_cartesian(clip = "off")
    return(p)
}

# generate three subplots (corresponding to A / B / C)
p1 <- plot_comparison(dt, "Stem_vs_Loop", "Odds ratio (Stem vs Loop)", "A")
p2 <- plot_comparison(dt, "Stem_vs_4FD",  "Odds ratio (Stem vs 4FD)",  "B")
p3 <- plot_comparison(dt, "Loop_vs_4FD",  "Odds ratio (Loop vs 4FD)",  "C")

# combine figures (three columns side by side), no overall title added
combined <- (p1 + p2 + p3) +
    plot_layout(ncol = 3)

# ----- save -----
pdf_file <- file.path(output_dir, "Figure_HKA_combined.pdf")
ggsave(pdf_file, combined, width = 15, height = 5, device = "pdf")
cat("PDF saved to", pdf_file, "\n")

png_file <- file.path(output_dir, "Figure_HKA_combined.png")
ggsave(png_file, combined, width = 15, height = 5, dpi = 300, device = "png")
cat("PNG saved to", png_file, "\n")

# print results table
print(dt[, .(Group, Comparison, OR, CI_lower, CI_upper)])
