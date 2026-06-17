#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

cat("=========================================\n")
cat("Result analysis start time:", format(Sys.time()), "\n")
cat("=========================================\n")

# Read species list
species_list <- fread("sp240_nonhuman.txt", header = FALSE)$V1
cat("Species count:", length(species_list), "\n")

# Read iMotif info
imotif_info <- fread("imotif_clean_final.bed", header = FALSE)
colnames(imotif_info) <- c("chr", "start", "end", "id", "score", "strand")
cat("iMotif count:", nrow(imotif_info), "\n")

# Check processed species
ar_files <- list.files("results/align_ratio", pattern = "*.imotif.ar.txt", full.names = TRUE)
cat("Processed species:", length(ar_files), "/240\n")

if(length(ar_files) < 240) {
    warning("Some species still unfinished, please wait for all tasks to complete")
}

# Build alignment rate matrix
cat("\nBuilding alignment rate matrix...\n")
ar_list <- list()

for(ar_file in ar_files) {
    species <- gsub(".imotif.ar.txt", "", basename(ar_file))
    ar_data <- fread(ar_file, header = FALSE)
    if(nrow(ar_data) > 0) {
        colnames(ar_data) <- c("id", "size", "covered", "sum", "mean0", "mean")
        ar_list[[species]] <- ar_data[, .(id, mean)]
    }
}

# 合并矩阵
ar_matrix <- Reduce(function(x, y) merge(x, y, by = "id", all = TRUE), ar_list)
colnames(ar_matrix) <- c("id", names(ar_list))
ar_matrix[is.na(ar_matrix)] <- 0

# Save matrix
fwrite(ar_matrix, "matrix/imotif_alignment_matrix.txt", sep = "\t")
cat("Matrix dimensions:", dim(ar_matrix), "\n")

# Calculate N1 and N2
cat("\nCalculating N1 and N2...\n")
mat <- as.matrix(ar_matrix[, -1])
rownames(mat) <- ar_matrix$id

N1 <- rowSums(mat >= 0.9, na.rm = TRUE)
N2 <- rowSums(mat <= 0.1, na.rm = TRUE)

# Group
result <- data.frame(
    id = ar_matrix$id,
    N1 = N1,
    N2 = N2,
    group = case_when(
        N1 >= 120 & N2 <= 25 ~ "C1",
        N1 >= 20 & N1 <= 50 & N2 <= 120 ~ "C2",
        N1 <= 50 & N2 >= 180 ~ "C3",
        TRUE ~ "Other"
    )
) %>% left_join(imotif_info, by = "id")

# Save results
fwrite(result, "imotif_conservation_results.tsv", sep = "\t")
cat("Results saved\n")

# Statistics
stats <- result %>% 
    group_by(group) %>% 
    summarise(
        count = n(),
        percentage = n() / nrow(.) * 100
    )
cat("\nGroup statistics:\n")
print(stats)

# Plot scatter plot
cat("\nPlotting scatter plot...\n")
p1 <- ggplot(result, aes(x = N2, y = N1)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_viridis_c() +
    geom_segment(x = 0, y = 120, xend = 25, yend = 120, linetype = 2, color = "red") +
    geom_segment(x = 25, y = 120, xend = 25, yend = 240, linetype = 2, color = "red") +
    geom_segment(x = 0, y = 20, xend = 120, yend = 20, linetype = 2, color = "blue") +
    geom_segment(x = 0, y = 50, xend = 120, yend = 50, linetype = 2, color = "blue") +
    geom_segment(x = 120, y = 20, xend = 120, yend = 50, linetype = 2, color = "blue") +
    geom_segment(x = 180, y = 0, xend = 180, yend = 50, linetype = 2, color = "green") +
    geom_segment(x = 180, y = 50, xend = 240, yend = 50, linetype = 2, color = "green") +
    coord_cartesian(xlim = c(0, 240), ylim = c(0, 240)) +
    labs(x = "N2 (≤10% alignment)", y = "N1 (≥90% alignment)",
         title = paste0("iMotif Conservation (n=", nrow(result), ")")) +
    theme_bw()

ggsave("plots/N1N2_density.png", p1, width = 8, height = 7, dpi = 300)

# 绘制分组柱状图
p2 <- ggplot(stats, aes(x = group, y = percentage, fill = group)) +
    geom_col(width = 0.6, color = "black") +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
    scale_fill_manual(values = c("C1" = "#DA635C", "C2" = "#5FB6B9", 
                                 "C3" = "#6387BB", "Other" = "#B0BDB0")) +
    labs(x = "Group", y = "Percentage (%)") +
    theme_bw() +
    theme(legend.position = "none")

ggsave("plots/group_barplot.png", p2, width = 6, height = 5, dpi = 300)

cat("\nFigures saved to plots/ directory\n")
cat("=========================================\n")
cat("Completion time:", format(Sys.time()), "\n")
cat("=========================================\n")
