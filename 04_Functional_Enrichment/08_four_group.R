# ============================================================================
# Four-group KEGG enrichment analysis and specific pathway comparison bar chart script (OUP style adapted)
# Input: four gene list files (hominid_genes.txt, homininae_genes.txt, hominini_genes.txt, humanSpecific_genes.txt)
# Output: bubble plots for each group (p<0.05 and p<0.1), full result CSV, and comparison bar charts for
#         Type II diabetes mellitus and Insulin resistance pathways
# Compliant with OUP illustration guidelines (enlarged font, colorblind-friendly palette, uncompressed TIFF,
# avoids light colors and transparency)
# ============================================================================

# Clean workspace
rm(list = ls())
gc()

# Set working directory (modify according to actual path)
setwd("D:/R/data/kegg_human")
cat("Working directory set to:", getwd(), "\n")

# Input file directory (gene list storage location)
input_dir <- file.path(getwd(), "input")
if (!dir.exists(input_dir)) {
  stop("Please ensure the input directory exists and place the four gene list files in it")
}
cat("Gene list directory:", input_dir, "\n")

# Create output directory
output_dir <- "kegg_results_four_groups"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# ============================================================================
# Step 1: Load required R packages
# ============================================================================
cat("\n=== Step 1: Loading required R packages ===\n")

required_packages <- c(
  "clusterProfiler", "org.Hs.eg.db",
  "ggplot2", "dplyr", "tidyr", "stringr",
  "patchwork", "grid", "ggrepel"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      if (!require("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pkg, character.only = TRUE)
  cat(pkg, "loaded\n")
}

# ============================================================================
# Step 2: Define colors and theme (compliant with OUP requirements)
# ============================================================================

# Four-group color scheme (colorblind-friendly, high contrast)
group_colors <- c(
  "Hominid" = "#1F77B4",       # Blue
  "Homininae" = "#FF7F0E",     # Orange
  "Hominini" = "#2CA02C",      # Green
  "HumanSpecific" = "#D62728"  # Red
)

# OUP theme function: font Arial, line width 0.5pt, text no smaller than 7pt
theme_oup <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      axis.text = element_text(size = rel(1.2), color = "black"),
      axis.title = element_text(size = rel(1.5), face = "plain"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = rel(1.2), face = "plain"),
      legend.text = element_text(size = rel(1)),
      legend.key.size = unit(0.6, "cm"),
      panel.grid.major = element_line(linewidth = 0.3, color = "#CCCCCC", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      plot.title = element_text(size = rel(1.8), face = "bold", hjust = 0.5, margin = margin(b = 15)),
      plot.margin = margin(15, 20, 15, 15),
      text = element_text(family = base_family)
    )
}

# Bubble plot function (modified color gradient to avoid light colors)
create_bubble_plot <- function(enrichment_df,
                               title = NULL,
                               color_palette = "#1F77B4",
                               top_n = 15,
                               p_cutoff = 0.05) {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)
  sig_df <- enrichment_df %>% filter(p.adjust < p_cutoff) %>% arrange(p.adjust)
  if (nrow(sig_df) == 0) {
    cat("  Warning: p.adjust <", p_cutoff, "no significant pathways\n")
    return(NULL)
  }
  plot_df <- sig_df %>%
    head(top_n) %>%
    mutate(
      GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"),
                                        function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      log10_padj = -log10(p.adjust),
      Description_short = ifelse(nchar(Description) > 50,
                                 paste0(substr(Description, 1, 47), "..."),
                                 Description),
      Description_short = factor(Description_short, levels = rev(Description_short))
    )
  plot_height <- max(5, nrow(plot_df) * 0.45 + 3)

  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_short)) +
    geom_point(aes(size = Count, fill = log10_padj),
               shape = 21, color = "black", stroke = 0.3) +
    scale_size_continuous(range = c(4, 12), name = "Gene count",
                          guide = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) +
    # Avoid pure white, use light gray as minimum value
    scale_fill_gradient(low = "#F7F7F7", high = color_palette,
                        name = expression(-log[10]("p.adjust")),
                        guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                               barwidth = unit(5, "cm"), barheight = unit(0.4, "cm"))) +
    labs(x = "Gene ratio", y = NULL, title = title) +
    theme_oup(base_size = 12) +
    theme(axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16, margin = margin(t = 10)),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15)))
  return(list(plot = p, height = plot_height))
}

# ============================================================================
# Step 3: Read gene lists
# ============================================================================
cat("\n=== Step 3: Reading gene lists ===\n")

read_gene_list <- function(filename) {
  file_path <- file.path(input_dir, filename)
  if (!file.exists(file_path)) {
    cat("Error: Gene file not found", file_path, "\n")
    return(NULL)
  }
  genes <- readLines(file_path)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)  # Remove version number
  return(genes)
}

groups <- list(
  Hominid = list(file = "hominid_genes.txt", name = "Hominid", color = group_colors["Hominid"]),
  Homininae = list(file = "homininae_genes.txt", name = "Homininae", color = group_colors["Homininae"]),
  Hominini = list(file = "hominini_genes.txt", name = "Hominini", color = group_colors["Hominini"]),
  HumanSpecific = list(file = "humanSpecific_genes.txt", name = "HumanSpecific", color = group_colors["HumanSpecific"])
)

gene_lists <- list()
for (g in names(groups)) {
  genes <- read_gene_list(groups[[g]]$file)
  if (!is.null(genes)) {
    gene_lists[[g]] <- genes
    cat(g, "gene count:", length(genes), "\n")
  } else {
    cat(g, "gene list missing, skipping\n")
  }
}

if (length(gene_lists) == 0) stop("No gene lists successfully read, please check files in input directory")

# ============================================================================
# Step 4: Gene ID conversion
# ============================================================================
cat("\n=== Step 4: Gene ID conversion ===\n")

convert_genes <- function(genes, list_name) {
  if (length(genes) == 0) return(NULL)
  gene_df <- tryCatch({
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat(list_name, "Gene conversion error:", e$message, "\n")
    return(NULL)
  })
  if (!is.null(gene_df) && nrow(gene_df) > 0) {
    cat(list_name, ": Successfully converted", nrow(gene_df), "/", length(genes), "genes\n")
    return(gene_df$ENTREZID)
  } else {
    cat(list_name, ": Gene conversion failed\n")
    return(NULL)
  }
}

entrez_lists <- list()
for (g in names(gene_lists)) {
  entrez <- convert_genes(gene_lists[[g]], g)
  if (!is.null(entrez)) entrez_lists[[g]] <- entrez
}

if (length(entrez_lists) == 0) stop("No genes successfully converted")

# ============================================================================
# Step 5: KEGG enrichment analysis
# ============================================================================
cat("\n=== Step 5: KEGG enrichment analysis ===\n")

run_kegg_analysis <- function(entrez_ids, list_name) {
  if (length(entrez_ids) < 5) {
    cat(list_name, ": Insufficient genes (<5), skipping KEGG analysis\n")
    return(NULL)
  }
  cat(list_name, ": Running KEGG enrichment analysis...\n")
  kegg_result <- tryCatch({
    enrichKEGG(gene = unique(entrez_ids),
               organism = "hsa",
               keyType = "kegg",
               pvalueCutoff = 0.2,
               pAdjustMethod = "BH",
               qvalueCutoff = 0.25,
               minGSSize = 5,
               maxGSSize = 500,
               use_internal_data = FALSE)
  }, error = function(e) {
    cat(list_name, "KEGG analysis error:", e$message, "\n")
    return(NULL)
  })
  if (is.null(kegg_result) || nrow(kegg_result) == 0) {
    cat(list_name, ": No enriched pathways\n")
    return(NULL)
  }
  df <- as.data.frame(kegg_result)
  cat(list_name, ": Enriched", nrow(df), "pathways\n")
  return(df)
}

kegg_results <- list()
for (g in names(entrez_lists)) {
  kegg_results[[g]] <- run_kegg_analysis(entrez_lists[[g]], g)
  if (!is.null(kegg_results[[g]])) {
    write.csv(kegg_results[[g]], file.path(output_dir, paste0(g, "_kegg_full.csv")), row.names = FALSE)
  }
}

# ============================================================================
# Step 6: Generate bubble plots (p<0.05 and p<0.1) and save as PDF and TIFF (uncompressed)
# ============================================================================
cat("\n=== Step 6: Generating bubble plots ===\n")

thresholds <- c(0.05, 0.10)

for (g in names(kegg_results)) {
  df <- kegg_results[[g]]
  if (is.null(df)) next
  n_genes <- length(entrez_lists[[g]])
  for (thresh in thresholds) {
    res <- create_bubble_plot(
      enrichment_df = df,
      title = paste0(g, " pG4 Genes\n(n = ", n_genes, " genes, p < ", thresh, ")"),
      color_palette = groups[[g]]$color,
      top_n = 15,
      p_cutoff = thresh
    )
    if (!is.null(res)) {
      file_suffix <- ifelse(thresh == 0.05, "main", "supp")
      pdf_file <- file.path(output_dir, paste0(g, "_bubble_", file_suffix, ".pdf"))
      tiff_file <- file.path(output_dir, paste0(g, "_bubble_", file_suffix, ".tiff"))
      # PDF using cairo_pdf for embedded fonts
      ggsave(pdf_file, res$plot, width = 10, height = res$height + 0.5, device = cairo_pdf, dpi = 300, limitsize = FALSE)
      # TIFF uncompressed, 600dpi
      ggsave(tiff_file, res$plot, width = 10, height = res$height + 0.5, device = "tiff", dpi = 600, compression = "none", limitsize = FALSE)
      cat(g, "p <", thresh, "bubble plot saved\n")
    }
  }
}

# ============================================================================
# Step 7: Specific pathway comparison bar chart (complete version)
# Includes: building comparison data, renaming groups, color mapping, adding -log10(p.adjust),
#           row arrangement, column arrangement, horizontal bar chart (three layouts)
# ============================================================================
# ============================================================================
# ============================================================================
# Step 7: Specific pathway comparison bar chart (final version)
# Function: build comparison data, draw row/column/horizontal bar charts
# Features: display only significance stars, narrower bars, vertically elongated, adjustable legend spacing and margins
# ============================================================================
cat("\n=== Step 7: Specific pathway comparison bar chart ===\n")

# ----- 0. Load required packages -----
library(ggplot2)
library(reshape2)   # For data reshaping

# ----- 1. Define target pathways (ID and description) -----
target_ids <- c("hsa04930", "hsa04931")
target_names <- c("hsa04930" = "Type II diabetes mellitus",
                  "hsa04931" = "Insulin resistance")

# ----- 2. Build comparison data frame (extract from kegg_results) -----
comparison_data <- data.frame()

for (g in names(kegg_results)) {
  df <- kegg_results[[g]]
  if (is.null(df)) next

  cat("\nProcessing group:", g, "(total", nrow(df), "pathways)\n")

  for (id in target_ids) {
    pw_row <- df[df$ID == id, ]
    if (nrow(pw_row) > 0) {
      # Ensure FoldEnrichment column exists
      if (!"FoldEnrichment" %in% colnames(pw_row)) {
        gene_ratio <- as.numeric(sapply(strsplit(pw_row$GeneRatio, "/"),
                                        function(x) as.numeric(x[1]) / as.numeric(x[2])))
        bg_ratio <- as.numeric(sapply(strsplit(pw_row$BgRatio, "/"),
                                      function(x) as.numeric(x[1]) / as.numeric(x[2])))
        pw_row$FoldEnrichment <- gene_ratio / bg_ratio
      }

      fold <- pw_row[[1, "FoldEnrichment"]]
      padj <- pw_row[[1, "p.adjust"]]
      cnt <- pw_row[[1, "Count"]]

      cat("  Found pathway:", target_names[id],
          "Fold enrichment =", fold, "p.adjust =", padj, "Gene count =", cnt, "\n")

      comparison_data <- rbind(comparison_data, data.frame(
        Group = g,
        Pathway = target_names[id],
        FoldEnrichment = fold,
        p_adjust = padj,
        GeneCount = cnt,
        stringsAsFactors = FALSE
      ))
    } else {
      cat("  Pathway not found:", target_names[id], "\n")
    }
  }
}

# ----- 3. Data cleaning and group renaming -----
comparison_data$Group <- as.character(comparison_data$Group)
comparison_data <- comparison_data[!is.na(comparison_data$FoldEnrichment), ]
comparison_data$Group[comparison_data$Group == "Hominid"] <- "Great ape"

# Custom group order (remove Hominini if no data)
group_order <- c("Great ape", "Homininae", "HumanSpecific")
comparison_data$Group <- factor(comparison_data$Group, levels = group_order)

# Pathway factor order
comparison_data$Pathway <- factor(comparison_data$Pathway,
                                  levels = c("Type II diabetes mellitus",
                                             "Insulin resistance"))

# ----- 4. Define colors -----
group_colors <- c(
  "Great ape"   = "#1f77b4",   # Blue
  "Homininae"   = "#9467bd",   # Purple
  "HumanSpecific" = "#ffd966"  # Yellow
)

# ----- 5. Add significance stars -----
comparison_data$sig <- ifelse(comparison_data$p_adjust < 0.001, "***",
                              ifelse(comparison_data$p_adjust < 0.01, "**",
                                     ifelse(comparison_data$p_adjust < 0.05, "*", "ns")))

# ----- 6. Custom theme function (adjustable legend spacing and margins) -----
# Parameter description:
#   base_size: base font size
#   legend_spacing_x: horizontal spacing between legend entries (unit: "lines", larger value = more spacing)
#   legend_margin: outer margin of the legend box, format margin(t, r, b, l) (unit: pt)
#   plot_margin: overall plot margin, format margin(t, r, b, l)
theme_oup <- function(base_size = 12,
                      legend_spacing_x = 5.5,      # Default 0.5 line height
                      legend_margin = margin(20, 20, 20, 20),  # 20pt on each side
                      plot_margin = margin(15, 20, 15, 35)) {  # top, right, bottom, left
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      # Horizontal spacing between legend entries
      legend.spacing.x = unit(legend_spacing_x, "lines"),
      # Legend box outer margin
      legend.margin = legend_margin,
      # Overall plot margin
      plot.margin = plot_margin
    )
}

# ========================
# Version A: Row arrangement (side-by-side)
# ========================
p_bar_row <- ggplot(comparison_data, aes(x = Group, y = FoldEnrichment, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Fold enrichment") +
  theme_oup(base_size = 12, legend_spacing_x = 0.5,
            legend_margin = margin(20, 20, 20, 20),
            plot_margin = margin(15, 20, 15, 15)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  geom_text(aes(label = sig, y = FoldEnrichment + 0.12 * max(FoldEnrichment, na.rm = TRUE)),
            position = position_dodge(width = 0.8),
            size = 4, color = "black", fontface = "bold")

# ========================
# Version B: Column arrangement (stacked vertically) -- vertically elongated
# ========================
p_bar_col <- ggplot(comparison_data, aes(x = Group, y = FoldEnrichment, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = group_colors) +
  labs(x = NULL, y = "Fold enrichment") +
  theme_oup(base_size = 12, legend_spacing_x = 0.5,
            legend_margin = margin(20, 20, 20, 20),
            plot_margin = margin(15, 20, 15, 15)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  geom_text(aes(label = sig, y = FoldEnrichment + 0.12 * max(FoldEnrichment, na.rm = TRUE)),
            position = position_dodge(width = 0.8),
            size = 4, color = "black", fontface = "bold")

# ========================
# Version C: Horizontal bar chart (coord_flip) + vertical arrangement
# ========================
p_bar_horiz <- ggplot(comparison_data, aes(x = Group, y = FoldEnrichment, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.6, color = "black", linewidth = 0.3) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = group_colors) +
  labs(y = "Fold enrichment", x = NULL) +
  theme_oup(base_size = 12, legend_spacing_x = 10.5,
            legend_margin = margin(50, 50, 50, 50),
            plot_margin = margin(15, 20, 15, 35)) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  geom_text(aes(label = sig, y = FoldEnrichment + 0.12 * max(FoldEnrichment, na.rm = TRUE)),
            position = position_dodge(width = 0.8),
            size = 4, color = "black", fontface = "bold", hjust = -0.1)

# ----- 7. Save plots (specify output directory, modify according to actual path) -----
output_dir <- getwd()  # Can be changed to specific path, e.g., "D:/R/data/kegg_human/kegg_results_four_groups"

# Row arrangement
ggsave(file.path(output_dir, "pathway_comparison_barplot_row.pdf"),
       p_bar_row, width = 10, height = 6, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_row.tiff"),
       p_bar_row, width = 10, height = 6, device = "tiff", dpi = 600, compression = "none")

# Column arrangement (vertically elongated)
ggsave(file.path(output_dir, "pathway_comparison_barplot_col.pdf"),
       p_bar_col, width = 6, height = 10, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_col.tiff"),
       p_bar_col, width = 6, height = 10, device = "tiff", dpi = 600, compression = "none")

# Horizontal bar chart (vertically elongated)
ggsave(file.path(output_dir, "pathway_comparison_barplot_horiz.pdf"),
       p_bar_horiz, width = 5, height = 8, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_horiz.tiff"),
       p_bar_horiz, width = 5, height = 8, device = "tiff", dpi = 600, compression = "none")

cat("\nSpecific pathway comparison bar charts (three layouts) saved to:", output_dir, "\n")

# ----- 8. Save fold enrichment matrix (optional) -----
if (nrow(comparison_data) > 0) {
  mat_wide <- reshape(comparison_data[, c("Group", "Pathway", "FoldEnrichment")],
                      idvar = "Group", timevar = "Pathway", direction = "wide")
  colnames(mat_wide) <- gsub("FoldEnrichment\\.", "", colnames(mat_wide))
  write.csv(mat_wide, file.path(output_dir, "pathway_fold_enrichment_matrix.csv"),
            row.names = FALSE)
  cat("Fold enrichment matrix saved\n")
}

cat("\n=== Step 7 complete ===\n")

cat("\n=== Step 7 complete ===\n")
# ============================================================================
# Complete
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Four-group KEGG enrichment analysis and comparison bar chart generation complete!\n")
cat("Results saved in:", normalizePath(output_dir), "\n")
cat("Note: TIFF images are uncompressed format, compliant with OUP requirements; bar charts have transparency removed for clarity.\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
