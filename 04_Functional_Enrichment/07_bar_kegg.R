# ============================================================================
# Bar chart version: KEGG enrichment analysis bar chart generation script (OUP guideline compliant, enlarged font)
# Improvements: take top15 significant pathways by p-value, then display sorted by GeneRatio (longest bar at top)
# Used for generating p<0.05 main plots and p<0.1 supplementary plots for hominid and humanSpecific gene lists
# ============================================================================

# Clean workspace
rm(list = ls())
gc()

# Set working directory (modify according to actual path)
setwd("D:/R/data/kegg_human")
cat("Working directory set to:", getwd(), "\n")

# ============================================================================
# Step 1: Load required R packages
# ============================================================================

cat("\n=== Step 1: Loading required R packages ===\n")

required_packages <- c(
  "clusterProfiler", "org.Hs.eg.db",
  "ggplot2", "dplyr", "tidyr", "stringr",
  "patchwork", "grid"
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
# Step 2: Create output directory
# ============================================================================

cat("\n=== Step 2: Creating output directory ===\n")

output_dir <- "bar_plots_dual_oup_largefont"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created directory:", output_dir, "\n")
}

# ============================================================================
# Step 3: Define OUP style theme (enlarged font version)
# ============================================================================

# Nature color scheme (colorblind-friendly)
nature_colors <- c(
  blue = "#1F77B4",
  orange = "#FF7F0E",
  green = "#2CA02C",
  red = "#D62728"
)

# OUP theme function: base font 12pt, line width 0.3-0.5pt
theme_oup <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis lines
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # Tick labels
      axis.text = element_text(size = rel(1.2), color = "black"),  # approx 14.4pt
      axis.title = element_text(size = rel(1.5), face = "plain"),  # approx 18pt
      # Legend
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = rel(1.2), face = "plain"), # approx 14.4pt
      legend.text = element_text(size = rel(1)),                     # 12pt
      legend.key.size = unit(0.6, "cm"),
      # Grid lines
      panel.grid.major = element_line(linewidth = 0.3, color = "#CCCCCC", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      # Panel border
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # Title
      plot.title = element_text(size = rel(1.8), face = "bold", hjust = 0.5, margin = margin(b = 15)), # approx 21.6pt
      # Plot margin
      plot.margin = margin(15, 20, 15, 15),
      text = element_text(family = base_family)
    )
}

# Enhanced bar chart function (supports custom p threshold and OUP theme)
# Improvement: first take top_n most significant pathways by p-value, then sort by GeneRatio ascending (longest bar at top)
create_bar_plot <- function(enrichment_df,
                            title = NULL,
                            color_palette = nature_colors["blue"],
                            top_n = 15,
                            p_cutoff = 0.05) {

  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }

  # Filter significant pathways (p.adjust < p_cutoff)
  sig_df <- enrichment_df %>%
    filter(p.adjust < p_cutoff) %>%
    arrange(p.adjust)   # Sort by p-value ascending (most significant first)

  if (nrow(sig_df) == 0) {
    cat("  Warning: p.adjust <", p_cutoff, "no significant pathways\n")
    return(NULL)
  }

  cat("  Significant pathways (p<", p_cutoff, "):", nrow(sig_df), "\n")

  # Take top_n pathways (most significant)
  plot_df <- sig_df %>%
    head(top_n) %>%
    mutate(
      GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"),
                                        function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      log10_padj = -log10(p.adjust),
      Description_short = ifelse(nchar(Description) > 50,
                                 paste0(substr(Description, 1, 47), "..."),
                                 Description)
    )

  # Key modification: sort by GeneRatio_num ascending so bar length increases from top to bottom
  plot_df <- plot_df %>%
    arrange(GeneRatio_num) %>%
    mutate(Description_short = factor(Description_short, levels = Description_short))

  # Dynamically calculate plot height (0.45 inches per row, plus title and legend space)
  plot_height <- max(5, nrow(plot_df) * 0.45 + 3)

  # Create bar chart (apply OUP theme)
  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_short)) +
    geom_col(aes(fill = log10_padj),
             width = 0.7, color = "black", linewidth = 0.3) +  # Bar border 0.3pt
    scale_fill_gradient(
      low = "white",
      high = color_palette,
      name = expression(-log[10]("p.adjust")),
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(5, "cm"),
        barheight = unit(0.4, "cm")
      )
    ) +
    labs(x = "Gene ratio", y = NULL, title = title) +
    theme_oup(base_size = 12) +  # Base font 12pt
    theme(
      axis.text.y = element_text(size = 14),          # 14pt
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))

  return(list(plot = p, height = plot_height))
}

# ============================================================================
# Step 4: Read and process gene lists (unchanged)
# ============================================================================

cat("\n=== Step 4: Reading and processing gene lists ===\n")

read_gene_list <- function(filename) {
  if (!file.exists(filename)) {
    cat("Error: Gene file not found", filename, "\n")
    return(NULL)
  }

  genes <- readLines(filename)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)  # Remove version number

  return(genes)
}

# Read gene lists
hominid_genes <- read_gene_list("hominid_genes.txt")
humanSpecific_genes <- read_gene_list("humanSpecific_genes.txt")

if (is.null(hominid_genes) && is.null(humanSpecific_genes)) {
  stop("Error: Both gene list files not found")
}

cat("Gene list statistics:\n")
if (!is.null(hominid_genes)) {
  cat("  Great ape shared genes:", length(hominid_genes), "\n")
}
if (!is.null(humanSpecific_genes)) {
  cat("  Human-specific genes:", length(humanSpecific_genes), "\n")
}

# ============================================================================
# Step 5: Gene ID conversion (unchanged)
# ============================================================================

cat("\n=== Step 5: Gene ID conversion ===\n")

convert_genes <- function(genes, list_name) {
  if (is.null(genes) || length(genes) == 0) {
    cat(list_name, ": No gene data\n")
    return(NULL)
  }

  gene_df <- tryCatch({
    bitr(genes,
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
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

hominid_entrez <- convert_genes(hominid_genes, "Great ape shared")
humanSpecific_entrez <- convert_genes(humanSpecific_genes, "Human-specific")

# ============================================================================
# Step 6: KEGG enrichment analysis (unchanged)
# ============================================================================

cat("\n=== Step 6: KEGG enrichment analysis ===\n")

run_kegg_analysis <- function(entrez_ids, list_name) {
  if (is.null(entrez_ids) || length(entrez_ids) < 5) {
    cat(list_name, ": Insufficient genes, skipping KEGG analysis\n")
    return(NULL)
  }

  cat(list_name, ": Running KEGG enrichment analysis...\n")

  kegg_result <- tryCatch({
    enrichKEGG(
      gene = unique(entrez_ids),
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.2,        # Relax threshold to capture more pathways
      pAdjustMethod = "BH",
      qvalueCutoff = 0.25,
      minGSSize = 5,
      maxGSSize = 500,
      use_internal_data = FALSE
    )
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

hominid_kegg <- run_kegg_analysis(hominid_entrez, "Great ape shared")
humanSpecific_kegg <- run_kegg_analysis(humanSpecific_entrez, "Human-specific")

# Save full results to CSV
if (!is.null(hominid_kegg)) {
  write.csv(hominid_kegg, file.path(output_dir, "hominid_kegg_full.csv"), row.names = FALSE)
}
if (!is.null(humanSpecific_kegg)) {
  write.csv(humanSpecific_kegg, file.path(output_dir, "humanSpecific_kegg_full.csv"), row.names = FALSE)
}

# ============================================================================
# Step 7: Generate bar charts (p<0.05 main plot and p<0.1 supplementary) and save as PDF (editable) and TIFF
# ============================================================================

cat("\n=== Step 7: Generating bar charts ===\n")

thresholds <- c(0.05, 0.10)

# Great ape shared bar chart
if (!is.null(hominid_kegg) && nrow(hominid_kegg) > 0) {
  for (thresh in thresholds) {
    res <- create_bar_plot(
      enrichment_df = hominid_kegg,
      title = paste0("Great Ape-Shared pG4 Genes\n(n = ",
                     length(hominid_entrez), " genes, p < ", thresh, ")"),
      color_palette = nature_colors["blue"],
      top_n = 15,
      p_cutoff = thresh
    )

    if (!is.null(res)) {
      file_suffix <- ifelse(thresh == 0.05, "main", "supp")

      # Save PDF (cairo_pdf ensures editable text)
      ggsave(
        filename = file.path(output_dir, paste0("hominid_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,                # Width increased to 10 inches
        height = res$height + 0.5, # Height slightly increased
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )

      # Save TIFF (600dpi, LZW compression, for printing)
      ggsave(
        filename = file.path(output_dir, paste0("hominid_bar_", file_suffix, ".tiff")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = "tiff",
        dpi = 600,
        compression = "lzw",
        limitsize = FALSE
      )

      cat("Great ape shared p <", thresh, "bar chart saved (PDF+TIFF)\n")
    }
  }
}

# Human-specific bar chart
if (!is.null(humanSpecific_kegg) && nrow(humanSpecific_kegg) > 0) {
  for (thresh in thresholds) {
    res <- create_bar_plot(
      enrichment_df = humanSpecific_kegg,
      title = paste0("Human-Specific pG4 Genes\n(n = ",
                     length(humanSpecific_entrez), " genes, p < ", thresh, ")"),
      color_palette = nature_colors["orange"],
      top_n = 15,
      p_cutoff = thresh
    )

    if (!is.null(res)) {
      file_suffix <- ifelse(thresh == 0.05, "main", "supp")

      ggsave(
        filename = file.path(output_dir, paste0("humanSpecific_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )

      ggsave(
        filename = file.path(output_dir, paste0("humanSpecific_bar_", file_suffix, ".tiff")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = "tiff",
        dpi = 600,
        compression = "lzw",
        limitsize = FALSE
      )

      cat("Human-specific p <", thresh, "bar chart saved (PDF+TIFF)\n")
    }
  }
}

# ============================================================================
# Step 8: Generate side-by-side comparison plot (p<0.05 main plot only)
# ============================================================================

if (!is.null(hominid_kegg) && !is.null(humanSpecific_kegg)) {
  res_hom <- create_bar_plot(hominid_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["blue"])
  res_hum <- create_bar_plot(humanSpecific_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["orange"])

  if (!is.null(res_hom) && !is.null(res_hum)) {
    # Remove individual titles, add unified top title
    p_hom <- res_hom$plot + labs(title = NULL) + theme(plot.margin = margin(15, 15, 15, 20))
    p_hum <- res_hum$plot + labs(title = NULL) + theme(plot.margin = margin(15, 20, 15, 15))

    combined_height <- max(res_hom$height, res_hum$height) + 1.0  # Reserve space for title

    # Combine using patchwork
    combined_plot <- (p_hom | p_hum) +
      plot_annotation(
        title = "KEGG Pathway Enrichment (p.adjust < 0.05)",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
      ) &
      theme(legend.position = "bottom")
    combined_plot <- combined_plot + theme(plot.margin = margin(r = 100))

    # Save PDF
    ggsave(
      filename = file.path(output_dir, "combined_bar_main.pdf"),
      plot = combined_plot,
      width = 18,
      height = combined_height,
      device = cairo_pdf,
      dpi = 300,
      limitsize = FALSE
    )

    # Save TIFF
    ggsave(
      filename = file.path(output_dir, "combined_bar_main.tiff"),
      plot = combined_plot,
      width = 18,
      height = combined_height,
      device = "tiff",
      dpi = 600,
      compression = "lzw",
      limitsize = FALSE
    )

    cat("Side-by-side comparison plot (main) saved (PDF+TIFF)\n")
  } else {
    cat("Side-by-side plot skipped: at least one dataset has no p<0.05 significant pathways\n")
  }
}

# ============================================================================
# Complete!
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Bar chart generation complete (OUP guideline compliant, enlarged font)!\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\nAll results saved in:", normalizePath(output_dir), "\n")
