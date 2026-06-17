# ============================================================================
# KEGG enrichment analysis bar chart generation script (dual threshold version, enlarged font, OUP guideline compliant)
# Template strand downstream vs non-template strand upstream genes
# Function: generate bar charts for p.adjust < 0.05 (main) and p.adjust < 0.1 (supplementary)
#           Automatically output pathway counts for each threshold, complete legend, adaptive height
#           All text significantly enlarged, suitable for scaling and editing in AI
#           Bars sorted by gene ratio ascending (longest bar at top)
# ============================================================================

# Clean workspace
rm(list = ls())
gc()

# Set working directory (modify according to actual path)
setwd("D:/R/data/kegg_human/tem_nontem")
cat("Working directory set to:", getwd(), "\n")

# ============================================================================
# Step 1: Load required R packages
# ============================================================================

required_packages <- c(
  "clusterProfiler", "org.Hs.eg.db",
  "ggplot2", "dplyr", "tidyr", "stringr",
  "patchwork", "grid"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  library(pkg, character.only = TRUE)
}

# ============================================================================
# Step 2: Create output directory
# ============================================================================

output_dir <- "bar_plots_dual_threshold_v3_oup_largefont"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Output directory:", output_dir, "\n")

# ============================================================================
# Step 3: Define OUP style theme (enlarged font version)
# ============================================================================

# Colorblind-friendly palette
nature_colors <- c(
  blue = "#1F77B4",
  orange = "#FF7F0E"
)

# OUP theme function (base font 14pt, proportional scaling)
theme_oup <- function(base_size = 14, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Axis lines
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # Tick labels
      axis.text = element_text(size = rel(1), color = "black"),      # 14pt
      axis.title = element_text(size = rel(1.2), face = "plain"),    # approx 16.8pt
      # Legend
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = rel(1.1), face = "plain"),  # approx 15.4pt
      legend.text = element_text(size = rel(1)),                      # 14pt
      legend.key.size = unit(0.7, "cm"),
      # Grid lines
      panel.grid.major = element_line(linewidth = 0.3, color = "#CCCCCC", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      # Panel border
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # Title
      plot.title = element_text(size = rel(1.4), face = "bold", hjust = 0.5, margin = margin(b = 15)), # approx 19.6pt
      # Plot margin
      plot.margin = margin(20, 60, 20, 60),
      text = element_text(family = base_family)
    )
}

# Bar chart plotting function (enlarged font version, sorted by GeneRatio ascending)
create_bar_plot <- function(enrichment_df,
                            title = "KEGG Pathway Enrichment",
                            color_palette = nature_colors["blue"],
                            top_n = 10,
                            p_cutoff = 0.05) {

  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)

  # Filter significant pathways at specified threshold
  sig_df <- enrichment_df %>%
    filter(p.adjust < p_cutoff) %>%
    arrange(p.adjust)

  if (nrow(sig_df) == 0) {
    cat("Warning: p.adjust <", p_cutoff, "no significant pathways\n")
    return(NULL)
  }

  cat("  p.adjust <", p_cutoff, ": found", nrow(sig_df), "significant pathways\n")

  # Take top_n pathways (most significant)
  plot_df <- sig_df %>%
    head(top_n) %>%
    mutate(
      GeneRatio_num = as.numeric(sapply(strsplit(GeneRatio, "/"),
                                        function(x) as.numeric(x[1]) / as.numeric(x[2]))),
      log10_padj = -log10(p.adjust),
      Description_short = ifelse(nchar(Description) > 60,
                                 paste0(substr(Description, 1, 57), "..."),
                                 Description)
    )

  # Key modification: sort by GeneRatio_num ascending so bar length increases from top to bottom (longest bar at top)
  plot_df <- plot_df %>%
    arrange(GeneRatio_num) %>%
    mutate(Description_short = factor(Description_short, levels = Description_short))

  # Dynamically calculate plot height (font enlarged, increase row height to 0.45 inches)
  plot_height <- max(5, nrow(plot_df) * 0.45 + 2.5)

  # Bar chart (using geom_col)
  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_short)) +
    geom_col(aes(fill = log10_padj),
             width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = color_palette,
                        name = expression(-log[10]("p.adjust")),
                        guide = guide_colorbar(barwidth = unit(5, "cm"),
                                               barheight = unit(0.4, "cm"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
    labs(x = "Gene ratio", y = NULL, title = title) +
    theme_oup(base_size = 14) +  # Base font 14pt
    theme(
      axis.text.y = element_text(size = 14),          # 14pt
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))  # Leave a small space on the right

  return(list(plot = p, height = plot_height))
}

# ============================================================================
# Step 4: Read gene lists
# ============================================================================

read_gene_list <- function(filename) {
  if (!file.exists(filename)) stop("File does not exist: ", filename)
  genes <- readLines(filename)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)
  return(genes)
}

template_genes <- read_gene_list("template_downstream_all_genes.txt")
nontemplate_genes <- read_gene_list("nontemplate_upstream_all_genes.txt")

cat("Template strand downstream gene count:", length(template_genes), "\n")
cat("Non-template strand upstream gene count:", length(nontemplate_genes), "\n")

# ============================================================================
# Step 5: Gene ID conversion
# ============================================================================

convert_to_entrez <- function(genes, list_name) {
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(list_name, "successfully converted", nrow(gene_df), "/", length(genes), "genes\n")
  return(gene_df$ENTREZID)
}

template_entrez <- convert_to_entrez(template_genes, "Template downstream")
nontemplate_entrez <- convert_to_entrez(nontemplate_genes, "Non-template upstream")

# ============================================================================
# Step 6: KEGG enrichment analysis (using relaxed thresholds to capture all possible pathways)
# ============================================================================

run_kegg <- function(entrez_ids, list_name) {
  if (length(entrez_ids) < 5) {
    cat(list_name, "insufficient genes (<5), skipping KEGG analysis\n")
    return(NULL)
  }

  kegg <- enrichKEGG(
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

  if (is.null(kegg) || nrow(kegg) == 0) {
    cat(list_name, "no enriched pathways\n")
    return(NULL)
  }

  df <- as.data.frame(kegg)
  sig_05 <- sum(df$p.adjust < 0.05)
  sig_10 <- sum(df$p.adjust < 0.10)
  cat(list_name, "enriched", nrow(df), "pathways (p.adjust < 0.05:", sig_05, "; <0.10:", sig_10, ")\n")
  return(df)
}

template_kegg <- run_kegg(template_entrez, "Template downstream")
nontemplate_kegg <- run_kegg(nontemplate_entrez, "Non-template upstream")

# ============================================================================
# Step 7: Save enrichment results to CSV
# ============================================================================

if (!is.null(template_kegg)) {
  write.csv(template_kegg, file.path(output_dir, "template_downstream_kegg_full.csv"), row.names = FALSE)
}
if (!is.null(nontemplate_kegg)) {
  write.csv(nontemplate_kegg, file.path(output_dir, "nontemplate_upstream_kegg_full.csv"), row.names = FALSE)
}

# ============================================================================
# Step 8: Generate bar charts (two threshold versions, adaptive height) and save PDF (editable) and TIFF
# ============================================================================

thresholds <- c(0.05, 0.10)

# Template downstream
if (!is.null(template_kegg)) {
  for (thresh in thresholds) {
    res <- create_bar_plot(
      template_kegg,
      title = paste0("Template Strand Downstream Genes\n(p.adjust < ", thresh, ")"),
      color_palette = nature_colors["blue"],
      top_n = 10,
      p_cutoff = thresh
    )
    if (!is.null(res)) {
      file_suffix <- ifelse(thresh == 0.05, "main", "supp")

      # PDF (cairo_pdf ensures editable text)
      ggsave(
        filename = file.path(output_dir, paste0("template_downstream_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,                    # Width increased to 10 inches
        height = res$height + 0.5,
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )

      # TIFF (600 dpi, LZW compression, for printing)
      ggsave(
        filename = file.path(output_dir, paste0("template_downstream_bar_", file_suffix, ".tiff")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = "tiff",
        dpi = 600,
        compression = "lzw",
        limitsize = FALSE
      )

      cat("Template downstream p <", thresh, "bar chart saved (height", round(res$height, 1), "inches)\n")
    }
  }
}

# Non-template upstream
if (!is.null(nontemplate_kegg)) {
  for (thresh in thresholds) {
    res <- create_bar_plot(
      nontemplate_kegg,
      title = paste0("Non-template Strand Upstream Genes\n(p.adjust < ", thresh, ")"),
      color_palette = nature_colors["orange"],
      top_n = 10,
      p_cutoff = thresh
    )
    if (!is.null(res)) {
      file_suffix <- ifelse(thresh == 0.05, "main", "supp")

      ggsave(
        filename = file.path(output_dir, paste0("nontemplate_upstream_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )

      ggsave(
        filename = file.path(output_dir, paste0("nontemplate_upstream_bar_", file_suffix, ".tiff")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = "tiff",
        dpi = 600,
        compression = "lzw",
        limitsize = FALSE
      )

      cat("Non-template upstream p <", thresh, "bar chart saved (height", round(res$height, 1), "inches)\n")
    }
  }
}

# ============================================================================
# Step 9: Side-by-side comparison plot (main plot p < 0.05 only, adaptive height, combined using patchwork)
# ============================================================================

if (!is.null(template_kegg) && !is.null(nontemplate_kegg)) {
  res_temp <- create_bar_plot(template_kegg, p_cutoff = 0.05, top_n = 10, color_palette = nature_colors["blue"])
  res_nontemp <- create_bar_plot(nontemplate_kegg, p_cutoff = 0.05, top_n = 10, color_palette = nature_colors["orange"])

  if (!is.null(res_temp) && !is.null(res_nontemp)) {
    # Remove individual titles, add unified top title
    p_temp <- res_temp$plot + labs(title = NULL) + theme(plot.margin = margin(15, 15, 15, 20))
    p_nontemp <- res_nontemp$plot + labs(title = NULL) + theme(plot.margin = margin(15, 20, 15, 15))

    combined_height <- max(res_temp$height, res_nontemp$height) + 1.0  # Reserve space for title

    # Combine using patchwork
    combined_plot <- (p_temp | p_nontemp) +
      plot_annotation(
        title = "KEGG Pathway Enrichment (p.adjust < 0.05)",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
      ) &
      theme(legend.position = "bottom")

    # Increase right margin to prevent legend from being cropped
    combined_plot <- combined_plot + theme(plot.margin = margin(r = 40))

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

    cat("Side-by-side comparison plot (main) saved, height", round(combined_height, 1), "inches\n")
  }
}

# ============================================================================
cat("\nAll tasks complete! Results saved in:", normalizePath(output_dir), "\n")
