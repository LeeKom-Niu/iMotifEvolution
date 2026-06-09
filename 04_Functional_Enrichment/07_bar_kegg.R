
rm(list = ls())
gc()

setwd("D:/R/data/kegg_human")
cat("工作目录设置为:", getwd(), "\n")


cat("\n=== 第一步：加载必要的R包 ===\n")

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
  cat(pkg, "已加载\n")
}


cat("\n=== 第二步：创建输出目录 ===\n")

output_dir <- "bar_plots_dual_oup_largefont"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建目录:", output_dir, "\n")
}


nature_colors <- c(
  blue = "
  orange = "
  green = "
  red = "
)

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
      panel.grid.major = element_line(linewidth = 0.3, color = "
      panel.grid.minor = element_blank(),
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      plot.title = element_text(size = rel(1.8), face = "bold", hjust = 0.5, margin = margin(b = 15)),
      plot.margin = margin(15, 20, 15, 15),
      text = element_text(family = base_family)
    )
}

create_bar_plot <- function(enrichment_df, 
                            title = NULL,
                            color_palette = nature_colors["blue"],
                            top_n = 15,
                            p_cutoff = 0.05) {
  
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }
  
  sig_df <- enrichment_df %>% 
    filter(p.adjust < p_cutoff) %>%
    arrange(p.adjust)
  
  if (nrow(sig_df) == 0) {
    cat("  警告：p.adjust <", p_cutoff, "无显著通路\n")
    return(NULL)
  }
  
  cat("  显著通路数 (p<", p_cutoff, "):", nrow(sig_df), "\n")
  
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
  
  plot_df <- plot_df %>%
    arrange(GeneRatio_num) %>%
    mutate(Description_short = factor(Description_short, levels = Description_short))
  
  plot_height <- max(5, nrow(plot_df) * 0.45 + 3)
  
  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_short)) +
    geom_col(aes(fill = log10_padj), 
             width = 0.7, color = "black", linewidth = 0.3) +
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
    theme_oup(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))
  
  return(list(plot = p, height = plot_height))
}


cat("\n=== 第四步：读取和处理基因列表 ===\n")

read_gene_list <- function(filename) {
  if (!file.exists(filename)) {
    cat("错误：找不到基因文件", filename, "\n")
    return(NULL)
  }
  
  genes <- readLines(filename)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)
  
  return(genes)
}

hominid_genes <- read_gene_list("hominid_genes.txt")
humanSpecific_genes <- read_gene_list("humanSpecific_genes.txt")

if (is.null(hominid_genes) && is.null(humanSpecific_genes)) {
  stop("错误：两个基因列表文件都未找到")
}

cat("基因列表统计:\n")
if (!is.null(hominid_genes)) {
  cat("  全猿共享基因:", length(hominid_genes), "个\n")
}
if (!is.null(humanSpecific_genes)) {
  cat("  人类特有基因:", length(humanSpecific_genes), "个\n")
}


cat("\n=== 第五步：基因ID转换 ===\n")

convert_genes <- function(genes, list_name) {
  if (is.null(genes) || length(genes) == 0) {
    cat(list_name, ": 没有基因数据\n")
    return(NULL)
  }
  
  gene_df <- tryCatch({
    bitr(genes, 
         fromType = "SYMBOL",
         toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    cat(list_name, "基因转换错误:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(gene_df) && nrow(gene_df) > 0) {
    cat(list_name, ": 成功转换", nrow(gene_df), "/", length(genes), "个基因\n")
    return(gene_df$ENTREZID)
  } else {
    cat(list_name, ": 基因转换失败\n")
    return(NULL)
  }
}

hominid_entrez <- convert_genes(hominid_genes, "全猿共享")
humanSpecific_entrez <- convert_genes(humanSpecific_genes, "人类特有")


cat("\n=== 第六步：KEGG富集分析 ===\n")

run_kegg_analysis <- function(entrez_ids, list_name) {
  if (is.null(entrez_ids) || length(entrez_ids) < 5) {
    cat(list_name, ": 基因数不足，跳过KEGG分析\n")
    return(NULL)
  }
  
  cat(list_name, ": 运行KEGG富集分析...\n")
  
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = unique(entrez_ids),
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.2,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.25,
      minGSSize = 5,
      maxGSSize = 500,
      use_internal_data = FALSE
    )
  }, error = function(e) {
    cat(list_name, "KEGG分析错误:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(kegg_result) || nrow(kegg_result) == 0) {
    cat(list_name, ": 没有富集通路\n")
    return(NULL)
  }
  
  df <- as.data.frame(kegg_result)
  cat(list_name, ": 富集到", nrow(df), "个通路\n")
  return(df)
}

hominid_kegg <- run_kegg_analysis(hominid_entrez, "全猿共享")
humanSpecific_kegg <- run_kegg_analysis(humanSpecific_entrez, "人类特有")

if (!is.null(hominid_kegg)) {
  write.csv(hominid_kegg, file.path(output_dir, "hominid_kegg_full.csv"), row.names = FALSE)
}
if (!is.null(humanSpecific_kegg)) {
  write.csv(humanSpecific_kegg, file.path(output_dir, "humanSpecific_kegg_full.csv"), row.names = FALSE)
}


cat("\n=== 第七步：生成条形图 ===\n")

thresholds <- c(0.05, 0.10)

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
      
      ggsave(
        filename = file.path(output_dir, paste0("hominid_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,
        height = res$height + 0.5,
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )
      
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
      
      cat("全猿共享 p <", thresh, "条形图已保存（PDF+TIFF）\n")
    }
  }
}

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
      
      cat("人类特有 p <", thresh, "条形图已保存（PDF+TIFF）\n")
    }
  }
}


if (!is.null(hominid_kegg) && !is.null(humanSpecific_kegg)) {
  res_hom <- create_bar_plot(hominid_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["blue"])
  res_hum <- create_bar_plot(humanSpecific_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["orange"])
  
  if (!is.null(res_hom) && !is.null(res_hum)) {
    p_hom <- res_hom$plot + labs(title = NULL) + theme(plot.margin = margin(15, 15, 15, 20))
    p_hum <- res_hum$plot + labs(title = NULL) + theme(plot.margin = margin(15, 20, 15, 15))
    
    combined_height <- max(res_hom$height, res_hum$height) + 1.0
    
    combined_plot <- (p_hom | p_hum) +
      plot_annotation(
        title = "KEGG Pathway Enrichment (p.adjust < 0.05)",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
      ) &
      theme(legend.position = "bottom")
    combined_plot <- combined_plot + theme(plot.margin = margin(r = 100))
    
    ggsave(
      filename = file.path(output_dir, "combined_bar_main.pdf"),
      plot = combined_plot,
      width = 18,
      height = combined_height,
      device = cairo_pdf,
      dpi = 300,
      limitsize = FALSE
    )
    
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
    
    cat("并排对比图（主图）已保存（PDF+TIFF）\n")
  } else {
    cat("并排图跳过：至少一个数据集无p<0.05显著通路\n")
  }
}


cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("条形图生成完成（OUP插图指南合规版，字体放大）！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n所有结果保存在:", normalizePath(output_dir), "\n")