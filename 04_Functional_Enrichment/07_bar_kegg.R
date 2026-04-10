# ============================================================================
# 条形图版：KEGG富集分析条形图生成脚本（OUP插图指南合规版，字体放大版）
# 改进点：按p值取top15显著通路后，按GeneRatio排序显示，最长条形在顶部
# 用于生成hominid和humanSpecific基因列表的p<0.05主图和p<0.1附图
# ============================================================================

# 清理工作空间
rm(list = ls())
gc()

# 设置工作目录（请根据实际情况修改）
setwd("D:/R/data/kegg_human")
cat("工作目录设置为:", getwd(), "\n")

# ============================================================================
# 第一步：加载必要的R包
# ============================================================================

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

# ============================================================================
# 第二步：创建输出目录
# ============================================================================

cat("\n=== 第二步：创建输出目录 ===\n")

output_dir <- "bar_plots_dual_oup_largefont"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建目录:", output_dir, "\n")
}

# ============================================================================
# 第三步：定义OUP风格主题（字体放大版）
# ============================================================================

# Nature配色方案（色盲友好）
nature_colors <- c(
  blue = "#1F77B4",
  orange = "#FF7F0E", 
  green = "#2CA02C",
  red = "#D62728"
)

# OUP主题函数：基础字体12pt，线条粗细0.3-0.5pt
theme_oup <- function(base_size = 12, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 轴线
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # 刻度标签
      axis.text = element_text(size = rel(1.2), color = "black"),  # 约14.4pt
      axis.title = element_text(size = rel(1.5), face = "plain"),  # 约18pt
      # 图例
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = rel(1.2), face = "plain"), # 约14.4pt
      legend.text = element_text(size = rel(1)),                     # 12pt
      legend.key.size = unit(0.6, "cm"),
      # 网格线
      panel.grid.major = element_line(linewidth = 0.3, color = "#CCCCCC", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      # 面板边框
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # 标题
      plot.title = element_text(size = rel(1.8), face = "bold", hjust = 0.5, margin = margin(b = 15)), # 约21.6pt
      # 图形边距
      plot.margin = margin(15, 20, 15, 15),
      text = element_text(family = base_family)
    )
}

# 增强的条形图函数（支持自定义p阈值和OUP主题）
# 改进点：先按p值取前top_n个最显著通路，然后按GeneRatio升序排列，使最长条形位于顶部
create_bar_plot <- function(enrichment_df, 
                            title = NULL,
                            color_palette = nature_colors["blue"],
                            top_n = 15,
                            p_cutoff = 0.05) {
  
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
    return(NULL)
  }
  
  # 筛选显著通路（p.adjust < p_cutoff）
  sig_df <- enrichment_df %>% 
    filter(p.adjust < p_cutoff) %>%
    arrange(p.adjust)   # 按p值升序（最显著在前）
  
  if (nrow(sig_df) == 0) {
    cat("  警告：p.adjust <", p_cutoff, "无显著通路\n")
    return(NULL)
  }
  
  cat("  显著通路数 (p<", p_cutoff, "):", nrow(sig_df), "\n")
  
  # 取前top_n个通路（最显著的通路）
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
  
  # 关键修改：按GeneRatio_num升序排序，使条形长度从上到下递增
  plot_df <- plot_df %>%
    arrange(GeneRatio_num) %>%
    mutate(Description_short = factor(Description_short, levels = Description_short))
  
  # 动态计算图形高度（每行0.45英寸，加上标题和图例空间）
  plot_height <- max(5, nrow(plot_df) * 0.45 + 3)
  
  # 创建条形图（应用OUP主题）
  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description_short)) +
    geom_col(aes(fill = log10_padj), 
             width = 0.7, color = "black", linewidth = 0.3) +  # 条形边框0.3pt
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
    theme_oup(base_size = 12) +  # 基础字体12pt
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
# 第四步：读取和处理基因列表（保持不变）
# ============================================================================

cat("\n=== 第四步：读取和处理基因列表 ===\n")

read_gene_list <- function(filename) {
  if (!file.exists(filename)) {
    cat("错误：找不到基因文件", filename, "\n")
    return(NULL)
  }
  
  genes <- readLines(filename)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)  # 去除版本号
  
  return(genes)
}

# 读取基因列表
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

# ============================================================================
# 第五步：基因ID转换（保持不变）
# ============================================================================

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

# ============================================================================
# 第六步：KEGG富集分析（保持不变）
# ============================================================================

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
      pvalueCutoff = 0.2,        # 放宽阈值以捕获更多通路
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

# 保存完整结果到CSV
if (!is.null(hominid_kegg)) {
  write.csv(hominid_kegg, file.path(output_dir, "hominid_kegg_full.csv"), row.names = FALSE)
}
if (!is.null(humanSpecific_kegg)) {
  write.csv(humanSpecific_kegg, file.path(output_dir, "humanSpecific_kegg_full.csv"), row.names = FALSE)
}

# ============================================================================
# 第七步：生成条形图（p<0.05主图和p<0.1附图）并保存为PDF（可编辑）和TIFF
# ============================================================================

cat("\n=== 第七步：生成条形图 ===\n")

thresholds <- c(0.05, 0.10)

# 全猿共享条形图
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
      
      # 保存PDF（cairo_pdf确保文字可编辑）
      ggsave(
        filename = file.path(output_dir, paste0("hominid_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,                # 宽度增至10英寸
        height = res$height + 0.5, # 高度稍增
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )
      
      # 保存TIFF（600dpi，LZW压缩，用于印刷）
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

# 人类特有条形图
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

# ============================================================================
# 第八步：生成并排对比图（仅p<0.05主图）
# ============================================================================

if (!is.null(hominid_kegg) && !is.null(humanSpecific_kegg)) {
  res_hom <- create_bar_plot(hominid_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["blue"])
  res_hum <- create_bar_plot(humanSpecific_kegg, p_cutoff = 0.05, top_n = 15, color_palette = nature_colors["orange"])
  
  if (!is.null(res_hom) && !is.null(res_hum)) {
    # 移除各自标题，添加统一顶部标题
    p_hom <- res_hom$plot + labs(title = NULL) + theme(plot.margin = margin(15, 15, 15, 20))
    p_hum <- res_hum$plot + labs(title = NULL) + theme(plot.margin = margin(15, 20, 15, 15))
    
    combined_height <- max(res_hom$height, res_hum$height) + 1.0  # 为标题留空间
    
    # 使用patchwork组合
    combined_plot <- (p_hom | p_hum) +
      plot_annotation(
        title = "KEGG Pathway Enrichment (p.adjust < 0.05)",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
      ) &
      theme(legend.position = "bottom")
    combined_plot <- combined_plot + theme(plot.margin = margin(r = 100))
    
    # 保存PDF
    ggsave(
      filename = file.path(output_dir, "combined_bar_main.pdf"),
      plot = combined_plot,
      width = 18,
      height = combined_height,
      device = cairo_pdf,
      dpi = 300,
      limitsize = FALSE
    )
    
    # 保存TIFF
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

# ============================================================================
# 完成！
# ============================================================================

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("条形图生成完成（OUP插图指南合规版，字体放大）！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\n所有结果保存在:", normalizePath(output_dir), "\n")