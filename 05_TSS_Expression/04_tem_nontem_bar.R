# ============================================================================
# KEGG富集分析条形图生成脚本（双阈值版，字体放大版，OUP插图指南合规）
# 模板链下游 vs 非模板链上游基因
# 功能：生成 p.adjust < 0.05（正图）和 p.adjust < 0.1（附图）的条形图
#       自动输出每个阈值下的通路数量，图例完整，高度自适应
#       所有文字显著放大，适合在 AI 中缩放编辑
#       条形按基因比率升序排列（最长条形在顶部）
# ============================================================================

# 清理工作空间
rm(list = ls())
gc()

# 设置工作目录（请根据实际情况修改）
setwd("D:/R/data/kegg_human/tem_nontem")
cat("工作目录设置为:", getwd(), "\n")

# ============================================================================
# 第一步：加载必要的R包
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
# 第二步：创建输出目录
# ============================================================================

output_dir <- "bar_plots_dual_threshold_v3_oup_largefont"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("输出目录:", output_dir, "\n")

# ============================================================================
# 第三步：定义 OUP 风格主题（字体放大版）
# ============================================================================

# 色盲友好配色
nature_colors <- c(
  blue = "#1F77B4",
  orange = "#FF7F0E"
)

# OUP 主题函数（基础字体 14pt，按比例放大）
theme_oup <- function(base_size = 14, base_family = "Arial") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # 轴线
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # 刻度标签
      axis.text = element_text(size = rel(1), color = "black"),      # 14pt
      axis.title = element_text(size = rel(1.2), face = "plain"),    # 约 16.8pt
      # 图例
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = rel(1.1), face = "plain"),  # 约 15.4pt
      legend.text = element_text(size = rel(1)),                      # 14pt
      legend.key.size = unit(0.7, "cm"),
      # 网格线
      panel.grid.major = element_line(linewidth = 0.3, color = "#CCCCCC", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      # 面板边框
      panel.border = element_rect(linewidth = 0.5, fill = NA),
      # 标题
      plot.title = element_text(size = rel(1.4), face = "bold", hjust = 0.5, margin = margin(b = 15)), # 约 19.6pt
      # 图形边距
      plot.margin = margin(20, 60, 20, 60),
      text = element_text(family = base_family)
    )
}

# 条形图绘制函数（字体放大版，按GeneRatio升序排列）
create_bar_plot <- function(enrichment_df, 
                            title = "KEGG Pathway Enrichment",
                            color_palette = nature_colors["blue"],
                            top_n = 10,
                            p_cutoff = 0.05) {
  
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)
  
  # 筛选指定阈值下的显著通路
  sig_df <- enrichment_df %>% 
    filter(p.adjust < p_cutoff) %>%
    arrange(p.adjust)
  
  if (nrow(sig_df) == 0) {
    cat("警告：p.adjust <", p_cutoff, "无显著通路\n")
    return(NULL)
  }
  
  cat("  p.adjust <", p_cutoff, "：找到", nrow(sig_df), "个显著通路\n")
  
  # 取前 top_n 个通路（最显著的）
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
  
  # 关键修改：按GeneRatio_num升序排序，使条形长度从上到下递增（最长条形在顶部）
  plot_df <- plot_df %>%
    arrange(GeneRatio_num) %>%
    mutate(Description_short = factor(Description_short, levels = Description_short))
  
  # 动态计算图形高度（字体变大，每行高度增加至 0.45 英寸）
  plot_height <- max(5, nrow(plot_df) * 0.45 + 2.5)
  
  # 条形图（使用geom_col）
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
    theme_oup(base_size = 14) +  # 基础字体 14pt
    theme(
      axis.text.y = element_text(size = 14),          # 14pt
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))  # 右侧留出少量空间
  
  return(list(plot = p, height = plot_height))
}

# ============================================================================
# 第四步：读取基因列表
# ============================================================================

read_gene_list <- function(filename) {
  if (!file.exists(filename)) stop("文件不存在: ", filename)
  genes <- readLines(filename)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)
  return(genes)
}

template_genes <- read_gene_list("template_downstream_all_genes.txt")
nontemplate_genes <- read_gene_list("nontemplate_upstream_all_genes.txt")

cat("模板链下游基因数:", length(template_genes), "\n")
cat("非模板链上游基因数:", length(nontemplate_genes), "\n")

# ============================================================================
# 第五步：基因ID转换
# ============================================================================

convert_to_entrez <- function(genes, list_name) {
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(list_name, "成功转换", nrow(gene_df), "/", length(genes), "个基因\n")
  return(gene_df$ENTREZID)
}

template_entrez <- convert_to_entrez(template_genes, "模板链下游")
nontemplate_entrez <- convert_to_entrez(nontemplate_genes, "非模板链上游")

# ============================================================================
# 第六步：KEGG富集分析（使用宽松阈值捕获所有可能通路）
# ============================================================================

run_kegg <- function(entrez_ids, list_name) {
  if (length(entrez_ids) < 5) {
    cat(list_name, "基因数不足5，跳过KEGG分析\n")
    return(NULL)
  }
  
  kegg <- enrichKEGG(
    gene = unique(entrez_ids),
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.2,        # 放宽阈值，捕获更多通路
    pAdjustMethod = "BH",
    qvalueCutoff = 0.25,
    minGSSize = 5,
    maxGSSize = 500,
    use_internal_data = FALSE
  )
  
  if (is.null(kegg) || nrow(kegg) == 0) {
    cat(list_name, "未富集到任何通路\n")
    return(NULL)
  }
  
  df <- as.data.frame(kegg)
  sig_05 <- sum(df$p.adjust < 0.05)
  sig_10 <- sum(df$p.adjust < 0.10)
  cat(list_name, "富集到", nrow(df), "个通路，其中 p.adjust < 0.05:", sig_05, "个；<0.10:", sig_10, "个\n")
  return(df)
}

template_kegg <- run_kegg(template_entrez, "模板链下游")
nontemplate_kegg <- run_kegg(nontemplate_entrez, "非模板链上游")

# ============================================================================
# 第七步：保存富集结果到CSV
# ============================================================================

if (!is.null(template_kegg)) {
  write.csv(template_kegg, file.path(output_dir, "template_downstream_kegg_full.csv"), row.names = FALSE)
}
if (!is.null(nontemplate_kegg)) {
  write.csv(nontemplate_kegg, file.path(output_dir, "nontemplate_upstream_kegg_full.csv"), row.names = FALSE)
}

# ============================================================================
# 第八步：生成条形图（两个阈值版本，使用自适应高度）并保存 PDF（可编辑）和 TIFF
# ============================================================================

thresholds <- c(0.05, 0.10)

# 模板链下游
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
      
      # PDF (cairo_pdf 确保文字可编辑)
      ggsave(
        filename = file.path(output_dir, paste0("template_downstream_bar_", file_suffix, ".pdf")),
        plot = res$plot,
        width = 10,                    # 宽度增至10英寸
        height = res$height + 0.5,
        device = cairo_pdf,
        dpi = 300,
        limitsize = FALSE
      )
      
      # TIFF (600 dpi, LZW 压缩，用于印刷)
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
      
      cat("模板链下游 p <", thresh, "条形图已保存（高度", round(res$height, 1), "英寸）\n")
    }
  }
}

# 非模板链上游
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
      
      cat("非模板链上游 p <", thresh, "条形图已保存（高度", round(res$height, 1), "英寸）\n")
    }
  }
}

# ============================================================================
# 第九步：并排对比图（仅对主图 p < 0.05，自适应高度，使用 patchwork 组合）
# ============================================================================

if (!is.null(template_kegg) && !is.null(nontemplate_kegg)) {
  res_temp <- create_bar_plot(template_kegg, p_cutoff = 0.05, top_n = 10, color_palette = nature_colors["blue"])
  res_nontemp <- create_bar_plot(nontemplate_kegg, p_cutoff = 0.05, top_n = 10, color_palette = nature_colors["orange"])
  
  if (!is.null(res_temp) && !is.null(res_nontemp)) {
    # 移除各自标题，添加统一顶部标题
    p_temp <- res_temp$plot + labs(title = NULL) + theme(plot.margin = margin(15, 15, 15, 20))
    p_nontemp <- res_nontemp$plot + labs(title = NULL) + theme(plot.margin = margin(15, 20, 15, 15))
    
    combined_height <- max(res_temp$height, res_nontemp$height) + 1.0  # 为标题留空间
    
    # 使用 patchwork 组合
    combined_plot <- (p_temp | p_nontemp) +
      plot_annotation(
        title = "KEGG Pathway Enrichment (p.adjust < 0.05)",
        theme = theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20))
        )
      ) &
      theme(legend.position = "bottom")
    
    # 增加右边距防止图例被裁剪
    combined_plot <- combined_plot + theme(plot.margin = margin(r = 40))
    
    # 保存 PDF
    ggsave(
      filename = file.path(output_dir, "combined_bar_main.pdf"),
      plot = combined_plot,
      width = 18,
      height = combined_height,
      device = cairo_pdf,
      dpi = 300,
      limitsize = FALSE
    )
    
    # 保存 TIFF
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
    
    cat("并排对比图（主图）已保存，高度", round(combined_height, 1), "英寸\n")
  }
}

# ============================================================================
cat("\n所有任务完成！结果保存在:", normalizePath(output_dir), "\n")