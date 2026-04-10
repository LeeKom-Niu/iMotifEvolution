# ============================================================================
# 四组 KEGG 富集分析及特定通路对比条形图脚本（OUP风格适配版）
# 输入：四个基因列表文件（hominid_genes.txt, homininae_genes.txt, hominini_genes.txt, humanSpecific_genes.txt）
# 输出：各组气泡图（p<0.05 和 p<0.1）、完整结果 CSV、以及 Type II diabetes mellitus 和 Insulin resistance 通路的对比条形图
# 符合 OUP 插图指南（字体放大，色盲友好配色，TIFF无压缩，避免浅色和透明度）
# ============================================================================

# 清理工作空间
rm(list = ls())
gc()

# 设置工作目录（请根据实际情况修改）
setwd("D:/R/data/kegg_human")
cat("工作目录设置为:", getwd(), "\n")

# 输入文件目录（基因列表存放位置）
input_dir <- file.path(getwd(), "input")
if (!dir.exists(input_dir)) {
  stop("请确保 input 目录存在，并将四个基因列表文件放入该目录")
}
cat("基因列表目录:", input_dir, "\n")

# 创建输出目录
output_dir <- "kegg_results_four_groups"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录:", output_dir, "\n")
}

# ============================================================================
# 第一步：加载必要的R包
# ============================================================================
cat("\n=== 第一步：加载必要的R包 ===\n")

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
  cat(pkg, "已加载\n")
}

# ============================================================================
# 第二步：定义配色和主题（符合OUP要求）
# ============================================================================

# 四组配色（色盲友好，高对比度）
group_colors <- c(
  "Hominid" = "#1F77B4",       # 蓝
  "Homininae" = "#FF7F0E",     # 橙
  "Hominini" = "#2CA02C",      # 绿
  "HumanSpecific" = "#D62728"  # 红
)

# OUP主题函数：字体Arial，线条粗细0.5pt，文字不小于7pt
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

# 气泡图函数（修改颜色梯度避免浅色）
create_bubble_plot <- function(enrichment_df, 
                               title = NULL,
                               color_palette = "#1F77B4",
                               top_n = 15,
                               p_cutoff = 0.05) {
  if (is.null(enrichment_df) || nrow(enrichment_df) == 0) return(NULL)
  sig_df <- enrichment_df %>% filter(p.adjust < p_cutoff) %>% arrange(p.adjust)
  if (nrow(sig_df) == 0) {
    cat("  警告：p.adjust <", p_cutoff, "无显著通路\n")
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
    # 避免纯白色，使用浅灰色作为最低值
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
# 第三步：读取基因列表
# ============================================================================
cat("\n=== 第三步：读取基因列表 ===\n")

read_gene_list <- function(filename) {
  file_path <- file.path(input_dir, filename)
  if (!file.exists(file_path)) {
    cat("错误：找不到基因文件", file_path, "\n")
    return(NULL)
  }
  genes <- readLines(file_path)
  genes <- genes[genes != ""]
  genes <- trimws(genes)
  genes <- gsub("\\.[0-9]+$", "", genes)  # 去除版本号
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
    cat(g, "基因数:", length(genes), "\n")
  } else {
    cat(g, "基因列表缺失，跳过\n")
  }
}

if (length(gene_lists) == 0) stop("没有成功读取任何基因列表，请检查 input 目录下的文件")

# ============================================================================
# 第四步：基因ID转换
# ============================================================================
cat("\n=== 第四步：基因ID转换 ===\n")

convert_genes <- function(genes, list_name) {
  if (length(genes) == 0) return(NULL)
  gene_df <- tryCatch({
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
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

entrez_lists <- list()
for (g in names(gene_lists)) {
  entrez <- convert_genes(gene_lists[[g]], g)
  if (!is.null(entrez)) entrez_lists[[g]] <- entrez
}

if (length(entrez_lists) == 0) stop("没有成功转换任何基因")

# ============================================================================
# 第五步：KEGG富集分析
# ============================================================================
cat("\n=== 第五步：KEGG富集分析 ===\n")

run_kegg_analysis <- function(entrez_ids, list_name) {
  if (length(entrez_ids) < 5) {
    cat(list_name, ": 基因数不足5，跳过KEGG分析\n")
    return(NULL)
  }
  cat(list_name, ": 运行KEGG富集分析...\n")
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

kegg_results <- list()
for (g in names(entrez_lists)) {
  kegg_results[[g]] <- run_kegg_analysis(entrez_lists[[g]], g)
  if (!is.null(kegg_results[[g]])) {
    write.csv(kegg_results[[g]], file.path(output_dir, paste0(g, "_kegg_full.csv")), row.names = FALSE)
  }
}

# ============================================================================
# 第六步：生成气泡图（p<0.05 和 p<0.1）并保存为PDF和TIFF（无压缩）
# ============================================================================
cat("\n=== 第六步：生成气泡图 ===\n")

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
      # PDF使用cairo_pdf嵌入字体
      ggsave(pdf_file, res$plot, width = 10, height = res$height + 0.5, device = cairo_pdf, dpi = 300, limitsize = FALSE)
      # TIFF无压缩，600dpi
      ggsave(tiff_file, res$plot, width = 10, height = res$height + 0.5, device = "tiff", dpi = 600, compression = "none", limitsize = FALSE)
      cat(g, "p <", thresh, "气泡图已保存\n")
    }
  }
}

# ============================================================================
# 第七步：特定通路对比条形图（完整版）
# 包括：构建对比数据、重命名分组、颜色映射、添加 -log10(p.adjust)、
#       行排列、列排列、横向条形图三种布局
# ============================================================================
# ============================================================================
# ============================================================================
# 第七步：特定通路对比条形图（最终版）
# 功能：构建对比数据、绘制行排列/列排列/横向条形图
# 特点：仅显示显著性星号、柱子宽度收窄、纵向细长、可调节图注间隙与页边距
# ============================================================================
cat("\n=== 第七步：特定通路对比条形图 ===\n")

# ----- 0. 加载所需包 -----
library(ggplot2)
library(reshape2)   # 用于数据重塑

# ----- 1. 定义目标通路（ID 和描述）-----
target_ids <- c("hsa04930", "hsa04931")
target_names <- c("hsa04930" = "Type II diabetes mellitus", 
                  "hsa04931" = "Insulin resistance")

# ----- 2. 构建对比数据框（从 kegg_results 提取）-----
comparison_data <- data.frame()

for (g in names(kegg_results)) {
  df <- kegg_results[[g]]
  if (is.null(df)) next
  
  cat("\n处理组别:", g, "（共有", nrow(df), "个通路）\n")
  
  for (id in target_ids) {
    pw_row <- df[df$ID == id, ]
    if (nrow(pw_row) > 0) {
      # 确保 FoldEnrichment 列存在
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
      
      cat("  找到通路:", target_names[id], 
          "富集倍数 =", fold, "p.adjust =", padj, "基因数 =", cnt, "\n")
      
      comparison_data <- rbind(comparison_data, data.frame(
        Group = g,
        Pathway = target_names[id],
        FoldEnrichment = fold,
        p_adjust = padj,
        GeneCount = cnt,
        stringsAsFactors = FALSE
      ))
    } else {
      cat("  未找到通路:", target_names[id], "\n")
    }
  }
}

# ----- 3. 数据清洗与分组重命名 -----
comparison_data$Group <- as.character(comparison_data$Group)
comparison_data <- comparison_data[!is.na(comparison_data$FoldEnrichment), ]
comparison_data$Group[comparison_data$Group == "Hominid"] <- "Great ape"

# 自定义分组顺序（移除无数据的 Hominini）
group_order <- c("Great ape", "Homininae", "HumanSpecific")
comparison_data$Group <- factor(comparison_data$Group, levels = group_order)

# 通路因子顺序
comparison_data$Pathway <- factor(comparison_data$Pathway, 
                                  levels = c("Type II diabetes mellitus", 
                                             "Insulin resistance"))

# ----- 4. 定义颜色 -----
group_colors <- c(
  "Great ape"   = "#1f77b4",   # 蓝色
  "Homininae"   = "#9467bd",   # 紫色
  "HumanSpecific" = "#ffd966"  # 黄色
)

# ----- 5. 添加显著性星号 -----
comparison_data$sig <- ifelse(comparison_data$p_adjust < 0.001, "***",
                              ifelse(comparison_data$p_adjust < 0.01, "**",
                                     ifelse(comparison_data$p_adjust < 0.05, "*", "ns")))

# ----- 6. 自定义主题函数（可调节图注间隙与页边距）-----
# 参数说明：
#   base_size: 基础字体大小
#   legend_spacing_x: 图例条目之间的水平间距（单位："lines"，数值越大间隔越大）
#   legend_margin: 图例框的外边距，格式 margin(t, r, b, l)（单位：pt）
#   plot_margin: 整个图形的页边距，格式 margin(t, r, b, l)
theme_oup <- function(base_size = 12,
                      legend_spacing_x = 5.5,      # 默认0.5行高度
                      legend_margin = margin(20, 20, 20, 20),  # 上下左右各20pt
                      plot_margin = margin(15, 20, 15, 35)) {  # 上右下左
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.3),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      # 图注条目水平间距
      legend.spacing.x = unit(legend_spacing_x, "lines"),
      # 图注框外边距
      legend.margin = legend_margin,
      # 整个图形页边距
      plot.margin = plot_margin
    )
}

# ========================
# 版本 A：行排列（左右并排）
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
# 版本 B：列排列（上下堆叠）—— 纵向细长
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
# 版本 C：横向条形图（coord_flip）+ 纵向排列
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

# ----- 7. 保存图形（指定输出目录，请根据实际路径修改）-----
output_dir <- getwd()  # 可改为具体路径，如 "D:/R/data/kegg_human/kegg_results_four_groups"

# 行排列
ggsave(file.path(output_dir, "pathway_comparison_barplot_row.pdf"), 
       p_bar_row, width = 10, height = 6, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_row.tiff"), 
       p_bar_row, width = 10, height = 6, device = "tiff", dpi = 600, compression = "none")

# 列排列（纵向细长）
ggsave(file.path(output_dir, "pathway_comparison_barplot_col.pdf"), 
       p_bar_col, width = 6, height = 10, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_col.tiff"), 
       p_bar_col, width = 6, height = 10, device = "tiff", dpi = 600, compression = "none")

# 横向条形图（纵向细长）
ggsave(file.path(output_dir, "pathway_comparison_barplot_horiz.pdf"), 
       p_bar_horiz, width = 5, height = 8, device = cairo_pdf, dpi = 300)
ggsave(file.path(output_dir, "pathway_comparison_barplot_horiz.tiff"), 
       p_bar_horiz, width = 5, height = 8, device = "tiff", dpi = 600, compression = "none")

cat("\n特定通路对比条形图（三种布局）已保存至:", output_dir, "\n")

# ----- 8. 保存富集倍数矩阵（可选）-----
if (nrow(comparison_data) > 0) {
  mat_wide <- reshape(comparison_data[, c("Group", "Pathway", "FoldEnrichment")],
                      idvar = "Group", timevar = "Pathway", direction = "wide")
  colnames(mat_wide) <- gsub("FoldEnrichment\\.", "", colnames(mat_wide))
  write.csv(mat_wide, file.path(output_dir, "pathway_fold_enrichment_matrix.csv"), 
            row.names = FALSE)
  cat("富集倍数矩阵已保存\n")
}

cat("\n=== 第七步完成 ===\n")

cat("\n=== 第七步完成 ===\n")
# ============================================================================
# 完成
# ============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("四组 KEGG 富集分析及对比条形图生成完成！\n")
cat("结果保存在:", normalizePath(output_dir), "\n")
cat("注：TIFF图片为无压缩格式，符合OUP要求；条形图已移除透明度以确保清晰度。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")