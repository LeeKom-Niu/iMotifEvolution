# ==============================================
# 最终版 · 7张图 · 纵向 + 无网格线 + 组织标题 + 精简年龄标签 + 图例间距加大
# - 无单元格边框（border_color = NA）
# - 图形宽度 10 英寸，高度 18/14 英寸
# - 纵轴仅显示基因名，组织名在图注中说明
# - 单独组织热图顶部显示组织名称
# - x轴年龄标签仅显示 21,30,40,50,60,70
# - 图例与热图主体间距增加 1.5 cm
# ==============================================
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

dir.create("plots_7figs_FINAL", showWarnings = FALSE)

gene_list <- c(
  "AKT1", "CACNA1A", "CPT1B", "CREB3L1", "GFPT1", "HK2",
  "INSR", "IRS2", "MAPK9", "MLX", "MLXIPL", "PDPK1",
  "PKM", "PPARGC1B", "PPP1CA", "PPP1CC", "PRKAA2", "PRKAG1",
  "PRKCD", "PRKCZ", "PTEN", "PTPN11", "PYGB", "RPS6KA2",
  "SOCS2", "STAT3"
)

tissue1 <- "Small Intestine - Terminal Ileum"
tissue2 <- "Colon - Transverse"
age_min <- 21
age_max <- 70
target_ages <- c(21, 30, 40, 50, 60, 70)

read_gene <- function(gene) {
  f <- paste0("aging_expression/", gene, "_acrossTissues.csv")
  df <- read.csv(f, na.strings = c("", NA), stringsAsFactors = F, fileEncoding = "UTF-8-BOM")
  colnames(df) <- c("Age", "Tissue", "Zscore")
  df %>%
    filter(!is.na(Zscore), Age >= age_min, Age <= age_max) %>%
    filter(Tissue %in% c(tissue1, tissue2)) %>%
    mutate(Age = as.integer(Age), 
           RowLabel = paste0(gene, " | ", Tissue))
}

cat("📥 读取 26 个基因表达...\n")
df_list <- lapply(gene_list, read_gene)
data <- bind_rows(df_list)

mat <- data %>%
  select(RowLabel, Age, Zscore) %>%
  pivot_wider(names_from = Age, values_from = Zscore) %>%
  column_to_rownames("RowLabel") %>%
  as.matrix()
mat[mat > 2] <- 2
mat[mat < -2] <- -2

calc_cor <- function(row_vals) {
  ages <- as.integer(colnames(mat))
  keep <- !is.na(row_vals)
  if(sum(keep) < 3) return(NA)
  cor(ages[keep], row_vals[keep], method = "pearson")
}
cor_vals <- apply(mat, 1, calc_cor)

mat_intestine <- mat[grepl(tissue1, rownames(mat)), , drop = FALSE]
mat_colon     <- mat[grepl(tissue2, rownames(mat)), , drop = FALSE]
cor_intestine <- cor_vals[rownames(mat_intestine)]
cor_colon     <- cor_vals[rownames(mat_colon)]

fixed_order <- c()
for (g in gene_list) {
  gi <- paste0(g, " | ", tissue1)
  gc <- paste0(g, " | ", tissue2)
  if (gi %in% rownames(mat)) fixed_order <- c(fixed_order, gi)
  if (gc %in% rownames(mat)) fixed_order <- c(fixed_order, gc)
}
mat_fixed <- mat[fixed_order, , drop = FALSE]
mat_all_sorted <- mat_fixed[order(cor_vals[fixed_order], decreasing = TRUE), , drop = FALSE]
mat_int_sorted <- mat_intestine[order(cor_intestine, decreasing = TRUE), , drop = FALSE]
mat_col_sorted <- mat_colon[order(cor_colon, decreasing = TRUE), , drop = FALSE]

get_heatmap_colors <- function() colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
legend_breaks <- seq(-2, 2, by = 1)

make_age_labels <- function(mat) {
  ages <- as.integer(colnames(mat))
  ifelse(ages %in% target_ages, as.character(ages), "")
}

# ===================== 绘图函数（增加 legend_spacing 参数） =====================
plot_no_cluster_oup <- function(mat, name, height_inch, width_inch = 10, title = NA) {
  mat_display <- mat
  rownames(mat_display) <- gsub(paste0(" \\| ", tissue1, "|", tissue2), "", rownames(mat_display))
  age_labels <- make_age_labels(mat)
  colors <- get_heatmap_colors()
  
  pdf(file = paste0("plots_7figs_FINAL/", name, ".pdf"),
      width = width_inch, height = height_inch,
      family = "Helvetica", useDingbats = FALSE)
  pheatmap(mat_display,
    color = colors, breaks = seq(-2, 2, length.out = 101),
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = TRUE,
    labels_col = age_labels,
    fontsize_row = 9, fontsize_col = 11, fontfamily = "Helvetica",
    treeheight_row = 0, treeheight_col = 0,
    border_color = NA, na_col = "white",
    legend = TRUE, legend_title = "Z-score",
    legend_breaks = legend_breaks, legend_labels = as.character(legend_breaks),
    legend_spacing = 1.5,          # 增加图例与主图间距（单位：cm）
    main = title)
  dev.off()
  
  png(paste0("plots_7figs_FINAL/", name, ".png"),
      width = width_inch, height = height_inch, units = "in", res = 300)
  pheatmap(mat_display,
    color = colors, breaks = seq(-2, 2, length.out = 101),
    cluster_rows = FALSE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = TRUE,
    labels_col = age_labels,
    fontsize_row = 9, fontsize_col = 11, fontfamily = "Helvetica",
    treeheight_row = 0, treeheight_col = 0,
    border_color = NA, na_col = "white",
    legend = TRUE, legend_title = "Z-score",
    legend_breaks = legend_breaks, legend_labels = as.character(legend_breaks),
    legend_spacing = 1.5,
    main = title)
  dev.off()
}

plot_cluster_oup <- function(mat, name, height_inch, width_inch = 10, tree_height = 50, title = NA) {
  mat_display <- mat
  rownames(mat_display) <- gsub(paste0(" \\| ", tissue1, "|", tissue2), "", rownames(mat_display))
  age_labels <- make_age_labels(mat)
  colors <- get_heatmap_colors()
  
  pdf(file = paste0("plots_7figs_FINAL/", name, ".pdf"),
      width = width_inch, height = height_inch,
      family = "Helvetica", useDingbats = FALSE)
  pheatmap(mat_display,
    color = colors, breaks = seq(-2, 2, length.out = 101),
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = TRUE,
    labels_col = age_labels,
    fontsize_row = 9, fontsize_col = 11, fontfamily = "Helvetica",
    treeheight_row = tree_height, treeheight_col = 0,
    border_color = NA, na_col = "white",
    legend = TRUE, legend_title = "Z-score",
    legend_breaks = legend_breaks, legend_labels = as.character(legend_breaks),
    legend_spacing = 1.5,
    main = title)
  dev.off()
  
  png(paste0("plots_7figs_FINAL/", name, ".png"),
      width = width_inch, height = height_inch, units = "in", res = 300)
  pheatmap(mat_display,
    color = colors, breaks = seq(-2, 2, length.out = 101),
    cluster_rows = TRUE, cluster_cols = FALSE,
    show_rownames = TRUE, show_colnames = TRUE,
    labels_col = age_labels,
    fontsize_row = 9, fontsize_col = 11, fontfamily = "Helvetica",
    treeheight_row = tree_height, treeheight_col = 0,
    border_color = NA, na_col = "white",
    legend = TRUE, legend_title = "Z-score",
    legend_breaks = legend_breaks, legend_labels = as.character(legend_breaks),
    legend_spacing = 1.5,
    main = title)
  dev.off()
}

# ------------------- 生成7张图 -------------------
cat("1/7 总图 不排序（无聚类）\n")
plot_no_cluster_oup(mat_fixed,      "1_ALL_fixed_NOclust",         18, title = NA)

cat("2/7 总图 排序（无聚类）\n")
plot_no_cluster_oup(mat_all_sorted, "2_ALL_sorted_NOclust",        18, title = NA)

cat("3/7 小肠 排序（无聚类）\n")
plot_no_cluster_oup(mat_int_sorted, "3_Intestine_sorted_NOclust",  14, title = "Small Intestine")

cat("4/7 结肠 排序（无聚类）\n")
plot_no_cluster_oup(mat_col_sorted, "4_Colon_sorted_NOclust",      14, title = "Colon")

cat("5/7 总图 固定顺序（聚类）\n")
plot_cluster_oup(mat_fixed,         "5_ALL_fixed_CLUST",           18, tree_height = 60, title = NA)

cat("6/7 小肠 排序（聚类）\n")
plot_cluster_oup(mat_int_sorted,    "6_Intestine_sorted_CLUST",    14, tree_height = 50, title = "Small Intestine")

cat("7/7 结肠 排序（聚类）\n")
plot_cluster_oup(mat_col_sorted,    "7_Colon_sorted_CLUST",        14, tree_height = 50, title = "Colon")

cat("\n🎉 完成！图例与热图主体间距已增大（legend_spacing = 1.5 cm）。\n")
