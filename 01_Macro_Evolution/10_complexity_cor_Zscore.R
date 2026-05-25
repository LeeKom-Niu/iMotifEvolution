# ============================================================
# 脚本1：i-Motif Z-score 偏相关分析（主分析）
# 输出目录：./results_zscore/
# 三张图：无控制 → 双控制 → 三控制
# 已修改：同时输出 PDF + PNG 格式
# ============================================================

library(ppcor)
library(ggplot2)
library(ggpubr)
library(dplyr)

dir.create("./results_zscore", showWarnings = FALSE, recursive = TRUE)

df <- read.csv("merged_analysis.csv", stringsAsFactors = FALSE)

# 进化顺序
complexity_order <- c(
  "Archaea", "Bacteria", "Protozoa", "Plant", "Fungi",
  "Invertebrate", "Vertebrate Other", "Mammalian"
)

df$complexity <- as.numeric(factor(df$Classification_9cat, levels = complexity_order))

required_cols <- c("robust_z_score", "complexity", "Genome_size_Mb",
                   "Gene_density_per_Mb", "GC_content")
df_clean <- df[complete.cases(df[, required_cols]), ]

cat(sprintf("原始数据: %d 行\n", nrow(df)))
cat(sprintf("清洗后: %d 行\n", nrow(df_clean)))

if (nrow(df_clean) == 0) stop("数据清洗后没有剩余行！请检查列名。")

df_clean$log_genome_size <- log10(df_clean$Genome_size_Mb)

# ---- 统计分析 ----
cat("\n========== Z-score 分析结果 ==========\n")

cor_simple <- cor.test(df_clean$robust_z_score, df_clean$complexity,
                       method = "spearman", exact = FALSE)
cat(sprintf("简单相关 rho = %.3f, P = %s\n", cor_simple$estimate,
           format.pval(cor_simple$p.value, digits = 3, eps = 1e-300)))

cor_partial <- pcor.test(
  x = df_clean$robust_z_score,
  y = df_clean$complexity,
  z = df_clean[, c("log_genome_size", "Gene_density_per_Mb")],
  method = "spearman"
)
cat(sprintf("偏相关（控制Size+Density） rho = %.3f, P = %s\n",
           cor_partial$estimate,
           format.pval(cor_partial$p.value, digits = 3, eps = 1e-300)))

cor_partial_strict <- pcor.test(
  x = df_clean$robust_z_score,
  y = df_clean$complexity,
  z = df_clean[, c("log_genome_size", "Gene_density_per_Mb", "GC_content")],
  method = "spearman"
)
cat(sprintf("偏相关（控制Size+Density+GC） rho = %.3f, P = %s\n",
           cor_partial_strict$estimate,
           format.pval(cor_partial_strict$p.value, digits = 3, eps = 1e-300)))

lm1 <- lm(robust_z_score ~ complexity, data = df_clean)
lm2 <- lm(robust_z_score ~ complexity + log_genome_size + Gene_density_per_Mb, data = df_clean)
lm3 <- lm(robust_z_score ~ complexity + log_genome_size + Gene_density_per_Mb + GC_content, data = df_clean)

cat(sprintf("模型1 R² = %.4f, complexity P = %s\n",
           summary(lm1)$r.squared,
           format.pval(coef(summary(lm1))[2,4], digits = 3, eps = 1e-300)))
cat(sprintf("模型2 R² = %.4f, complexity P = %s\n",
           summary(lm2)$r.squared,
           format.pval(coef(summary(lm2))[2,4], digits = 3, eps = 1e-300)))
cat(sprintf("模型3 R² = %.4f, complexity P = %s\n",
           summary(lm3)$r.squared,
           format.pval(coef(summary(lm3))[2,4], digits = 3, eps = 1e-300)))

# ---- 可视化准备 ----
class_colors <- c(
  "Archaea" = "#E41A1C",
  "Bacteria" = "#377EB8",
  "Protozoa" = "#984EA3",
  "Plant" = "#4DAF4A",
  "Fungi" = "#FF7F00",
  "Invertebrate" = "#A65628",
  "Vertebrate Other" = "#F781BF",
  "Mammalian" = "#999999"
)

df_clean$Classification_ordered <- factor(df_clean$Classification_9cat,
                                          levels = names(class_colors))

df_clean$complexity_label <- factor(df_clean$complexity, levels = 1:8,
  labels = names(class_colors))

# 计算残差
df_clean$residual_Z_partial <- resid(
  lm(robust_z_score ~ log_genome_size + Gene_density_per_Mb, data = df_clean))
df_clean$residual_complexity_partial <- resid(
  lm(complexity ~ log_genome_size + Gene_density_per_Mb, data = df_clean))

df_clean$residual_Z_strict <- resid(
  lm(robust_z_score ~ log_genome_size + Gene_density_per_Mb + GC_content, data = df_clean))
df_clean$residual_complexity_strict <- resid(
  lm(complexity ~ log_genome_size + Gene_density_per_Mb + GC_content, data = df_clean))

# P值辅助函数
fmt_p <- function(p) {
  ifelse(p < 1e-300, "P < 1e-300",
         paste0("P = ", format(p, scientific = TRUE, digits = 3)))
}

# 通用主题
theme_paper <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# ---- 统一保存图片函数：同时输出 PDF + PNG ----
save_plot <- function(plot_obj, filename, width=8, height=7, dpi=300) {
  # 保存PDF
  ggsave(paste0("./results_zscore/", filename, ".pdf"), 
         plot = plot_obj, width = width, height = height, device = "pdf")
  # 保存PNG（高清300dpi）
  ggsave(paste0("./results_zscore/", filename, ".png"), 
         plot = plot_obj, width = width, height = height, 
         dpi = dpi, device = "png")
}

# ---- 图1：无控制散点图 ----
subtitle_raw <- paste0(
  "Spearman rho = ", round(cor_simple$estimate, 3),
  ", ", fmt_p(cor_simple$p.value),
  "\nNo controlling variables"
)

p1 <- ggplot(df_clean, aes(x = complexity, y = robust_z_score)) +
  geom_jitter(aes(color = Classification_ordered), size = 2.5, alpha = 0.7, width = 0.2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
  scale_color_manual(values = class_colors, name = "Clade",
                     breaks = names(class_colors)) +
  scale_x_continuous(breaks = 1:8,
                     labels = names(class_colors)) +
  labs(x = "Biological Complexity (Evolutionary Order)",
       y = "i-Motif Robust Z-score",
       title = "i-Motif Enrichment (Z-score) vs. Biological Complexity",
       subtitle = subtitle_raw) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
save_plot(p1, "Figure_Raw_Correlation", width=8, height=7)
cat("图1（无控制）已保存 PDF + PNG\n")

# ---- 图2：双控制散点图 ----
subtitle_partial <- paste0(
  "Partial Spearman rho = ", round(cor_partial$estimate, 3),
  ", ", fmt_p(cor_partial$p.value),
  "\nControlling for: log10(genome size) + gene density"
)

p2 <- ggplot(df_clean, aes(x = residual_complexity_partial, y = residual_Z_partial)) +
  geom_point(aes(color = Classification_ordered), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
  scale_color_manual(values = class_colors, name = "Clade",
                     breaks = names(class_colors)) +
  labs(x = "Biological Complexity\n(adjusted for genome size & gene density)",
       y = "i-Motif Robust Z-score\n(adjusted for genome size & gene density)",
       title = "i-Motif Enrichment (Z-score) vs. Biological Complexity",
       subtitle = subtitle_partial) +
  theme_paper +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
save_plot(p2, "Figure_PartialCorrelation_Double", width=8, height=7)
cat("图2（双控制）已保存 PDF + PNG\n")

# ---- 图3：三控制散点图 ----
subtitle_strict <- paste0(
  "Partial Spearman rho = ", round(cor_partial_strict$estimate, 3),
  ", ", fmt_p(cor_partial_strict$p.value),
  "\nControlling for: log10(genome size) + gene density + GC content"
)

p3 <- ggplot(df_clean, aes(x = residual_complexity_strict, y = residual_Z_strict)) +
  geom_point(aes(color = Classification_ordered), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.2) +
  scale_color_manual(values = class_colors, name = "Clade",
                     breaks = names(class_colors)) +
  labs(x = "Biological Complexity\n(adjusted for genome size, gene density & GC%)",
       y = "i-Motif Robust Z-score\n(adjusted for genome size, gene density & GC%)",
       title = "i-Motif Enrichment (Z-score) vs. Biological Complexity",
       subtitle = subtitle_strict) +
  theme_paper +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
save_plot(p3, "Figure_PartialCorrelation_Triple", width=8, height=7)
cat("图3（三控制）已保存 PDF + PNG\n")

# ---- 图4：箱线图 ----
p4 <- ggplot(df_clean, aes(x = complexity_label, y = robust_z_score, fill = complexity_label)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.3, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.1, color = "grey30") +
  scale_fill_manual(values = class_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = "dotted", color = "grey70", linewidth = 0.5) +
  labs(x = "Biological Complexity (Evolutionary Order)",
       y = "i-Motif Robust Z-score",
       title = "Distribution of i-Motif Enrichment (Z-score) across Evolutionary Grades") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 15))
save_plot(p4, "Figure_Zscore_by_Complexity", width=12, height=6)
cat("图4（箱线图）已保存 PDF + PNG\n")

# ---- 图5：分面回归图 ----
p5 <- ggplot(df_clean, aes(x = log_genome_size, y = robust_z_score)) +
  geom_point(aes(color = complexity_label), size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  facet_wrap(~ complexity_label, ncol = 4, scales = "free_x") +
  scale_color_manual(values = class_colors, guide = "none") +
  labs(x = "log10(Genome Size / Mb)", y = "i-Motif Robust Z-score",
       title = "Z-score vs. Genome Size (stratified by evolutionary grade)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 15), strip.text = element_text(size = 10))
save_plot(p5, "Figure_Zscore_vs_GenomeSize_by_Complexity", width=14, height=10)
cat("图5（分面回归图）已保存 PDF + PNG\n")

cat("\n========== Z-score 分析完成 ==========\n")
