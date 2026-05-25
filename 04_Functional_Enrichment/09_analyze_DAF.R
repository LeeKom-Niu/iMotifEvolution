library(ggplot2)

# 读取 iM 位点数据
imotif <- read.table("imotif_af_info.tsv", header = FALSE, col.names = c("chr", "pos", "AF"))
# 读取背景数据
bg <- read.table("background_af_info.tsv", header = FALSE, col.names = c("chr", "pos", "AF"))

# 过滤缺失或无效 AF
imotif <- imotif[!is.na(imotif$AF) & imotif$AF >= 0 & imotif$AF <= 1, ]
bg     <- bg[!is.na(bg$AF) & bg$AF >= 0 & bg$AF <= 1, ]

# 计算 minor allele frequency (MAF)
imotif$MAF <- pmin(imotif$AF, 1 - imotif$AF)
bg$MAF     <- pmin(bg$AF, 1 - bg$AF)

cat("iM 位点数:", nrow(imotif), "\n")
cat("背景位点数:", nrow(bg), "\n")
cat("iM 位点平均 MAF:", mean(imotif$MAF), "\n")
cat("背景平均 MAF:", mean(bg$MAF), "\n")

# KS 检验
ks <- ks.test(imotif$MAF, bg$MAF)
cat("KS test p-value:", ks$p.value, "\n")

# 绘制叠加直方图
p <- ggplot() +
  # 背景先画，再用 iM 画在上面
  geom_histogram(data = bg, aes(x = MAF, y = after_stat(density), fill = "Genomic background"),
                 alpha = 0.4, binwidth = 0.05, boundary = 0) +
  geom_histogram(data = imotif, aes(x = MAF, y = after_stat(density), fill = "Human-specific iMs"),
                 alpha = 0.6, binwidth = 0.05, boundary = 0) +
  scale_fill_manual(name = "",   # 去掉图例标题（或改成 name="Group"）
                    values = c("Human-specific iMs" = "red",
                               "Genomic background" = "grey70")) +
  labs(x = "Minor Allele Frequency (MAF)", y = "Density",
       title = "Frequency spectrum of human-specific i-Motif loci",
       subtitle = paste0("KS test p = ", format.pval(ks$p.value, digits = 3),
                         " | iM mean MAF = ", round(mean(imotif$MAF), 3))) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")   # 图例放到底部更清晰

ggsave("daf_distribution.pdf", p, width = 8, height = 5)
ggsave("daf_distribution.png", p, width = 8, height = 5, dpi = 300)
