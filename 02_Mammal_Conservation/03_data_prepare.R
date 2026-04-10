library(tidyverse)

# 读取物种列表
species_list <- readLines("sp240_nonhuman.txt")

# 初始化一个空的数据框，行名为 iMotif 名称
# 先读取第一个物种获取行名
first <- read.table("results/align_ratio/Acinonyx_jubatus.imotif.ar.txt", 
                    col.names = c("name","size","covered","sum","mean0","mean"))
imotif_names <- first$name

# 构建矩阵
mat <- matrix(NA, nrow = length(imotif_names), ncol = length(species_list))
colnames(mat) <- species_list
rownames(mat) <- imotif_names

# 填充每个物种的 mean 值
for (i in seq_along(species_list)) {
  sp <- species_list[i]
  file_path <- paste0("results/align_ratio/", sp, ".imotif.ar.txt")
  if (file.exists(file_path)) {
    dat <- read.table(file_path, col.names = c("name","size","covered","sum","mean0","mean"))
    # 确保顺序一致
    mat[dat$name, i] <- dat$mean
  } else {
    warning(paste("Missing file for", sp))
  }
}

# 保存矩阵
write.table(mat, "imotif_coverage_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE)
