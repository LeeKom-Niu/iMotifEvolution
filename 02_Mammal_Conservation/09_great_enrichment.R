# 06_great_enrichment.R
# 功能：对C1、C2、C3组的iMotif进行GREAT富集分析（GO:BP）
# 输入文件：imotif_clean.bed4, imotif_group.tsv
# 输出文件：great_C1_results.tsv, great_C2_results.tsv, great_C3_results.tsv

# 加载必要的包
library(rGREAT)
library(GenomicRanges)
library(dplyr)

# 1. 读取iMotif坐标文件（BED4格式）
imotif_bed <- read.table("imotif_clean.bed4", header = FALSE,
                         col.names = c("chr", "start", "end", "name"),
                         stringsAsFactors = FALSE)

# 2. 读取分组信息
group_df <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)
# 确保分组数据与iMotif坐标对应（通过name列匹配）
imotif_bed <- imotif_bed %>%
  left_join(group_df[, c("iMotif", "group")], by = c("name" = "iMotif"))

# 检查是否有缺失group的情况
if (any(is.na(imotif_bed$group))) {
  warning("部分iMotif无分组信息，这些区域将被排除在分组富集之外，但背景包含所有iMotif。")
}

# 3. 按分组提取前三列坐标
C1_bed <- imotif_bed %>% filter(group == "C1") %>% select(chr, start, end)
C2_bed <- imotif_bed %>% filter(group == "C2") %>% select(chr, start, end)
C3_bed <- imotif_bed %>% filter(group == "C3") %>% select(chr, start, end)
bg_bed <- imotif_bed %>% select(chr, start, end)   # 所有iMotif作为背景

# 4. 转换为GRanges对象
make_gr <- function(df) {
  makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
}
C1_gr <- make_gr(C1_bed)
C2_gr <- make_gr(C2_bed)
C3_gr <- make_gr(C3_bed)
bg_gr <- make_gr(bg_bed)

# 5. 运行GREAT
# 设置基因组为hg38，基因集为GO:BP
res_C1 <- great(C1_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C2 <- great(C2_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C3 <- great(C3_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)

# 6. 提取富集表，并进行过滤（使用实际列名 p_adjust_hyper 和 p_adjust）
extract_filtered <- function(res) {
  tb <- getEnrichmentTable(res)
  # 打印列名以帮助调试
  cat("当前结果表的列名：\n")
  print(colnames(tb))
  
  # 使用 p_adjust_hyper（超几何检验FDR）和 p_adjust（二项检验FDR）过滤
  if ("p_adjust_hyper" %in% colnames(tb) & "p_adjust" %in% colnames(tb)) {
    tb_filtered <- tb %>%
      filter(p_adjust_hyper < 0.05 & p_adjust < 0.05) %>%
      arrange(p_adjust_hyper, p_adjust)
  } else {
    stop("无法识别校正后列名，请根据打印的列名手动修改过滤条件。")
  }
  return(tb_filtered)
}

tb_C1 <- extract_filtered(res_C1)
tb_C2 <- extract_filtered(res_C2)
tb_C3 <- extract_filtered(res_C3)

# 7. 保存结果到文件
write.table(tb_C1, "great_C1_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C2, "great_C2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C3, "great_C3_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 打印统计信息
cat("C1组富集到的GO term数量:", nrow(tb_C1), "\n")
cat("C2组富集到的GO term数量:", nrow(tb_C2), "\n")
cat("C3组富集到的GO term数量:", nrow(tb_C3), "\n")
