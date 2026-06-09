
library(rGREAT)
library(GenomicRanges)
library(dplyr)
imotif_bed <- read.table("imotif_clean.bed4", header = FALSE,
                         col.names = c("chr", "start", "end", "name"),
                         stringsAsFactors = FALSE)
group_df <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)
imotif_bed <- imotif_bed %>%
  left_join(group_df[, c("iMotif", "group")], by = c("name" = "iMotif"))
if (any(is.na(imotif_bed$group))) {
  warning("部分iMotif无分组信息，这些区域将被排除在分组富集之外，但背景包含所有iMotif。")
}
C1_bed <- imotif_bed %>% filter(group == "C1") %>% select(chr, start, end)
C2_bed <- imotif_bed %>% filter(group == "C2") %>% select(chr, start, end)
C3_bed <- imotif_bed %>% filter(group == "C3") %>% select(chr, start, end)
bg_bed <- imotif_bed %>% select(chr, start, end)
make_gr <- function(df) {
  makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
}
C1_gr <- make_gr(C1_bed)
C2_gr <- make_gr(C2_bed)
C3_gr <- make_gr(C3_bed)
bg_gr <- make_gr(bg_bed)
res_C1 <- great(C1_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C2 <- great(C2_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C3 <- great(C3_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
extract_filtered <- function(res) {
  tb <- getEnrichmentTable(res)
  cat("当前结果表的列名：\n")
  print(colnames(tb))
  
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
write.table(tb_C1, "great_C1_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C2, "great_C2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C3, "great_C3_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
cat("C1组富集到的GO term数量:", nrow(tb_C1), "\n")
cat("C2组富集到的GO term数量:", nrow(tb_C2), "\n")
cat("C3组富集到的GO term数量:", nrow(tb_C3), "\n")
