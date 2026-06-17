# 06_great_enrichment.R
# Function: Perform GREAT enrichment analysis (GO:BP) on C1, C2, C3 iMotif groups
# Input files: imotif_clean.bed4, imotif_group.tsv
# Output files: great_C1_results.tsv, great_C2_results.tsv, great_C3_results.tsv

# Load required packages
library(rGREAT)
library(GenomicRanges)
library(dplyr)

# 1. Read iMotif coordinate file (BED4 format)
imotif_bed <- read.table("imotif_clean.bed4", header = FALSE,
                         col.names = c("chr", "start", "end", "name"),
                         stringsAsFactors = FALSE)

# 2. Read group information
group_df <- read.table("imotif_group.tsv", header = TRUE, stringsAsFactors = FALSE)
# Ensure group data matches iMotif coordinates (via name column)
imotif_bed <- imotif_bed %>%
  left_join(group_df[, c("iMotif", "group")], by = c("name" = "iMotif"))

# 检查是否有缺失group的情况
if (any(is.na(imotif_bed$group))) {
  warning("Some iMotifs lack group info; these will be excluded from group enrichment, but background includes all iMotifs.")
}

# 3. Extract first three columns by group
C1_bed <- imotif_bed %>% filter(group == "C1") %>% select(chr, start, end)
C2_bed <- imotif_bed %>% filter(group == "C2") %>% select(chr, start, end)
C3_bed <- imotif_bed %>% filter(group == "C3") %>% select(chr, start, end)
bg_bed <- imotif_bed %>% select(chr, start, end)   # All iMotifs as background

# 4. Convert to GRanges objects
make_gr <- function(df) {
  makeGRangesFromDataFrame(df, seqnames.field = "chr", start.field = "start", end.field = "end")
}
C1_gr <- make_gr(C1_bed)
C2_gr <- make_gr(C2_bed)
C3_gr <- make_gr(C3_bed)
bg_gr <- make_gr(bg_bed)

# 5. Run GREAT
# Set genome to hg38, gene set to GO:BP
res_C1 <- great(C1_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C2 <- great(C2_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)
res_C3 <- great(C3_gr, gene_sets = "GO:BP", tss_source = "hg38", background = bg_gr)

# 6. Extract enrichment table and filter (using actual column names p_adjust_hyper and p_adjust)
extract_filtered <- function(res) {
  tb <- getEnrichmentTable(res)
  # Print column names for debugging
  cat("Current enrichment table column names:\n")
  print(colnames(tb))
  
  # Filter using p_adjust_hyper (hypergeometric test FDR) and p_adjust (binomial test FDR)
  if ("p_adjust_hyper" %in% colnames(tb) & "p_adjust" %in% colnames(tb)) {
    tb_filtered <- tb %>%
      filter(p_adjust_hyper < 0.05 & p_adjust < 0.05) %>%
      arrange(p_adjust_hyper, p_adjust)
  } else {
    stop("Cannot recognize adjusted column names. Please modify filter criteria based on printed column names.")
  }
  return(tb_filtered)
}

tb_C1 <- extract_filtered(res_C1)
tb_C2 <- extract_filtered(res_C2)
tb_C3 <- extract_filtered(res_C3)

# 7. Save results to files
write.table(tb_C1, "great_C1_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C2, "great_C2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tb_C3, "great_C3_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# 打印统计信息
cat("GO terms enriched in C1 group:", nrow(tb_C1), "\n")
cat("GO terms enriched in C2 group:", nrow(tb_C2), "\n")
cat("GO terms enriched in C3 group:", nrow(tb_C3), "\n")
