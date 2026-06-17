library(dplyr)

mat <- as.matrix(read.table("imotif_coverage_matrix.txt", header = TRUE, row.names = 1))
N1 <- rowSums(mat >= 0.9, na.rm = TRUE)
N2 <- rowSums(mat <= 0.1, na.rm = TRUE)
N3 <- rowSums(mat >= 0.5, na.rm = TRUE)

group <- case_when(
  N1 >= 120 & N2 <= 25 ~ "C1",
  N1 >= 20 & N1 <= 50 & N2 <= 120 ~ "C2",
  N1 <= 50 & N2 >= 180 ~ "C3",
  TRUE ~ "Other"
)

group_df <- data.frame(iMotif = rownames(mat), N1, N2, N3, group)
write.table(group_df, "imotif_group.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
