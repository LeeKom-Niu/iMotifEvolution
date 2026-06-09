library(tidyverse)
species_list <- readLines("sp240_nonhuman.txt")
first <- read.table("results/align_ratio/Acinonyx_jubatus.imotif.ar.txt", 
                    col.names = c("name","size","covered","sum","mean0","mean"))
imotif_names <- first$name
mat <- matrix(NA, nrow = length(imotif_names), ncol = length(species_list))
colnames(mat) <- species_list
rownames(mat) <- imotif_names
for (i in seq_along(species_list)) {
  sp <- species_list[i]
  file_path <- paste0("results/align_ratio/", sp, ".imotif.ar.txt")
  if (file.exists(file_path)) {
    dat <- read.table(file_path, col.names = c("name","size","covered","sum","mean0","mean"))
    mat[dat$name, i] <- dat$mean
  } else {
    warning(paste("Missing file for", sp))
  }
}
write.table(mat, "imotif_coverage_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE)
