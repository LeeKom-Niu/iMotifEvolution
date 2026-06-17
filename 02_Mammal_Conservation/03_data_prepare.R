library(tidyverse)

# 读取物种列表
species_list <- readLines("sp240_nonhuman.txt")

# Initialize an empty data frame, row names are iMotif names
# First read the first species to get row names
first <- read.table("results/align_ratio/Acinonyx_jubatus.imotif.ar.txt", 
                    col.names = c("name","size","covered","sum","mean0","mean"))
imotif_names <- first$name

# Build matrix
mat <- matrix(NA, nrow = length(imotif_names), ncol = length(species_list))
colnames(mat) <- species_list
rownames(mat) <- imotif_names

# Fill mean values for each species
for (i in seq_along(species_list)) {
  sp <- species_list[i]
  file_path <- paste0("results/align_ratio/", sp, ".imotif.ar.txt")
  if (file.exists(file_path)) {
    dat <- read.table(file_path, col.names = c("name","size","covered","sum","mean0","mean"))
    # Ensure consistent order
    mat[dat$name, i] <- dat$mean
  } else {
    warning(paste("Missing file for", sp))
  }
}

# Save matrix
write.table(mat, "imotif_coverage_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE)
