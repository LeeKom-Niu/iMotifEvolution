#!/usr/bin/env Rscript
# 02_fisher_test_hka.R
# read counts for the promoter region, run Fisher test

library(data.table)

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
# Root directory for input/output files.
# On the original cluster: Sys.setenv(IMOTIF_BASE = "/datapool/home/2023200496/niulk/my_project")
BASE_DIR <- Sys.getenv("IMOTIF_BASE", unset = ".")

workdir <- file.path(BASE_DIR, "19_modified_analysis/09_HKA_PhyloFit")
output_dir <- file.path(workdir, "output")

region <- "promoter"
counts_file <- file.path(output_dir, paste0("counts_", region, ".tsv"))
if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
}

dt <- fread(counts_file, header = TRUE)
if (nrow(dt) == 0) {
    stop("Counts file is empty: ", counts_file)
}

P4 <- dt[1, P_4fd]
D4 <- dt[1, D_4fd]

all_results <- data.table()

for (i in 1:nrow(dt)) {
    group <- dt[i, Group]
    P_stem <- dt[i, P_stem]
    D_stem <- dt[i, D_stem]
    P_loop <- dt[i, P_loop]
    D_loop <- dt[i, D_loop]

    if (any(c(P_stem, D_stem, P_loop, D_loop) == 0)) {
        cat("Skipping", group, region, "due to zero counts\n")
        next
    }

    # stem vs loop
    mat1 <- matrix(c(P_stem, D_stem, P_loop, D_loop), nrow = 2, byrow = TRUE)
    test1 <- fisher.test(mat1)
    all_results <- rbind(all_results, data.table(
        Region = region,
        Group = group,
        Comparison = "Stem_vs_Loop",
        OR = test1$estimate,
        p_value = test1$p.value,
        CI_lower = test1$conf.int[1],
        CI_upper = test1$conf.int[2]
    ))

    # stem vs 4FD
    mat2 <- matrix(c(P_stem, D_stem, P4, D4), nrow = 2, byrow = TRUE)
    test2 <- fisher.test(mat2)
    all_results <- rbind(all_results, data.table(
        Region = region,
        Group = group,
        Comparison = "Stem_vs_4FD",
        OR = test2$estimate,
        p_value = test2$p.value,
        CI_lower = test2$conf.int[1],
        CI_upper = test2$conf.int[2]
    ))

    # loop vs 4FD
    mat3 <- matrix(c(P_loop, D_loop, P4, D4), nrow = 2, byrow = TRUE)
    test3 <- fisher.test(mat3)
    all_results <- rbind(all_results, data.table(
        Region = region,
        Group = group,
        Comparison = "Loop_vs_4FD",
        OR = test3$estimate,
        p_value = test3$p.value,
        CI_lower = test3$conf.int[1],
        CI_upper = test3$conf.int[2]
    ))
}

output_file <- file.path(output_dir, "fisher_results_promoter.tsv")
fwrite(all_results, output_file, sep = "\t")
cat("Fisher results saved to", output_file, "\n")
print(all_results)
