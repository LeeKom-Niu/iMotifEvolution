#!/bin/bash
# run_HKA_promoter.sh
set -e

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/home/2023200496/niulk/my_project
BASE_DIR="${IMOTIF_BASE:-.}"

WORKDIR="${BASE_DIR}/19_modified_analysis/09_HKA_PhyloFit"
cd "$WORKDIR"

echo "Step 1: generate counts_promoter.tsv"
bash 01_prepare_hka_counts.sh

echo "Step 2: Fisher test"
Rscript 02_fisher_test_hka.R

echo "Step 3: plot HKA promoter figure"
Rscript 03_plot_hka_combined.R

echo "Done: output/Figure_HKA_combined.pdf / .png"
