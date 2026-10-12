#!/bin/bash
# run_PhyloFit_promoter.sh
set -e

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
#   /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/home/2023200496/niulk/my_project
BASE_DIR="${IMOTIF_BASE:-.}"

WORKDIR="${BASE_DIR}/19_modified_analysis/09_HKA_PhyloFit"
cd "$WORKDIR"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate base
export MPLBACKEND=Agg

echo "Step 1: generate MAF, PhyloFit models, and write branch_lengths.tsv"
sbatch --wait 04_run_phylofit_models.sh

echo "Step 2: plot PhyloFit promoter figure"
python 05_plot_phylofit_promoter.py

echo "Done: Figure_phylofit_promoter.pdf / .png"
