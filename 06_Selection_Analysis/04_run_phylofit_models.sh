#!/bin/bash
#SBATCH --job-name=phylofit_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=30-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=phylofit_batch_%j.out
#SBATCH --error=phylofit_batch_%j.err

# 04_run_phylofit_models.sh
# Purpose: for each group × region × structure, generate MAF and fit a PhyloFit model,
#       then extract the hs1 branch length from the model file and write it to branch_lengths.tsv

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
#   /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/home/2023200496/niulk/my_project
BASE_DIR="${IMOTIF_BASE:-.}"

# load conda environment (adjust to your actual path)
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate base  # or the name of the environment where you installed cactus and phast

# cd into the working directory
cd "${BASE_DIR}/19_modified_analysis/09_HKA_PhyloFit" || exit 1

# create output directory
mkdir -p maf_output phylofit_models

# define species list (with .2)
SPECIES_LIST="hs1,GCA_028858775.2,GCA_029289425.2,GCA_029281585.2,GCA_028885655.2,GCA_028885625.2"

# define the tree for each group (simplify names, without .2)
declare -A TREES
TREES[hominid]="((hs1,(GCA_028858775,GCA_029289425)),GCA_029281585,(GCA_028885655,GCA_028885625));"
TREES[homininae]="((hs1,(GCA_028858775,GCA_029289425)),GCA_029281585);"
TREES[hominini]="(hs1,(GCA_028858775,GCA_029289425));"

# ensure filter_maf_first.py is in the current directory
if [ ! -f "filter_maf_first.py" ]; then
    echo "ERROR: filter_maf_first.py not found in current directory."
    exit 1
fi

# branch length results table (header), for 05_plot_phylofit_promoter.py to read
echo -e "group\tregion\tpart\ths1_length" > branch_lengths.tsv

# loop over all groups, regions and structures
for group in hominid homininae hominini; do
    tree="${TREES[$group]}"
    for region in promoter; do
        for part in stems loops; do
            bed="hal_input/${group}_${region}_${part}.bed"
            maf="maf_output/${group}_${region}_${part}.maf"
            mod="phylofit_models/${group}_${region}_${part}.mod"
            if [ -f "$bed" ]; then
                echo "Processing $bed -> $maf"
                # generate MAF (using --onlyOrthologs and --onlySequenceNames)
                hal2maf hal_data/7-t2t-apes-2023v2.hal "$maf.tmp" \
                    --refGenome hs1 \
                    --refTargets "$bed" \
                    --targetGenomes "$SPECIES_LIST" \
                    --noAncestors \
                    --onlyOrthologs 
                # filter MAF (deduplicate and simplify names)
                python filter_maf_first.py "$maf.tmp" "$maf"
                rm "$maf.tmp"
                # run PhyloFit
                echo "Running PhyloFit on $maf -> $mod"
                phyloFit --msa-format MAF --tree "$tree" --min-informative 50 "$maf" > "$mod"
                # extract hs1 branch length from the TREE line of the model file, write to the results table
                hs1_len=$(grep "^TREE:" "$mod" | grep -oP "hs1:\K[0-9.]+" | head -1)
                if [ -z "$hs1_len" ]; then
                    hs1_len="NA"
                fi
                echo -e "${group}\t${region}\t${part}\t${hs1_len}" >> branch_lengths.tsv
                echo "$group $region $part hs1_length: $hs1_len"
            else
                echo "Warning: $bed not found, skipping."
            fi
        done
    done
done

echo "All jobs completed. Branch lengths written to branch_lengths.tsv"
