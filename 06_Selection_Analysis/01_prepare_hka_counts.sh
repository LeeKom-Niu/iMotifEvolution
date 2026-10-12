#!/bin/bash
# 01_prepare_hka_counts.sh
# Purpose: extract stems/loops in promoter regions, count divergence, SNP and 4FD counts
# Output: counts_promoter.tsv (counts for all groups)

set -e

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/19_modified_analysis/09_HKA_PhyloFit
#   /datapool/home/2023200496/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/functionalOutputs/Homo_sapiens/promoter_regions.bed.gz
# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/home/2023200496/niulk/my_project
BASE_DIR="${IMOTIF_BASE:-.}"

# -------- configuration --------
WORKDIR="${BASE_DIR}/19_modified_analysis/09_HKA_PhyloFit"
INPUT_DIR="${WORKDIR}/input"
OUTPUT_DIR="${WORKDIR}/output"

STEM_BED="${INPUT_DIR}/stems.bed"
LOOP_BED="${INPUT_DIR}/loops.bed"
DIVERGENCE_BED="${INPUT_DIR}/all_divergence_genomewide.bed"
SNP_BED="${INPUT_DIR}/hprc.snps.t2t.bed"
PROMOTER_GZ="${BASE_DIR}/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/functionalOutputs/Homo_sapiens/promoter_regions.bed.gz"
FOURFD_ALL="${INPUT_DIR}/4fd_all.bed"

mkdir -p "${OUTPUT_DIR}"

# -------- Step 1: convert divergence site coordinates --------
if [ ! -f "${OUTPUT_DIR}/all_divergence_pos.bed" ]; then
    echo "Converting divergence BED to standard integers..."
    awk -v OFS='\t' '{print $1, int($2), int($3)}' "${DIVERGENCE_BED}" > "${OUTPUT_DIR}/all_divergence_pos.bed"
fi
DIVERGENCE_POS="${OUTPUT_DIR}/all_divergence_pos.bed"

# -------- Step 2: generate 4FD sites (aligned regions only) and deduplicate --------
if [ ! -f "${OUTPUT_DIR}/4fd_in_divergence.bed" ]; then
    echo "Generating 4FD sites within divergence regions..."
    bedtools intersect -a "${FOURFD_ALL}" -b "${DIVERGENCE_POS}" -u | sort -k1,1 -k2,2n -k3,3n | uniq > "${OUTPUT_DIR}/4fd_in_divergence.bed"
else
    echo "4FD in divergence file already exists. Skipping."
fi
FOURFD="${OUTPUT_DIR}/4fd_in_divergence.bed"

# -------- Step 3: count SNP and divergence on 4FD --------
echo "Counting P_4fd and D_4fd..."
P_4fd=$(bedtools intersect -a "${SNP_BED}" -b "${FOURFD}" | wc -l)
D_4fd=$(bedtools intersect -a "${DIVERGENCE_POS}" -b "${FOURFD}" | wc -l)
echo -e "4FD\t${P_4fd}\t${D_4fd}" > "${OUTPUT_DIR}/4fd_counts.txt"

# -------- Step 4: process promoter regions for each group --------
echo -e "Group\tP_stem\tD_stem\tP_loop\tD_loop\tP_4fd\tD_4fd" > "${OUTPUT_DIR}/counts_promoter.tsv"

# hardcode groups directly to avoid variable pollution
for group in hominid homininae hominini; do
    echo "Processing group: ${group}"
    G4_BED="${INPUT_DIR}/groups_0based/human.${group}G4s_0based.bed"
    if [ ! -f "${G4_BED}" ]; then
        echo "Warning: ${G4_BED} not found. Skipping ${group}."
        continue
    fi
    
    # extract G4s within promoters
    PROM_G4="${OUTPUT_DIR}/${group}_promoter_G4s.bed"
    bedtools intersect -a "${G4_BED}" -b <(zcat "${PROMOTER_GZ}") -wa > "${PROM_G4}"
    
    # extract the stems and loops corresponding to these G4s
    PROM_STEM="${OUTPUT_DIR}/${group}_promoter_stems.bed"
    PROM_LOOP="${OUTPUT_DIR}/${group}_promoter_loops.bed"
    bedtools intersect -a "${STEM_BED}" -b "${PROM_G4}" -wa > "${PROM_STEM}"
    bedtools intersect -a "${LOOP_BED}" -b "${PROM_G4}" -wa > "${PROM_LOOP}"
    
    # count
    D_stem=$(bedtools intersect -a "${DIVERGENCE_POS}" -b "${PROM_STEM}" | wc -l)
    P_stem=$(bedtools intersect -a "${SNP_BED}" -b "${PROM_STEM}" | wc -l)
    D_loop=$(bedtools intersect -a "${DIVERGENCE_POS}" -b "${PROM_LOOP}" | wc -l)
    P_loop=$(bedtools intersect -a "${SNP_BED}" -b "${PROM_LOOP}" | wc -l)
    
    echo -e "${group}\t${P_stem}\t${D_stem}\t${P_loop}\t${D_loop}\t${P_4fd}\t${D_4fd}" >> "${OUTPUT_DIR}/counts_promoter.tsv"
done

echo "All done! Counts saved to ${OUTPUT_DIR}/counts_promoter.tsv"
echo "Now run Rscript 02_fisher_test_hka.R"
