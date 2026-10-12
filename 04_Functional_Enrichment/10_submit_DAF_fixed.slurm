#!/bin/bash
#SBATCH --job-name=daf_imotif_fixed
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=daf_fixed_%j.out
#SBATCH --error=daf_fixed_%j.err

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
#   /datapool/home/2023200496/miniconda3/envs/bcftools

# Conda lives under $HOME (not under IMOTIF_BASE).
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate $HOME/miniconda3/envs/bcftools

VCF="ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
BED_IMOTIF="human_specific.bed3"
BACKGROUND_POS="fixed_background_100k.pos"

echo "=== Extract human-specific iM site frequency ==="
bcftools query -R ${BED_IMOTIF} -f '%CHROM\t%POS\t%AF\n' ${VCF} > imotif_af_info.tsv
echo "iM site extraction done, line count: $(wc -l < imotif_af_info.tsv)"
head -3 imotif_af_info.tsv

# Generate or reuse the fixed background site set
if [ ! -f ${BACKGROUND_POS} ]; then
    echo "=== Generate fixed background site set (100k) ==="
    bcftools query -f '%CHROM\t%POS\n' ${VCF} | shuf -n 100000 > ${BACKGROUND_POS}
else
    echo "=== Reusing existing fixed background site set: ${BACKGROUND_POS} ==="
fi

echo "=== Extract background site frequency ==="
bcftools query -R ${BACKGROUND_POS} -f '%CHROM\t%POS\t%AF\n' ${VCF} > background_af_info.tsv
echo "Background site extraction done, line count: $(wc -l < background_af_info.tsv)"

# Run R analysis (using improved script)
echo "=== Run R analysis ==="
conda activate R
Rscript analyze_DAF_fixed.R
echo "=== All done ==="
