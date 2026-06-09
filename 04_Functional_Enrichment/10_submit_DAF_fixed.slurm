#!/bin/bash
#SBATCH --job-name=daf_imotif_fixed
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=daf_fixed_%j.out
#SBATCH --error=daf_fixed_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate /datapool/home/2023200496/miniconda3/envs/bcftools

VCF="ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"
BED_IMOTIF="human_specific.bed3"
BACKGROUND_POS="fixed_background_100k.pos"

echo "=== 提取人类特异 iM 位点频率 ==="
bcftools query -R ${BED_IMOTIF} -f '%CHROM\t%POS\t%AF\n' ${VCF} > imotif_af_info.tsv
echo "iM 位点提取完成，行数：$(wc -l < imotif_af_info.tsv)"
head -3 imotif_af_info.tsv

# 生成或复用固定背景位点集
if [ ! -f ${BACKGROUND_POS} ]; then
    echo "=== 生成固定背景位点集 (100k) ==="
    bcftools query -f '%CHROM\t%POS\n' ${VCF} | shuf -n 100000 > ${BACKGROUND_POS}
else
    echo "=== 使用已有固定背景位点集: ${BACKGROUND_POS} ==="
fi

echo "=== 提取背景位点频率 ==="
bcftools query -R ${BACKGROUND_POS} -f '%CHROM\t%POS\t%AF\n' ${VCF} > background_af_info.tsv
echo "背景位点提取完成，行数：$(wc -l < background_af_info.tsv)"

# 运行 R 分析（使用改进脚本）
echo "=== 运行 R 分析 ==="
conda activate R
Rscript analyze_DAF_fixed.R
echo "=== 全部完成 ==="
