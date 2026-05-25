#!/bin/bash
#SBATCH --job-name=daf_imotif
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=daf_%j.out
#SBATCH --error=daf_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate /datapool/home/2023200496/miniconda3/envs/bcftools

VCF="ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz"

# 使用修复后的纯坐标BED文件（无chr前缀）
BED_IMOTIF="human_specific.bed3"

echo "=== 提取人类特异 iM 位点频率 ==="
bcftools query -R ${BED_IMOTIF} -f '%CHROM\t%POS\t%AF\n' ${VCF} > imotif_af_info.tsv
echo "iM 位点提取完成，行数：$(wc -l < imotif_af_info.tsv)"
head -3 imotif_af_info.tsv

# 背景位点（同样用无chr前缀的坐标）
echo "=== 提取背景位点 ==="
bcftools query -f '%CHROM\t%POS\n' ${VCF} | shuf -n 100000 > random_100k_sites.pos
bcftools query -R random_100k_sites.pos -f '%CHROM\t%POS\t%AF\n' ${VCF} > background_af_info.tsv
echo "背景位点提取完成，行数：$(wc -l < background_af_info.tsv)"

# R分析
echo "=== 运行 R 分析 ==="
conda activate R
Rscript analyze_DAF.R
echo "=== 全部完成 ==="
