#!/bin/bash
#SBATCH --job-name=imotif_conservation
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --time=30-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=conservation_%j.out
#SBATCH --error=conservation_%j.err

echo "========================================="
echo "开始时间: $(date)"
echo "工作目录: $(pwd)"
echo "========================================="

# 激活环境
source ~/miniconda3/bin/activate cactus

# 创建目录结构
mkdir -p logs
mkdir -p results/align_base
mkdir -p results/align_bedgraph
mkdir -p results/align_bw
mkdir -p results/align_ratio
mkdir -p results/split_base
mkdir -p split_genome
mkdir -p matrix
mkdir -p plots

echo "目录创建完成"

# 步骤1: 生成人类基因组单碱基窗口文件
echo "生成人类基因组单碱基窗口文件..."
if [ ! -f human_genome_windows.4col.bed ]; then
    bedtools makewindows -g hg38.chrom.sizes -w 1 > human_genome_windows.bed
    awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3}' human_genome_windows.bed > human_genome_windows.4col.bed
    echo "完成: $(wc -l < human_genome_windows.4col.bed) 个窗口"
fi

# ========== 修改点1: 使用作者的分块大小和4位后缀 ==========
echo "分割基因组窗口文件..."
total_lines=$(wc -l < human_genome_windows.4col.bed)
lines_per_file=25750000  # 作者使用的每块2575万行
echo "每块行数: $lines_per_file"

# 删除旧的分块（如果存在）
rm -rf split_genome/*

# 使用4位数字后缀 (0000-9999)
split -l $lines_per_file -d -a 4 human_genome_windows.4col.bed split_genome/genome_part.

# 生成分块列表
ls split_genome/ > split_genome_parts.txt
part_count=$(wc -l < split_genome_parts.txt)
echo "完成: $part_count 个分块"
# ======================================================

# 步骤3: 创建物种处理脚本
cat > process_species.sh << 'INNER_EOF'
#!/bin/bash
#SBATCH --job-name=imotif_species
#SBATCH --nodes=1
# ========== 修改点2: 调整每个物种任务的资源 ==========
#SBATCH --cpus-per-task=10      # 每个物种用4个CPU（根据集群情况调整）
#SBATCH --mem=40G               # 每个物种用16G内存
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
# ==================================================
#SBATCH --output=logs/species_%A_%a.out
#SBATCH --error=logs/species_%A_%a.err

source ~/miniconda3/bin/activate cactus
export PATH=~/miniconda3/envs/cactus/bin:$PATH
export LD_LIBRARY_PATH=~/miniconda3/envs/cactus/lib:$LD_LIBRARY_PATH

# 获取当前物种
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sp240_nonhuman.txt)
if [ -z "$species" ]; then
    echo "错误: 无法获取物种名称"
    exit 1
fi

echo "处理物种: $species (第 $SLURM_ARRAY_TASK_ID/240 个)"
log_file="logs/${species}.log"
echo "[$(date)] 开始处理 $species" > $log_file

# 临时目录
tmp_dir="results/split_base/${species}_tmp"
mkdir -p $tmp_dir

# 步骤1: halLiftover映射每个分块
cat split_genome_parts.txt | while read part; do
    part_base="split_genome/$part"
    part_output="$tmp_dir/${part%.bed}.aligned.bed"
    
    if [ -f "$part_base" ]; then
        echo "[$(date)] 处理分块: $part" >> $log_file
        halLiftover \
            241-mammalian-2020v2.hal \
            Homo_sapiens \
            "$part_base" \
            "$species" \
            "$part_output" 2>> $log_file
    fi
done

# 合并结果
echo "[$(date)] 合并分块结果..." >> $log_file
cat $tmp_dir/*.aligned.bed > "results/align_base/${species}.aligned.bed" 2>> $log_file
rm -rf $tmp_dir

# 检查映射结果
if [ -s "results/align_base/${species}.aligned.bed" ]; then
    mapped_lines=$(wc -l < "results/align_base/${species}.aligned.bed")
    echo "[$(date)] 映射成功: $mapped_lines 行" >> $log_file
    
    # 步骤2: 生成bedgraph
    echo "[$(date)] 生成bedgraph..." >> $log_file
    awk '{print $4}' "results/align_base/${species}.aligned.bed" | \
        awk -F '_' '{print $1"\t"$2"\t"$3}' | \
        bedtools sort | \
        bedtools merge > "results/align_base/${species}.aligned.merged.bed"
    
    bedtools genomecov -i "results/align_base/${species}.aligned.merged.bed" \
        -g hg38.chrom.sizes -bga > "results/align_bedgraph/${species}.bedgraph"
    
    # 步骤3: 转换为bigWig
    echo "[$(date)] 转换为bigWig..." >> $log_file
    bedGraphToBigWig "results/align_bedgraph/${species}.bedgraph" \
        hg38.chrom.sizes \
        "results/align_bw/${species}.bw" 2>> $log_file
    
    # 步骤4: 计算iMotif比对率
    echo "[$(date)] 计算iMotif比对率..." >> $log_file
    bigWigAverageOverBed \
        "results/align_bw/${species}.bw" \
        imotif_clean_final.bed \
        "results/align_ratio/${species}.imotif.ar.txt" 2>> $log_file
    
    if [ -f "results/align_ratio/${species}.imotif.ar.txt" ]; then
        mapped=$(wc -l < "results/align_ratio/${species}.imotif.ar.txt")
        echo "[$(date)] $species: 成功处理 $mapped 个iMotif" >> $log_file
        echo "$species: $mapped" >> "results/mapping_summary.txt"
    fi
else
    echo "[$(date)] $species: 无映射结果" >> $log_file
    echo "$species: 无映射结果" >> "results/mapping_failed.txt"
fi

echo "[$(date)] $species 处理完成" >> $log_file
INNER_EOF

chmod +x process_species.sh

# ========== 修改点3: 调整并发数 ==========
total_species=$(wc -l < sp240_nonhuman.txt)
concurrent_jobs=54  # 根据集群空闲CPU调整
echo "提交 $total_species 个并行任务，每次并发 $concurrent_jobs 个..."
sbatch --array=1-$total_species%$concurrent_jobs process_species.sh
# ======================================

echo "========================================="
echo "主脚本提交完成！"
echo "使用以下命令监控进度："
echo "  squeue -u $USER"
echo "  tail -f logs/species_*.log"
echo "========================================="
