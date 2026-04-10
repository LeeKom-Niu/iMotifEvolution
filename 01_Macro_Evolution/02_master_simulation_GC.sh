#!/bin/bash
# master_simulation_GC.sh - GC含量背景模拟实验主控脚本

BASE_DIR="/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"

echo "=== GC含量背景模拟实验 ==="
echo "开始时间: $(date)"
echo "工作目录: $BASE_DIR"
echo "当前节点: $(hostname)"

# 直接设置参数
GENOME_SIZE=10000000  # 10Mb
REPLICATES=30
GC_LEVELS=("0.10" "0.20" "0.30" "0.40" "0.50" "0.60" "0.70" "0.80" "0.90")

echo "基因组大小: $GENOME_SIZE bp"
echo "重复次数: $REPLICATES"
echo "GC含量梯度数量: ${#GC_LEVELS[@]}"
echo "GC含量梯度: ${GC_LEVELS[@]}"

# 创建目录结构
mkdir -p $BASE_DIR/scripts
mkdir -p $BASE_DIR/sequences
mkdir -p $BASE_DIR/results
mkdir -p $BASE_DIR/logs
mkdir -p $BASE_DIR/figures

# 为每个GC含量梯度提交作业
for i in "${!GC_LEVELS[@]}"; do
    gc_level="${GC_LEVELS[$i]}"
    gc_percent=$(python3 -c "print(int(float('$gc_level') * 100))")
    
    echo "提交GC含量 ${gc_level} (${gc_percent}%) 到SLURM队列..."
    
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=im_GC${gc_percent}
#SBATCH --cpus-per-task=60
#SBATCH --mem=250G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=$BASE_DIR/logs/GC${gc_percent}_%j.out
#SBATCH --error=$BASE_DIR/logs/GC${gc_percent}_%j.err

source /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh
conda activate im-seeker

cd $BASE_DIR
echo "=== 开始处理GC含量 ${gc_level} ==="
echo "运行节点: \$(hostname)"
echo "开始时间: \$(date)"
echo "重复次数: $REPLICATES"
echo "基因组大小: $GENOME_SIZE bp"

python3 $BASE_DIR/scripts/parallel_simulation_GC.py $gc_level $REPLICATES $GENOME_SIZE

echo "GC含量 ${gc_level} 完成于 \$(date)"
EOF

done

echo "所有梯度作业已提交到SLURM队列！"
echo "使用 'squeue -u \$USER' 查看作业状态"
echo "结束时间: $(date)"
