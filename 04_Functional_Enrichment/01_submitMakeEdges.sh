#!/bin/bash

#SBATCH --job-name=makeEdges
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=120G
#SBATCH --time=100-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=logs/makeEdges_%j.out
#SBATCH --error=logs/makeEdges_%j.err

# Parameters
# $1 = 染色体编号 (如: 1, 2, ..., 22, X, Y)

CHROM=$1
PROJECT_DIR="/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
ORGANIZED_DIR="${PROJECT_DIR}/organized_files/organized_consistent"
OUTPUT_DIR="${PROJECT_DIR}/output"
SCRIPT_DIR="${PROJECT_DIR}"

# Create log directory
mkdir -p logs

echo "=== Starting processing chromosome hsa${CHROM} ==="
echo "Time: $(date)"
echo "Project directory: ${PROJECT_DIR}"
echo "Input directory: ${ORGANIZED_DIR}/hsa${CHROM}"
echo "Output directory: ${OUTPUT_DIR}/hsa${CHROM}"

# 检查输入目录
if [ ! -d "${ORGANIZED_DIR}/hsa${CHROM}" ]; then
    echo "错误: 输入目录不存在: ${ORGANIZED_DIR}/hsa${CHROM}"
    exit 1
fi

# 检查 makeEdges.py 脚本
if [ ! -f "${SCRIPT_DIR}/makeEdges.py" ]; then
    echo "Error: Cannot find makeEdges.py script: ${SCRIPT_DIR}/makeEdges.py"
    exit 1
fi

# 创建输出目录
mkdir -p "${OUTPUT_DIR}/hsa${CHROM}"

echo "Cleaning old output files..."
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" ]; then
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs"
fi

if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"
fi

# 获取文件列表
files=($(ls "${ORGANIZED_DIR}/hsa${CHROM}"/*.rmredundant.df 2>/dev/null))
file_count=${#files[@]}

if [ $file_count -eq 0 ]; then
    echo "Warning: No .rmredundant.df files found"
    exit 0
fi

echo "Found ${file_count} files to process"

# Process files
echo "Starting file processing..."
processed=0
failed=0

for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        echo "  Processing: ${filename}"
        
        # 使用 Python 处理文件
        python3 "${SCRIPT_DIR}/makeEdges.py" \
            "$file" \
            "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" \
            >> "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" 2>&1
        
        if [ $? -eq 0 ]; then
            ((processed++))
        else
            echo "    Failed: ${filename}"
            ((failed++))
        fi
    fi
done

# 排序去重 alignedUnique.nds 文件
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    echo "Sorting and deduplicating alignedUnique.nds file..."
    mv "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" \
       "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak"
    
    sort -u "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak" \
        > "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"
    
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak"
fi

# 统计结果
edges=0
nodes=0

if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" ]; then
    edges=$(wc -l < "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs")
fi

if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    nodes=$(wc -l < "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds")
fi

echo ""
echo "=== Processing Complete ==="
echo "Time: $(date)"
echo "Processed files: ${processed}/${file_count} (failed: ${failed})"
echo "Edges generated: ${edges}"
echo "Alignment-unique nodes: ${nodes}"
echo "输出文件:"
echo "  ${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs"
echo "  ${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"

# 保存统计信息
echo "${CHROM},${file_count},${processed},${failed},${edges},${nodes},$(date)" \
    >> "${PROJECT_DIR}/makeEdges_stats_test.csv"

echo "=== hsa${CHROM} Processing Complete ==="
