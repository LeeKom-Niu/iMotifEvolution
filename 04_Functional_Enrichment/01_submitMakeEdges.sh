
CHROM=$1
PROJECT_DIR="/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
ORGANIZED_DIR="${PROJECT_DIR}/organized_files/organized_consistent"
OUTPUT_DIR="${PROJECT_DIR}/output"
SCRIPT_DIR="${PROJECT_DIR}"
mkdir -p logs
echo "=== 开始处理染色体 hsa${CHROM} ==="
echo "时间: $(date)"
echo "项目目录: ${PROJECT_DIR}"
echo "输入目录: ${ORGANIZED_DIR}/hsa${CHROM}"
echo "输出目录: ${OUTPUT_DIR}/hsa${CHROM}"
if [ ! -d "${ORGANIZED_DIR}/hsa${CHROM}" ]; then
    echo "错误: 输入目录不存在: ${ORGANIZED_DIR}/hsa${CHROM}"
    exit 1
fi
if [ ! -f "${SCRIPT_DIR}/makeEdges.py" ]; then
    echo "错误: 找不到 makeEdges.py 脚本: ${SCRIPT_DIR}/makeEdges.py"
    exit 1
fi
mkdir -p "${OUTPUT_DIR}/hsa${CHROM}"
echo "清理旧输出文件..."
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" ]; then
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs"
fi
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"
fi
files=($(ls "${ORGANIZED_DIR}/hsa${CHROM}"/*.rmredundant.df 2>/dev/null))
file_count=${
if [ $file_count -eq 0 ]; then
    echo "警告: 没有找到 .rmredundant.df 文件"
    exit 0
fi
echo "找到 ${file_count} 个文件需要处理"
echo "开始并行处理文件..."
processed=0
failed=0
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        echo "  处理: ${filename}"
        
        python3 "${SCRIPT_DIR}/makeEdges.py" \
            "$file" \
            "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" \
            >> "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" 2>&1
        
        if [ $? -eq 0 ]; then
            ((processed++))
        else
            echo "    ✗ 处理失败: ${filename}"
            ((failed++))
        fi
    fi
done
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    echo "排序去重 alignedUnique.nds 文件..."
    mv "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" \
       "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak"
    
    sort -u "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak" \
        > "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"
    
    rm -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds.bak"
fi
edges=0
nodes=0
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs" ]; then
    edges=$(wc -l < "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs")
fi
if [ -f "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds" ]; then
    nodes=$(wc -l < "${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds")
fi
echo ""
echo "=== 处理完成 ==="
echo "时间: $(date)"
echo "处理文件: ${processed}/${file_count} (失败: ${failed})"
echo "生成边: ${edges} 条"
echo "生成对齐唯一节点: ${nodes} 个"
echo "输出文件:"
echo "  ${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}.egs"
echo "  ${OUTPUT_DIR}/hsa${CHROM}/hsa${CHROM}_alignedUnique.nds"
echo "${CHROM},${file_count},${processed},${failed},${edges},${nodes},$(date)" \
    >> "${PROJECT_DIR}/makeEdges_stats_test.csv"
echo "=== hsa${CHROM} 处理完成 ==="
