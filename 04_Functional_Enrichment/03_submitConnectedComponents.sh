
PROJECT_DIR="/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
SCRIPT_DIR="${PROJECT_DIR}/scripts/connectedComponents"
OUTPUT_DIR="${PROJECT_DIR}/output"
PYTHON_SCRIPT="${PROJECT_DIR}/makeConnectedComponents.py"
echo "=== 开始 Connected Components 分析 ==="
echo "时间: $(date)"
echo "项目目录: ${PROJECT_DIR}"
echo "输出目录: ${OUTPUT_DIR}"
echo "Python脚本: ${PYTHON_SCRIPT}"
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    echo "错误: 找不到 makeConnectedComponents.py 脚本"
    echo "尝试从原仓库复制..."
    cp "${PROJECT_DIR}/GreatApeT2T-G4s-main/src/mapseaAndPostWorkflow/makeConnectedComponents.py" "${PYTHON_SCRIPT}"
    
    if [ ! -f "${PYTHON_SCRIPT}" ]; then
        echo "错误: 无法复制 makeConnectedComponents.py"
        exit 1
    fi
fi
python3 -c "import networkx as nx; print('networkx 版本:', nx.__version__)" || {
    echo "错误: networkx 库未安装"
    echo "安装命令: pip install networkx pandas"
    exit 1
}
mkdir -p ../logs
echo "开始并行处理所有染色体..."
parallel -v -j $SLURM_NTASKS '
    chrom={}
    egs_file="'${OUTPUT_DIR}'/hsa${chrom}/hsa${chrom}.egs"
    graph_file="'${OUTPUT_DIR}'/hsa${chrom}/hsa${chrom}.graph"
    
    if [ -f "${egs_file}" ]; then
        edges=$(wc -l < "${egs_file}")
        if [ ${edges} -gt 0 ]; then
            echo "处理 hsa${chrom} (${edges} 条边)..."
            python3 "'${PYTHON_SCRIPT}'" "${egs_file}" > "${graph_file}"
            if [ $? -eq 0 ]; then
                components=$(wc -l < "${graph_file}")
                echo "  hsa${chrom}: 生成 ${components} 个连通组件"
            else
                echo "  hsa${chrom}: 处理失败"
            fi
        else
            echo "跳过 hsa${chrom} (无边数据)"
            touch "${graph_file}"
        fi
    else
        echo "跳过 hsa${chrom} (.egs 文件不存在)"
    fi
' ::: {1..22} X Y
echo ""
echo "=== Connected Components 处理完成 ==="
echo "时间: $(date)"
echo ""
echo "=== 结果统计 ==="
total_components=0
for chrom in {1..22} X Y; do
    graph_file="${OUTPUT_DIR}/hsa${chrom}/hsa${chrom}.graph"
    if [ -f "${graph_file}" ]; then
        components=$(wc -l < "${graph_file}" 2>/dev/null || echo "0")
        if [ $components -gt 0 ]; then
            echo "hsa${chrom}: ${components} 个连通组件"
            total_components=$((total_components + components))
        fi
    fi
done
echo ""
echo "总计: ${total_components} 个连通组件"
echo "输出文件位置: ${OUTPUT_DIR}/hsa*/hsa*.graph"
echo "chromosome,connected_components,timestamp" > "${PROJECT_DIR}/connectedComponents_stats.csv"
for chrom in {1..22} X Y; do
    graph_file="${OUTPUT_DIR}/hsa${chrom}/hsa${chrom}.graph"
    if [ -f "${graph_file}" ]; then
        components=$(wc -l < "${graph_file}" 2>/dev/null || echo "0")
        echo "${chrom},${components},$(date)" >> "${PROJECT_DIR}/connectedComponents_stats.csv"
    fi
done
