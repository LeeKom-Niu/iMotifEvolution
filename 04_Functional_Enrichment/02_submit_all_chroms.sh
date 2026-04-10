#!/bin/bash

# 批量提交所有染色体的 makeEdges 任务

PROJECT_DIR="/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
SCRIPT_DIR="${PROJECT_DIR}/scripts/makeEdges"
LOG_DIR="${PROJECT_DIR}/logs"

# 创建日志目录
mkdir -p ${LOG_DIR}

# 清理旧的统计文件
echo "chromosome,total_files,processed,failed,edges,nodes,timestamp" > ${PROJECT_DIR}/makeEdges_stats.csv

echo "=== 开始批量提交 makeEdges 任务 ==="
echo "时间: $(date)"
echo ""

# 提交所有染色体的任务
job_ids=""

for chrom in {1..22} X Y; do
    # 检查是否有对应的输入文件
    if [ -d "${PROJECT_DIR}/organized_files/organized_consistent/hsa${chrom}" ]; then
        file_count=$(ls "${PROJECT_DIR}/organized_files/organized_consistent/hsa${chrom}"/*.rmredundant.df 2>/dev/null | wc -l)
        
        if [ $file_count -gt 0 ]; then
            echo "提交染色体 hsa${chrom} (${file_count} 个文件)..."
            
            # 提交作业
            job_id=$(sbatch \
                --job-name="makeEdges_hsa${chrom}" \
                --nodes=1 \
                --ntasks=1 \
                --cpus-per-task=10 \
                --mem=40G \
                --time=100-00:00:00 \
                --partition=life-zhanghk \
                --output="${LOG_DIR}/makeEdges_hsa${chrom}_%j.out" \
                --error="${LOG_DIR}/makeEdges_hsa${chrom}_%j.err" \
                ${SCRIPT_DIR}/submitMakeEdges.sh ${chrom} | awk '{print $4}')
            
            if [ -n "$job_id" ]; then
                echo "  作业ID: ${job_id}"
                job_ids="${job_ids} ${job_id}"
            else
                echo "  提交失败"
            fi
        else
            echo "跳过染色体 hsa${chrom} (无文件)"
        fi
    else
        echo "跳过染色体 hsa${chrom} (目录不存在)"
    fi
done

echo ""
echo "=== 提交完成 ==="
echo "时间: $(date)"
echo "已提交作业ID: ${job_ids}"
echo ""

# 创建监控脚本
cat > ${PROJECT_DIR}/monitor_makeEdges.sh << MONITOR_EOF
#!/bin/bash
echo "=== makeEdges 作业监控 ==="
echo "时间: \$(date)"
echo ""
echo "正在运行的作业:"
squeue -u \$USER -o "%.18i %.9P %.35j %.8u %.8T %.10M %.9l %.6D %R" | grep makeEdges
echo ""
echo "已完成作业统计:"
if [ -f "${PROJECT_DIR}/makeEdges_stats.csv" ]; then
    tail -n +2 "${PROJECT_DIR}/makeEdges_stats.csv" | while IFS=',' read -r chrom total processed failed edges nodes timestamp; do
        echo "hsa\${chrom}: 处理 \${processed}/\${total}, 边: \${edges}, 节点: \${nodes}"
    done
fi
MONITOR_EOF

chmod +x ${PROJECT_DIR}/monitor_makeEdges.sh

echo "监控脚本已创建: ${PROJECT_DIR}/monitor_makeEdges.sh"
echo "运行: ./monitor_makeEdges.sh 查看进度"
