#!/bin/bash
# MAP-SEA并行运行脚本（输出文件名与示例一致）

cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea/scripts

# 路径设置
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results"
TEMP_BASE="${BASE_DIR}/temp"
MAF_DIR="${BASE_DIR}/maffiles"
MAPSEA_SCRIPT="${BASE_DIR}/mapsea-main/src/mapsea.py"
REFINER_SCRIPT="${BASE_DIR}/mapsea-main/src/refiner.py"

# 资源设置
CPUS_PER_TASK=20
MEM_PER_TASK=250G
TIME_PER_TASK=300-00:00:00
MAX_PARALLEL=1

# 加载环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# 物种名称映射
declare -A SPECIES_NAME_MAP=(
    ["human"]="Homo_sapiens"
    ["chimp"]="Pan_troglodytes"
    ["bonobo"]="Pan_paniscus"
    ["gorilla"]="Gorilla_gorilla"
    ["bornean"]="Pongo_pygmaeus"
    ["sumatran"]="Pongo_abelii"
)

# 创建输出目录
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BASE_DIR}/scripts/logs"
mkdir -p "${TEMP_BASE}"

# 创建任务列表
TASK_FILE="${OUTPUT_DIR}/tasks.txt"
echo "# 任务列表" > "$TASK_FILE"
echo "# 格式: 物种1,物种2,chr1,chr2,MAF文件,输出前缀" >> "$TASK_FILE"

TASK_COUNT=0
echo "扫描MAF文件..."

# 扫描新的命名格式：chr{数字}_human_vs_chr{数字}_{物种}.maf.gz
for MAF_FILE in "$MAF_DIR"/chr*_human_vs_*.maf.gz; do
    if [ -f "$MAF_FILE" ]; then
        filename=$(basename "$MAF_FILE")
        
        # 新格式: chr10_human_vs_chr8_bonobo.maf.gz
        if [[ $filename =~ ^chr([0-9XY]+)_human_vs_chr([0-9XY]+)_([a-z]+)\.maf\.gz$ ]]; then
            chr1="${BASH_REMATCH[1]}"
            chr2="${BASH_REMATCH[2]}"
            species2="${BASH_REMATCH[3]}"
            species1="human"
            
            # 输出前缀：chr10_human_vs_chr8_bonobo（与示例一致）
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"
            
            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "添加任务: ${output_prefix}"
            
        # 旧格式: human_chr10_bonobo_chr8.maf.gz（兼容性）
        elif [[ $filename =~ ^human_chr([0-9XY]+)_([a-z]+)_chr([0-9XY]+)\.maf\.gz$ ]]; then
            chr1="${BASH_REMATCH[1]}"
            species2="${BASH_REMATCH[2]}"
            chr2="${BASH_REMATCH[3]}"
            species1="human"
            
            # 输出前缀：chr10_human_vs_chr8_bonobo（转换为示例格式）
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"
            
            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "添加任务: ${output_prefix} (从旧格式转换)"
        fi
    fi
done

echo "任务数: $TASK_COUNT"

if [ $TASK_COUNT -eq 0 ]; then
    echo "错误: 未找到任务"
    echo "提示: 确保MAF文件已重命名为 chr{数字}_human_vs_chr{数字}_{物种}.maf.gz 格式"
    exit 1
fi

# 创建任务脚本
TASK_SCRIPT="${BASE_DIR}/scripts/task_final.sh"
cat > "$TASK_SCRIPT" << 'TASK_EOF'
#!/bin/bash

#SBATCH --job-name=mapsea_final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=logs/mapsea_final_%A_%a.out
#SBATCH --error=logs/mapsea_final_%A_%a.err

TASK_ID=${SLURM_ARRAY_TASK_ID}

cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea

# 路径
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results"
TEMP_BASE="${BASE_DIR}/temp"
MAF_DIR="${BASE_DIR}/maffiles"
MAPSEA_SCRIPT="${BASE_DIR}/mapsea-main/src/mapsea.py"
REFINER_SCRIPT="${BASE_DIR}/mapsea-main/src/refiner.py"

# 加载环境
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# 物种映射
declare -A SPECIES_NAME_MAP=(
    ["human"]="Homo_sapiens"
    ["chimp"]="Pan_troglodytes"
    ["bonobo"]="Pan_paniscus"
    ["gorilla"]="Gorilla_gorilla"
    ["bornean"]="Pongo_pygmaeus"
    ["sumatran"]="Pongo_abelii"
)

# 读取任务
TASK_FILE="${OUTPUT_DIR}/tasks.txt"
TASK_INFO=$(sed -n "$((TASK_ID+2))p" "$TASK_FILE")

if [ -z "$TASK_INFO" ]; then
    exit 0
fi

IFS=',' read -r species1 species2 chr1 chr2 maf_filename output_prefix <<< "$TASK_INFO"

# 检查是否已处理（使用新的输出前缀）
OUTPUT_DAT="${OUTPUT_DIR}/${output_prefix}.dat"
OUTPUT_DF="${OUTPUT_DIR}/${output_prefix}.rmredundant.df"

if [ -f "${OUTPUT_DF}" ]; then
    echo "任务已完成: ${output_prefix}"
    exit 0
fi

# 检查文件
MAF_FILE="${MAF_DIR}/${maf_filename}"
if [ ! -f "$MAF_FILE" ]; then
    echo "错误: MAF文件不存在 - $MAF_FILE"
    exit 1
fi

species1_sci="${SPECIES_NAME_MAP[$species1]}"
species2_sci="${SPECIES_NAME_MAP[$species2]}"
BED_FILE1="${BED_FOLDER}/${species1_sci}/chr${chr1}.bed"
BED_FILE2="${BED_FOLDER}/${species2_sci}/chr${chr2}.bed"

if [ ! -f "$BED_FILE1" ] || [ ! -f "$BED_FILE2" ]; then
    echo "错误: BED文件不存在"
    echo "  $BED_FILE1: $(if [ -f "$BED_FILE1" ]; then echo '存在'; else echo '不存在'; fi)"
    echo "  $BED_FILE2: $(if [ -f "$BED_FILE2" ]; then echo '存在'; else echo '不存在'; fi)"
    exit 1
fi

# 临时目录
TEMP_DIR="${TEMP_BASE}/task_final_${TASK_ID}"
mkdir -p "$TEMP_DIR"

echo "处理任务 $TASK_ID: ${output_prefix}"
echo "MAF文件: $maf_filename"
echo "输出文件: ${output_prefix}.dat / ${output_prefix}.rmredundant.df"
echo "BED文件:"
echo "  $species1 ($species1_sci): $BED_FILE1"
echo "  $species2 ($species2_sci): $BED_FILE2"

# 运行mapsea
python "$MAPSEA_SCRIPT" \
    -m "$MAF_FILE" \
    -b "$BED_FOLDER" \
    -o "${OUTPUT_DAT}" \
    -t "$TEMP_DIR" \
    -r 1.0 \
    -d "$SPECIES_DICT" \
    -c 20 > "${TEMP_DIR}/mapsea.log" 2>&1

# 检查是否成功生成.dat文件
if [ ! -f "${OUTPUT_DAT}" ]; then
    echo "mapsea失败: 未生成.dat文件"
    echo "日志最后50行:"
    tail -50 "${TEMP_DIR}/mapsea.log"
    exit 1
fi

# 检查.dat文件是否有元数据头部
if ! head -1 "${OUTPUT_DAT}" | grep -q "^## {METADATA}"; then
    echo "警告: 生成的.dat文件没有元数据头部，修复中..."
    tmp_file="${OUTPUT_DAT}.tmp"
    
    cat > "$tmp_file" << METADATA_FIX
## {METADATA}
## INPUT FILE: $MAF_FILE
## OUTPUT FILE: ${OUTPUT_DAT}
## HSA MAP: None
## SPECIES DICTIONARY: $SPECIES_DICT
## INTERSECTION RATIO (f): 1.00
## 
## {TYPE}
## Absent in .BED file: 
##   FGAP: No sequence present
##   FNotA: Sequence is non-annotated
##   GAP.xx: Sequence has xx% gaps
## Present in .BED file: 
##   partMAF: Sequence partly in alignment
##   fullMAF: Sequence fully in alignment
## 
## {SPECIESID}
## 1: Homo_sapiens
## 2: Pan_troglodytes
## 3: Pan_paniscus
## 4: Gorilla_gorilla
## 5: Pongo_pygmaeus
## 6: Pongo_abelii
## 
## {STRUCTURE}
## #BLOCK ID
## NUMBER
## SPECIESID@CHR:STARTQUERY_LENGTHTYPESCORESTRANDALIGNMENT_LENGTHSEQUENCE

METADATA_FIX
    
    cat "${OUTPUT_DAT}" >> "$tmp_file"
    mv "$tmp_file" "${OUTPUT_DAT}"
    echo "元数据头部修复完成"
fi

echo "mapsea完成: ${OUTPUT_DAT}"

# 运行refiner
python "$REFINER_SCRIPT" \
    -d "${OUTPUT_DAT}" \
    -f 3 \
    -o "${OUTPUT_DF}" \
    -m \
    -c 20 > "${TEMP_DIR}/refiner.log" 2>&1

if [ $? -ne 0 ] || [ ! -f "${OUTPUT_DF}" ]; then
    echo "refiner失败"
    echo "日志最后30行:"
    tail -30 "${TEMP_DIR}/refiner.log"
    exit 1
fi

echo "refiner完成: ${OUTPUT_DF}"

# 清理临时目录
rm -rf "$TEMP_DIR"

echo "任务 $TASK_ID 完成: ${output_prefix}"
exit 0
TASK_EOF

chmod +x "$TASK_SCRIPT"

# 提交作业
JOB_CMD="sbatch"
JOB_CMD="$JOB_CMD --array=1-${TASK_COUNT}"
JOB_CMD="$JOB_CMD --cpus-per-task=${CPUS_PER_TASK}"
JOB_CMD="$JOB_CMD --mem=${MEM_PER_TASK}"
JOB_CMD="$JOB_CMD --time=${TIME_PER_TASK}"
JOB_CMD="$JOB_CMD --partition=life-zhanghk"
JOB_CMD="$JOB_CMD --output=${BASE_DIR}/scripts/logs/mapsea_final_%A_%a.out"
JOB_CMD="$JOB_CMD --error=${BASE_DIR}/scripts/logs/mapsea_final_%A_%a.err"

if [ $TASK_COUNT -gt $MAX_PARALLEL ]; then
    JOB_CMD="$JOB_CMD --array=1-${TASK_COUNT}%${MAX_PARALLEL}"
fi

JOB_CMD="$JOB_CMD $TASK_SCRIPT"

echo "提交命令: $JOB_CMD"
JOB_ID=$($JOB_CMD | awk '{print $4}')

if [ -n "$JOB_ID" ]; then
    echo "作业ID: $JOB_ID"
    echo "任务数: $TASK_COUNT"
    echo "最大并行: $MAX_PARALLEL"
    echo "资源: ${CPUS_PER_TASK}核心, ${MEM_PER_TASK}内存, ${TIME_PER_TASK}时间"
    echo "输出文件格式: chr{数字}_{物种1}_vs_chr{数字}_{物种2}.{dat|rmredundant.df}"
fi

# 创建检查脚本
cat > "${BASE_DIR}/scripts/check_final.sh" << 'CHECK_EOF'
#!/bin/bash
cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea
COMPLETED=$(find results -name "*.rmredundant.df" 2>/dev/null | wc -l)
TOTAL=$(cat results/tasks.txt 2>/dev/null | wc -l)
if [ $TOTAL -ge 2 ]; then
    TOTAL=$((TOTAL - 2))
    echo "完成: $COMPLETED/$TOTAL"
    echo -e "\n最近完成的任务:"
    find results -name "*.rmredundant.df" -type f -exec ls -lt {} \; 2>/dev/null | head -5 | awk '{print $6" "$7" "$8" "$9}'
else
    echo "完成: $COMPLETED/0 (任务文件不存在)"
fi
CHECK_EOF

chmod +x "${BASE_DIR}/scripts/check_final.sh"
echo "检查脚本: ${BASE_DIR}/scripts/check_final.sh"
