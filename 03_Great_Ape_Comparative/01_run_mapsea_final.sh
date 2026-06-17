#!/bin/bash
# MAP-SEA parallel execution script (output file names match example)

cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea/scripts

# Path settings
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results"
TEMP_BASE="${BASE_DIR}/temp"
MAF_DIR="${BASE_DIR}/maffiles"
MAPSEA_SCRIPT="${BASE_DIR}/mapsea-main/src/mapsea.py"
REFINER_SCRIPT="${BASE_DIR}/mapsea-main/src/refiner.py"

# Resource settings
CPUS_PER_TASK=20
MEM_PER_TASK=250G
TIME_PER_TASK=300-00:00:00
MAX_PARALLEL=1

# Load environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# Species name mapping
declare -A SPECIES_NAME_MAP=(
    ["human"]="Homo_sapiens"
    ["chimp"]="Pan_troglodytes"
    ["bonobo"]="Pan_paniscus"
    ["gorilla"]="Gorilla_gorilla"
    ["bornean"]="Pongo_pygmaeus"
    ["sumatran"]="Pongo_abelii"
)

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BASE_DIR}/scripts/logs"
mkdir -p "${TEMP_BASE}"

# Create task list
TASK_FILE="${OUTPUT_DIR}/tasks.txt"
echo "# Task list" > "$TASK_FILE"
echo "# Format: species1,species2,chr1,chr2,MAF file,output prefix" >> "$TASK_FILE"

TASK_COUNT=0
echo "Scanning MAF files..."

# Scan new naming format: chr{N}_human_vs_chr{N}_{species}.maf.gz
for MAF_FILE in "$MAF_DIR"/chr*_human_vs_*.maf.gz; do
    if [ -f "$MAF_FILE" ]; then
        filename=$(basename "$MAF_FILE")
        
        # New format: chr10_human_vs_chr8_bonobo.maf.gz
        if [[ $filename =~ ^chr([0-9XY]+)_human_vs_chr([0-9XY]+)_([a-z]+)\.maf\.gz$ ]]; then
            chr1="${BASH_REMATCH[1]}"
            chr2="${BASH_REMATCH[2]}"
            species2="${BASH_REMATCH[3]}"
            species1="human"

            # Output prefix: chr10_human_vs_chr8_bonobo (matches example)
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"
            
            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "Adding task: ${output_prefix}"

        # Old format: human_chr10_bonobo_chr8.maf.gz (compatibility)
        elif [[ $filename =~ ^human_chr([0-9XY]+)_([a-z]+)_chr([0-9XY]+)\.maf\.gz$ ]]; then
            chr1="${BASH_REMATCH[1]}"
            species2="${BASH_REMATCH[2]}"
            chr2="${BASH_REMATCH[3]}"
            species1="human"
            
            # Output prefix: chr10_human_vs_chr8_bonobo (converted to example format)
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"
            
            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "Adding task: ${output_prefix} (converted from old format)"
        fi
    fi
done

echo "任务数: $TASK_COUNT"

if [ $TASK_COUNT -eq 0 ]; then
    echo "Error: No tasks found"
    echo "Hint: Ensure MAF files are renamed to chr{N}_human_vs_chr{N}_{species}.maf.gz format"
    exit 1
fi

# Create task script
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

# Paths
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results"
TEMP_BASE="${BASE_DIR}/temp"
MAF_DIR="${BASE_DIR}/maffiles"
MAPSEA_SCRIPT="${BASE_DIR}/mapsea-main/src/mapsea.py"
REFINER_SCRIPT="${BASE_DIR}/mapsea-main/src/refiner.py"

# Load environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# Species mapping
declare -A SPECIES_NAME_MAP=(
    ["human"]="Homo_sapiens"
    ["chimp"]="Pan_troglodytes"
    ["bonobo"]="Pan_paniscus"
    ["gorilla"]="Gorilla_gorilla"
    ["bornean"]="Pongo_pygmaeus"
    ["sumatran"]="Pongo_abelii"
)

# Read task
TASK_FILE="${OUTPUT_DIR}/tasks.txt"
TASK_INFO=$(sed -n "$((TASK_ID+2))p" "$TASK_FILE")

if [ -z "$TASK_INFO" ]; then
    exit 0
fi

IFS=',' read -r species1 species2 chr1 chr2 maf_filename output_prefix <<< "$TASK_INFO"

# Check if already processed (using new output prefix)
OUTPUT_DAT="${OUTPUT_DIR}/${output_prefix}.dat"
OUTPUT_DF="${OUTPUT_DIR}/${output_prefix}.rmredundant.df"

if [ -f "${OUTPUT_DF}" ]; then
    echo "Task already completed: ${output_prefix}"
    exit 0
fi

# Check files
MAF_FILE="${MAF_DIR}/${maf_filename}"
if [ ! -f "$MAF_FILE" ]; then
    echo "Error: MAF file not found - $MAF_FILE"
    exit 1
fi

species1_sci="${SPECIES_NAME_MAP[$species1]}"
species2_sci="${SPECIES_NAME_MAP[$species2]}"
BED_FILE1="${BED_FOLDER}/${species1_sci}/chr${chr1}.bed"
BED_FILE2="${BED_FOLDER}/${species2_sci}/chr${chr2}.bed"

if [ ! -f "$BED_FILE1" ] || [ ! -f "$BED_FILE2" ]; then
    echo "Error: BED file not found"
    echo "  $BED_FILE1: $(if [ -f "$BED_FILE1" ]; then echo 'Exists'; else echo 'Not found'; fi)"
    echo "  $BED_FILE2: $(if [ -f "$BED_FILE2" ]; then echo 'Exists'; else echo 'Not found'; fi)"
    exit 1
fi

# Temporary directory
TEMP_DIR="${TEMP_BASE}/task_final_${TASK_ID}"
mkdir -p "$TEMP_DIR"

echo "Processing task $TASK_ID: ${output_prefix}"
echo "MAF file: $maf_filename"
echo "Output files: ${output_prefix}.dat / ${output_prefix}.rmredundant.df"
echo "BED files:"
echo "  $species1 ($species1_sci): $BED_FILE1"
echo "  $species2 ($species2_sci): $BED_FILE2"

# Run mapsea
python "$MAPSEA_SCRIPT" \
    -m "$MAF_FILE" \
    -b "$BED_FOLDER" \
    -o "${OUTPUT_DAT}" \
    -t "$TEMP_DIR" \
    -r 1.0 \
    -d "$SPECIES_DICT" \
    -c 20 > "${TEMP_DIR}/mapsea.log" 2>&1

# Check if .dat file was generated successfully
if [ ! -f "${OUTPUT_DAT}" ]; then
    echo "mapsea failed: .dat file not generated"
    echo "Last 50 lines of log:"
    tail -50 "${TEMP_DIR}/mapsea.log"
    exit 1
fi

# Check if .dat file has metadata header
if ! head -1 "${OUTPUT_DAT}" | grep -q "^## {METADATA}"; then
    echo "Warning: Generated .dat file has no metadata header, fixing..."
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
    echo "Metadata header fix completed"
fi

echo "mapsea completed: ${OUTPUT_DAT}"

# Run refiner
python "$REFINER_SCRIPT" \
    -d "${OUTPUT_DAT}" \
    -f 3 \
    -o "${OUTPUT_DF}" \
    -m \
    -c 20 > "${TEMP_DIR}/refiner.log" 2>&1

if [ $? -ne 0 ] || [ ! -f "${OUTPUT_DF}" ]; then
    echo "refiner failed"
    echo "Last 30 lines of log:"
    tail -30 "${TEMP_DIR}/refiner.log"
    exit 1
fi

echo "refiner completed: ${OUTPUT_DF}"

# Clean up temporary directory
rm -rf "$TEMP_DIR"

echo "Task $TASK_ID completed: ${output_prefix}"
exit 0
TASK_EOF

chmod +x "$TASK_SCRIPT"

# Submit job
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

echo "Submit command: $JOB_CMD"
JOB_ID=$($JOB_CMD | awk '{print $4}')

if [ -n "$JOB_ID" ]; then
    echo "Job ID: $JOB_ID"
    echo "Task count: $TASK_COUNT"
    echo "Max parallel: $MAX_PARALLEL"
    echo "Resources: ${CPUS_PER_TASK} cores, ${MEM_PER_TASK} memory, ${TIME_PER_TASK} time"
    echo "Output file format: chr{N}_{species1}_vs_chr{N}_{species2}.{dat|rmredundant.df}"
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
echo "Check script: ${BASE_DIR}/scripts/check_final.sh"
