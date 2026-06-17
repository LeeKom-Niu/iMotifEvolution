#!/bin/bash
# Process ape_vs_ape directory parallel script (output file names match example)

cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea_ape/scripts

# Path settings
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea_ape"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results_ape_vs_ape"
TEMP_BASE="${BASE_DIR}/temp_ape"
MAF_DIR="${BASE_DIR}/maffiles/ape_vs_ape"
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

# Create output directories first
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${BASE_DIR}/scripts/logs"
mkdir -p "${TEMP_BASE}"

# Species name mapping
declare -A SPECIES_NAME_MAP=(
    ["human"]="Homo_sapiens"
    ["chimp"]="Pan_troglodytes"
    ["bonobo"]="Pan_paniscus"
    ["gorilla"]="Gorilla_gorilla"
    ["bornean"]="Pongo_pygmaeus"
    ["sumatran"]="Pongo_abelii"
)

# Create task list
TASK_FILE="${OUTPUT_DIR}/tasks.txt"
echo "# ape_vs_ape task list" > "$TASK_FILE"
echo "# Format: species1,species2,chr1,chr2,MAF file,output prefix" >> "$TASK_FILE"

TASK_COUNT=0
echo "Scanning APE comparison MAF files..."

# Scan new naming format: chr{N}_{species1}_vs_chr{N}_{species2}.maf.gz
for MAF_FILE in "$MAF_DIR"/chr*_*_vs_*.maf.gz; do
    if [ -f "$MAF_FILE" ]; then
        filename=$(basename "$MAF_FILE")

        # New format: chr10_bonobo_vs_chr10_bornean.maf.gz
        if [[ $filename =~ ^chr([0-9XY]+)_([a-z]+)_vs_chr([0-9XY]+)_([a-z]+)\.maf\.gz$ ]]; then
            chr1="${BASH_REMATCH[1]}"
            species1="${BASH_REMATCH[2]}"
            chr2="${BASH_REMATCH[3]}"
            species2="${BASH_REMATCH[4]}"

            # Output prefix: chr10_bonobo_vs_chr10_bornean (matches example)
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"

            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "Adding task: ${output_prefix}"

        # Old format: bonobo_chr10_bornean_chr10.maf.gz (compatibility)
        elif [[ $filename =~ ^([a-z]+)_chr([0-9XY]+)_([a-z]+)_chr([0-9XY]+)\.maf\.gz$ ]]; then
            species1="${BASH_REMATCH[1]}"
            chr1="${BASH_REMATCH[2]}"
            species2="${BASH_REMATCH[3]}"
            chr2="${BASH_REMATCH[4]}"

            # Output prefix: chr10_bonobo_vs_chr10_bornean (converted to example format)
            output_prefix="chr${chr1}_${species1}_vs_chr${chr2}_${species2}"

            TASK_COUNT=$((TASK_COUNT + 1))
            echo "${species1},${species2},${chr1},${chr2},${filename},${output_prefix}" >> "$TASK_FILE"
            echo "Adding task: ${output_prefix} (converted from old format)"
        fi
    fi
done

echo "ape_vs_ape task count: $TASK_COUNT"

if [ $TASK_COUNT -eq 0 ]; then
    echo "Error: No tasks found"
    echo "Hint: Ensure APE MAF files are renamed to chr{N}_{species1}_vs_chr{N}_{species2}.maf.gz format"
    exit 1
fi

# Create task script
TASK_SCRIPT="${BASE_DIR}/scripts/task_ape_final.sh"
cat > "$TASK_SCRIPT" << 'TASK_EOF'
#!/bin/bash

#SBATCH --job-name=mapsea_ape_final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=250G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=logs/mapsea_ape_final_%A_%a.out
#SBATCH --error=logs/mapsea_ape_final_%A_%a.err

TASK_ID=${SLURM_ARRAY_TASK_ID}

cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea_ape

# Paths
BASE_DIR="/datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea_ape"
BED_FOLDER="${BASE_DIR}/bedfiles"
SPECIES_DICT="${BASE_DIR}/species_dict.json"
OUTPUT_DIR="${BASE_DIR}/results_ape_vs_ape"
TEMP_BASE="${BASE_DIR}/temp_ape"
MAF_DIR="${BASE_DIR}/maffiles/ape_vs_ape"
MAPSEA_SCRIPT="${BASE_DIR}/mapsea-main/src/mapsea.py"
REFINER_SCRIPT="${BASE_DIR}/mapsea-main/src/refiner.py"

# Load environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# Species name to scientific name mapping (for BED files)
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

# Get scientific names for BED files
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
TEMP_DIR="${TEMP_BASE}/task_ape_final_${TASK_ID}"
mkdir -p "$TEMP_DIR"

echo "Processing APE comparison task $TASK_ID: ${output_prefix}"
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

echo "APE task $TASK_ID completed: ${output_prefix}"
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
JOB_CMD="$JOB_CMD --output=${BASE_DIR}/scripts/logs/mapsea_ape_final_%A_%a.out"
JOB_CMD="$JOB_CMD --error=${BASE_DIR}/scripts/logs/mapsea_ape_final_%A_%a.err"

if [ $TASK_COUNT -gt $MAX_PARALLEL ]; then
    JOB_CMD="$JOB_CMD --array=1-${TASK_COUNT}%${MAX_PARALLEL}"
fi

JOB_CMD="$JOB_CMD $TASK_SCRIPT"

echo "Submit command: $JOB_CMD"
JOB_ID=$($JOB_CMD | awk '{print $4}')

if [ -n "$JOB_ID" ]; then
    echo "Job ID: $JOB_ID"
    echo "ape_vs_ape tasks: $TASK_COUNT"
    echo "Max parallel: $MAX_PARALLEL"
    echo "Resources: ${CPUS_PER_TASK} cores, ${MEM_PER_TASK} memory, ${TIME_PER_TASK} time"
    echo "Output directory: $OUTPUT_DIR"
    echo "Output file format: chr{N}_{species1}_vs_chr{N}_{species2}.{dat|rmredundant.df}"
fi

# Create check script
cat > "${BASE_DIR}/scripts/check_ape_final.sh" << 'CHECK_EOF'
#!/bin/bash
cd /datapool/home/2023200496/niulk/my_project/18_mammalian/03_map_sea_ape
COMPLETED=$(find results_ape_vs_ape -name "*.rmredundant.df" 2>/dev/null | wc -l)
TOTAL=$(cat results_ape_vs_ape/tasks.txt 2>/dev/null | wc -l)
if [ $TOTAL -ge 2 ]; then
    TOTAL=$((TOTAL - 2))
    echo "ape_vs_ape completed: $COMPLETED/$TOTAL"
    echo -e "\nRecently completed tasks:"
    find results_ape_vs_ape -name "*.rmredundant.df" -type f -exec ls -lt {} \; 2>/dev/null | head -5 | awk '{print $6" "$7" "$8" "$9}'
else
    echo "ape_vs_ape completed: $COMPLETED/0 (task file not found)"
fi
CHECK_EOF

chmod +x "${BASE_DIR}/scripts/check_ape_final.sh"
echo "Check script: ${BASE_DIR}/scripts/check_ape_final.sh"
