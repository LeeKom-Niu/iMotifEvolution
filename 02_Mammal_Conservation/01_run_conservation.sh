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
echo "Start time: $(date)"
echo "工作目录: $(pwd)"
echo "========================================="

# Activate environment
source ~/miniconda3/bin/activate cactus

# Create directory structure
mkdir -p logs
mkdir -p results/align_base
mkdir -p results/align_bedgraph
mkdir -p results/align_bw
mkdir -p results/align_ratio
mkdir -p results/split_base
mkdir -p split_genome
mkdir -p matrix
mkdir -p plots

echo "Directory created"

# Step 1: Generate human genome single-base window file
echo "Generating human genome single-base window file..."
if [ ! -f human_genome_windows.4col.bed ]; then
    bedtools makewindows -g hg38.chrom.sizes -w 1 > human_genome_windows.bed
    awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3}' human_genome_windows.bed > human_genome_windows.4col.bed
    echo "Completed: $(wc -l < human_genome_windows.4col.bed) windows"
fi

# ========== Modification 1: Use author's chunk size and 4-digit suffix ==========
echo "Splitting genome window file..."
total_lines=$(wc -l < human_genome_windows.4col.bed)
lines_per_file=25750000  # Author uses 25.75M lines per chunk
echo "Lines per file: $lines_per_file"

# Remove old chunks (if they exist)
rm -rf split_genome/*

# Use 4-digit numeric suffix (0000-9999)
split -l $lines_per_file -d -a 4 human_genome_windows.4col.bed split_genome/genome_part.

# Generate chunk list
ls split_genome/ > split_genome_parts.txt
part_count=$(wc -l < split_genome_parts.txt)
echo "Completed: $part_count chunks"
# ======================================================

# Step 3: Create species processing script
cat > process_species.sh << 'INNER_EOF'
#!/bin/bash
#SBATCH --job-name=imotif_species
#SBATCH --nodes=1
# ========== Modification 2: Adjust resources per species task ==========
#SBATCH --cpus-per-task=10      # 4 CPUs per species (adjust based on cluster)
#SBATCH --mem=40G               # 16G memory per species
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
# ==================================================
#SBATCH --output=logs/species_%A_%a.out
#SBATCH --error=logs/species_%A_%a.err

source ~/miniconda3/bin/activate cactus
export PATH=~/miniconda3/envs/cactus/bin:$PATH
export LD_LIBRARY_PATH=~/miniconda3/envs/cactus/lib:$LD_LIBRARY_PATH

# Get current species
species=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sp240_nonhuman.txt)
if [ -z "$species" ]; then
    echo "Error: Cannot get species name"
    exit 1
fi

echo "Processing species: $species (job $SLURM_ARRAY_TASK_ID/240)"
log_file="logs/${species}.log"
echo "[$(date)] Starting processing $species" > $log_file

# Temporary directory
tmp_dir="results/split_base/${species}_tmp"
mkdir -p $tmp_dir

# Step 1: halLiftover mapping for each chunk
cat split_genome_parts.txt | while read part; do
    part_base="split_genome/$part"
    part_output="$tmp_dir/${part%.bed}.aligned.bed"
    
    if [ -f "$part_base" ]; then
        echo "[$(date)] Processing chunk: $part" >> $log_file
        halLiftover \
            241-mammalian-2020v2.hal \
            Homo_sapiens \
            "$part_base" \
            "$species" \
            "$part_output" 2>> $log_file
    fi
done

# Merge results
echo "[$(date)] Merging chunk results..." >> $log_file
cat $tmp_dir/*.aligned.bed > "results/align_base/${species}.aligned.bed" 2>> $log_file
rm -rf $tmp_dir

# Check mapping results
if [ -s "results/align_base/${species}.aligned.bed" ]; then
    mapped_lines=$(wc -l < "results/align_base/${species}.aligned.bed")
    echo "[$(date)] Mapping successful: $mapped_lines lines" >> $log_file
    
    # Step 2: Generate bedgraph
    echo "[$(date)] Generating bedgraph..." >> $log_file
    awk '{print $4}' "results/align_base/${species}.aligned.bed" | \
        awk -F '_' '{print $1"\t"$2"\t"$3}' | \
        bedtools sort | \
        bedtools merge > "results/align_base/${species}.aligned.merged.bed"
    
    bedtools genomecov -i "results/align_base/${species}.aligned.merged.bed" \
        -g hg38.chrom.sizes -bga > "results/align_bedgraph/${species}.bedgraph"
    
    # Step 3: Convert to bigWig
    echo "[$(date)] Converting to bigWig..." >> $log_file
    bedGraphToBigWig "results/align_bedgraph/${species}.bedgraph" \
        hg38.chrom.sizes \
        "results/align_bw/${species}.bw" 2>> $log_file
    
    # Step 4: Calculate iMotif alignment rate
    echo "[$(date)] Calculating iMotif alignment rate..." >> $log_file
    bigWigAverageOverBed \
        "results/align_bw/${species}.bw" \
        imotif_clean_final.bed \
        "results/align_ratio/${species}.imotif.ar.txt" 2>> $log_file
    
    if [ -f "results/align_ratio/${species}.imotif.ar.txt" ]; then
        mapped=$(wc -l < "results/align_ratio/${species}.imotif.ar.txt")
        echo "[$(date)] $species: Successfully processed $mapped iMotifs" >> $log_file
        echo "$species: $mapped" >> "results/mapping_summary.txt"
    fi
else
    echo "[$(date)] $species: No mapping results" >> $log_file
    echo "$species: No mapping results" >> "results/mapping_failed.txt"
fi

echo "[$(date)] $species processing completed" >> $log_file
INNER_EOF

chmod +x process_species.sh

# ========== Modification 3: Adjust concurrency ==========
total_species=$(wc -l < sp240_nonhuman.txt)
concurrent_jobs=54  # Adjust based on cluster available CPUs
echo "Submitting $total_species parallel tasks, concurrency $concurrent_jobs..."
sbatch --array=1-$total_species%$concurrent_jobs process_species.sh
# ======================================

echo "========================================="
echo "Master script submitted!"
echo "Monitor progress with:"
echo "  squeue -u $USER"
echo "  tail -f logs/species_*.log"
echo "========================================="
