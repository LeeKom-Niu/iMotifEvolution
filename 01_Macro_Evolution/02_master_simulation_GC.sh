#!/bin/bash
# master_simulation_GC.sh - Master script for GC content background simulation

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC
#   /datapool/home/2023200496/miniconda3/etc/profile.d/conda.sh

# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/home/2023200496/niulk/my_project
# For local use: set IMOTIF_BASE to your own data directory.
BASE_DIR="${IMOTIF_BASE:-.}"

# In this script BASE_DIR means the simulation study directory (a subdirectory of
# the project root); keep that original meaning intact for all downstream uses.
BASE_DIR="${BASE_DIR}/1555_simulation_study_extended_GC"

echo "=== GC Content Background Simulation ==="
echo "Start time: $(date)"
echo "Working directory: $BASE_DIR"
echo "Current node: $(hostname)"

# Set parameters directly
GENOME_SIZE=10000000  # 10Mb
REPLICATES=30
GC_LEVELS=("0.10" "0.20" "0.30" "0.40" "0.50" "0.60" "0.70" "0.80" "0.90")

echo "Genome size: $GENOME_SIZE bp"
echo "Number of replicates: $REPLICATES"
echo "GC level count: ${#GC_LEVELS[@]}"
echo "GC levels: ${GC_LEVELS[@]}"

# Create directory structure
mkdir -p $BASE_DIR/scripts
mkdir -p $BASE_DIR/sequences
mkdir -p $BASE_DIR/results
mkdir -p $BASE_DIR/logs
mkdir -p $BASE_DIR/figures

# Submit a job for each GC level
for i in "${!GC_LEVELS[@]}"; do
    gc_level="${GC_LEVELS[$i]}"
    gc_percent=$(python3 -c "print(int(float('$gc_level') * 100))")
    
    echo "Submitting GC level ${gc_level} (${gc_percent}%) to SLURM queue..."
    
    sbatch << EOF
#!/bin/bash
#SBATCH --job-name=im_GC${gc_percent}
#SBATCH --cpus-per-task=60
#SBATCH --mem=250G
#SBATCH --time=300-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=$BASE_DIR/logs/GC${gc_percent}_%j.out
#SBATCH --error=$BASE_DIR/logs/GC${gc_percent}_%j.err

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate im-seeker

cd $BASE_DIR
echo "=== Starting processing GC content ${gc_level} ==="
echo "Run node: \$(hostname)"
echo "Start time: \$(date)"
echo "Number of replicates: $REPLICATES"
echo "Genome size: $GENOME_SIZE bp"

python3 $BASE_DIR/scripts/parallel_simulation_GC.py $gc_level $REPLICATES $GENOME_SIZE

echo "GC content ${gc_level} completed at \$(date)"
EOF

done

echo "All gradient jobs submitted to SLURM queue!"
echo "Use 'squeue -u \$USER' to check job status"
echo "End time: $(date)"
