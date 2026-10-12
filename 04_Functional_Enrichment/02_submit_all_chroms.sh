#!/bin/bash

# Batch submit makeEdges tasks for all chromosomes

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset

# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/life-zhanghk/niulk/my_project
# For local use: set IMOTIF_BASE to your own data directory.
BASE_DIR="${IMOTIF_BASE:-.}"

PROJECT_DIR="${BASE_DIR}/18_mammalian/03_upset"
SCRIPT_DIR="${PROJECT_DIR}/scripts/makeEdges"
LOG_DIR="${PROJECT_DIR}/logs"

# Create log directory
mkdir -p ${LOG_DIR}

# Clean old stats file
echo "chromosome,total_files,processed,failed,edges,nodes,timestamp" > ${PROJECT_DIR}/makeEdges_stats.csv

echo "=== Starting batch submission of makeEdges tasks ==="
echo "Time: $(date)"
echo ""

# Submit jobs for all chromosomes
job_ids=""

for chrom in {1..22} X Y; do
    # Check whether the corresponding input files exist
    if [ -d "${PROJECT_DIR}/organized_files/organized_consistent/hsa${chrom}" ]; then
        file_count=$(ls "${PROJECT_DIR}/organized_files/organized_consistent/hsa${chrom}"/*.rmredundant.df 2>/dev/null | wc -l)
        
        if [ $file_count -gt 0 ]; then
            echo "Submitting chromosome hsa${chrom} (${file_count} files)..."
            
            # Submit job
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
                echo "  Job ID: ${job_id}"
                job_ids="${job_ids} ${job_id}"
            else
                echo "  Submission failed"
            fi
        else
            echo "Skipping chromosome hsa${chrom} (no files)"
        fi
    else
        echo "Skipping chromosome hsa${chrom} (directory not found)"
    fi
done

echo ""
echo "=== Submission Complete ==="
echo "Time: $(date)"
echo "Submitted job IDs: ${job_ids}"
echo ""

# Create monitor script
cat > ${PROJECT_DIR}/monitor_makeEdges.sh << MONITOR_EOF
#!/bin/bash
echo "=== makeEdges job monitor ==="
echo "Time: \$(date)"
echo ""
echo "Running jobs:"
squeue -u \$USER -o "%.18i %.9P %.35j %.8u %.8T %.10M %.9l %.6D %R" | grep makeEdges
echo ""
echo "Completed job count:"
if [ -f "${PROJECT_DIR}/makeEdges_stats.csv" ]; then
    tail -n +2 "${PROJECT_DIR}/makeEdges_stats.csv" | while IFS=',' read -r chrom total processed failed edges nodes timestamp; do
        echo "hsa\${chrom}: processed \${processed}/\${total}, edges: \${edges}, nodes: \${nodes}"
    done
fi
MONITOR_EOF

chmod +x ${PROJECT_DIR}/monitor_makeEdges.sh

echo "Monitor script created: ${PROJECT_DIR}/monitor_makeEdges.sh"
echo "Run: ./monitor_makeEdges.sh to view progress"
