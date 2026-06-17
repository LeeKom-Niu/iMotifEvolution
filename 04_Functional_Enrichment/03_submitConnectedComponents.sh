#!/bin/bash

#SBATCH --job-name=connectedComponents
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G  # Total memory = 25 tasks * 16G
#SBATCH --time=50-00:00:00
#SBATCH --partition=life-zhanghk-fat
#SBATCH --output=../logs/connectedComponents_%j.out
#SBATCH --error=../logs/connectedComponents_%j.err

PROJECT_DIR="/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
SCRIPT_DIR="${PROJECT_DIR}/scripts/connectedComponents"
OUTPUT_DIR="${PROJECT_DIR}/output"
PYTHON_SCRIPT="${PROJECT_DIR}/makeConnectedComponents.py"

echo "=== Starting Connected Components Analysis ==="
echo "Time: $(date)"
echo "Project directory: ${PROJECT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Python script: ${PYTHON_SCRIPT}"

# Check the required script
if [ ! -f "${PYTHON_SCRIPT}" ]; then
    echo "Error: makeConnectedComponents.py script not found"
    echo "Attempting to copy from original repository..."
    cp "${PROJECT_DIR}/GreatApeT2T-G4s-main/src/mapseaAndPostWorkflow/makeConnectedComponents.py" "${PYTHON_SCRIPT}"
    
    if [ ! -f "${PYTHON_SCRIPT}" ]; then
        echo "Error: Cannot copy makeConnectedComponents.py"
        exit 1
    fi
fi

# Check networkx library
python3 -c "import networkx as nx; print('networkx version:', nx.__version__)" || {
    echo "Error: networkx library not installed"
    echo "Install command: pip install networkx pandas"
    exit 1
}

# Create log directory
mkdir -p ../logs

echo "Starting parallel processing of all chromosomes..."

# Use GNU Parallel to process all chromosomes in parallel
# Note: directly calls the python script instead of submitting sub-jobs via sbatch
parallel -v -j $SLURM_NTASKS '
    chrom={}
    egs_file="'${OUTPUT_DIR}'/hsa${chrom}/hsa${chrom}.egs"
    graph_file="'${OUTPUT_DIR}'/hsa${chrom}/hsa${chrom}.graph"
    
    if [ -f "${egs_file}" ]; then
        edges=$(wc -l < "${egs_file}")
        if [ ${edges} -gt 0 ]; then
            echo "Processing hsa${chrom} (${edges} edges)..."
            python3 "'${PYTHON_SCRIPT}'" "${egs_file}" > "${graph_file}"
            if [ $? -eq 0 ]; then
                components=$(wc -l < "${graph_file}")
                echo "  hsa${chrom}: Generated ${components} connected components"
            else
                echo "  hsa${chrom}: Processing failed"
            fi
        else
            echo "Skipping hsa${chrom} (no edge data)"
            touch "${graph_file}"  # Create empty file
        fi
    else
        echo "Skipping hsa${chrom} (.egs file does not exist)"
    fi
' ::: {1..22} X Y

echo ""
echo "=== Connected Components Processing Complete ==="
echo "Time: $(date)"

# Results statistics
echo ""
echo "=== Results Statistics ==="
total_components=0
for chrom in {1..22} X Y; do
    graph_file="${OUTPUT_DIR}/hsa${chrom}/hsa${chrom}.graph"
    if [ -f "${graph_file}" ]; then
        components=$(wc -l < "${graph_file}" 2>/dev/null || echo "0")
        if [ $components -gt 0 ]; then
            echo "hsa${chrom}: ${components} connected components"
            total_components=$((total_components + components))
        fi
    fi
done

echo ""
echo "Total: ${total_components} connected components"
echo "Output files location: ${OUTPUT_DIR}/hsa*/hsa*.graph"

# Save statistics
echo "chromosome,connected_components,timestamp" > "${PROJECT_DIR}/connectedComponents_stats.csv"
for chrom in {1..22} X Y; do
    graph_file="${OUTPUT_DIR}/hsa${chrom}/hsa${chrom}.graph"
    if [ -f "${graph_file}" ]; then
        components=$(wc -l < "${graph_file}" 2>/dev/null || echo "0")
        echo "${chrom},${components},$(date)" >> "${PROJECT_DIR}/connectedComponents_stats.csv"
    fi
done
