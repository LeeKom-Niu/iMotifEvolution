#!/bin/bash
# ====================================================================
# T2T iM predicted peaks vs promoter vs ATAC analysis (optimized version)
# Output directory: results_v2/
# Features: TSS ± 1kb promoter, filter score ≥ 0.5
# ====================================================================

# Original server paths (this code was developed on an HPC cluster):
#   /datapool/life-zhanghk/niulk/my_project/19_modified_analysis/06_ATAC/T2T_ATAC
# Root directory for input/output files.
# On the original cluster: export IMOTIF_BASE=/datapool/life-zhanghk/niulk/my_project
BASE_DIR="${IMOTIF_BASE:-.}"

cd "${BASE_DIR}/19_modified_analysis/06_ATAC/T2T_ATAC" || exit
mkdir -p results_v2

# 1. Generate TSS ±1kb promoter (if not generated yet)
if [ ! -f promoter_TSSpm1kb.bed ]; then
    echo "Generating TSS ±1kb promoter regions..."
    awk -v OFS='\t' '{
        if ($6 == "+") {
            tss = $3;
            new_start = tss - 1000;
            new_end = tss + 1000;
        } else if ($6 == "-") {
            tss = $2;
            new_start = tss - 1000;
            new_end = tss + 1000;
        }
        if (new_start < 0) new_start = 0;
        print $1, new_start, new_end, $4, $5, $6;
    }' promoter_regions.bed > promoter_TSSpm1kb.bed
fi

# 2. Filter iM peaks (keep score >= 0.5, score is column 4)
echo "Filtering iM peaks: keep score >= 0.5 (column 4)..."
awk '$4 >= 0.5' Homo_sapiens_all.bed > Homo_sapiens_filtered.bed
echo "Line count after filtering: $(wc -l < Homo_sapiens_filtered.bed)"

IM_BED="Homo_sapiens_filtered.bed"
PROMOTER_BED="promoter_TSSpm1kb.bed"
ATAC_BED="GSE274771_ATAC_WT_peaks.narrowPeak"

echo "========================================"
echo "T2T analysis (TSS±1kb, score≥0.5)"
echo "========================================"

echo "=== 1. Count lines in raw files ==="
wc -l ${IM_BED} ${PROMOTER_BED} ${ATAC_BED}

echo "=== 2. Compute overlap between iM and promoters (iM promoter peaks) ==="
bedtools intersect -a ${IM_BED} -b ${PROMOTER_BED} -u > results_v2/im_promoter.bed
IM_PROMOTER=$(wc -l < results_v2/im_promoter.bed)
echo "iM promoter peak count: ${IM_PROMOTER}"

echo "=== 3. Compute overlap between iM promoter peaks and ATAC ==="
bedtools intersect -a results_v2/im_promoter.bed -b ${ATAC_BED} -u > results_v2/im_promoter_atac.bed
IM_PROMOTER_ATAC=$(wc -l < results_v2/im_promoter_atac.bed)
echo "iM promoter ∩ ATAC: ${IM_PROMOTER_ATAC}"

echo "=== 4. Count ATAC total ==="
ATAC_TOTAL=$(wc -l < ${ATAC_BED})
echo "ATAC peak total: ${ATAC_TOTAL}"

ONLY_IM=$((IM_PROMOTER - IM_PROMOTER_ATAC))
ONLY_ATAC=$((ATAC_TOTAL - IM_PROMOTER_ATAC))

echo "iM promoter only: ${ONLY_IM}"
echo "ATAC only: ${ONLY_ATAC}"
echo "Overlap ratio: $(echo "scale=2; ${IM_PROMOTER_ATAC}*100/${IM_PROMOTER}" | bc)%"

cat > results_v2/statistics.txt << EOF
T2T analysis statistics (TSS ±1kb promoter, score >= 0.5)
===============================================
iM predicted peak total (after filtering): $(wc -l < ${IM_BED})
Promoter region count (TSS±1kb): $(wc -l < ${PROMOTER_BED})
ATAC peak total: ${ATAC_TOTAL}
iM promoter peak count: ${IM_PROMOTER}
iM promoter ∩ ATAC: ${IM_PROMOTER_ATAC}
iM promoter only: ${ONLY_IM}
ATAC only: ${ONLY_ATAC}
Overlap ratio: $(echo "scale=2; ${IM_PROMOTER_ATAC}*100/${IM_PROMOTER}" | bc)%
EOF

echo "Statistics saved to results_v2/statistics.txt"

# ========== Plot Venn diagram ==========
echo "=== 5. Plot Venn diagram ==="
cat > results_v2/plot_venn.py << 'EOF'
import matplotlib.pyplot as plt
import sys

if len(sys.argv) != 4:
    print("Usage: python plot_venn.py A_total B_total A_and_B")
    sys.exit(1)

A_total = int(sys.argv[1])
B_total = int(sys.argv[2])
A_and_B = int(sys.argv[3])

only_A = A_total - A_and_B
only_B = B_total - A_and_B

print(f"iM promoter peak total: {A_total}")
print(f"ATAC peak total: {B_total}")
print(f"Intersection: {A_and_B}")
print(f"iM promoter only: {only_A}")
print(f"ATAC only: {only_B}")

fig, ax = plt.subplots(figsize=(6, 5))
r = 1.0
cx1, cy1 = -0.5, 0.0
cx2, cy2 = 0.5, 0.0

circle1 = plt.Circle((cx1, cy1), r, color='#E69F00', alpha=0.4)
circle2 = plt.Circle((cx2, cy2), r, color='#56B4E9', alpha=0.4)
ax.add_patch(circle1)
ax.add_patch(circle2)

ax.text(cx1, cy1 - r - 0.2, 'iM in Promoters', ha='center', va='top', fontsize=10)
ax.text(cx2, cy2 - r - 0.2, 'ATAC (GSE274771)', ha='center', va='top', fontsize=10)

ax.text(cx1 - 0.6, cy1, f'{only_A}', ha='center', va='center', fontsize=12, fontweight='bold')
ax.text(cx2 + 0.6, cy2, f'{only_B}', ha='center', va='center', fontsize=12, fontweight='bold')
ax.text((cx1+cx2)/2, cy1 - 0.1, f'{A_and_B}', ha='center', va='center', fontsize=12, fontweight='bold')

ax.set_xlim(-2, 2)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect('equal')
ax.axis('off')
plt.title('iM Promoter vs ATAC (T2T, score≥0.5)', fontsize=12)
plt.tight_layout()
plt.savefig("venn_im_promoter_vs_atac_T2T_v2.pdf", dpi=300)
plt.savefig("venn_im_promoter_vs_atac_T2T_v2.png", dpi=300)
print("Venn diagram saved")
EOF

cd results_v2
python3 plot_venn.py ${IM_PROMOTER} ${ATAC_TOTAL} ${IM_PROMOTER_ATAC}

echo "Analysis complete! Results are in the results_v2/ directory."
