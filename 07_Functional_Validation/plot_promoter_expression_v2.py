#!/usr/bin/env python3
"""
plot_promoter_expression_v2.py
Pure violin plot + jittered scatter (remove outliers, raw TPM, no boxplot)
Show the expression distribution of genes associated with unique iM and unique G4 promoters
Add significance markers (generated dynamically by the Mann-Whitney U test)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

# Set random seed so the jitter is reproducible
np.random.seed(42)

# Read data
gene_tpm = pd.read_csv("HEK293T_gene_TPM.tsv", sep='\t')
gene_tpm_exp = gene_tpm[
    (gene_tpm['TPM_rep1'] > 0) |
    (gene_tpm['TPM_rep2'] > 0) |
    (gene_tpm['TPM_rep3'] > 0)
].copy()

# Read overlap files (with gene IDs)
im_overlap = pd.read_csv("iM_unique_promoter_overlap.bed", sep='\t', header=None)
im_genes = im_overlap.iloc[:, 9].unique()
g4_overlap = pd.read_csv("G4_unique_promoter_overlap.bed", sep='\t', header=None)
g4_genes = g4_overlap.iloc[:, 9].unique()

# Get expression values
im_expr = gene_tpm_exp[gene_tpm_exp['gene_id'].isin(im_genes)]['TPM_mean'].values
g4_expr = gene_tpm_exp[gene_tpm_exp['gene_id'].isin(g4_genes)]['TPM_mean'].values
all_expr = gene_tpm_exp['TPM_mean'].values

# Remove outliers (IQR method)
def remove_outliers(x):
    if len(x) == 0:
        return x
    q25, q75 = np.percentile(x, [25, 75])
    iqr = q75 - q25
    lower = q25 - 1.5 * iqr
    upper = q75 + 1.5 * iqr
    return x[(x >= lower) & (x <= upper)]

im_clean = remove_outliers(im_expr)
g4_clean = remove_outliers(g4_expr)
all_clean = remove_outliers(all_expr)

data_to_plot = [im_clean, g4_clean, all_clean]
labels = ['iM unique near TSS', 'G4 unique near TSS', 'All Genes']
colors = ['#E69F00', '#56B4E9', '#B3B3B3']

fig, ax = plt.subplots(figsize=(6, 5))

# Pure violin plot (no internal markers)
vp = ax.violinplot(data_to_plot, positions=[1,2,3], widths=0.6,
                   showmeans=False, showmedians=False, showextrema=False)
for i, body in enumerate(vp['bodies']):
    body.set_facecolor(colors[i])
    body.set_alpha(0.6)

# Jittered scatter
for i, data in enumerate(data_to_plot):
    x_pos = i + 1
    jitter = np.random.normal(0, 0.08, size=len(data))
    ax.scatter(x_pos + jitter, data, s=4, alpha=0.3, color=colors[i], edgecolor='none')

# ========== Significance marker ==========
# Mann-Whitney U test (two-sided): iM unique near TSS vs All Genes
u_stat, p_value = scipy.stats.mannwhitneyu(im_clean, all_clean, alternative='two-sided')

if p_value < 0.001:
    stars = '***'
elif p_value < 0.01:
    stars = '**'
elif p_value < 0.05:
    stars = '*'
else:
    stars = 'ns'

print("Mann-Whitney U test (iM unique near TSS vs All Genes)")
print(f"  n(iM) = {len(im_clean)}, n(All) = {len(all_clean)}")
print(f"  U = {u_stat:.1f}, p = {p_value:.6g}  ->  {stars}")

y_max = max([max(data) for data in data_to_plot if len(data) > 0])
y_line = y_max * 1.03          # line position slightly lowered
y_text = y_line * 1.02         # asterisk position lowered accordingly

# iM-only vs All Genes (positions 1 vs 3)
x1, x2 = 1, 3
ax.plot([x1, x1, x2, x2], [y_line, y_line*1.02, y_line*1.02, y_line], 'k-', linewidth=1)
ax.text((x1+x2)/2, y_text, stars, ha='center', va='bottom', fontsize=14, fontweight='bold')
# ====================================

ax.set_xticks([1,2,3])
ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel('Gene Expression (TPM)', fontsize=12)
ax.set_title('HEK293T:Expression near TSS', fontsize=12)
ax.set_ylim(bottom=0)
ax.yaxis.grid(True, linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("promoter_v2_unique.pdf", dpi=300)
plt.savefig("promoter_v2_unique.png", dpi=300)
print("Version 2 (unique) saved as promoter_v2_unique.pdf/.png")
