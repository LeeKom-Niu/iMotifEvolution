#!/usr/bin/env python3
"""
Cross-species promoter i-Motif presence and gene expression association analysis (raw TPM bar chart version)
Grouping rules:
  With i-Motif: has i-motif overlap (max_imotif_score non-empty)
  Without i-Motif: no i-motif overlap
Output:
  - One main directory per threshold, containing subdirectories for each species and corresponding
    single-species bar charts and combined plot
"""

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from matplotlib import rcParams

# ==================== Configuration ====================
BASE_DIR = os.getcwd()
GENE_BED_DIR = os.path.join(BASE_DIR, "gene_bed_files")
IMOTIF_BED_DIR = os.path.join(BASE_DIR, "imotif_bed")
TPM_DIR = os.path.join(BASE_DIR, "TPM_matrices")
BASE_OUTPUT_DIR = os.path.join(BASE_DIR, "imotif_analysis_results_barplot_raw_trim_high")

BEDTOOLS = "/datapool/home/2023200496/envs/bedtools/bin/bedtools"

SPECIES = ["sumatran", "gorilla", "human", "chimp", "bonobo"]
# Plotting order (can be modified as needed)
DISPLAY_ORDER = ["human", "bonobo", "chimp", "sumatran", "gorilla"]

# Species -> i-motif BED file mapping
SPECIES_TO_IMOTIF_BED = {
    "sumatran": "Pongo_abelii_all.bed",
    "gorilla": "Gorilla_gorilla_all.bed",
    "human": "Homo_sapiens_all.bed",
    "chimp": "Pan_troglodytes_all.bed",
    "bonobo": "Pan_paniscus_all.bed",
}

# Error bar type: 'sem' or 'std'
ERRORBAR_TYPE = 'sem'   # Standard error

# ========== New configuration: list of high-expression trimming percentiles to try ==========
TRIM_PERCENTILES = [95]   # Exclude values above this percentile
# =====================================================

def set_oup_style():
    plt.style.use('default')
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 11,
        'axes.linewidth': 0.5,
        'grid.linewidth': 0.3,
        'lines.linewidth': 1.0,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def print_file_head(file_path, n=5):
    if not os.path.exists(file_path):
        print(f"    File does not exist: {file_path}")
        return
    with open(file_path) as f:
        lines = [next(f) for _ in range(n)]
    for i, line in enumerate(lines):
        print(f"      Line {i+1}: {line.strip()}")

def create_promoter_bed(gene_bed_path, output_bed):
    df = pd.read_csv(gene_bed_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
    promoters = []
    for _, row in df.iterrows():
        chrom = row['chr']
        gene_id = row['gene_id']
        strand = row['strand']
        if strand == '+':
            prom_start = max(0, row['start'] - 1000)
            prom_end = row['start']
        else:
            prom_start = row['end']
            prom_end = row['end'] + 1000
        promoters.append([chrom, prom_start, prom_end, gene_id, 0, strand])
    promoter_df = pd.DataFrame(promoters)
    promoter_df.to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"  -> Generated promoter regions: {len(promoter_df)} regions -> {output_bed}")

def run_bedtools_intersect(promoter_bed, imotif_bed, output_intersect):
    cmd = [BEDTOOLS, "intersect", "-a", promoter_bed, "-b", imotif_bed,
           "-wa", "-wb", "-s"]
    try:
        with open(output_intersect, 'w') as fout:
            subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"  X bedtools command failed: {e.stderr.decode()}")
        return False
    if os.path.getsize(output_intersect) == 0:
        print(f"  W Overlap result is empty, no i-motif in promoter regions")
        return False
    else:
        print(f"  -> Overlap result non-empty: {output_intersect}")
        return True

def get_gene_max_imotif_score(intersect_file):
    df_intersect = pd.read_csv(intersect_file, sep='\t', header=None)
    if df_intersect.shape[1] != 12:
        print(f"   Warning: Overlap file has {df_intersect.shape[1]} columns, expected 12, attempting to continue...")
    gene_col = 3
    score_col = 10
    if df_intersect.shape[1] <= score_col:
        print(f"   Error: Overlap file has insufficient columns, cannot extract score")
        return pd.DataFrame(columns=['gene_id', 'max_imotif_score'])
    df_max = df_intersect.groupby(gene_col)[score_col].max().reset_index()
    df_max.columns = ['gene_id', 'max_imotif_score']
    return df_max

def load_tpm_matrix(species):
    tpm_file = os.path.join(TPM_DIR, f"{species}_TPM_matrix.csv")
    if not os.path.exists(tpm_file):
        raise FileNotFoundError(f"TPM matrix not found: {tpm_file}")
    tpm_df = pd.read_csv(tpm_file, index_col=0)
    tpm_df['mean_TPM'] = tpm_df.mean(axis=1)
    tpm_df['log2_TPM'] = np.log2(tpm_df['mean_TPM'] + 1)
    tpm_df = tpm_df.reset_index().rename(columns={'index': 'GeneID'})
    return tpm_df

def analyze_species(species, output_dir, trim_percentile, all_species_data):
    """
    Analyze a single species, save results to all_species_data list, and generate bar chart for that species (raw TPM)
    trim_percentile: percentile for trimming high expression outliers (e.g., 95), or None for no trimming
    """
    print(f"\n{'='*60}")
    print(f"Analyzing species: {species} (trim >{trim_percentile}% quantile)")
    print(f"{'='*60}")

    species_out = os.path.join(output_dir, species)
    ensure_dir(species_out)

    # 1. Gene BED
    gene_bed = os.path.join(GENE_BED_DIR, f"{species}.genes.bed")
    if not os.path.exists(gene_bed):
        print(f"  X Gene BED file does not exist: {gene_bed}, skipping")
        return

    # 2. i-motif BED
    imotif_bed_name = SPECIES_TO_IMOTIF_BED.get(species)
    if not imotif_bed_name:
        print(f"  X No i-motif BED mapping found for {species}, skipping")
        return
    imotif_bed = os.path.join(IMOTIF_BED_DIR, imotif_bed_name)
    if not os.path.exists(imotif_bed):
        print(f"  X i-motif BED file does not exist: {imotif_bed}, skipping")
        return

    print("\n  Gene BED first 3 lines:")
    print_file_head(gene_bed, 3)
    print("\n  i-motif BED first 3 lines:")
    print_file_head(imotif_bed, 3)

    # 3. Generate promoter BED
    promoter_bed = os.path.join(species_out, f"{species}_promoter.bed")
    create_promoter_bed(gene_bed, promoter_bed)

    # 4. bedtools intersect
    intersect_file = os.path.join(species_out, f"{species}_promoter_imotif_intersect.bed")
    has_overlap = run_bedtools_intersect(promoter_bed, imotif_bed, intersect_file)
    if not has_overlap:
        print(f"  -> No overlap found, skipping further analysis for {species}")
        return

    # 5. Extract max i-motif score per gene
    max_score_df = get_gene_max_imotif_score(intersect_file)
    print(f"  -> Found {len(max_score_df)} genes with i-motif overlap")

    # 6. Load TPM matrix
    tpm_df = load_tpm_matrix(species)
    print(f"  -> TPM matrix gene count: {len(tpm_df)}")

    # 7. Merge
    merged = tpm_df.merge(max_score_df, left_on='GeneID', right_on='gene_id', how='left')
    merged = merged.drop(columns=['gene_id'])

    # Group: with i-motif vs without i-motif
    merged['group'] = merged['max_imotif_score'].apply(
        lambda x: 'With i-Motif' if pd.notna(x) else 'Without i-Motif'
    )

    # Filter out genes with zero expression
    merged_filtered = merged[merged['mean_TPM'] > 0].copy()
    print(f"  -> Genes after filtering (expression>0): {len(merged_filtered)}")

    # Trim high expression outliers (if trim_percentile is not None)
    if trim_percentile is not None:
        all_tpm = merged_filtered['mean_TPM']
        high_thresh = np.percentile(all_tpm, trim_percentile)
        n_before = len(merged_filtered)
        merged_filtered = merged_filtered[merged_filtered['mean_TPM'] <= high_thresh].copy()
        n_after = len(merged_filtered)
        print(f"  -> After trimming >{trim_percentile}% quantile (threshold {high_thresh:.2f}), gene count from {n_before} to {n_after}")

    cat_counts = merged_filtered['group'].value_counts()
    print("\n  Group statistics:")
    for cat in ['With i-Motif', 'Without i-Motif']:
        cnt = cat_counts.get(cat, 0)
        print(f"    {cat}: {cnt} genes")

    # Extract raw TPM per group
    with_data = merged_filtered[merged_filtered['group'] == 'With i-Motif']
    without_data = merged_filtered[merged_filtered['group'] == 'Without i-Motif']

    if len(with_data) == 0 or len(without_data) == 0:
        print("  -> One group is empty, cannot perform statistical test and plot, skipping")
        # Still save data file
        data_out = os.path.join(species_out, f"{species}_imotif_analysis_data.csv")
        merged_filtered[['GeneID', 'group', 'max_imotif_score', 'mean_TPM', 'log2_TPM']].to_csv(data_out, index=False)
        print(f"  -> Data saved to: {data_out}")
        return

    # Statistical test (Mann-Whitney U) and raw TPM fold change
    tpm_with = with_data['mean_TPM']
    tpm_without = without_data['mean_TPM']
    p_tpm = stats.mannwhitneyu(tpm_with, tpm_without).pvalue
    median_with = tpm_with.median()
    median_without = tpm_without.median()
    fold_change = median_with / median_without if median_without > 0 else np.inf

    print("\n Mann-Whitney U test p-value (raw TPM):")
    print(f"  Raw TPM: With vs Without: p = {p_tpm:.2e}")
    print(f"\n Raw TPM medians:")
    print(f"   With i-Motif:  {median_with:.2f}")
    print(f"   Without i-Motif: {median_without:.2f}")
    print(f"  Fold change (With/Without): {fold_change:.2f}")

    # Calculate mean and standard error (or standard deviation) for raw TPM
    mean_with = tpm_with.mean()
    mean_without = tpm_without.mean()
    if ERRORBAR_TYPE == 'sem':
        error_with = tpm_with.sem()
        error_without = tpm_without.sem()
    else:  # std
        error_with = tpm_with.std()
        error_without = tpm_without.std()

    # Save data for combined plot
    all_species_data.append({
        'species': species,
        'mean_with': mean_with,
        'mean_without': mean_without,
        'error_with': error_with,
        'error_without': error_without,
        'n_with': len(tpm_with),
        'n_without': len(tpm_without),
        'p_tpm': p_tpm,
        'fold_change': fold_change
    })

    # ========== Draw single-species bar chart (raw TPM) ==========
    colors = {'With i-Motif': '#E64B35', 'Without i-Motif': '#4DBBD5'}
    categories = ['With i-Motif', 'Without i-Motif']
    means = [mean_with, mean_without]
    errors = [error_with, error_without]
    ns = [len(tpm_with), len(tpm_without)]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # Draw bar chart
    bars = ax.bar(categories, means, yerr=errors, capsize=5,
                  color=[colors[cat] for cat in categories],
                  edgecolor='black', linewidth=1, alpha=0.8,
                  error_kw={'linewidth': 1.5, 'ecolor': 'black'})

    # Add significance marker and fold change label
    y_max = max(means + errors) * 1.1  # Reserve space for annotation
    if p_tpm < 0.001: sig = '***'
    elif p_tpm < 0.01: sig = '**'
    elif p_tpm < 0.05: sig = '*'
    else: sig = 'ns'

    if p_tpm < 0.05:
        # Draw a line connecting the two bars
        x1, x2 = 0, 1
        y = y_max
        ax.plot([x1, x1, x2, x2], [y, y+0.05*y_max, y+0.05*y_max, y], 'k-', linewidth=0.8)
        ax.text((x1+x2)/2, y+0.08*y_max, sig, ha='center', va='bottom', fontsize=12)
        # Show fold change
        ax.text((x1+x2)/2, y_max*0.9, f'Fold = {fold_change:.2f}', ha='center', va='bottom', fontsize=10, style='italic')
    else:
        # Show only fold change and p-value
        ax.text(0.5, y_max*0.9, f'Fold = {fold_change:.2f}\np = {p_tpm:.2e}',
                ha='center', va='bottom', fontsize=9, style='italic')

    # Set labels and title
    ax.set_ylabel('Mean TPM (raw)', fontsize=14)
    ax.set_title(f'{species.capitalize()}', loc='left', fontweight='bold', fontsize=16)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([f'{cat}\n(n={n})' for cat, n in zip(categories, ns)], fontsize=12)
    ax.set_ylim(0, y_max * 1.05)  # Leave some top space

    # Optional: add grid lines
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    # Save figure
    pdf_path = os.path.join(species_out, f"{species}_barplot_raw.pdf")
    plt.savefig(pdf_path, format='pdf')
    plt.savefig(os.path.join(species_out, f"{species}_barplot_raw.tiff"), format='tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(os.path.join(species_out, f"{species}_barplot_raw.png"), format='png', dpi=150)
    plt.close()
    print(f"  -> Bar chart saved to: {species_out}/")

    # Save data file
    data_out = os.path.join(species_out, f"{species}_imotif_analysis_data.csv")
    merged_filtered[['GeneID', 'group', 'max_imotif_score', 'mean_TPM', 'log2_TPM']].to_csv(data_out, index=False)
    print(f"  -> Data saved to: {data_out}")

    print(f"\nOK {species} analysis complete!")

def plot_combined_barplot(all_species_data, output_dir, trim_percentile):
    """
    Generate combined grouped bar chart for all species (single plot, x-axis = species, landscape layout)
    Sort by DISPLAY_ORDER, optimize annotation positions
    trim_percentile: used for title description
    """
    # Reorder data according to DISPLAY_ORDER
    data_dict = {data['species']: data for data in all_species_data}
    sorted_data = [data_dict[sp] for sp in DISPLAY_ORDER if sp in data_dict]
    n = len(sorted_data)
    if n == 0:
        print("No valid data, cannot generate combined plot.")
        return

    # Landscape layout: width proportional to species count, fixed small height for landscape shape
    fig_width = max(8, n * 1.2)   # ~1.2 inches per species
    fig_height = 5                # Fixed height for landscape layout
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    # Color configuration
    colors = {'With i-Motif': '#E64B35', 'Without i-Motif': '#4DBBD5'}

    # Prepare data
    species_names = [data['species'].capitalize() for data in sorted_data]
    means_with = [data['mean_with'] for data in sorted_data]
    means_without = [data['mean_without'] for data in sorted_data]
    errors_with = [data['error_with'] for data in sorted_data]
    errors_without = [data['error_without'] for data in sorted_data]
    p_values = [data['p_tpm'] for data in sorted_data]
    folds = [data['fold_change'] for data in sorted_data]

    # Set grouped bar positions
    x = np.arange(n)  # Species positions
    width = 0.35      # Bar width
    # Draw bars (With i-Motif on left, Without i-Motif on right)
    bars_with = ax.bar(x - width/2, means_with, width, yerr=errors_with, capsize=5,
                       color=colors['With i-Motif'], edgecolor='black', linewidth=1, alpha=0.8,
                       label='With i-Motif', error_kw={'linewidth': 1.5, 'ecolor': 'black'})
    bars_without = ax.bar(x + width/2, means_without, width, yerr=errors_without, capsize=5,
                          color=colors['Without i-Motif'], edgecolor='black', linewidth=1, alpha=0.8,
                          label='Without i-Motif', error_kw={'linewidth': 1.5, 'ecolor': 'black'})

    # Calculate maximum height for each bar group (for placing markers)
    max_heights = [max(means_with[i] + errors_with[i], means_without[i] + errors_without[i]) for i in range(n)]
    # Reserve space for markers: 1.3x max height, but at least 20% above max
    y_max_plot = max(max_heights) * 1.3
    if y_max_plot < 1e-5:
        y_max_plot = 1.0

    # Add baseline (y=0 line)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, zorder=0)

    # Add significance markers and fold change labels
    for i in range(n):
        # Current bar group maximum (considering error bars)
        top_with = means_with[i] + errors_with[i]
        top_without = means_without[i] + errors_without[i]
        max_top = max(top_with, top_without)
        # Marker start y-coordinate: bar group max + spacing (dynamic, 5% of max height)
        y_start = max_top + 0.05 * y_max_plot
        # If marker exceeds y-axis limit, adjust y-axis limit
        if y_start > y_max_plot:
            y_max_plot = y_start + 0.05 * y_max_plot

        # Draw horizontal line connecting the two bars
        x_left = x[i] - width/2 + width
        x_right = x[i] + width/2
        ax.plot([x_left, x_right], [y_start, y_start], 'k-', linewidth=0.8)

        # Determine star rating based on p-value
        p = p_values[i]
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'
        else:
            sig = 'ns'

        # Annotate stars (above horizontal line)
        ax.text(x[i], y_start + 0.02 * y_max_plot, sig, ha='center', va='bottom', fontsize=10)

        # Annotate fold change (below horizontal line, inside bar group)
        fold = folds[i]
        # Select the taller of the two bars, place fold text inside
        higher_bar_top = max(means_with[i], means_without[i])
        if higher_bar_top > 0:
            y_fold = higher_bar_top + 0.02 * y_max_plot
            # Ensure it does not exceed the horizontal line
            if y_fold >= y_start:
                y_fold = y_start - 0.03 * y_max_plot
        else:
            y_fold = y_start - 0.03 * y_max_plot
        ax.text(x[i], y_fold, f'Fold={fold:.2f}', ha='center', va='bottom', fontsize=8, style='italic')

    # Set axis
    ax.set_xticks(x)
    ax.set_xticklabels(species_names, rotation=0, ha='center', fontsize=11)
    ax.set_ylabel('Mean TPM', fontsize=14)
    if trim_percentile is not None:
        title = f'Expression comparison (genes with/without i-Motif in promoters)\n(high expression trimmed >{trim_percentile}% quantile)'
    else:
        title = 'Expression comparison (genes with/without i-Motif in promoters)\n(no trimming)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', frameon=True, fontsize=10)
    ax.set_ylim(0, y_max_plot)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    # Save figure
    outfile = os.path.join(output_dir, "combined_barplot_raw.pdf")
    plt.savefig(outfile, format='pdf')
    plt.savefig(os.path.join(output_dir, "combined_barplot_raw.png"), format='png', dpi=150)
    plt.close()
    print(f"\nCombined plot (landscape) saved: {outfile}")

def main():
    set_oup_style()
    ensure_dir(BASE_OUTPUT_DIR)

    # Loop over each trimming threshold
    for trim_percentile in TRIM_PERCENTILES:
        print(f"\n{'='*60}")
        print(f"Starting processing for trim threshold: >{trim_percentile}% quantile")
        print(f"{'='*60}")

        # Create output directory for this threshold
        output_dir = os.path.join(BASE_OUTPUT_DIR, f"trim_high_{trim_percentile}")
        ensure_dir(output_dir)

        all_species_data = []  # Collect all species data for this threshold

        for sp in SPECIES:
            try:
                analyze_species(sp, output_dir, trim_percentile, all_species_data)
            except Exception as e:
                print(f"  X Error analyzing {sp}: {e}")
                import traceback
                traceback.print_exc()
                continue

        # Generate combined plot for this threshold
        if all_species_data:
            plot_combined_barplot(all_species_data, output_dir, trim_percentile)
        else:
            print("No species successfully analyzed, cannot generate combined plot.")

        print(f"\nThreshold {trim_percentile}% processing complete, results in: {output_dir}")

    print("\n" + "="*60)
    print("All thresholds processed! Overall results directory: " + BASE_OUTPUT_DIR)
    print("="*60)

if __name__ == "__main__":
    main()
