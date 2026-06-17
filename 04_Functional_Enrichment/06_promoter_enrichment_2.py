#!/usr/bin/env python3
"""
Calculate and plot pG4 enrichment in promoter regions (Figure 5A human promoter part)
Corrected version: y-axis goes from old to young (top to bottom): Great ape -> Homininae -> Hominini -> Human-specific
Compliant with OUP illustration guidelines: font >=7pt, line width 0.25-1pt, colorblind-friendly palette, editable PDF text.
Removed title and y-axis label, baseline retains only dashed line without text, y-axis names abbreviated (removed "IMs").
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
import os
import gzip

# ===== OUP illustration guideline style settings (uniformly enlarged font) =====
def set_oup_style():
    """Set plotting style compliant with OUP guidelines (font uniformly enlarged to 16pt)"""
    rcParams.update({
        # Font settings
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 14,                     # Base font (slightly larger)
        # Axes
        'axes.labelsize': 16,                 # Axis label (Fold Enrichment)
        'axes.titlesize': 16,                 # Title (removed)
        'axes.linewidth': 0.5,                 # Axis line width
        'axes.edgecolor': 'black',
        'axes.labelpad': 10,
        # Ticks
        'xtick.labelsize': 16,                 # X-axis tick labels
        'ytick.labelsize': 16,                 # Y-axis tick labels (group names)
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        # Legend
        'legend.fontsize': 12,                  # Legend text
        'legend.frameon': False,
        # Lines
        'lines.linewidth': 1.0,                # Horizontal line width
        # Figure dimensions
        'figure.dpi': 300,
        'savefig.dpi': 600,
        'figure.figsize': (11, 6.5),            # Slightly enlarged to accommodate large font
        # PDF text editable
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

set_oup_style()

# ===== Configuration parameters =====
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
INPUT_DIR = os.path.join(BASE_DIR, "functionalOutputs/Homo_sapiens")
OUTPUT_DIR = os.path.join(BASE_DIR, "enrichment_plots_oup")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Color mapping (retained)
COLOR_MAP = {
    "Great ape": "#004d71",          # Dark blue
    "Homininae": "#8759a1",          # Purple
    "Hominini": "#f75a78",            # Pink
    "Human-specific": "#ffa600"       # Orange
}

# Group mapping (note: removed IMs)
GROUP_MAP = {
    "hominid": "Great ape",
    "homininae": "Homininae",
    "hominini": "Hominini",
    "humanSpecific": "Human-specific"
}

# For internal processing, keep full name to color mapping, but use abbreviations for display
DISPLAY_GROUP = GROUP_MAP

# ===== Basic genome information =====
GENOME_LENGTH = 3117275501

def calculate_promoter_length():
    """Calculate total promoter region length"""
    promoter_file = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/functionalOutputs/Homo_sapiens/promoter_regions.bed.gz"

    if not os.path.exists(promoter_file):
        print(f"Warning: Promoter file does not exist {promoter_file}")
        print("Using value from author's code: 18,842,577 bp")
        return 18842577

    total_length = 0
    with gzip.open(promoter_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                start = int(parts[1])
                end = int(parts[2])
                total_length += (end - start)

    print(f"Total promoter length: {total_length:,} bp")
    return total_length

PROMOTER_LENGTH = calculate_promoter_length()

# ===== Total pG4 counts per group =====
GROUP_TOTALS = {
    "hominid": 314942,
    "homininae": 124431,
    "hominini": 34964,
    "humanSpecific": 104483
}

# ===== Main functions =====
def count_pg4_in_promoters():
    """Count pG4s in promoters for each group, output using abbreviated names"""
    print("=== Counting pG4s in promoters ===")

    counts = {}

    for group_short, group_display in GROUP_MAP.items():
        input_file = os.path.join(INPUT_DIR, f"allhsaG.intersected.betn.human_promoter.{group_short}G4s.bed.gz")

        if not os.path.exists(input_file):
            print(f"Warning: File does not exist {input_file}")
            counts[group_display] = 0
            continue

        pg4_ids = set()
        with gzip.open(input_file, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pg4_id = f"{parts[0]}:{parts[1]}-{parts[2]}:{parts[5]}"
                    pg4_ids.add(pg4_id)

        counts[group_display] = len(pg4_ids)
        print(f"  {group_display}: {len(pg4_ids):,} unique pG4s")

    return counts

def calculate_fold_enrichment(promoter_counts):
    """
    Calculate fold enrichment
    """
    print("\n=== Calculating fold enrichment ===")

    results = []

    for group_display, promoter_count in promoter_counts.items():
        group_short = [k for k, v in GROUP_MAP.items() if v == group_display][0]
        total_pg4 = GROUP_TOTALS[group_short]

        pG4_ratio = promoter_count / total_pg4
        promoter_ratio = PROMOTER_LENGTH / GENOME_LENGTH
        fold_enrichment = pG4_ratio / promoter_ratio

        percentage_in_promoter = (promoter_count / total_pg4) * 100
        percentage_genome = (PROMOTER_LENGTH / GENOME_LENGTH) * 100

        # Use abbreviated name as group name
        results.append({
            'Group': group_display,
            'Short_Name': group_short,
            'Promoter_Count': promoter_count,
            'Total_pG4': total_pg4,
            '%_in_Promoter': percentage_in_promoter,
            '%_Genome_Promoter': percentage_genome,
            'Fold_Enrichment': fold_enrichment,
            'Color': COLOR_MAP[group_display]  # Color directly from abbreviated name
        })

        print(f"  {group_display}:")
        print(f"    pG4 in promoter: {promoter_count:,} / {total_pg4:,} = {percentage_in_promoter:.2f}%")
        print(f"    Promoter / Genome: {PROMOTER_LENGTH:,} / {GENOME_LENGTH:,} = {percentage_genome:.4f}%")
        print(f"    Fold enrichment: {fold_enrichment:.2f}")

    return pd.DataFrame(results)

def plot_enrichment(results_df):
    """
    Plot horizontal lollipop chart
    Y-axis order: top to bottom as Great ape -> Homininae -> Hominini -> Human-specific
    """
    print("\n=== Plotting enrichment (horizontal lollipop chart) ===")

    # Sort by evolutionary order
    order = ["Great ape", "Homininae", "Hominini", "Human-specific"]
    results_df['Order'] = results_df['Group'].map({g: i for i, g in enumerate(order)})
    results_df = results_df.sort_values('Order', ascending=True)

    fig, ax = plt.subplots(figsize=(11, 6.5))

    y_pos = np.arange(len(results_df))

    # Draw horizontal lines
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax.hlines(
            y=i,
            xmin=0,
            xmax=row['Fold_Enrichment'],
            color=row['Color'],
            linewidth=1.0,
            alpha=0.7
        )

    # Draw points
    ax.scatter(
        results_df['Fold_Enrichment'],
        y_pos,
        c=results_df['Color'],
        s=200,
        edgecolor='black',
        linewidth=0.5,
        zorder=5
    )

    # Add count labels (slightly larger font)
    for i, (_, row) in enumerate(results_df.iterrows()):
        label_text = f"{row['Promoter_Count']:,} ({row['Fold_Enrichment']:.2f}x)"
        ax.text(
            row['Fold_Enrichment'] + 0.2,
            i,
            label_text,
            va='center',
            ha='left',
            fontsize=12,          # Value label 12pt
            fontweight='normal',
            color=row['Color']
        )

    # Set Y-axis tick labels (abbreviated names)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Group'], fontsize=16)  # Consistent with tick label size
    ax.set_ylabel('')                 # No y-axis label
    ax.invert_yaxis()

    # Set X-axis (font controlled by axes.labelsize, explicitly specified here for consistency)
    ax.set_xlabel('Fold Enrichment', fontsize=16)
    ax.set_title('')

    # Baseline dashed line
    ax.axvline(x=1, color='black', linestyle='--', linewidth=0.5, alpha=0.5)

    # Grid lines
    ax.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.3)
    ax.set_axisbelow(True)

    # X-axis range
    x_max = results_df['Fold_Enrichment'].max() * 1.45
    ax.set_xlim([0, x_max])

    # Legend (using abbreviated names)
    legend_elements = []
    for _, row in results_df.iterrows():
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=row['Color'], markersize=12,
                   label=f"{row['Group']}\n{row['Promoter_Count']:,} in promoters")
        )
    ax.legend(handles=legend_elements, loc='lower right', fontsize=12, frameon=False)

    plt.tight_layout()

    # Save files
    pdf_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.pdf")
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    print(f"  PDF saved: {pdf_path}")

    tiff_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.tiff")
    plt.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    print(f"  TIFF saved: {tiff_path}")

    png_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"  PNG saved: {png_path}")

    svg_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    print(f"  SVG saved: {svg_path}")

    plt.close()

    return fig, ax

def compare_with_paper(results_df):
    """Compare with original Figure 5A (unchanged)"""
    print("\n=== Comparison with original Figure 5A ===")

    paper_values = {
        "Great ape": {"fold": 6.71, "count": "~?"},
        "Homininae": {"fold": 5.77, "count": "~?"},
        "Hominini": {"fold": 4.67, "count": "~?"},
        "Human-specific": {"fold": 4.11, "count": 2753}
    }

    print(f"{'Group':<20} {'Your Fold':<12} {'Paper Fold':<12} {'Difference':<12} {'Your Count':<12} {'Paper Count':<12}")
    print("-" * 85)

    for _, row in results_df.iterrows():
        group = row['Group']
        your_fold = row['Fold_Enrichment']
        your_count = row['Promoter_Count']

        if group in paper_values:
            paper_fold = paper_values[group]['fold']
            paper_count = paper_values[group]['count']
            diff_fold = your_fold - paper_fold
            diff_percent = (diff_fold / paper_fold) * 100
            print(f"{group:<20} {your_fold:<12.2f} {paper_fold:<12.2f} {diff_percent:<12.1f}% {your_count:<12} {paper_count:<12}")
        else:
            print(f"{group:<20} {your_fold:<12.2f} {'N/A':<12} {'N/A':<12} {your_count:<12} {'N/A':<12}")

def save_results(results_df):
    """Save calculation results (unchanged)"""
    print("\n=== Saving results ===")

    output_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_results.csv")
    results_df.to_csv(output_csv, index=False)
    print(f"  Detailed results: {output_csv}")

    summary_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_summary.csv")
    summary = results_df[['Group', 'Promoter_Count', 'Total_pG4', '%_in_Promoter', 'Fold_Enrichment']].copy()
    summary.to_csv(summary_csv, index=False)
    print(f"  Summary results: {summary_csv}")

    print("\nEnrichment analysis summary:")
    print("-" * 80)
    print(f"{'Group':<20} {'Promoter':<12} {'Total':<12} {'% in Promoter':<15} {'Fold Enrichment':<15}")
    print("-" * 80)

    for _, row in results_df.iterrows():
        print(f"{row['Group']:<20} {row['Promoter_Count']:<12,} {row['Total_pG4']:<12,} {row['%_in_Promoter']:<15.2f} {row['Fold_Enrichment']:<15.2f}")

def main():
    """Main function"""
    print("=" * 80)
    print("pG4 Promoter Enrichment Analysis (Human) - Horizontal Lollipop Chart (font 16pt, abbreviated y-axis)")
    print("=" * 80)

    promoter_counts = count_pg4_in_promoters()
    results_df = calculate_fold_enrichment(promoter_counts)
    plot_enrichment(results_df)
    compare_with_paper(results_df)
    save_results(results_df)

    print("\n" + "=" * 80)
    print(f"Analysis complete! Results saved in: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
