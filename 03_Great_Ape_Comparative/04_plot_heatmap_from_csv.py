#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Heatmap plotting script (OUP adapted, vertical narrow version, colorbar flush with bottom)
Read i-Motif density data from CSV file, generate editable PDF vector graphics.
Colorbar placed horizontally at the bottom, flush with the heatmap; figure width narrowed to emphasize vertical layout.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
import os
import sys

# ==================== OUP Style Settings ====================
def set_oup_style():
    """
    Set plot style compliant with OUP illustration guidelines:
    - Font: Arial/Helvetica, TrueType embedding (PDF editable)
    - Font size: base 12pt, other elements scaled appropriately
    - Lines: moderate thickness
    - Colormap: viridis (scientifically standard)
    """
    # Font embedding settings (ensures PDF text is editable)
    rcParams['pdf.fonttype'] = 42        # TrueType font
    rcParams['ps.fonttype'] = 42         # PostScript TrueType

    # Font family
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['font.size'] = 12           # Base font size

    # Figure size: width reduced to 5.5 inches for vertical layout; height 11 inches to accommodate all chromosomes
    rcParams['figure.figsize'] = (5.5, 11)
    rcParams['figure.dpi'] = 300

    # Lines and markers
    rcParams['lines.linewidth'] = 1.0
    rcParams['lines.markersize'] = 4

    # Axes
    rcParams['axes.linewidth'] = 0.8
    rcParams['axes.labelpad'] = 6
    rcParams['axes.titlepad'] = 15
    rcParams['axes.labelsize'] = 14      # Axis label font size

    # Ticks
    rcParams['xtick.major.width'] = 0.8
    rcParams['ytick.major.width'] = 0.8
    rcParams['xtick.minor.width'] = 0.6
    rcParams['ytick.minor.width'] = 0.6
    rcParams['xtick.labelsize'] = 10     # X-axis tick label size (slightly reduced for narrow width)
    rcParams['ytick.labelsize'] = 11     # Y-axis tick label size

    # Legend
    rcParams['legend.fontsize'] = 11
    rcParams['legend.frameon'] = False

# ==================== Data Loading ====================
def load_heatmap_data(file_path):
    """Load heatmap data from CSV file, return DataFrame"""
    try:
        df = pd.read_csv(file_path, index_col=0)
        print(f"Data file loaded successfully: {file_path}")
        print(f"Data shape: {df.shape}")
        print(f"Row index: {list(df.index)}")
        print(f"Column names: {list(df.columns)}")
        print("\nData preview:")
        print(df.head())
        df = df.astype(float)
        return df
    except Exception as e:
        print(f"Error loading data file: {e}")
        return None

# ==================== Data Statistical Analysis ====================
def analyze_data(data_df):
    """Print data statistics"""
    print("\n" + "=" * 50)
    print("Data Statistical Analysis:")
    print("=" * 50)
    print(f"Data shape: {data_df.shape}")
    print(f"Rows (chromosomes): {len(data_df.index)}")
    print(f"Columns (species): {len(data_df.columns)}")

    # Missing values
    missing_counts = data_df.isnull().sum()
    print("\nMissing value statistics:")
    for species, count in missing_counts.items():
        if count > 0:
            print(f"{species}: {count} missing values")

    # Global statistics
    print("\nDensity value statistics (/Mb):")
    print(f"Global minimum: {data_df.min().min():.6f}")
    print(f"Global maximum: {data_df.max().max():.6f}")
    print(f"Global mean: {data_df.mean().mean():.6f}")
    print(f"Global median: {data_df.stack().median():.6f}")
    print(f"Global standard deviation: {data_df.stack().std():.6f}")

    # Species statistics
    print("\nPer-species statistics (/Mb):")
    for species in data_df.columns:
        species_data = data_df[species].dropna()
        if len(species_data) > 0:
            print(f"{species}: mean={species_data.mean():.6f}, max={species_data.max():.6f}")

    # Chromosome ranking
    print("\nChromosome mean density ranking (top 5):")
    mean_by_chromosome = data_df.mean(axis=1)
    top_chromosomes = mean_by_chromosome.sort_values(ascending=False).head(5)
    for idx, (chr_name, density) in enumerate(top_chromosomes.items(), 1):
        print(f"  {idx}. Chr {chr_name}: {density:.6f}")

# ==================== Heatmap Plotting ====================
def plot_heatmap_from_data(data_df, output_dir, filename_prefix="iMotif_density_heatmap", annotate=True):
    """
    Plot heatmap and save as editable PDF
    annotate: whether to display values in cells
    Colorbar placed horizontally at the bottom, flush with heatmap (pad=0.03)
    Figure width narrowed to emphasize vertical layout
    """
    set_oup_style()
    os.makedirs(output_dir, exist_ok=True)

    # Prepare data
    data = data_df.copy()

    # Create figure
    fig, ax = plt.subplots(figsize=(5.5, 11))  # Narrow vertical dimensions

    # Colormap
    cmap = 'viridis'

    # Annotation format
    if annotate:
        annot = True
        annot_kws = {
            'size': 9,                 # Cell number font size (slightly reduced for narrow cells)
            'color': 'white',
            'fontweight': 'bold'
        }
        fmt = '.3f'
    else:
        annot = False
        annot_kws = None
        fmt = None

    # Plot heatmap
    heatmap = sns.heatmap(
        data,
        cmap=cmap,
        linewidths=0.5,
        linecolor='white',
        square=False,                   # Not forced to square; naturally stretched vertically
        annot=annot,
        fmt=fmt,
        annot_kws=annot_kws,
        cbar_kws={
            'orientation': 'horizontal',   # Horizontal colorbar
            'location': 'bottom',          # Place at bottom
            'label': 'i-Motif density (/Mb)',
            'shrink': 0.7,                 # Colorbar length ratio (slightly narrower than heatmap width)
            'pad': 0.03,                   # Key modification: spacing greatly reduced, flush with heatmap
            'ticks': plt.MaxNLocator(5)
        },
        ax=ax
    )

    # Colorbar label font
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('i-Motif density (/Mb)', fontsize=14)
    cbar.ax.tick_params(labelsize=11)

    # Y-axis labels (chromosomes)
    y_labels = []
    for label in data_df.index:
        if label in ['2a', '2b']:
            y_labels.append(f'Chr {label}')
        else:
            y_labels.append(f'Chr {label}')
    ax.set_yticklabels(y_labels, rotation=0, fontsize=11)

    # X-axis labels (species) mapping
    species_names = {
        'bonobo': 'Bonobo',
        'chimp': 'Chimpanzee',
        'human': 'Human',
        'gorilla': 'Gorilla',
        'sumatran': 'S. orangutan',
        'bornean': 'B. orangutan'
    }
    x_labels = []
    for col in data_df.columns:
        if col in species_names:
            x_labels.append(species_names[col])
        else:
            x_labels.append(col)
    # Keep horizontal, font size 10
    ax.set_xticklabels(x_labels, rotation=0, ha='center', fontsize=10)

    # Adjust layout
    plt.tight_layout()

    # Save as editable PDF
    suffix = "_with_values" if annotate else "_color_only"
    pdf_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', format='pdf')

    # Save PNG preview
    png_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight', format='png')

    print(f"Editable PDF saved to: {pdf_path}")
    print(f"PNG preview saved to: {png_path}")

    plt.close(fig)
    return heatmap

# ==================== Main Program ====================
def main():
    print("=" * 60)
    print("i-Motif Density Heatmap Generator (vertical narrow, colorbar flush with bottom)")
    print("=" * 60)

    # File path
    data_file = "heatmap_data_fixed.csv"
    output_dir = "./heatmap_output_narrow"

    if not os.path.exists(data_file):
        print(f"Error: Data file '{data_file}' not found!")
        sys.exit(1)

    heatmap_df = load_heatmap_data(data_file)
    if heatmap_df is None:
        sys.exit(1)

    analyze_data(heatmap_df)

    os.makedirs(output_dir, exist_ok=True)

    print("\nPlotting heatmap with values...")
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=True)

    print("\nPlotting heatmap without values...")
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=False)

    print("\n" + "=" * 60)
    print("Heatmap generation complete!")
    print(f"Output directory: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()
