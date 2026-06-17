#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_human_TSS_oup.py - Plot human TSS region i-Motif enrichment curve (compliant with OUP illustration guidelines)
Fix: increased margins, font embedding settings, ensure PDF displays completely in AI.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# ========== OUP illustration guideline style settings ==========
plt.rcParams.update({
    # Font (base size 8pt, all text >=7pt)
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,

    # Line widths (strictly controlled within 0.25-1 pt)
    'lines.linewidth': 0.8,
    'axes.linewidth': 0.5,
    'grid.linewidth': 0.3,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,

    # PDF font editable
    'pdf.fonttype': 42,
    'ps.fonttype': 42,

    # Other
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.linestyle': ':',
    'grid.alpha': 0.3,
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

# ========== Colorblind-friendly palette ==========
COLOR_TEMPLATE = '#1F77B4'      # Blue
COLOR_NONTEMPLATE = '#D62728'   # Red

# ========== Parameter settings ==========
DATA_DIR = "primate_TSS_TES_enrichment_results"
SPECIES = "Human"
REGION = "TSS"
SMOOTH_SIGMA = 60
WINDOW_SIZE = 1000
OUTPUT_DIR = "human_TSS_figure"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_human_tss_data():
    """Load human TSS data"""
    file_path = os.path.join(DATA_DIR, SPECIES, f"{REGION}_results.tsv")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Data file not found: {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    if 'position' not in df.columns:
        raise ValueError("Missing 'position' column in data file")
    df = df[(df['position'] >= -WINDOW_SIZE) & (df['position'] <= WINDOW_SIZE)].copy()
    return df

def apply_smoothing(df, sigma):
    """Apply Gaussian smoothing to enrichment columns"""
    if sigma > 0:
        df['template_enrich_smooth'] = gaussian_filter1d(df['template_enrich'], sigma=sigma)
        df['non_template_enrich_smooth'] = gaussian_filter1d(df['non_template_enrich'], sigma=sigma)
    else:
        df['template_enrich_smooth'] = df['template_enrich']
        df['non_template_enrich_smooth'] = df['non_template_enrich']
    return df

def plot_human_tss(df):
    """Plot human TSS enrichment curve"""
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    # Draw two curves
    ax.plot(df['position'], df['template_enrich_smooth'],
            label='Template strand', color=COLOR_TEMPLATE, linewidth=0.8)
    ax.plot(df['position'], df['non_template_enrich_smooth'],
            label='Non-template strand', color=COLOR_NONTEMPLATE, linewidth=0.8)

    # Mark TSS position
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.5, alpha=0.7, label='TSS')

    # Axis range
    ax.set_xlim(-WINDOW_SIZE, WINDOW_SIZE)
    ax.set_xticks(np.arange(-1000, 1001, 200))

    # Calculate y-axis range, avoid invalid data
    y_min = min(df['template_enrich_smooth'].min(), df['non_template_enrich_smooth'].min())
    y_max = max(df['template_enrich_smooth'].max(), df['non_template_enrich_smooth'].max())
    if np.isnan(y_min) or np.isnan(y_max):
        raise ValueError("Smoothed data contains NaN, please check raw data")
    y_range = y_max - y_min
    ax.set_ylim(max(0, y_min - 0.1*y_range), y_max + 0.1*y_range)

    # Labels
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Normalized enrichment")
    ax.set_title("Human i-Motif enrichment around TSS", fontweight='normal')

    # Legend
    ax.legend(loc='upper right', frameon=False)

    # Grid lines
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)

    plt.tight_layout()
    return fig

def save_figure(fig, base_name):
    """Save as PDF (AI-editable) and TIFF (print)"""
    pdf_path = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
    tiff_path = os.path.join(OUTPUT_DIR, f"{base_name}.tiff")

    # PDF: increase margins, ensure content is not cropped
    fig.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
    print(f"OK PDF saved: {pdf_path}")

    # TIFF: 600dpi, LZW compression
    fig.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight', pad_inches=0.1,
                pil_kwargs={'compression': 'tiff_lzw'})
    print(f"OK TIFF saved: {tiff_path}")

def main():
    print("="*60)
    print("Plotting human TSS region i-Motif enrichment curve (compliant with OUP guidelines)")
    print("="*60)

    try:
        # 1. Load data
        print(f"Loading {SPECIES} {REGION} data...")
        df_raw = load_human_tss_data()
        print(f"  Data points: {len(df_raw)}")
        print(f"  Position range: {df_raw['position'].min()} ~ {df_raw['position'].max()} bp")

        # 2. Apply smoothing
        if SMOOTH_SIGMA > 0:
            print(f"Applying Gaussian smoothing, sigma = {SMOOTH_SIGMA}")
            df = apply_smoothing(df_raw, SMOOTH_SIGMA)
        else:
            df = apply_smoothing(df_raw, 0)

        # Check if smoothed data is valid
        if df['template_enrich_smooth'].isnull().all() or df['non_template_enrich_smooth'].isnull().all():
            raise ValueError("Smoothed data is all NaN")

        # 3. Plot
        print("Generating figure...")
        fig = plot_human_tss(df)

        # 4. Save
        save_figure(fig, "human_TSS_enrichment_oup")
        plt.close(fig)

        print("\nComplete!")
        print(f"Output directory: {os.path.abspath(OUTPUT_DIR)}")

    except Exception as e:
        print(f"\nError: {e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())
