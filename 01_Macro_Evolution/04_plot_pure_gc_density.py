#!/usr/bin/env python3
"""
plot_pure_gc_density.py - Pure GC content vs i-motif density relationship plot
Generate Nature journal style core figures
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

# Set Nature journal pure style
def set_nature_pure_style():
    """Set Nature journal pure plot style"""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.minor.width': 0.6,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'figure.autolayout': False,
        'axes.grid': True,
        'grid.alpha': 0.2,
        'grid.linewidth': 0.5
    })

def load_data():
    """Load data"""
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    combined_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"

    if not combined_file.exists():
        raise FileNotFoundError(f"Data file not found: {combined_file}")

    df = pd.read_csv(combined_file)

    # Calculate summary statistics
    summary = df.groupby('GC_content').agg({
        'Density_IM_per_Mb': ['mean', 'std', 'count', 'sem'],
        'Actual_C_content': 'mean',
        'Actual_G_content': 'mean'
    }).round(4)

    summary.columns = ['Density_mean', 'Density_std', 'Replicate_count', 'Density_sem',
                      'Actual_C_mean', 'Actual_G_mean']
    summary = summary.reset_index()

    # Calculate confidence intervals (95% CI)
    summary['ci_lower'] = summary['Density_mean'] - 1.96 * summary['Density_sem']
    summary['ci_upper'] = summary['Density_mean'] + 1.96 * summary['Density_sem']

    return df, summary

def create_confidence_interval_plot(df, summary):
    """Create scatter plot with confidence interval"""
    fig, ax = plt.subplots(figsize=(4.0, 3.0))  # Nature standard size

    # Use viridis colormap
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary)))

    # Plot all replicate data points (with slight jitter)
    for idx, (_, row) in enumerate(summary.iterrows()):
        gc_content = row['GC_content']
        subset = df[df['GC_content'] == gc_content]

        # Add slight horizontal jitter to avoid overlap
        jitter = np.random.normal(0, 0.0015, len(subset))
        ax.scatter(gc_content + jitter, subset['Density_IM_per_Mb'],
                  alpha=0.5, s=25, color=colors[idx], edgecolors='white', linewidth=0.3,
                  zorder=2)

    # Plot mean curve
    ax.plot(summary['GC_content'], summary['Density_mean'],
           color='#2E4057', linewidth=2, marker='o', markersize=6,
           markerfacecolor='white', markeredgecolor='#2E4057', markeredgewidth=1.5,
           zorder=3)

    # Plot 95% confidence interval
    ax.fill_between(summary['GC_content'],
                   summary['ci_lower'], summary['ci_upper'],
                   alpha=0.15, color='#2E4057', zorder=1)

    # Set axis labels
    ax.set_xlabel('GC content', fontsize=10)
    ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10)

    # Set grid
    ax.grid(True, alpha=0.2, linewidth=0.5, linestyle='-', zorder=0)

    # Set ticks
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xticklabels(['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])

    # Set axis range
    ax.set_xlim(0.08, 0.92)

    # Add confidence interval annotation
    ax.text(0.02, 0.98, 'Mean ± 95% CI', transform=ax.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, pad=0.3))

    plt.tight_layout()
    return fig

def create_fitted_curve_plot(df, summary):
    """Create plot with fitted curve"""
    fig, ax = plt.subplots(figsize=(4.0, 3.0))  # Nature standard size

    # Use viridis colormap
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary)))

    # Plot all replicate data points
    for idx, (_, row) in enumerate(summary.iterrows()):
        gc_content = row['GC_content']
        subset = df[df['GC_content'] == gc_content]

        jitter = np.random.normal(0, 0.0015, len(subset))
        ax.scatter(gc_content + jitter, subset['Density_IM_per_Mb'],
                  alpha=0.4, s=25, color=colors[idx], edgecolors='white', linewidth=0.3,
                  zorder=2)

    # Plot mean curve
    ax.plot(summary['GC_content'], summary['Density_mean'],
           color='#2E4057', linewidth=2, marker='o', markersize=6,
           markerfacecolor='white', markeredgecolor='#2E4057', markeredgewidth=1.5,
           zorder=3, label='Mean')

    # Try multiple fits, select the best
    x = summary['GC_content'].values
    y = summary['Density_mean'].values

    # 1. Quadratic polynomial fit
    try:
        coeff_quad = np.polyfit(x, y, 2)
        poly_quad = np.poly1d(coeff_quad)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit_quad = poly_quad(x_fit)
        r2_quad = 1 - np.sum((y - poly_quad(x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_quad = -np.inf

    # 2. Cubic polynomial fit
    try:
        coeff_cubic = np.polyfit(x, y, 3)
        poly_cubic = np.poly1d(coeff_cubic)
        y_fit_cubic = poly_cubic(x_fit)
        r2_cubic = 1 - np.sum((y - poly_cubic(x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_cubic = -np.inf

    # 3. Exponential fit
    try:
        # y = a * exp(b*x)
        log_y = np.log(y)
        coeff_exp = np.polyfit(x, log_y, 1)
        a = np.exp(coeff_exp[1])
        b = coeff_exp[0]
        y_fit_exp = a * np.exp(b * x_fit)
        r2_exp = 1 - np.sum((y - a*np.exp(b*x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_exp = -np.inf

    # Select fit with highest R²
    fits = {
        'Quadratic': (r2_quad, y_fit_quad, '#D62728'),
        'Cubic': (r2_cubic, y_fit_cubic, '#FF7F0E'),
        'Exponential': (r2_exp, y_fit_exp, '#2CA02C')
    }

    best_fit_name = max(fits, key=lambda k: fits[k][0])
    best_r2, best_y_fit, best_color = fits[best_fit_name]

    # Plot best fit curve
    ax.plot(x_fit, best_y_fit, '--', color=best_color, linewidth=2,
           alpha=0.8, zorder=4, label=f'{best_fit_name} fit')

    # Add R² value
    if best_r2 > 0:
        ax.text(0.98, 0.98, f'R² = {best_r2:.3f}\n({best_fit_name})',
                transform=ax.transAxes, fontsize=8, verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, pad=0.3))

    # Set axis labels
    ax.set_xlabel('GC content', fontsize=10)
    ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10)

    # Set grid
    ax.grid(True, alpha=0.2, linewidth=0.5, linestyle='-', zorder=0)

    # Set ticks
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xticklabels(['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])

    # Set axis range
    ax.set_xlim(0.08, 0.92)

    # Add legend
    ax.legend(loc='upper left', frameon=True, fancybox=True,
              framealpha=0.8, edgecolor='none', fontsize=8)

    plt.tight_layout()
    return fig, best_fit_name, best_r2

def create_correlation_heatmap(df):
    """Create correlation heatmap (optional)"""
    # Calculate Spearman and Pearson correlations
    corr_cols = ['Actual_GC_content', 'Actual_C_content', 'Actual_G_content', 'Density_IM_per_Mb']

    # Pearson correlation
    pearson_corr = df[corr_cols].corr(method='pearson')

    # Spearman correlation
    spearman_corr = df[corr_cols].corr(method='spearman')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.8))

    # Rename columns for display
    display_names = ['GC', 'C', 'G', 'Density']
    pearson_corr.index = display_names
    pearson_corr.columns = display_names
    spearman_corr.index = display_names
    spearman_corr.columns = display_names

    # Pearson heatmap
    im1 = ax1.imshow(pearson_corr.values, cmap='coolwarm', vmin=-1, vmax=1)
    for i in range(len(pearson_corr)):
        for j in range(len(pearson_corr)):
            ax1.text(j, i, f'{pearson_corr.iloc[i, j]:.2f}',
                    ha='center', va='center',
                    color='white' if abs(pearson_corr.iloc[i, j]) > 0.5 else 'black',
                    fontsize=8)

    ax1.set_xticks(range(len(display_names)))
    ax1.set_yticks(range(len(display_names)))
    ax1.set_xticklabels(display_names)
    ax1.set_yticklabels(display_names)
    ax1.set_title('Pearson Correlation', fontsize=10, pad=10)
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    # Spearman heatmap
    im2 = ax2.imshow(spearman_corr.values, cmap='coolwarm', vmin=-1, vmax=1)
    for i in range(len(spearman_corr)):
        for j in range(len(spearman_corr)):
            ax2.text(j, i, f'{spearman_corr.iloc[i, j]:.2f}',
                    ha='center', va='center',
                    color='white' if abs(spearman_corr.iloc[i, j]) > 0.5 else 'black',
                    fontsize=8)

    ax2.set_xticks(range(len(display_names)))
    ax2.set_yticks(range(len(display_names)))
    ax2.set_xticklabels(display_names)
    ax2.set_yticklabels(display_names)
    ax2.set_title('Spearman Correlation', fontsize=10, pad=10)
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)

    plt.tight_layout()
    return fig

def save_figures(figures_dict, base_dir):
    """Save all figures"""
    output_dir = Path(base_dir) / "figures_pure"
    output_dir.mkdir(exist_ok=True)

    for name, fig in figures_dict.items():
        fig.savefig(output_dir / f"{name}.png", dpi=600, bbox_inches='tight')
        fig.savefig(output_dir / f"{name}.pdf", bbox_inches='tight')
        fig.savefig(output_dir / f"{name}.svg", bbox_inches='tight')
        print(f"Saved: {name}.png/pdf/svg")

    print(f"\nAll pure figures saved to: {output_dir}")

def main():
    """Main function"""
    # Set Nature pure style
    set_nature_pure_style()

    print("=== Pure GC Content vs i-Motif Density Relationship Plot ===")

    # Load data
    try:
        df, summary = load_data()
        print(f"Data loaded successfully, {len(df)} simulations total")
        print(f"GC content gradients: {len(summary)} levels (10%-90%)")

        # Display basic statistics
        print(f"\nDensity Summary Statistics:")
        print(f"  Highest density: {summary['Density_mean'].max():.1f} IM/Mb @ GC{summary.loc[summary['Density_mean'].idxmax(), 'GC_content']*100:.0f}%")
        print(f"  Lowest density: {summary['Density_mean'].min():.1f} IM/Mb @ GC{summary.loc[summary['Density_mean'].idxmin(), 'GC_content']*100:.0f}%")
        print(f"  Mean CV: {((summary['Density_std'] / summary['Density_mean']) * 100).mean():.1f}%")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please run combine_results_GC.py first to merge data")
        return

    # Create figures
    print("\nCreating figures...")

    # 1. Scatter plot with confidence interval
    fig_ci = create_confidence_interval_plot(df, summary)

    # 2. Plot with fitted curve
    fig_fit, fit_name, fit_r2 = create_fitted_curve_plot(df, summary)
    print(f"Best fit: {fit_name}, R² = {fit_r2:.3f}")

    # 3. Optional: Correlation heatmap
    try:
        fig_corr = create_correlation_heatmap(df)
        print("Correlation heatmap created")
        include_corr = True
    except:
        print("Skipping correlation heatmap (possible data issue)")
        include_corr = False

    # Save figures
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    figures = {
        'gc_density_confidence_interval': fig_ci,
        'gc_density_fitted_curve': fig_fit,
    }

    if include_corr:
        figures['correlation_heatmap'] = fig_corr

    save_figures(figures, base_dir)

    # Calculate and display correlations
    print("\n=== Correlation Analysis ===")

    # Pearson correlation (linear relationship)
    pearson_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'], method='pearson')
    pearson_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'], method='pearson')

    # Spearman correlation (monotonic relationship)
    spearman_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'], method='spearman')
    spearman_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'], method='spearman')

    print(f"Pearson correlation:")
    print(f"  GC content vs density: {pearson_gc:.3f}")
    print(f"  C content vs density: {pearson_c:.3f}")

    print(f"\nSpearman correlation:")
    print(f"  GC content vs density: {spearman_gc:.3f}")
    print(f"  C content vs density: {spearman_c:.3f}")

    # Determine which correlation to use
    if abs(spearman_gc - pearson_gc) > 0.1:
        print("\nSuggest using Spearman correlation (data may not be normally distributed)")
    else:
        print("\nPearson and Spearman correlations are close, either can be used")

    # Display figures
    plt.show()

    print("\n=== Complete ===")
    print("Pure figures generated, suitable for publication.")

if __name__ == "__main__":
    main()
