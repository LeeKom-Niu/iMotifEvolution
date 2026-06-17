#!/usr/bin/env python3
"""
figure2_zscore_distribution_BD.py - Generate robust Z-score kernel density estimation distribution plots (BD version)
Modification: Remove Viral category from the first plot (four-category becomes three-category)
Contains:
B. Three-category Z-score distribution (Bacteria, Archaea, Eukaryota)
D. Eukaryotic category Z-score distribution (Protozoa first)
Keep original style, all categories with shading
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, norm, gaussian_kde
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.signal import savgol_filter
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import json

# Set clean style
def set_clean_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'font.size': 9,
        'pdf.fonttype': 42,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })

# Color scheme (consistent with original, but without viral color)
three_category_colors = {
    'Bacteria': '#2CA02C',     # Green
    'Archaea': '#1F77B4',      # Blue
    'Eukaryota': '#D62728',    # Red
}

# Viral color (reserved, but not used in the first plot)
viral_color = '#17BECF'        # Cyan

nine_category_colors = {
    'Bacteria': '#2CA02C',     # Green
    'Archaea': '#1F77B4',      # Blue
    'Viral': '#17BECF',        # Cyan
    'Fungi': '#9467BD',        # Purple
    'Plant': '#8C564B',        # Brown
    'Invertebrate': '#E377C2', # Pink
    'Protozoa': '#BCBD22',     # Yellow-green
    'Vertebrate Other': '#7F7F7F',  # Gray
    'Mammalian': '#D62728',    # Red (mammals stay red)
}

class ZScoreDistributionBD:
    def __init__(self, simulation_file, real_genome_file, z_threshold=1.96):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.z_threshold = z_threshold

        # Three categories (without viral)
        self.three_categories = ['Bacteria', 'Archaea', 'Eukaryota']

        # Nine categories (kept as before)
        self.nine_categories = [
            'Bacteria', 'Archaea', 'Viral',
            'Fungi', 'Plant', 'Invertebrate',
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]

        # Eukaryotic subgroups (Protozoa first, keep original style)
        self.eukaryotic_categories = [
            'Protozoa',     # Protozoa placed first
            'Fungi',
            'Plant',
            'Invertebrate',
            'Vertebrate Other',
            'Mammalian'     # Mammals keep red
        ]

        # Data
        self.simulation_data = None
        self.real_data = None
        self.real_data_3cat = None  # Three-category data (modified: 3cat instead of 4cat)
        self.real_data_9cat = None  # Nine-category data
        self.background_model = None

    def load_and_preprocess_data(self):
        """Load and preprocess data"""
        print("=== Loading Data ===")

        # Load simulation data (theoretical background)
        self.simulation_data = pd.read_csv(self.simulation_file)
        print(f"Theoretical background data: {len(self.simulation_data)} rows")

        # Load real genome data
        print(f"Loading real genome data: {self.real_genome_file}")
        self.real_data = pd.read_csv(self.real_genome_file, sep='\t')

        # Standardize column names
        column_mapping = {}
        for col in self.real_data.columns:
            col_lower = col.lower()
            if 'genome' in col_lower or 'species' in col_lower:
                column_mapping['Genome'] = col
            elif 'class' in col_lower or 'domain' in col_lower:
                column_mapping['Classification'] = col
            elif 'gc' in col_lower:
                column_mapping['GC_content'] = col
            elif 'density' in col_lower:
                column_mapping['Genomic_density'] = col

        for new_name, old_name in column_mapping.items():
            if old_name in self.real_data.columns:
                self.real_data = self.real_data.rename(columns={old_name: new_name})

        # Clean data
        if 'Classification' in self.real_data.columns:
            self.real_data = self.real_data[self.real_data['Classification'] != 'Classification']
            self.real_data['Classification'] = self.real_data['Classification'].str.strip()

            # Replace Vertebrate Mammalian with Mammalian
            self.real_data['Classification'] = self.real_data['Classification'].replace(
                'Vertebrate Mammalian', 'Mammalian'
            )

        # Convert data types
        self.real_data['GC_content'] = pd.to_numeric(self.real_data['GC_content'], errors='coerce')
        self.real_data['Genomic_density'] = pd.to_numeric(self.real_data['Genomic_density'], errors='coerce')

        # Convert GC content from percentage to decimal
        if self.real_data['GC_content'].max() > 1:
            self.real_data['GC_content'] = self.real_data['GC_content'] / 100.0

        # Create three-category data (modified: 3cat instead of 4cat, removing Viral)
        self.real_data_3cat = self.real_data.copy()

        # Define eukaryotic subgroups
        eukaryotic_subcategories = [
            'Fungi', 'Plant', 'Invertebrate',
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]

        # Map eukaryotic subgroups to Eukaryota
        eukaryotic_mask = self.real_data_3cat['Classification'].isin(eukaryotic_subcategories)
        self.real_data_3cat.loc[eukaryotic_mask, 'Classification_3cat'] = 'Eukaryota'

        # Keep other categories unchanged (only Bacteria and Archaea, excluding Viral)
        other_categories = ['Bacteria', 'Archaea']  # Modified: removed Viral
        for cat in other_categories:
            mask = self.real_data_3cat['Classification'] == cat
            self.real_data_3cat.loc[mask, 'Classification_3cat'] = cat

        # Note: Viral data is excluded from three-category classification

        # Create nine-category data (kept as before, includes Viral)
        self.real_data_9cat = self.real_data[self.real_data['Classification'].isin(self.nine_categories)].copy()
        self.real_data_9cat['Classification_9cat'] = self.real_data_9cat['Classification']

        print(f"Loaded {len(self.real_data):,} real genomes")
        print(f"Three-category data (no viral): {len(self.real_data_3cat.dropna(subset=['Classification_3cat'])):,} genomes")
        print(f"Nine-category data (includes viral): {len(self.real_data_9cat):,} genomes")

        # Count by category
        if 'Classification_3cat' in self.real_data_3cat.columns:
            three_cat_counts = self.real_data_3cat['Classification_3cat'].value_counts()
            print(f"\nThree-category statistics:")
            for cat, count in three_cat_counts.items():
                print(f"  {cat}: {count:,}")

        if 'Classification_9cat' in self.real_data_9cat.columns:
            nine_cat_counts = self.real_data_9cat['Classification_9cat'].value_counts()
            print(f"\nNine-category statistics:")
            for cat, count in nine_cat_counts.items():
                if cat in self.eukaryotic_categories:  # Only print eukaryotic categories
                    print(f"  {cat}: {count:,}")

        return True

    def build_background_model(self):
        """Build background model from simulation data"""
        print("\n=== Building Background Model ===")

        # Find column names
        possible_gc_columns = ['GC_content', 'GC', 'gc_content', 'gc']
        possible_density_columns = ['Density_IM_per_Mb', 'density_IM_per_Mb', 'density', 'Density', 'median_density']

        gc_col = None
        density_col = None

        for col in possible_gc_columns:
            if col in self.simulation_data.columns:
                gc_col = col
                break

        for col in possible_density_columns:
            if col in self.simulation_data.columns:
                density_col = col
                break

        if not gc_col or not density_col:
            # Try to find columns by pattern matching
            for col in self.simulation_data.columns:
                if 'gc' in col.lower():
                    gc_col = col
                if 'density' in col.lower() or 'median' in col.lower():
                    density_col = col

        if not gc_col or not density_col:
            raise ValueError(f"Cannot find GC and density columns in simulation data. Available columns: {list(self.simulation_data.columns)}")

        print(f"Using GC column: '{gc_col}', Density column: '{density_col}'")

        # Extract data
        background_data = self.simulation_data[[gc_col, density_col]].copy()
        background_data = background_data.rename(columns={gc_col: 'GC_content', density_col: 'density'})

        # Convert GC content
        if background_data['GC_content'].max() > 1:
            background_data['GC_content'] = background_data['GC_content'] / 100.0

        # Remove NaN values
        background_data = background_data.dropna()

        # Sort by GC content
        background_data = background_data.sort_values('GC_content')

        # Extract values
        gc_values = background_data['GC_content'].values
        density_values = background_data['density'].values

        # Sort by GC content
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]

        # Create background model
        self.background_model = interp1d(
            gc_values, density_values,
            kind='linear',
            bounds_error=False,
            fill_value=(density_values[0], density_values[-1])
        )

        return background_data

    def calculate_deviations(self):
        """Calculate robust Z-score"""
        print(f"\n=== Calculating Robust Z-score ===")

        # Calculate deviation for three-category data
        self.real_data_3cat['expected_density'] = self.background_model(self.real_data_3cat['GC_content'])
        deviations_3cat = self.real_data_3cat['Genomic_density'] - self.real_data_3cat['expected_density']

        # Robust Z-score (using overall distribution of three-category data)
        median_deviation = deviations_3cat.median()
        mad = (deviations_3cat - median_deviation).abs().median()

        if mad > 0:
            self.real_data_3cat['robust_z_score'] = (deviations_3cat - median_deviation) / (1.4826 * mad)
            print(f"Robust Z-score calculation: median_deviation={median_deviation:.2f}, MAD={mad:.2f}")

            # Apply same Z-score calculation to nine-category data
            self.real_data_9cat['expected_density'] = self.background_model(self.real_data_9cat['GC_content'])
            deviations_9cat = self.real_data_9cat['Genomic_density'] - self.real_data_9cat['expected_density']
            self.real_data_9cat['robust_z_score'] = (deviations_9cat - median_deviation) / (1.4826 * mad)
        else:
            self.real_data_3cat['robust_z_score'] = 0
            self.real_data_9cat['robust_z_score'] = 0

        return True

    def create_figure_2_BD(self):
        """Create figure 2: Robust Z-score kernel density estimation distribution (BD version) - one row, two panels"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

        # 2B: Three-category Z-score distribution (modified: four-category changed to three-category)
        self._plot_zscore_by_3categories(axes[0])

        # 2D: Eukaryotic category Z-score distribution (Protozoa first)
        self._plot_zscore_by_eukaryotic_categories_protozoa_first(axes[1])

        #plt.suptitle('Figure : Robust Z-score Kernel Density Estimation Distribution', fontsize=14, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig

    def _plot_zscore_by_3categories(self, ax):
        """Plot three-category Z-score distribution (without viral)"""
        # Collect data
        categories = []
        zscore_data = []

        for category in self.three_categories:
            cat_data = self.real_data_3cat[self.real_data_3cat['Classification_3cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 10:
                categories.append(category)
                zscore_data.append(z_scores)

        if len(categories) >= 2:
            # Create KDE curves
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)

            # First draw filled areas (background)
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)

            # Then draw curves (foreground)
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2.5, label=category)

            # Add threshold lines (enriched and depleted annotations)
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6,
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6,
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')

            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            #ax.set_title('A. Z-score Distribution by Three Categories')

            # Legend: only show categories and threshold lines
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            #ax.set_title('A. Z-score Distribution by Three Categories')

    def _plot_zscore_by_eukaryotic_categories_protozoa_first(self, ax):
        """Plot eukaryotic category Z-score distribution (Protozoa first, keep original style)"""
        # Use custom order: Protozoa first
        categories_ordered = self.eukaryotic_categories.copy()

        # Collect data
        categories = []
        zscore_data = []

        for category in categories_ordered:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) >= 3:  # Lower threshold for eukaryotic categories
                categories.append(category)
                zscore_data.append(z_scores)

        if len(categories) >= 2:
            # Create KDE curves
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)

            # First draw filled areas (background)
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)

            # Then draw curves (foreground)
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2, label=category)

            # Add threshold lines (enriched and depleted annotations)
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6,
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6,
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')

            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            #ax.set_title('B. Z-score Distribution by Eukaryotic Categories')

            # Legend: only show categories and threshold lines (two columns)
            ax.legend(fontsize=8, loc='upper right', ncol=2)
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            #ax.set_title('B. Z-score Distribution by Eukaryotic Categories')

    def perform_statistical_analysis(self):
        """Perform statistical analysis"""
        results = {}

        print("\n=== Statistical Analysis ===")

        # 1. Three-category Z-score statistics (modified: 3cat instead of 4cat)
        print("\n1. Three-category Z-score statistics (no viral):")
        for category in self.three_categories:
            cat_data = self.real_data_3cat[self.real_data_3cat['Classification_3cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 0:
                mean_val = z_scores.mean()
                std_val = z_scores.std()
                median_val = z_scores.median()
                n_val = len(z_scores)
                results[f'{category}_zscore'] = {
                    'mean': float(mean_val),
                    'std': float(std_val),
                    'median': float(median_val),
                    'n': n_val
                }
                print(f"   {category}: mean={mean_val:.2f}, std={std_val:.2f}, median={median_val:.2f}, n={n_val:,}")

        # 2. Eukaryotic Z-score statistics (in display order)
        print("\n2. Eukaryotic Z-score statistics (Protozoa first):")
        for category in self.eukaryotic_categories:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 0:
                mean_val = z_scores.mean()
                median_val = z_scores.median()
                n_val = len(z_scores)
                results[f'{category}_zscore'] = {
                    'mean': float(mean_val),
                    'median': float(median_val),
                    'n': n_val
                }
                print(f"   {category}: mean={mean_val:.2f}, median={median_val:.2f}, n={n_val:,}")

        return results

    def export_results(self, output_dir):
        """Export results"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== Exporting Results to {output_dir} ===")

        # Save three-category Z-score data (modified: 3cat instead of 4cat)
        results_file_3cat = output_dir / "zscore_data_3categories.csv"
        self.real_data_3cat[['Genome', 'Classification_3cat', 'GC_content',
                           'Genomic_density', 'expected_density', 'robust_z_score']].to_csv(
                               results_file_3cat, index=False)
        print(f"Three-category Z-score data (no viral): {results_file_3cat}")

        # Save eukaryotic Z-score data (in display order)
        eukaryotic_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'].isin(self.eukaryotic_categories)]
        # Sort by display order
        eukaryotic_data['Classification_9cat'] = pd.Categorical(
            eukaryotic_data['Classification_9cat'],
            categories=self.eukaryotic_categories,
            ordered=True
        )
        eukaryotic_data = eukaryotic_data.sort_values('Classification_9cat')

        results_file_eukaryotic = output_dir / "zscore_data_eukaryotic_categories.csv"
        eukaryotic_data[['Genome', 'Classification_9cat', 'GC_content',
                       'Genomic_density', 'expected_density', 'robust_z_score']].to_csv(
                           results_file_eukaryotic, index=False)
        print(f"Eukaryotic Z-score data (Protozoa first): {results_file_eukaryotic}")

        # Save statistical results
        stats_results = self.perform_statistical_analysis()
        if stats_results:
            stats_file = output_dir / "zscore_statistics.json"
            with open(stats_file, 'w') as f:
                json.dump(stats_results, f, indent=2, ensure_ascii=False)
            print(f"Z-score statistics: {stats_file}")

        return output_dir

def main():
    """Main function"""
    set_clean_style()

    # File paths
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"

    print("="*80)
    print("Robust Z-score Kernel Density Estimation Distribution Analysis (BD version)")
    print(f"Z threshold: {1.96} (95% confidence interval)")
    print("Main modification: Removed Viral category from the first plot (four-category)")
    print("Now showing:")
    print("  Panel A: Three-category Z-score distribution (Bacteria, Archaea, Eukaryota)")
    print("  Panel B: Eukaryotic Z-score distribution (Protozoa first)")
    print("="*80)

    # Initialize analyzer
    analyzer = ZScoreDistributionBD(simulation_file, real_genome_file, z_threshold=1.96)

    try:
        # Execute analysis
        analyzer.load_and_preprocess_data()
        analyzer.build_background_model()
        analyzer.calculate_deviations()

        # Create BD version figure
        print("\n=== Creating BD Version Figure ===")
        print("Creating Figure 2 (BD version): Robust Z-score KDE distribution...")
        fig = analyzer.create_figure_2_BD()

        # Export results
        output_dir = Path(base_dir) / "figure2_BD_analysis"
        analyzer.export_results(output_dir)

        # Save figure
        print("\n=== Saving Figure ===")
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.png", dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.pdf", bbox_inches='tight')
        print(f"Figure 2 (BD version): {output_dir}/figure2_zscore_distribution_BD.png/pdf")

        # Display statistical summary
        print("\n" + "="*80)
        print("Statistical Summary")
        print("="*80)

        stats = analyzer.perform_statistical_analysis()

        print("\nThree-category Z-score summary (no viral):")
        for category in analyzer.three_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
                print(f"  {category}:")
                print(f"    Mean: {cat_stats['mean']:.2f}, Median: {cat_stats['median']:.2f}")
                print(f"    Std: {cat_stats['std']:.2f}, N: {cat_stats['n']:,}")

        print("\nEukaryotic Z-score summary (in display order):")
        for category in analyzer.eukaryotic_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
                print(f"  {category}:")
                print(f"    Mean: {cat_stats['mean']:.2f}, Median: {cat_stats['median']:.2f}")
                print(f"    N: {cat_stats['n']:,}")

        print("\n" + "="*80)
        print("Analysis complete! BD version figure and results saved.")
        print("="*80)

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
