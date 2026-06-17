#!/usr/bin/env python3
"""
figure3_deviation_composition.py - Create Figure 3: Percentage composition of deviation categories across biological groups
Modifications:
1. For Bacteria with very small enrichment proportion, mark with red font above the bar
2. Bar chart directly to top, no need to mark n values
3. Sort from low to high
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
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'figure.dpi': 150,
        'savefig.dpi': 150,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })

# Color scheme
four_category_colors = {
    'Bacteria': '#2CA02C',     # Green
    'Archaea': '#1F77B4',      # Blue
    'Eukaryota': '#D62728',    # Red
}

nine_category_colors = {
    'Bacteria': '#2CA02C',     # Green
    'Archaea': '#1F77B4',      # Blue
    'Fungi': '#9467BD',        # Purple
    'Plant': '#8C564B',        # Brown
    'Invertebrate': '#E377C2', # Pink
    'Protozoa': '#BCBD22',     # Yellow-green
    'Vertebrate Other': '#7F7F7F',  # Gray
    'Mammalian': '#D62728',    # Red
}

# Deviation category colors
deviation_colors = {
    'Enriched': '#D62728',     # Red
    'Normal': '#7F7F7F',       # Gray
    'Depleted': '#1F77B4'      # Blue
}

class GCBackgroundComparison:
    def __init__(self, simulation_file, real_genome_file, z_threshold=1.96):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.z_threshold = z_threshold
        self.confidence_level = (1 - 2 * (1 - norm.cdf(z_threshold))) * 100

        # Four categories (without viral)
        self.four_categories = ['Bacteria', 'Archaea', 'Eukaryota']

        # Nine categories (updated to Mammalian, without viral)
        self.nine_categories = [
            'Bacteria', 'Archaea',
            'Fungi', 'Plant', 'Invertebrate',
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]

        # Eukaryotic subgroups (for mapping to Eukaryota)
        self.eukaryotic_subcategories = [
            'Fungi', 'Plant', 'Invertebrate',
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]

        # Data
        self.simulation_data = None
        self.real_data = None
        self.real_data_4cat = None  # Four-category data
        self.real_data_9cat = None  # Nine-category data
        self.background_model = None
        self.background_stats = None

        # For visualization
        self.smooth_background_x = None
        self.smooth_background_y = None

    def load_and_preprocess_data(self):
        """Load and preprocess data"""
        print("=== Loading Data ===")

        # Load simulation data (theoretical background)
        self.simulation_data = pd.read_csv(self.simulation_file)
        print(f"Theoretical background data: {len(self.simulation_data)} rows")
        print(f"Simulation data columns: {list(self.simulation_data.columns)}")

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
            print("Converted GC content from percentage to decimal")

        # Create four-category data (without viral)
        self.real_data_4cat = self.real_data.copy()

        # Map eukaryotic subgroups to Eukaryota
        eukaryotic_mask = self.real_data_4cat['Classification'].isin(self.eukaryotic_subcategories)
        self.real_data_4cat.loc[eukaryotic_mask, 'Classification_4cat'] = 'Eukaryota'

        # Keep other categories unchanged, but exclude viral
        other_categories = ['Bacteria', 'Archaea']
        for cat in other_categories:
            mask = self.real_data_4cat['Classification'] == cat
            self.real_data_4cat.loc[mask, 'Classification_4cat'] = cat

        # Filter out viral and other unclassified
        self.real_data_4cat = self.real_data_4cat[self.real_data_4cat['Classification_4cat'].isin(self.four_categories)].copy()

        # Create nine-category data (without viral)
        self.real_data_9cat = self.real_data[self.real_data['Classification'].isin(self.nine_categories)].copy()
        self.real_data_9cat['Classification_9cat'] = self.real_data_9cat['Classification']

        print(f"Loaded {len(self.real_data):,} real genomes")
        print(f"Four-category data (no viral): {len(self.real_data_4cat):,} genomes")
        print(f"Nine-category data (no viral): {len(self.real_data_9cat):,} genomes")

        # Print counts per category
        print("\nFour-category counts (no viral):")
        for category in self.four_categories:
            count = (self.real_data_4cat['Classification_4cat'] == category).sum()
            if count > 0:
                print(f"  {category}: {count:,} genomes")

        print("\nNine-category counts (no viral):")
        for category in self.nine_categories:
            count = (self.real_data_9cat['Classification_9cat'] == category).sum()
            if count > 0:
                print(f"  {category}: {count:,} genomes")

        return True

    def build_background_model(self):
        """Build smooth background model from simulation data"""
        print("\n=== Building Background Model ===")

        # Check if simulation data has necessary columns
        print(f"Simulation data columns: {list(self.simulation_data.columns)}")

        # Try different column names
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
            print("Converted simulation GC content from percentage to decimal")

        # Remove NaN values
        background_data = background_data.dropna()

        # Sort by GC content
        background_data = background_data.sort_values('GC_content')

        print(f"Background data points: {len(background_data)}")
        print(f"GC range: {background_data['GC_content'].min():.3f} - {background_data['GC_content'].max():.3f}")
        print(f"Density range: {background_data['density'].min():.3f} - {background_data['density'].max():.3f}")

        # Extract values
        gc_values = background_data['GC_content'].values
        density_values = background_data['density'].values

        # Sort by GC content
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]

        # Create smooth curve with more points
        self.smooth_background_x = np.linspace(gc_values.min(), gc_values.max(), 500)

        # Smooth curve using Savitzky-Golay filter
        if len(gc_values) > 10:
            # Apply Savitzky-Golay filter
            window_length = min(15, len(gc_values) - 1)
            if window_length % 2 == 0:
                window_length -= 1  # Ensure odd

            smoothed_density = savgol_filter(density_values, window_length=window_length, polyorder=3)

            # Use UnivariateSpline for very smooth interpolation
            try:
                # Use small smoothing factor
                spline = UnivariateSpline(gc_values, smoothed_density, s=0.001, k=3)
                self.smooth_background_y = spline(self.smooth_background_x)
                print("Created very smooth background model using UnivariateSpline + Savitzky-Golay")
            except:
                # Fallback to cubic spline
                cs = interp1d(gc_values, smoothed_density, kind='cubic', fill_value='extrapolate')
                self.smooth_background_y = cs(self.smooth_background_x)
                print("Created smooth background model using cubic spline")
        else:
            # Direct cubic spline
            cs = interp1d(gc_values, density_values, kind='cubic', fill_value='extrapolate')
            self.smooth_background_y = cs(self.smooth_background_x)
            print("Created smooth background model using cubic spline")

        # Store for prediction (use linear interpolation to avoid overfitting)
        self.background_model = interp1d(
            gc_values, density_values,
            kind='linear',
            bounds_error=False,
            fill_value=(density_values[0], density_values[-1])
        )

        # Save background statistics
        self.background_stats = background_data

        print(f"Background model built successfully, GC range: {gc_values.min():.3f} - {gc_values.max():.3f}")

        return background_data

    def calculate_deviations(self):
        """Calculate deviation metrics"""
        print(f"\n=== Calculating Deviation Metrics (Z threshold={self.z_threshold}, {self.confidence_level:.1f}% confidence interval) ===")

        # Calculate deviation for four-category data
        self.real_data_4cat['expected_density'] = self.background_model(self.real_data_4cat['GC_content'])
        deviations_4cat = self.real_data_4cat['Genomic_density'] - self.real_data_4cat['expected_density']

        # Robust Z-score (using overall distribution of four-category data)
        median_deviation = deviations_4cat.median()
        mad = (deviations_4cat - median_deviation).abs().median()

        if mad > 0:
            self.real_data_4cat['robust_z_score'] = (deviations_4cat - median_deviation) / (1.4826 * mad)
            print(f"Robust Z-score calculation: median_deviation={median_deviation:.2f}, MAD={mad:.2f}")

            # Apply same Z-score calculation to nine-category data
            self.real_data_9cat['expected_density'] = self.background_model(self.real_data_9cat['GC_content'])
            deviations_9cat = self.real_data_9cat['Genomic_density'] - self.real_data_9cat['expected_density']
            self.real_data_9cat['robust_z_score'] = (deviations_9cat - median_deviation) / (1.4826 * mad)
        else:
            self.real_data_4cat['robust_z_score'] = 0
            self.real_data_9cat['robust_z_score'] = 0

        # Classification (four-category)
        conditions = [
            self.real_data_4cat['robust_z_score'] > self.z_threshold,
            self.real_data_4cat['robust_z_score'] < -self.z_threshold,
            (self.real_data_4cat['robust_z_score'] >= -self.z_threshold) &
            (self.real_data_4cat['robust_z_score'] <= self.z_threshold)
        ]
        choices = ['Enriched', 'Depleted', 'Normal']
        self.real_data_4cat['deviation_category'] = np.select(conditions, choices, default='Normal')

        # Classification (nine-category)
        conditions_9cat = [
            self.real_data_9cat['robust_z_score'] > self.z_threshold,
            self.real_data_9cat['robust_z_score'] < -self.z_threshold,
            (self.real_data_9cat['robust_z_score'] >= -self.z_threshold) &
            (self.real_data_9cat['robust_z_score'] <= self.z_threshold)
        ]
        self.real_data_9cat['deviation_category'] = np.select(conditions_9cat, choices, default='Normal')

        # Statistics
        print("\nFour-category deviation statistics (no viral):")
        category_counts = self.real_data_4cat['deviation_category'].value_counts()
        for category, count in category_counts.items():
            percentage = count / len(self.real_data_4cat) * 100
            print(f"  {category}: {count:,} genomes ({percentage:.1f}%)")

        print("\nNine-category deviation statistics (no viral):")
        category_counts_9 = self.real_data_9cat['deviation_category'].value_counts()
        for category, count in category_counts_9.items():
            percentage = count / len(self.real_data_9cat) * 100
            print(f"  {category}: {count:,} genomes ({percentage:.1f}%)")

        return True

    def create_figure_3_deviation_composition(self):
        """Create figure 3: Percentage composition of deviation categories across biological groups (AC version)"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # 3A: Four-category deviation composition (stacked bar chart)
        self._plot_deviation_composition_4cat_stacked(axes[0])

        # 3B: Nine-category deviation composition (stacked bar chart)
        self._plot_deviation_composition_9cat_stacked(axes[1])

        plt.suptitle('Figure 3: Percentage Composition of Deviation Categories (without Viral)', fontsize=14, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig

    def _plot_deviation_composition_4cat_stacked(self, ax):
        """Plot four-category deviation composition (stacked bar chart, mark enrichment percentage)"""
        # Calculate deviation composition per category
        composition_data = {}

        for category in self.four_categories:
            cat_data = self.real_data_4cat[self.real_data_4cat['Classification_4cat'] == category]
            if len(cat_data) > 0:
                total = len(cat_data)
                composition = {}
                for dev_cat in ['Enriched', 'Normal', 'Depleted']:
                    count = (cat_data['deviation_category'] == dev_cat).sum()
                    composition[dev_cat] = count / total * 100
                composition_data[category] = composition

        if composition_data:
            categories = list(composition_data.keys())
            enriched_pcts = [composition_data[c]['Enriched'] for c in categories]
            normal_pcts = [composition_data[c]['Normal'] for c in categories]
            depleted_pcts = [composition_data[c]['Depleted'] for c in categories]

            x = np.arange(len(categories))
            bottom = np.zeros(len(categories))

            # Plot stacked bar chart
            enriched_bars = ax.bar(x, enriched_pcts, label='Enriched', color=deviation_colors['Enriched'],
                                  alpha=0.8, bottom=bottom)
            bottom += enriched_pcts

            normal_bars = ax.bar(x, normal_pcts, label='Normal', color=deviation_colors['Normal'],
                                alpha=0.8, bottom=bottom)
            bottom += normal_pcts

            depleted_bars = ax.bar(x, depleted_pcts, label='Depleted', color=deviation_colors['Depleted'],
                                  alpha=0.8, bottom=bottom)

            # Mark enrichment percentage on the enriched portion
            for i, (bar, pct) in enumerate(zip(enriched_bars, enriched_pcts)):
                height = bar.get_height()
                if height > 0:
                    # For very small enrichment proportion (e.g., Bacteria: 0.9%), mark above the bar in red
                    if pct < 5:  # Less than 5%
                        # Mark above the bar top, using red font
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height + 1,
                               f'{pct:.0f}%', ha='center', va='bottom', fontsize=9,
                               fontweight='bold', color='red')
                    else:
                        # Mark inside the bar, using white font
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height/2,
                               f'{pct:.0f}%', ha='center', va='center', fontsize=9,
                               fontweight='bold', color='white')

            ax.set_xlabel('Biological Category')
            ax.set_ylabel('Percentage (%)')
            ax.set_title('A. Four-Category Deviation Composition (without Viral)')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.legend(fontsize=8, loc='upper left')

            # Set y-axis range so bars go to the top
            ax.set_ylim(0, 105)  # Leave 5% space for top labels

            ax.grid(True, alpha=0.3, axis='y', linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('A. Four-Category Deviation Composition')

    def _plot_deviation_composition_9cat_stacked(self, ax):
        """Plot nine-category deviation composition (stacked bar chart, mark enrichment percentage, sorted by enrichment from low to high)"""
        # Calculate deviation composition per category
        composition_data = {}

        for category in self.nine_categories:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            if len(cat_data) > 0:
                total = len(cat_data)
                composition = {}
                for dev_cat in ['Enriched', 'Normal', 'Depleted']:
                    count = (cat_data['deviation_category'] == dev_cat).sum()
                    composition[dev_cat] = count / total * 100
                composition_data[category] = composition

        if composition_data:
            # Sort by enrichment proportion from low to high
            categories_sorted = sorted(composition_data.keys(),
                                      key=lambda x: composition_data[x]['Enriched'],
                                      reverse=False)  # reverse=False means low to high

            categories = categories_sorted
            enriched_pcts = [composition_data[c]['Enriched'] for c in categories]
            normal_pcts = [composition_data[c]['Normal'] for c in categories]
            depleted_pcts = [composition_data[c]['Depleted'] for c in categories]

            x = np.arange(len(categories))
            bottom = np.zeros(len(categories))

            # Plot stacked bar chart
            enriched_bars = ax.bar(x, enriched_pcts, label='Enriched', color=deviation_colors['Enriched'],
                                  alpha=0.8, bottom=bottom)
            bottom += enriched_pcts

            normal_bars = ax.bar(x, normal_pcts, label='Normal', color=deviation_colors['Normal'],
                                alpha=0.8, bottom=bottom)
            bottom += normal_pcts

            depleted_bars = ax.bar(x, depleted_pcts, label='Depleted', color=deviation_colors['Depleted'],
                                  alpha=0.8, bottom=bottom)

            # Mark enrichment percentage on the enriched portion
            for i, (bar, pct) in enumerate(zip(enriched_bars, enriched_pcts)):
                height = bar.get_height()
                if height > 0:
                    # For very small enrichment proportion, mark above the bar in red
                    if pct < 5:  # Less than 5%
                        # Mark above the bar top, using red font
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height + 1,
                               f'{pct:.0f}%', ha='center', va='bottom', fontsize=9,
                               fontweight='bold', color='red')
                    else:
                        # Mark inside the bar, using white font
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height/2,
                               f'{pct:.0f}%', ha='center', va='center', fontsize=9,
                               fontweight='bold', color='white')

            ax.set_xlabel('Biological Category')
            ax.set_ylabel('Percentage (%)')
            ax.set_title('B. Nine-Category Deviation Composition (Sorted by Enrichment from Low to High, without Viral)')
            ax.set_xticks(x)
            ax.set_xticklabels(categories, rotation=45, ha='right')
            ax.legend(fontsize=8, loc='upper left')

            # Set y-axis range so bars go to the top
            ax.set_ylim(0, 105)  # Leave 5% space for top labels

            ax.grid(True, alpha=0.3, axis='y', linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('B. Nine-Category Deviation Composition')

def main():
    """Main function"""
    set_clean_style()

    # File paths
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"

    print("="*80)
    print("Creating Figure 3: Percentage Composition of Deviation Categories Across Biological Groups")
    print(f"Z threshold: {1.96} (95% confidence interval)")
    print("="*80)

    # Initialize analyzer
    analyzer = GCBackgroundComparison(simulation_file, real_genome_file, z_threshold=1.96)

    try:
        # Execute analysis
        analyzer.load_and_preprocess_data()
        analyzer.build_background_model()
        analyzer.calculate_deviations()

        # Create Figure 3
        print("\n=== Creating Figure 3 ===")
        print("Creating Figure 3: Deviation composition across biological groups (AC version)...")
        fig3 = analyzer.create_figure_3_deviation_composition()

        # Export results
        output_dir = Path(base_dir) / "figure3_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save figure
        print("\n=== Saving Figure ===")
        fig3.savefig(output_dir / "figure3_deviation_composition_AC.png", dpi=300, bbox_inches='tight')
        fig3.savefig(output_dir / "figure3_deviation_composition_AC.pdf", bbox_inches='tight')
        print(f"Figure 3 (AC version): {output_dir}/figure3_deviation_composition_AC.png/pdf")

        # Display figure
        plt.show()

        print("\n" + "="*80)
        print("Figure 3 created successfully! Figure saved.")
        print("="*80)

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
