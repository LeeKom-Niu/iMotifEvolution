#!/usr/bin/env python3
"""
figure1_two_versions_final.py - Final version: dual-version comparison figure
Version 1: Full information version (with n values)
Version 2: Publication clean version (concise, no n values)
Theoretical background curve description more accurate
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set font - use available fonts on Linux
def set_font_settings():
    """Set font, compatible with Linux environment"""
    import matplotlib
    # Get available fonts
    available_fonts = matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')

    # Prefer DejaVu Sans, a good font usually available on Linux
    preferred_fonts = ['DejaVu Sans', 'Liberation Sans', 'FreeSans', 'Nimbus Sans']

    selected_font = None
    for font_name in preferred_fonts:
        # Check if these fonts are available
        for font_path in available_fonts:
            if font_name.lower() in font_path.lower():
                selected_font = font_name
                break
        if selected_font:
            break

    if not selected_font:
        # If preferred fonts not found, use first available sans-serif font
        for font_path in available_fonts:
            font_prop = matplotlib.font_manager.FontProperties(fname=font_path)
            font_name = font_prop.get_name()
            if 'sans' in font_name.lower():
                selected_font = font_name
                break

    if not selected_font:
        selected_font = 'sans-serif'  # Use system default

    print(f"Using font: {selected_font}")

    # Set Nature style, using found font
    plt.rcParams.update({
        # Font settings
        'font.family': 'sans-serif',
        'font.sans-serif': [selected_font],
        'font.size': 8,
        'pdf.fonttype': 42,

        # Axes
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'axes.linewidth': 0.6,
        'axes.unicode_minus': False,
        'axes.labelweight': 'normal',

        # Ticks
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,

        # Legend
        'legend.fontsize': 7,
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '#CCCCCC',
        'legend.fancybox': False,

        # Lines
        'lines.linewidth': 1.8,
        'lines.markersize': 4.0,
        'lines.markeredgewidth': 0.5,

        # Scatter
        'scatter.marker': 'o',
        'scatter.edgecolors': 'white',

        # Figure
        'figure.dpi': 300,
        'figure.constrained_layout.use': True,
        'savefig.dpi': 600,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,

        # Grid
        'grid.alpha': 0.2,
        'grid.linestyle': ':',
        'grid.linewidth': 0.5,
    })

# Optimized color scheme
two_category_colors = {
    'Prokaryote': '#1F77B4',     # Dark blue - Prokaryotes
    'Eukaryota': '#D62728',      # Red - Eukaryotes
}

# Eukaryotic subgroup colors
eukaryotic_subgroup_colors = {
    'Protozoa': '#BCBD22',       # Olive green - Protozoa
    'Fungi': '#9467BD',          # Purple - Fungi
    'Plant': '#8C564B',          # Brown - Plants
    'Invertebrate': '#E377C2',   # Pink - Invertebrates
    'Vertebrate Other': '#7F7F7F',  # Gray - Other Vertebrates
    'Mammalian': '#D62728',      # Red - Mammals
}

# Theoretical background curve color - use more neutral color
BACKGROUND_COLOR = '#2CA02C'      # Soft green, distinct from data points
BACKGROUND_FILL_COLOR = '#C7E9C0' # Light green fill

class FinalDensityPlotter:
    def __init__(self, simulation_file, real_genome_file):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.background_data = None
        self.real_data = None
        self.smooth_background = None

        # Classification definitions
        self.two_categories = ['Prokaryote', 'Eukaryota']
        # Eukaryotic groups
        self.eukaryotic_subgroups = [
            'Protozoa', 'Fungi', 'Plant',
            'Invertebrate', 'Vertebrate Other', 'Mammalian'
        ]

        # Statistics
        self.stats = {}

        # Adaptive background curve range
        self.bg_gc_min = None
        self.bg_gc_max = None

    def load_and_preprocess_data(self):
        """Load and preprocess all data"""
        print("=== Loading Data ===")

        # 1. Load real genome data
        print(f"1. Loading real genome data: {self.real_genome_file}")
        try:
            real_data = pd.read_csv(self.real_genome_file, sep='\t')
            print(f"   Loaded {len(real_data)} rows successfully")
        except Exception as e:
            print(f"   Failed to load: {e}")
            return False

        # Standardize column names
        column_mapping = {}
        for col in real_data.columns:
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
            if old_name in real_data.columns:
                real_data = real_data.rename(columns={old_name: new_name})

        # Clean data
        if 'Classification' in real_data.columns:
            real_data = real_data[real_data['Classification'] != 'Classification']
            real_data['Classification'] = real_data['Classification'].str.strip()
            real_data['Classification'] = real_data['Classification'].replace({
                'Vertebrate Mammalian': 'Mammalian',
                'Vertebrate_Other': 'Vertebrate Other'
            })

        # Convert data types
        real_data['GC_content'] = pd.to_numeric(real_data['GC_content'], errors='coerce')
        real_data['Genomic_density'] = pd.to_numeric(real_data['Genomic_density'], errors='coerce')

        # Convert GC content from percentage to decimal
        if real_data['GC_content'].max() > 1:
            real_data['GC_content'] = real_data['GC_content'] / 100.0
            print("   Converted GC content from percentage to decimal")

        # Remove invalid data
        original_len = len(real_data)
        real_data = real_data.dropna(subset=['GC_content', 'Genomic_density', 'Classification'])
        print(f"   Data cleaning: {original_len} -> {len(real_data)} rows")

        # Filter out prokaryotes with GC content = 0
        prokaryotic_mask = real_data['Classification'].isin(['Bacteria', 'Archaea'])
        gc_zero_mask = (real_data['GC_content'] == 0) & prokaryotic_mask
        if gc_zero_mask.any():
            print(f"   Removed {gc_zero_mask.sum()} prokaryotic data points with GC=0")
            real_data = real_data[~gc_zero_mask]

        # Create classification columns
        self._create_classification_columns(real_data)

        self.real_data = real_data

        print(f"\n2. Data statistics:")
        print(f"   - Total genomes: {len(real_data):,}")
        print(f"   - GC range: {real_data['GC_content'].min():.3f} - {real_data['GC_content'].max():.3f}")
        print(f"   - Density range: {real_data['Genomic_density'].min():.3f} - {real_data['Genomic_density'].max():.3f}")

        # Calculate adaptive background curve range
        self._calculate_adaptive_background_range()

        return True

    def _create_classification_columns(self, data):
        """Create classification columns"""
        # Two categories (prokaryotes and eukaryotes)
        data['Classification_2cat'] = np.nan

        # Prokaryotes: Bacteria + Archaea
        prokaryotic_mask = data['Classification'].isin(['Bacteria', 'Archaea'])
        data.loc[prokaryotic_mask, 'Classification_2cat'] = 'Prokaryote'

        # Eukaryotes
        eukaryotic_categories = ['Protozoa', 'Fungi', 'Plant',
                               'Invertebrate', 'Vertebrate Other', 'Mammalian']
        eukaryotic_mask = data['Classification'].isin(eukaryotic_categories)
        data.loc[eukaryotic_mask, 'Classification_2cat'] = 'Eukaryota'

        # Eukaryotic subgroups
        data['Classification_eukaryote'] = np.nan
        data.loc[eukaryotic_mask, 'Classification_eukaryote'] = data.loc[eukaryotic_mask, 'Classification']

        # Count per category
        for cat_type in ['2cat', 'eukaryote']:
            col = f'Classification_{cat_type}'
            if col in data.columns:
                counts = data[col].value_counts()
                self.stats[cat_type] = counts.to_dict()
                print(f"   - {cat_type} categories: {len(counts)} groups")
                for cat, count in counts.items():
                    print(f"     {cat}: {count:,}")

    def _calculate_adaptive_background_range(self):
        """Calculate adaptive background curve range based on real data"""
        if self.real_data is None:
            return

        # Get GC range of all data
        all_gc = self.real_data['GC_content'].values
        gc_min = np.percentile(all_gc, 2)  # 2nd percentile as lower bound
        gc_max = np.percentile(all_gc, 98)  # 98th percentile as upper bound

        # Ensure reasonable range
        gc_min = max(0.0, gc_min - 0.05)  # Expand slightly
        gc_max = min(0.9, gc_max + 0.05)  # Limit within 0.9

        self.bg_gc_min = gc_min
        self.bg_gc_max = gc_max

        print(f"\n3. Adaptive background curve range:")
        print(f"   - GC range: {self.bg_gc_min:.3f} - {self.bg_gc_max:.3f}")

    def load_theoretical_background(self):
        """Load random simulation background curve"""
        print("\n=== Loading Random Simulation Background Curve ===")

        print(f"1. Loading simulation data: {self.simulation_file}")
        try:
            sim_data = pd.read_csv(self.simulation_file)

            # Find GC and density columns
            gc_col = None
            density_col = None

            for col in sim_data.columns:
                col_lower = col.lower()
                if 'gc' in col_lower:
                    gc_col = col
                if 'density' in col_lower or ('im' in col_lower and 'mb' in col_lower):
                    density_col = col

            if not gc_col or not density_col:
                # If not found, use first two columns
                gc_col = sim_data.columns[0]
                density_col = sim_data.columns[1]
                print(f"   Using default columns: GC={gc_col}, Density={density_col}")

            bg_data = sim_data[[gc_col, density_col]].copy()
            bg_data.columns = ['GC_content', 'density']
            print(f"   Created background curve from simulation data ({len(bg_data)} rows)")
        except Exception as e:
            print(f"   Failed to load simulation data: {e}")
            return False

        # Ensure GC content is decimal
        if bg_data['GC_content'].max() > 1:
            bg_data['GC_content'] = bg_data['GC_content'] / 100.0
            print("   Converted GC content from percentage to decimal")

        # Remove duplicate GC values (take mean)
        print(f"   Before processing: {len(bg_data)} rows")
        bg_data = bg_data.groupby('GC_content', as_index=False)['density'].mean()
        print(f"   After removing duplicates: {len(bg_data)} rows")

        # Ensure data is sorted
        bg_data = bg_data.sort_values('GC_content').dropna().reset_index(drop=True)

        print(f"\n2. Random simulation background data:")
        print(f"   - Data points: {len(bg_data)}")
        print(f"   - GC range: {bg_data['GC_content'].min():.3f} - {bg_data['GC_content'].max():.3f}")
        print(f"   - Density range: {bg_data['density'].min():.3f} - {bg_data['density'].max():.3f}")

        self.background_data = bg_data
        return True

    def create_adaptive_smooth_curve(self):
        """Create adaptive smooth random simulation background curve"""
        print("\n=== Creating Adaptive Smooth Random Simulation Background Curve ===")

        if self.background_data is None:
            print("   No background data available")
            return False

        # Extract data
        gc_values = self.background_data['GC_content'].values
        density_values = self.background_data['density'].values

        # Ensure strictly increasing with no duplicates
        unique_indices = np.unique(gc_values, return_index=True)[1]
        gc_values = gc_values[unique_indices]
        density_values = density_values[unique_indices]

        # Ensure sorted
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]

        print(f"1. Raw data:")
        print(f"   - Points: {len(gc_values)}")
        print(f"   - GC range: {gc_values.min():.3f} - {gc_values.max():.3f}")

        # Filter data based on adaptive range
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            # Expand range slightly for better interpolation
            extended_min = max(0.0, self.bg_gc_min - 0.05)
            extended_max = min(0.9, self.bg_gc_max + 0.05)

            mask = (gc_values >= extended_min) & (gc_values <= extended_max)
            if mask.any():
                gc_values = gc_values[mask]
                density_values = density_values[mask]
                print(f"   - After adaptive filtering: {len(gc_values)} points")
                print(f"   - Filter range: {extended_min:.3f} - {extended_max:.3f}")
            else:
                print(f"   No data within adaptive range, using raw data")

        # Ensure enough data points
        if len(gc_values) < 3:
            print(f"   Insufficient data points ({len(gc_values)}), using original range")
            gc_values = self.background_data['GC_content'].values
            density_values = self.background_data['density'].values
            gc_values = gc_values[:]
            density_values = density_values[:]

        # Extend data to adaptive range boundaries
        if self.bg_gc_min is not None and gc_values.min() > self.bg_gc_min:
            # Add point at low end
            low_gc = np.array([self.bg_gc_min])
            low_density = np.interp(low_gc, gc_values, density_values)
            gc_values = np.concatenate([low_gc, gc_values])
            density_values = np.concatenate([low_density, density_values])
            print(f"   Extended low GC end to {self.bg_gc_min:.3f}")

        if self.bg_gc_max is not None and gc_values.max() < self.bg_gc_max:
            # Add point at high end
            high_gc = np.array([self.bg_gc_max])
            high_density = np.interp(high_gc, gc_values, density_values)
            gc_values = np.concatenate([gc_values, high_gc])
            density_values = np.concatenate([density_values, high_density])
            print(f"   Extended high GC end to {self.bg_gc_max:.3f}")

        # Create adaptive sampling points
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            smooth_gc = np.linspace(self.bg_gc_min, self.bg_gc_max, 300)
        else:
            smooth_gc = np.linspace(gc_values.min(), gc_values.max(), 300)

        print(f"2. Creating adaptive smooth curve:")
        print(f"   - Sampling points: {len(smooth_gc)}")
        print(f"   - Final GC range: {smooth_gc.min():.3f} - {smooth_gc.max():.3f}")

        try:
            # Use cubic spline interpolation
            interp_func = interp1d(gc_values, density_values, kind='cubic',
                                 fill_value='extrapolate', bounds_error=False)
            smooth_density = interp_func(smooth_gc)
            print("   Using cubic spline interpolation")
        except Exception as e:
            print(f"   Cubic spline failed: {e}")
            try:
                # Try quadratic spline
                interp_func = interp1d(gc_values, density_values, kind='quadratic',
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
                print("   Using quadratic spline interpolation")
            except:
                # Use linear interpolation
                interp_func = interp1d(gc_values, density_values, kind='linear',
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
                print("   Using linear interpolation")

        # Ensure non-negative values
        smooth_density = np.maximum(smooth_density, 0)

        # Apply slight smoothing
        smooth_density = gaussian_filter1d(smooth_density, sigma=1.5)

        print(f"3. Smooth curve statistics:")
        print(f"   - Density range: {smooth_density.min():.3f} - {smooth_density.max():.3f}")
        print(f"   - Mean density: {smooth_density.mean():.3f}")

        # Save smooth curve
        self.smooth_background = pd.DataFrame({
            'GC_content': smooth_gc,
            'density': smooth_density
        })

        return True

    def create_complete_version_figure(self):
        """Create full information version figure (with n values)"""
        print("\n=== Creating Full Information Version Figure ===")

        # Create figure: two panels side by side
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

        # 1. Prokaryote vs Eukaryote comparison
        print("1. Creating prokaryote vs eukaryote comparison...")
        self._plot_prok_euk_comparison(axes[0], show_n=True)

        # 2. Eukaryotic subgroup comparison
        print("2. Creating eukaryotic subgroup comparison...")
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=True)

        print("Full information version Figure created successfully")
        return fig

    def create_clean_version_figure(self):
        """Create publication clean version figure (no n values)"""
        print("\n=== Creating Publication Clean Version Figure ===")

        # Create figure: two panels side by side
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)

        # 1. Prokaryote vs Eukaryote comparison
        print("1. Creating prokaryote vs eukaryote comparison...")
        self._plot_prok_euk_comparison(axes[0], show_n=False)

        # 2. Eukaryotic subgroup comparison
        print("2. Creating eukaryotic subgroup comparison...")
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=False)

        print("Publication clean version Figure created successfully")
        return fig

    def _plot_prok_euk_comparison(self, ax, show_n=True):
        """Plot prokaryote vs eukaryote comparison"""
        # Plot random simulation background curve (adaptive range)
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values

            # Use semi-transparent fill for subtle background
            ax.fill_between(bg_gc, 0, bg_density,
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')

            # Plot background curve
            ax.plot(bg_gc, bg_density,
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)

        # Plot scatter data
        if self.real_data is not None and 'Classification_2cat' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_2cat'])

            for category in self.two_categories:
                cat_data = plot_data[plot_data['Classification_2cat'] == category]
                if len(cat_data) > 0:
                    color = two_category_colors[category]

                    # Create legend label
                    label = f"{category}"
                    if show_n and '2cat' in self.stats and category in self.stats['2cat']:
                        label += f" (n={self.stats['2cat'][category]:,})"

                    # Use slightly larger points for better display
                    ax.scatter(
                        cat_data['GC_content'],
                        cat_data['Genomic_density'],
                        color=color,
                        alpha=0.75,
                        s=12,
                        label=label,
                        edgecolors='white',
                        linewidth=0.3,
                        zorder=10
                    )

        # Set axis labels
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')

        # Adaptive axis limits
        self._set_adaptive_axis_limits(ax, '2cat')

        # Set ticks
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)

        # Add grid
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)

        # Add legend
        self._add_legend_to_axis(ax, show_n=show_n)

    def _plot_eukaryote_subgroup_comparison(self, ax, show_n=True):
        """Plot eukaryotic subgroup comparison"""
        # Plot random simulation background curve (adaptive range)
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values

            ax.fill_between(bg_gc, 0, bg_density,
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')

            ax.plot(bg_gc, bg_density,
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)

        # Plot scatter data
        if self.real_data is not None and 'Classification_eukaryote' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_eukaryote'])

            for category in self.eukaryotic_subgroups:
                cat_data = plot_data[plot_data['Classification_eukaryote'] == category]
                if len(cat_data) > 0:
                    color = eukaryotic_subgroup_colors[category]

                    # Create legend label
                    label = f"{category}"
                    if show_n and 'eukaryote' in self.stats and category in self.stats['eukaryote']:
                        label += f" (n={self.stats['eukaryote'][category]:,})"

                    ax.scatter(
                        cat_data['GC_content'],
                        cat_data['Genomic_density'],
                        color=color,
                        alpha=0.8,
                        s=14,  # Eukaryotic data points slightly larger for differentiation
                        label=label,
                        edgecolors='white',
                        linewidth=0.3,
                        zorder=10
                    )

        # Set axis labels
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')

        # Adaptive axis limits
        self._set_adaptive_axis_limits(ax, 'eukaryote')

        # Set ticks
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)

        # Add grid
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)

        # Add legend - corrected position: now in upper left
        self._add_eukaryotic_legend_to_axis(ax, show_n=show_n)

    def _set_adaptive_axis_limits(self, ax, cat_type):
        """Adaptively set axis limits"""
        if self.real_data is None:
            return

        col = f'Classification_{cat_type}'
        if col not in self.real_data.columns:
            return

        plot_data = self.real_data.dropna(subset=[col])
        if len(plot_data) == 0:
            return

        # X-axis range: based on data distribution
        x_data = plot_data['GC_content'].values

        # Use 2nd and 98th percentiles, expand slightly
        x_min = np.percentile(x_data, 2) - 0.05
        x_max = np.percentile(x_data, 98) + 0.05

        # Ensure reasonable range
        x_min = max(-0.02, x_min)
        x_max = min(0.92, x_max)

        ax.set_xlim(x_min, x_max)

        # Y-axis range: based on data and background curve
        y_data = plot_data['Genomic_density'].values

        # Use 99th percentile as reference
        y_99 = np.percentile(y_data, 99)

        # If background curve exists, consider its maximum
        if self.smooth_background is not None:
            # Only consider background curve within current X range
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values

            # Filter background data within X range
            mask = (bg_gc >= x_min) & (bg_gc <= x_max)
            if mask.any():
                bg_max_in_range = bg_density[mask].max()
            else:
                bg_max_in_range = bg_density.max()

            y_max = max(y_99, bg_max_in_range) * 1.15
        else:
            y_max = y_99 * 1.15

        # Ensure minimum is 0
        ax.set_ylim(-0.02 * y_max, y_max)

        print(f"   - {cat_type} axis range: X={x_min:.3f}-{x_max:.3f}, Y=0-{y_max:.1f}")

    def _add_legend_to_axis(self, ax, show_n=True):
        """Add legend to axis (for prokaryote vs eukaryote)"""
        handles, labels = ax.get_legend_handles_labels()

        if len(handles) == 0:
            return

        # Separate background and data legends
        bg_handles = []
        bg_labels = []
        data_handles = []
        data_labels = []

        for handle, label in zip(handles, labels):
            if 'background' in label.lower():
                bg_handles.append(handle)
                bg_labels.append(label)
            else:
                data_handles.append(handle)
                data_labels.append(label)

        # Add data legend
        if data_handles:
            # Prokaryote vs Eukaryote use single column, upper left
            leg1 = ax.legend(data_handles, data_labels,
                           loc='upper left', fontsize=7.5,
                           handletextpad=0.5,
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=1)

            # Add background legend
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR,
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])

                # Add background legend in upper right
                ax.legend([bg_patch], [bg_labels[0]],
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)

                # Restore first legend
                ax.add_artist(leg1)

    def _add_eukaryotic_legend_to_axis(self, ax, show_n=True):
        """Add legend to axis (for eukaryotic subgroups) - corrected position: upper left"""
        handles, labels = ax.get_legend_handles_labels()

        if len(handles) == 0:
            return

        # Separate background and data legends
        bg_handles = []
        bg_labels = []
        data_handles = []
        data_labels = []

        for handle, label in zip(handles, labels):
            if 'background' in label.lower():
                bg_handles.append(handle)
                bg_labels.append(label)
            else:
                data_handles.append(handle)
                data_labels.append(label)

        # Add data legend
        if data_handles:
            # Eukaryotic subgroups use two columns, upper left
            leg1 = ax.legend(data_handles, data_labels,
                           loc='upper left', fontsize=7,
                           handletextpad=0.5,
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=2)

            # Add background legend
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR,
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])

                # Add background legend in upper right
                ax.legend([bg_patch], [bg_labels[0]],
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)

                # Restore first legend
                ax.add_artist(leg1)

    def save_both_versions(self, output_dir):
        """Save both versions of the figure"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== Saving Both Versions of Figures ===")

        saved_files = []

        # 1. Save full information version
        print("\n1. Saving full information version...")
        fig_complete = self.create_complete_version_figure()
        saved_complete = self._save_single_figure(fig_complete, output_dir, "figure1_complete_version")
        saved_files.extend(saved_complete)

        # 2. Save publication clean version
        print("\n2. Saving publication clean version...")
        fig_clean = self.create_clean_version_figure()
        saved_clean = self._save_single_figure(fig_clean, output_dir, "figure1_clean_version")
        saved_files.extend(saved_clean)

        plt.close('all')

        return saved_files

    def _save_single_figure(self, fig, output_dir, filename):
        """Save single figure in multiple formats"""
        formats = [
            ('tiff', 'tiff', {'pil_kwargs': {'compression': 'tiff_lzw'}}),
            ('pdf', 'pdf', {}),
            ('png', 'png', {}),
            ('svg', 'svg', {}),
        ]

        saved_files = []
        for fmt_name, fmt_ext, kwargs in formats:
            file_path = output_dir / f"{filename}.{fmt_ext}"
            try:
                fig.savefig(file_path, dpi=600, format=fmt_ext, **kwargs)
                print(f"   {fmt_name.upper()} format: {file_path}")
                saved_files.append(file_path)
            except Exception as e:
                print(f"   Failed to save {fmt_name.upper()}: {e}")

        return saved_files

    def export_statistics(self, output_dir):
        """Export statistics"""
        print("\n=== Exporting Statistics ===")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        stats_data = []

        # Collect all statistics
        for cat_type in ['2cat', 'eukaryote']:
            col = f'Classification_{cat_type}'
            if self.real_data is not None and col in self.real_data.columns:
                cat_data = self.real_data.dropna(subset=[col])

                for category in cat_data[col].unique():
                    sub_data = cat_data[cat_data[col] == category]

                    stats_row = {
                        'Category_Type': cat_type,
                        'Category': category,
                        'Count': len(sub_data),
                        'GC_Mean': sub_data['GC_content'].mean(),
                        'GC_Std': sub_data['GC_content'].std(),
                        'GC_Min': sub_data['GC_content'].min(),
                        'GC_Max': sub_data['GC_content'].max(),
                        'Density_Mean': sub_data['Genomic_density'].mean(),
                        'Density_Std': sub_data['Genomic_density'].std(),
                        'Density_Min': sub_data['Genomic_density'].min(),
                        'Density_Max': sub_data['Genomic_density'].max(),
                        'Density_Median': sub_data['Genomic_density'].median()
                    }
                    stats_data.append(stats_row)

        # Save to CSV
        if stats_data:
            stats_df = pd.DataFrame(stats_data)
            stats_path = output_dir / "figure1_statistics.csv"
            stats_df.to_csv(stats_path, index=False)
            print(f"Statistics: {stats_path}")

            # Print summary
            print("\nData summary:")
            for cat_type in ['2cat', 'eukaryote']:
                cat_stats = stats_df[stats_df['Category_Type'] == cat_type]
                if not cat_stats.empty:
                    print(f"\n{cat_type} classification:")
                    for _, row in cat_stats.iterrows():
                        print(f"  {row['Category']}: {row['Count']:,} genomes, "
                              f"GC={row['GC_Mean']:.3f}+-{row['GC_Std']:.3f} "
                              f"({row['GC_Min']:.3f}-{row['GC_Max']:.3f}), "
                              f"Density={row['Density_Mean']:.2f}+-{row['Density_Std']:.2f} IM/Mb")

def main():
    """Main function"""
    # First set font
    set_font_settings()

    # File paths
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"

    print("=" * 70)
    print("Figure 1: Two-Version Dual-Panel Comparison - Final")
    print("=" * 70)
    print("Version 1: Full information version (with n values)")
    print("Version 2: Publication clean version (no n values)")
    print("=" * 70)
    print("Figure content:")
    print("1. Prokaryotes vs Eukaryotes")
    print("2. Eukaryotic subgroups: Protozoa, Fungi, Plants, Invertebrates, Other Vertebrates, Mammals")
    print("=" * 70)
    print("Improvements:")
    print("1. Eukaryotic legend position: upper left")
    print("2. Background curve description: Random simulation background")
    print("3. Background curve color: green (better distinction from data points)")
    print("=" * 70)

    # Initialize plotter
    plotter = FinalDensityPlotter(simulation_file, real_genome_file)

    try:
        # 1. Load real data
        if not plotter.load_and_preprocess_data():
            print("Real data loading failed")
            return

        # 2. Load/create random simulation background curve
        if not plotter.load_theoretical_background():
            print("Random simulation background loading failed")
            return

        # 3. Create adaptive smooth random simulation background curve
        if not plotter.create_adaptive_smooth_curve():
            print("Adaptive smooth background curve creation failed")
            return

        # 4. Export statistics
        output_dir = Path(base_dir) / "figures" / "figure1_two_versions"
        plotter.export_statistics(output_dir)

        # 5. Save both versions of figures
        saved_files = plotter.save_both_versions(output_dir)

        print("\n" + "=" * 70)
        print("Both versions of Figure 1 generated successfully!")
        print("=" * 70)
        print("Main output files:")

        # Display by version
        print("\nFull information version:")
        for file_path in saved_files:
            if "complete_version" in str(file_path):
                print(f"  {file_path.suffix.upper()}: {file_path.name}")

        print("\nPublication clean version:")
        for file_path in saved_files:
            if "clean_version" in str(file_path):
                print(f"  {file_path.suffix.upper()}: {file_path.name}")

        print("=" * 70)

        # Display key information
        print("\nVersion comparison:")
        print("1. Full information version:")
        print("   - Includes sample counts per category (n=)")
        print("   - Suitable for internal reports and reviewer review")
        print("2. Publication clean version:")
        print("   - Concise, shows only category names")
        print("   - Suitable for final publication figures")
        print("3. Unified improvements:")
        print("   - Legend position: both in upper left (eukaryotic uses two columns)")
        print("   - Background description: Random simulation background")
        print("   - Color scheme: green background curve, good contrast with data points")
        print("=" * 70)

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
