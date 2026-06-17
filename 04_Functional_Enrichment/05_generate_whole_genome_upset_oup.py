#!/usr/bin/env python3
"""
Generate whole genome G4 sharing pattern Upset plot
Focus on whole genome analysis, maintaining the exact same style as the original script
Modifications:
  1. Output editable PDF (pdf.fonttype=42)
  2. Also output high-resolution PNG (300 dpi)
  3. New high-resolution TIFF (600 dpi, LZW compression) for publication
  4. Provides example comments for legend position adjustment
"""

import os
import sys
from itertools import combinations
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import upsetplot as upsplt
import warnings
from collections import defaultdict

# Set PDF fonts to editable (Type 42 TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42  # Also set PostScript

warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=UserWarning)

# Set project root directory
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"

# Set matplotlib parameters
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'

# Colorblind-friendly palette
CBPalette = {
    "Black": "#000000",
    "Orange": "#F4A637",
    "Light blue": "#B6DBFF",
    "Vermilion": "#DB5829",
    "Mid blue": "#7BB0DF",
    "Maroon": "#894B45",
    "Dark blue": "#1964B0",
    "Light purple": "#D2BBD7",
    "Light teal": "#00C992",
    "Purple": "#AE75A2",
    "Teal": "#008A69",
    "Dark purple": "#882D71",
    "Dark teal": "#386350",
    "Grey": "#DEDEDE",
    "Yellow": "#E9DC6D"
}

# Species name mapping dictionary
speciesnos = {
    1: ['human', 'Homo_sapiens', 'hs1', 'H. sapiens'],
    2: ['bonobo', 'Pan_paniscus', 'pan', 'P. paniscus'],
    3: ['chimp', 'Pan_troglodytes', 'pan', 'P. troglodytes'],
    4: ['gorilla', 'Gorilla_gorilla', 'gor', 'G. gorilla'],
    5: ['sorang', 'Pongo_abelii', 'pon', 'P. abelii'],
    6: ['borang', 'Pongo_pygmaeus', 'pon', 'P. pygmaeus']
}

def comma_formatter(x, pos):
    '''Format with thousands separator'''
    return '{:,.0f}'.format(x)

def splitIDsandarrange(df):
    '''Split IDs and rearrange column order'''
    dfOut = df[0].str.split('|', expand=True)
    dfOut.columns = ['species', 'chrom', 'start', 'strand', 'length', 'score']
    dfOut['species'] = dfOut['species'].astype(int)
    dfOut['start'] = dfOut['start'].astype(int)
    dfOut['length'] = dfOut['length'].astype(int)
    dfOut['end'] = dfOut['start'] + dfOut['length']
    dfOut.sort_values(by=['species', 'chrom', 'start'], inplace=True)
    dfOut['chrom'] = dfOut['chrom'].apply(lambda x: f'chr{x}')
    dfOut.drop(columns=['length'], inplace=True)
    dfOut.reset_index(drop=True, inplace=True)
    dfOut["dummy"] = '.'
    dfOut = dfOut[['species', 'chrom', 'start', 'end', 'dummy', 'score', 'strand']]
    return dfOut

def stackedBarUpset(upsetDatadf, species, alignedUnique):
    '''Adjust data for stacked bar Upset plot'''
    upsetDatadf.reset_index(inplace=True)

    for nos, target in enumerate(species):
        # Find rows specific to this species
        condition = ((upsetDatadf[target] == True) &
                     (upsetDatadf[[s for s in species if s != target]].eq(False).all(axis=1)))
        condition_indices = upsetDatadf[condition].index

        # Assign aligned and unaligned species-specific pG4s
        if len(condition_indices) > 0:
            # First alignedUnique[nos] are aligned
            upsetDatadf.loc[condition_indices[:alignedUnique[nos]], "forUnique01"] = "Aligned species-specific IMs"
            # The rest are unaligned
            upsetDatadf.loc[condition_indices[alignedUnique[nos]:], "forUnique01"] = "Unaligned species-specific IMs"

    # Mark shared pG4s
    upsetDatadf["forUnique01"] = upsetDatadf["forUnique01"].fillna("Shared IMs")
    upsetDatadf.set_index(species, inplace=True)
    return upsetDatadf

def analyze_intersection_patterns(presAbs_df, species_list, title="Whole Genome"):
    """Analyze and output intersection pattern statistics"""

    print(f"\n{'='*80}")
    print(f"Whole Genome Intersection Pattern Analysis")
    print('='*80)

    # Count each pattern
    pattern_counts = {}
    total_ims = len(presAbs_df)

    # Convert DataFrame to list for processing
    data_matrix = presAbs_df[species_list].values

    # Count all possible intersection patterns
    for i in range(len(data_matrix)):
        # Create pattern string, e.g., "111000" means first 3 species present, last 3 absent
        pattern_str = ''.join(['1' if x == 1 else '0' for x in data_matrix[i]])

        if pattern_str not in pattern_counts:
            pattern_counts[pattern_str] = 0
        pattern_counts[pattern_str] += 1

    # Sort by count
    sorted_patterns = sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)

    # Output console information
    print(f"\nTotal IMs: {total_ims:,}")
    print(f"\nDifferent intersection patterns: {len(pattern_counts)}")
    print("\nDetailed statistics:")
    print("-" * 100)
    print(f"{'Pattern':<15} {'Binary':<25} {'Species Combination':<40} {'Count':>10} {'Percentage':>10}")
    print("-" * 100)

    pattern_data = []
    for pattern_str, count in sorted_patterns:
        # Parse pattern string
        species_present = []
        pattern_binary = ' '.join(pattern_str[i:i+1] for i in range(0, len(pattern_str), 1))

        for i, species in enumerate(species_list):
            if pattern_str[i] == '1':
                species_present.append(species)

        # Generate readable species combination description
        if len(species_present) == 0:
            species_combo = "No species"
        elif len(species_present) == 1:
            species_combo = f"{species_present[0]} (specific)"
        elif len(species_present) == len(species_list):
            species_combo = "All species shared"
        else:
            species_combo = ' + '.join(species_present)

        percentage = (count / total_ims) * 100
        print(f"{pattern_str:<15} {pattern_binary:<25} {species_combo:<40} {count:>10,} {percentage:>9.2f}%")

        pattern_data.append({
            'pattern': pattern_str,
            'pattern_binary': pattern_binary,
            'species_count': len(species_present),
            'species_combination': species_combo,
            'count': count,
            'percentage': percentage
        })

    # Statistics by number of shared species
    print("\nStatistics by number of shared species:")
    print("-" * 60)

    species_count_stats = defaultdict(int)
    for pattern_str, count in pattern_counts.items():
        species_count = pattern_str.count('1')
        species_count_stats[species_count] += count

    for count in sorted(species_count_stats.keys()):
        total = species_count_stats[count]
        percentage = (total / total_ims) * 100
        if count == 1:
            print(f"{count} species (species-specific): {total:>15,} ({percentage:>6.2f}%)")
        elif count == len(species_list):
            print(f"{count} species (all species shared): {total:>12,} ({percentage:>6.2f}%)")
        else:
            print(f"{count} species: {total:>25,} ({percentage:>6.2f}%)")

    # Output summary statistics
    print("\nSummary statistics:")
    print("-" * 60)
    print(f"Total IMs: {total_ims:,}")
    print(f"Total species: {len(species_list)}")
    print(f"Total patterns: {len(pattern_counts)}")
    print(f"Average IMs per pattern: {total_ims/len(pattern_counts):,.1f}")

    return pd.DataFrame(pattern_data), total_ims

def generate_whole_genome_upset(df, alignedUniqueGQs, output_dir="plots/whole_genome", stats_dir="stats/whole_genome"):
    '''Generate whole genome Upset plot, output editable PDF, high-resolution PNG and TIFF'''

    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    # Species order - identical to original script
    order = ["P. paniscus", "P. troglodytes", "H. sapiens", "G. gorilla", "P. pygmaeus", "P. abelii"]

    print(f"\n{'='*100}")
    print(f"Processing: Whole Genome")
    print('='*100)

    # Generate presence-absence matrix
    print("Generating presence-absence matrix...")
    presabs = []
    grouped_df = df.groupby("ID")

    total_groups = len(df.groupby("ID"))
    print(f"Total IMs groups: {total_groups:,}")

    for idx, (name, group) in enumerate(grouped_df):
        if (idx + 1) % 100000 == 0:
            print(f"  Processing progress: {idx + 1:,}/{total_groups:,} ({((idx+1)/total_groups*100):.1f}%)")

        specieslist = [0] * 6
        presspecies = group["SPECIES"].unique()
        for index in presspecies:
            specieslist[index-1] = 1
        presabs.append(specieslist)

    presAbs_df = pd.DataFrame(presabs, columns=[i[3] for i in speciesnos.values()])
    presAbs_dfUpset = presAbs_df[order]
    speciesUpset = order

    # Analyze intersection patterns
    print("\nAnalyzing intersection patterns...")
    pattern_df, total_ims = analyze_intersection_patterns(presAbs_dfUpset, speciesUpset, "Whole Genome")

    # Save pattern statistics to CSV
    stats_filename = "whole_genome_intersection_patterns.csv"
    stats_path = os.path.join(stats_dir, stats_filename)
    pattern_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
    print(f"\nDetailed statistics saved to: {stats_path}")

    # Get aligned species-specific pG4s counts
    print("\nGetting aligned species-specific pG4s counts...")
    alignedUniqueUpset = []

    for s in range(1, 7):
        # Whole genome: select all records for this species
        specAlignedUniqueGQs = alignedUniqueGQs[alignedUniqueGQs['species'] == s]
        count = specAlignedUniqueGQs.shape[0]
        alignedUniqueUpset.append(count)
        print(f"  {speciesnos[s][3]}: {count:,} aligned species-specific pG4s")

    # Reorder according to order list - identical to original script
    human_element = alignedUniqueUpset.pop(0)
    alignedUniqueUpset.insert(2, human_element)
    sorang_element = alignedUniqueUpset.pop(4)
    alignedUniqueUpset.insert(5, sorang_element)
    alignedUniqueUpset = np.array(alignedUniqueUpset)

    # Generate Upset plot data
    print("\nGenerating Upset plot data...")
    presAbsMatrixNormUpsetDict = {}
    for column in presAbs_dfUpset.columns:
        indices = [i for i, value in enumerate(presAbs_dfUpset[column]) if value == 1]
        presAbsMatrixNormUpsetDict[column] = indices

    # Create Upset plot data
    upsetData = upsplt.from_contents(presAbsMatrixNormUpsetDict)
    upsetData["forUnique01"] = "Shared IMs"
    upsetData = stackedBarUpset(upsetData, speciesUpset, alignedUniqueUpset)

    # Create Upset plot - identical parameters to original script
    print("Creating Upset plot...")
    upset = upsplt.UpSet(
        upsetData,
        sort_by="cardinality",
        sort_categories_by="-input",
        facecolor=CBPalette["Black"],
        show_counts=False,
        totals_plot_elements=6,
        intersection_plot_elements=0
    )

    # Add stacked bar chart - identical parameters to original script
    upset.add_stacked_bars(
        by="forUnique01",
        colors=[CBPalette["Dark blue"], CBPalette["Vermilion"], CBPalette["Light blue"]],
        elements=10
    )

    # Draw plot - same dimensions as original script
    fig = plt.figure(figsize=(12, 8))
    plot_result = upset.plot(fig=fig)

    # Format stacked bar area - identical to original script
    plot_result["extra0"].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["extra0"].set_yticklabels(plot_result["extra0"].get_yticklabels(), fontsize=11)

    # Format matrix area - identical to original script
    plot_result["matrix"].set_yticklabels(
        ["S. orangutan", "B. orangutan", "Gorilla", "Human", "Chimpanzee", "Bonobo"],
        fontsize=12
    )

    # Format totals area - identical to original script
    plot_result["totals"].xaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["totals"].set_xlabel("\nTotal IMs in species", fontsize=16)
    plot_result["totals"].set_xticklabels(plot_result["totals"].get_xticklabels(), fontsize=7)

    # Set overall labels - identical to original script
    plt.ylabel("Number of IMs\n", fontsize=18)
    plt.grid(alpha=0.5, linestyle="--")

    # ========== Legend position adjustment notes ==========
    # upsetplot legends are usually located in the stacked bar subplot ("extra0").
    # You can get the legend handle and adjust its position via:
    #   legend = plot_result["extra0"].get_legend()
    #   if legend:
    #       legend.set_bbox_to_anchor((1.05, 1))  # Adjust to the right of the subplot
    # Or use more precise parameters: legend.set_bbox_to_anchor((x, y), loc='upper left')
    # For full control, you can also create a new legend directly:
    #   handles, labels = plot_result["extra0"].get_legend_handles_labels()
    #   plot_result["extra0"].legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
    # Note: modifications must be made before plt.savefig.
    # By default, upsetplot places the legend automatically. Uncomment the code below as needed.
    #
    # Example: Move legend outside the top-right corner of the figure
    legend = plot_result["extra0"].get_legend()
    if legend:
        legend.set_bbox_to_anchor((0.25, 1))
    # =====================================

    # Save editable PDF
    pdf_path = os.path.join(output_dir, "whole_genome.filtered_upset.pdf")
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    print(f"\nPDF saved (editable text): {pdf_path}")

    # Save SVG (vector format)
    svg_path = os.path.join(output_dir, "whole_genome.filtered_upset.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    print(f"SVG format saved: {svg_path}")

    # Save high-resolution PNG (300 dpi)
    png_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_300dpi.png")
    plt.savefig(png_highres_path, format='png', dpi=300, bbox_inches='tight')
    print(f"High-resolution PNG saved (300 dpi): {png_highres_path}")

    # Save high-resolution TIFF (600 dpi, LZW compression) - preferred for publication
    # Note: requires pillow library for TIFF format and LZW compression
    try:
        tiff_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_600dpi.tiff")
        plt.savefig(tiff_highres_path, format='tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        print(f"High-resolution TIFF saved (600 dpi, LZW): {tiff_highres_path}")
    except Exception as e:
        print(f"Warning: Could not save TIFF. Please ensure pillow is installed: pip install pillow")
        print(f"Error: {e}")

    plt.close()

    return pattern_df, total_ims

def main():
    """Main function - handles whole genome data specifically, maintaining identical logic to the original script"""

    print("=== Whole Genome G4 Sharing Pattern Upset Plot Generation Script ===")
    print(f"Project directory: {BASE_DIR}")

    # Create output directories - using same structure as original script
    stats_dir = os.path.join(BASE_DIR, "output/stats")
    plots_dir = os.path.join(BASE_DIR, "output/plots/upsetPlots")
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    # Load aligned species-specific pG4s dataset
    print("\n1. Loading aligned species-specific pG4s data...")
    aligned_unique_path = os.path.join(BASE_DIR, "output/datasets/alignedUniquehsaG.egs")

    if not os.path.exists(aligned_unique_path):
        print(f"Error: File does not exist - {aligned_unique_path}")
        sys.exit(1)

    try:
        alignedUniquedf = pd.read_csv(aligned_unique_path, sep='\t', header=None)
        alignedUniqueGQs = splitIDsandarrange(alignedUniquedf)
        print(f"  Successfully loaded: {len(alignedUniqueGQs)} records")
    except Exception as e:
        print(f"  Failed to load: {e}")
        sys.exit(1)

    # Load whole genome data
    print("\n2. Loading whole genome data...")
    whole_genome_path = os.path.join(BASE_DIR, "output/datasets/allhsaG.graph.df")

    if not os.path.exists(whole_genome_path):
        print(f"Error: Whole genome data file does not exist - {whole_genome_path}")
        sys.exit(1)

    try:
        # Read full data, consistent with original script
        df_whole = pd.read_csv(whole_genome_path, header=0, sep="\t", low_memory=False)
        print(f"  Successfully loaded: {len(df_whole)} rows")
    except Exception as e:
        print(f"  Failed to load: {e}")
        sys.exit(1)

    # Check data quality
    print("\n3. Checking data quality...")
    print(f"  Unique ID count: {df_whole['ID'].nunique():,}")
    print(f"  Species distribution:")
    species_counts = df_whole['SPECIES'].value_counts().sort_index()
    for species_id, count in species_counts.items():
        species_name = speciesnos.get(species_id, ['Unknown'])[3]
        print(f"    {species_name}: {count:,} records")

    # Generate whole genome Upset plot
    print("\n4. Generating whole genome Upset plot...")
    try:
        pattern_df, total_ims = generate_whole_genome_upset(
            df_whole, alignedUniqueGQs, plots_dir, stats_dir
        )

        # Output final summary
        print(f"\n{'='*100}")
        print("Whole genome analysis completed!")
        print('='*100)
        print(f"Total IMs: {total_ims:,}")
        print(f"Unique sharing patterns: {len(pattern_df)}")
        print(f"Average IMs per pattern: {pattern_df['count'].mean():,.1f}")

        print(f"\nOutput directories:")
        print(f"  Plots: {plots_dir}")
        print(f"  Statistics: {stats_dir}")
        print('='*100)

    except Exception as e:
        print(f"\nError generating Upset plot: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
