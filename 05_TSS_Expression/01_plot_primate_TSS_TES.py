#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

# Set parameters
window_size = 1000  # Upstream/downstream range
relative_positions = np.arange(-window_size, window_size + 1)
output_dir = "primate_TSS_TES_enrichment_results"
n_processes = 6  # One process per species

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Species mapping table (display name: file info)
species_info = {
    "Human": {
        "gene_file": "input/gene_bed/human.genes.bed",
        "imotif_file": "input/imotif_bed/Homo_sapiens_all.bed",
        "color": "#FF6B6B"  # Red
    },
    "Chimpanzee": {
        "gene_file": "input/gene_bed/chimp.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_troglodytes_all.bed",
        "color": "#4ECDC4"  # Cyan
    },
    "Bonobo": {
        "gene_file": "input/gene_bed/bonobo.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_paniscus_all.bed",
        "color": "#45B7D1"  # Blue
    },
    "Gorilla": {
        "gene_file": "input/gene_bed/gorilla.genes.bed",
        "imotif_file": "input/imotif_bed/Gorilla_gorilla_all.bed",
        "color": "#96CEB4"  # Green
    },
    "Sumatran Orangutan": {
        "gene_file": "input/gene_bed/sumatran.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_abelii_all.bed",
        "color": "#FECA57"  # Orange
    },
    "Bornean Orangutan": {
        "gene_file": "input/gene_bed/bornean.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_pygmaeus_all.bed",
        "color": "#FF9FF3"  # Pink
    }
}

def calculate_species_enrichment(species_name, species_data):
    """Calculate TSS/TES enrichment for a single species"""
    print(f"Processing {species_name}...")

    try:
        # Read gene data - BED format has 6 columns, we only need columns 1,2,3,6
        genes = pd.read_csv(species_data["gene_file"], sep="\t", header=None,
                          usecols=[0, 1, 2, 5],  # chromosome, start, end, strand
                          names=["chrom", "start", "end", "strand"])

        # Read i-motif data - BED format has 6 columns, we only need columns 1,2,3,6
        imotifs = pd.read_csv(species_data["imotif_file"], sep="\t", header=None,
                            usecols=[0, 1, 2, 5],  # chromosome, start, end, strand
                            names=["chrom", "start", "end", "strand"])

        # Filter valid strands
        genes = genes[genes["strand"].isin(["+", "-"])].copy()
        imotifs = imotifs[imotifs["strand"].isin(["+", "-"])].copy()

        print(f"  {species_name}: {len(genes)} genes, {len(imotifs)} i-motifs")

        # Group i-motifs by strand
        imotif_plus = imotifs[imotifs["strand"] == "+"]
        imotif_minus = imotifs[imotifs["strand"] == "-"]

        # Convert to numpy arrays for speed
        imotif_data = {
            '+': {
                'starts': imotif_plus["start"].values.astype('int32'),
                'ends': imotif_plus["end"].values.astype('int32')
            },
            '-': {
                'starts': imotif_minus["start"].values.astype('int32'),
                'ends': imotif_minus["end"].values.astype('int32')
            }
        }

        # Initialize result arrays
        results = {
            "TSS": {
                "template": np.zeros(len(relative_positions), dtype=np.int32),
                "non_template": np.zeros(len(relative_positions), dtype=np.int32)
            },
            "TES": {
                "template": np.zeros(len(relative_positions), dtype=np.int32),
                "non_template": np.zeros(len(relative_positions), dtype=np.int32)
            }
        }

        # Process each gene
        total_genes = len(genes)
        for idx, gene in genes.iterrows():
            strand = gene["strand"]

            # Determine TSS and TES positions
            if strand == "+":
                tss = gene["start"]
                tes = gene["end"]
            else:  # "-" strand
                tss = gene["end"]
                tes = gene["start"]

            # Determine template and non-template strands
            template_strand = "-" if strand == "+" else "+"

            # Process TSS region
            tmpl_data = imotif_data[template_strand]
            nontmpl_data = imotif_data[strand]

            # Template strand computation - TSS
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TSS"]["template"][i] += count

            # Non-template strand computation - TSS
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TSS"]["non_template"][i] += count

            # Template strand computation - TES
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TES"]["template"][i] += count

            # Non-template strand computation - TES
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TES"]["non_template"][i] += count

            # Progress display (every 10000 genes)
            if (idx + 1) % 10000 == 0:
                print(f"  {species_name}: processed {idx+1}/{total_genes} genes")

        # Normalization
        def normalize(arr):
            total = arr.sum()
            return arr / (total/len(relative_positions)) if total > 0 else arr.astype(np.float64)

        normalized_results = {
            region: {
                "template": normalize(results[region]["template"]),
                "non_template": normalize(results[region]["non_template"]),
                "template_raw": results[region]["template"].copy(),
                "non_template_raw": results[region]["non_template"].copy()
            } for region in results
        }

        # Save results for this species
        species_dir = os.path.join(output_dir, species_name.replace(" ", "_"))
        os.makedirs(species_dir, exist_ok=True)

        for region in ["TSS", "TES"]:
            df = pd.DataFrame({
                "position": relative_positions,
                "template_enrich": normalized_results[region]["template"],
                "non_template_enrich": normalized_results[region]["non_template"],
                "template_raw": normalized_results[region]["template_raw"],
                "non_template_raw": normalized_results[region]["non_template_raw"]
            })
            df.to_csv(os.path.join(species_dir, f"{region}_results.tsv"), sep="\t", index=False)

        return {
            "species": species_name,
            "data": normalized_results,
            "color": species_data["color"],
            "genes": len(genes),
            "imotifs": len(imotifs)
        }

    except Exception as e:
        print(f"Error processing {species_name}: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

def plot_species_comparison(all_results):
    """Generate comparison plots for all species"""

    # Set plot style
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")

    # 1. TSS template strand comparison
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    # TSS template strand
    ax = axes[0]
    for result in all_results:
        if result:
            ax.plot(relative_positions, result["data"]["TSS"]["template"],
                   label=result["species"], color=result["color"], linewidth=2, alpha=0.8)
    ax.set_title("Template Strand Enrichment around TSS", fontsize=14, fontweight='bold')
    ax.set_xlabel("Relative Position to TSS (bp)", fontsize=12)
    ax.set_ylabel("Normalized Enrichment", fontsize=12)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # TSS non-template strand
    ax = axes[1]
    for result in all_results:
        if result:
            ax.plot(relative_positions, result["data"]["TSS"]["non_template"],
                   label=result["species"], color=result["color"], linewidth=2, alpha=0.8)
    ax.set_title("Non-template Strand Enrichment around TSS", fontsize=14, fontweight='bold')
    ax.set_xlabel("Relative Position to TSS (bp)", fontsize=12)
    ax.set_ylabel("Normalized Enrichment", fontsize=12)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # TES template strand
    ax = axes[2]
    for result in all_results:
        if result:
            ax.plot(relative_positions, result["data"]["TES"]["template"],
                   label=result["species"], color=result["color"], linewidth=2, alpha=0.8)
    ax.set_title("Template Strand Enrichment around TES", fontsize=14, fontweight='bold')
    ax.set_xlabel("Relative Position to TES (bp)", fontsize=12)
    ax.set_ylabel("Normalized Enrichment", fontsize=12)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # TES non-template strand
    ax = axes[3]
    for result in all_results:
        if result:
            ax.plot(relative_positions, result["data"]["TES"]["non_template"],
                   label=result["species"], color=result["color"], linewidth=2, alpha=0.8)
    ax.set_title("Non-template Strand Enrichment around TES", fontsize=14, fontweight='bold')
    ax.set_xlabel("Relative Position to TES (bp)", fontsize=12)
    ax.set_ylabel("Normalized Enrichment", fontsize=12)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "all_species_comparison.png"), dpi=300, bbox_inches='tight')
    plt.show()

    # 2. Separate TSS and TES plots, each species individually
    for region in ["TSS", "TES"]:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()

        for idx, result in enumerate(all_results):
            if result and idx < len(axes):
                ax = axes[idx]

                # Template strand
                ax.plot(relative_positions, result["data"][region]["template"],
                       label="Template Strand", color='blue', linewidth=2, alpha=0.8)
                # Non-template strand
                ax.plot(relative_positions, result["data"][region]["non_template"],
                       label="Non-template Strand", color='red', linewidth=2, alpha=0.8)

                ax.set_title(f"{result['species']} - {region}", fontsize=12, fontweight='bold')
                ax.set_xlabel("Relative Position (bp)", fontsize=10)
                ax.set_ylabel("Enrichment", fontsize=10)
                ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)

                # Add statistics info
                stats_text = f"Genes: {result['genes']:,}\ni-motifs: {result['imotifs']:,}"
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                       fontsize=8, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.suptitle(f"i-motif Enrichment around {region} - Primate Species", fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{region}_by_species.png"), dpi=300, bbox_inches='tight')
        plt.show()

def calculate_combined_statistics(all_results):
    """Calculate combined statistics"""
    print("\n" + "="*60)
    print("PRIMATE SPECIES ANALYSIS SUMMARY")
    print("="*60)

    summary_data = []
    for result in all_results:
        if result:
            # Find max enrichment position
            tss_template_max = np.argmax(result["data"]["TSS"]["template"])
            tss_nontemplate_max = np.argmax(result["data"]["TSS"]["non_template"])
            tes_template_max = np.argmax(result["data"]["TES"]["template"])
            tes_nontemplate_max = np.argmax(result["data"]["TES"]["non_template"])

            summary_data.append({
                "Species": result["species"],
                "Genes": f"{result['genes']:,}",
                "i-motifs": f"{result['imotifs']:,}",
                "TSS Template Max": f"{relative_positions[tss_template_max]} bp",
                "TSS Non-template Max": f"{relative_positions[tss_nontemplate_max]} bp",
                "TES Template Max": f"{relative_positions[tes_template_max]} bp",
                "TES Non-template Max": f"{relative_positions[tes_nontemplate_max]} bp"
            })

    # Display table
    df_summary = pd.DataFrame(summary_data)
    print("\nDetailed Statistics:")
    print(df_summary.to_string(index=False))

    # Save as CSV
    df_summary.to_csv(os.path.join(output_dir, "summary_statistics.csv"), index=False)

    # Calculate average enrichment profiles
    print("\n" + "="*60)
    print("AVERAGE ENRICHMENT PROFILES")
    print("="*60)

    for region in ["TSS", "TES"]:
        for strand_type in ["template", "non_template"]:
            # Collect all species data
            all_curves = []
            for result in all_results:
                if result:
                    all_curves.append(result["data"][region][strand_type])

            if all_curves:
                avg_curve = np.mean(all_curves, axis=0)
                std_curve = np.std(all_curves, axis=0)

                # Save average curve
                avg_df = pd.DataFrame({
                    "position": relative_positions,
                    "average_enrichment": avg_curve,
                    "std_dev": std_curve
                })
                avg_df.to_csv(os.path.join(output_dir, f"{region}_{strand_type}_average.tsv"),
                            sep="\t", index=False)

                print(f"\n{region} - {strand_type}:")
                print(f"  Peak position: {relative_positions[np.argmax(avg_curve)]} bp")
                print(f"  Peak value: {avg_curve.max():.4f}")
                print(f"  Average enrichment: {avg_curve.mean():.4f} +/- {std_curve.mean():.4f}")

def main():
    """Main function"""
    print("="*60)
    print("PRIMATE i-MOTIF ENRICHMENT ANALYSIS")
    print("="*60)
    print(f"Species to analyze: {list(species_info.keys())}")
    print(f"Window size: +/-{window_size} bp")
    print(f"Output directory: {output_dir}")
    print("="*60)

    # Process all species in parallel
    print("\nProcessing species in parallel...")
    with Pool(processes=min(n_processes, len(species_info))) as pool:
        args = [(name, data) for name, data in species_info.items()]
        all_results = pool.starmap(calculate_species_enrichment, args)

    # Filter out failed results
    valid_results = [r for r in all_results if r is not None]

    if valid_results:
        print(f"\nSuccessfully processed {len(valid_results)} out of {len(species_info)} species")

        # Generate comparison plots
        print("\nGenerating plots...")
        plot_species_comparison(valid_results)

        # Calculate statistics
        calculate_combined_statistics(valid_results)

        # Generate final report
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        print(f"Results saved in: {os.path.abspath(output_dir)}")
        print("\nGenerated files:")
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if file.endswith(('.png', '.tsv', '.csv')):
                    print(f"  - {os.path.relpath(os.path.join(root, file), output_dir)}")
    else:
        print("\nERROR: No species were successfully processed!")
        return 1

    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
