
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')
window_size = 1000
relative_positions = np.arange(-window_size, window_size + 1)
output_dir = "primate_TSS_TES_enrichment_results"
n_processes = 6
os.makedirs(output_dir, exist_ok=True)
species_info = {
    "Human": {
        "gene_file": "input/gene_bed/human.genes.bed",
        "imotif_file": "input/imotif_bed/Homo_sapiens_all.bed",
        "color": "
    },
    "Chimpanzee": {
        "gene_file": "input/gene_bed/chimp.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_troglodytes_all.bed",
        "color": "
    },
    "Bonobo": {
        "gene_file": "input/gene_bed/bonobo.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_paniscus_all.bed",
        "color": "
    },
    "Gorilla": {
        "gene_file": "input/gene_bed/gorilla.genes.bed",
        "imotif_file": "input/imotif_bed/Gorilla_gorilla_all.bed",
        "color": "
    },
    "Sumatran Orangutan": {
        "gene_file": "input/gene_bed/sumatran.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_abelii_all.bed",
        "color": "
    },
    "Bornean Orangutan": {
        "gene_file": "input/gene_bed/bornean.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_pygmaeus_all.bed",
        "color": "
    }
}
def calculate_species_enrichment(species_name, species_data):
    """计算单个物种的TSS/TES富集度"""
    
    try:
        genes = pd.read_csv(species_data["gene_file"], sep="\t", header=None, 
                          usecols=[0, 1, 2, 5],
                          names=["chrom", "start", "end", "strand"])
        
        imotifs = pd.read_csv(species_data["imotif_file"], sep="\t", header=None,
                            usecols=[0, 1, 2, 5],
                            names=["chrom", "start", "end", "strand"])
        
        genes = genes[genes["strand"].isin(["+", "-"])].copy()
        imotifs = imotifs[imotifs["strand"].isin(["+", "-"])].copy()
        
        
        imotif_plus = imotifs[imotifs["strand"] == "+"]
        imotif_minus = imotifs[imotifs["strand"] == "-"]
        
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
        
        total_genes = len(genes)
        for idx, gene in genes.iterrows():
            strand = gene["strand"]
            
            if strand == "+":
                tss = gene["start"]
                tes = gene["end"]
            else:
                tss = gene["end"]
                tes = gene["start"]
            
            template_strand = "-" if strand == "+" else "+"
            
            tmpl_data = imotif_data[template_strand]
            nontmpl_data = imotif_data[strand]
            
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TSS"]["template"][i] += count
            
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TSS"]["non_template"][i] += count
            
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TES"]["template"][i] += count
            
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TES"]["non_template"][i] += count
            
            if (idx + 1) % 10000 == 0:
        
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
        import traceback
        traceback.print_exc()
        return None
def plot_species_comparison(all_results):
    """绘制所有物种的对比图"""
    
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
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
    
    for region in ["TSS", "TES"]:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()
        
        for idx, result in enumerate(all_results):
            if result and idx < len(axes):
                ax = axes[idx]
                
                ax.plot(relative_positions, result["data"][region]["template"],
                       label="Template Strand", color='blue', linewidth=2, alpha=0.8)
                ax.plot(relative_positions, result["data"][region]["non_template"],
                       label="Non-template Strand", color='red', linewidth=2, alpha=0.8)
                
                ax.set_title(f"{result['species']} - {region}", fontsize=12, fontweight='bold')
                ax.set_xlabel("Relative Position (bp)", fontsize=10)
                ax.set_ylabel("Enrichment", fontsize=10)
                ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                
                stats_text = f"Genes: {result['genes']:,}\ni-motifs: {result['imotifs']:,}"
                ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                       fontsize=8, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.suptitle(f"i-motif Enrichment around {region} - Primate Species", fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{region}_by_species.png"), dpi=300, bbox_inches='tight')
        plt.show()
def calculate_combined_statistics(all_results):
    """计算合并统计信息"""
    
    summary_data = []
    for result in all_results:
        if result:
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
    
    df_summary = pd.DataFrame(summary_data)
    
    df_summary.to_csv(os.path.join(output_dir, "summary_statistics.csv"), index=False)
    
    
    for region in ["TSS", "TES"]:
        for strand_type in ["template", "non_template"]:
            all_curves = []
            for result in all_results:
                if result:
                    all_curves.append(result["data"][region][strand_type])
            
            if all_curves:
                avg_curve = np.mean(all_curves, axis=0)
                std_curve = np.std(all_curves, axis=0)
                
                avg_df = pd.DataFrame({
                    "position": relative_positions,
                    "average_enrichment": avg_curve,
                    "std_dev": std_curve
                })
                avg_df.to_csv(os.path.join(output_dir, f"{region}_{strand_type}_average.tsv"), 
                            sep="\t", index=False)
                
def main():
    """主函数"""
    
    with Pool(processes=min(n_processes, len(species_info))) as pool:
        args = [(name, data) for name, data in species_info.items()]
        all_results = pool.starmap(calculate_species_enrichment, args)
    
    valid_results = [r for r in all_results if r is not None]
    
    if valid_results:
        
        plot_species_comparison(valid_results)
        
        calculate_combined_statistics(valid_results)
        
        for root, dirs, files in os.walk(output_dir):
            for file in files:
                if file.endswith(('.png', '.tsv', '.csv')):
    else:
        return 1
    
    return 0
if __name__ == "__main__":
    import sys
    sys.exit(main())
