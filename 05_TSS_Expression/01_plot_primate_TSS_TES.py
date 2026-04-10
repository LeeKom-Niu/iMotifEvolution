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

# 设置参数
window_size = 1000  # 上下游范围
relative_positions = np.arange(-window_size, window_size + 1)
output_dir = "primate_TSS_TES_enrichment_results"
n_processes = 6  # 每个物种一个进程

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

# 物种映射表 (显示名: 文件名)
species_info = {
    "Human": {
        "gene_file": "input/gene_bed/human.genes.bed",
        "imotif_file": "input/imotif_bed/Homo_sapiens_all.bed",
        "color": "#FF6B6B"  # 红色
    },
    "Chimpanzee": {
        "gene_file": "input/gene_bed/chimp.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_troglodytes_all.bed",
        "color": "#4ECDC4"  # 青色
    },
    "Bonobo": {
        "gene_file": "input/gene_bed/bonobo.genes.bed",
        "imotif_file": "input/imotif_bed/Pan_paniscus_all.bed",
        "color": "#45B7D1"  # 蓝色
    },
    "Gorilla": {
        "gene_file": "input/gene_bed/gorilla.genes.bed",
        "imotif_file": "input/imotif_bed/Gorilla_gorilla_all.bed",
        "color": "#96CEB4"  # 绿色
    },
    "Sumatran Orangutan": {
        "gene_file": "input/gene_bed/sumatran.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_abelii_all.bed",
        "color": "#FECA57"  # 橙色
    },
    "Bornean Orangutan": {
        "gene_file": "input/gene_bed/bornean.genes.bed",
        "imotif_file": "input/imotif_bed/Pongo_pygmaeus_all.bed",
        "color": "#FF9FF3"  # 粉色
    }
}

def calculate_species_enrichment(species_name, species_data):
    """计算单个物种的TSS/TES富集度"""
    print(f"Processing {species_name}...")
    
    try:
        # 读取基因数据 - BED格式有6列，我们只需要第1,2,3,6列
        genes = pd.read_csv(species_data["gene_file"], sep="\t", header=None, 
                          usecols=[0, 1, 2, 5],  # 染色体，起始，终止，链
                          names=["chrom", "start", "end", "strand"])
        
        # 读取i-motif数据 - BED格式有6列，我们只需要第1,2,3,6列
        imotifs = pd.read_csv(species_data["imotif_file"], sep="\t", header=None,
                            usecols=[0, 1, 2, 5],  # 染色体，起始，终止，链
                            names=["chrom", "start", "end", "strand"])
        
        # 过滤有效链
        genes = genes[genes["strand"].isin(["+", "-"])].copy()
        imotifs = imotifs[imotifs["strand"].isin(["+", "-"])].copy()
        
        print(f"  {species_name}: {len(genes)} genes, {len(imotifs)} i-motifs")
        
        # 按链分组i-motif
        imotif_plus = imotifs[imotifs["strand"] == "+"]
        imotif_minus = imotifs[imotifs["strand"] == "-"]
        
        # 转换为numpy数组提高速度
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
        
        # 初始化结果数组
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
        
        # 处理每个基因
        total_genes = len(genes)
        for idx, gene in genes.iterrows():
            strand = gene["strand"]
            
            # 确定TSS和TES位置
            if strand == "+":
                tss = gene["start"]
                tes = gene["end"]
            else:  # "-"链
                tss = gene["end"]
                tes = gene["start"]
            
            # 确定模板链和非模板链
            template_strand = "-" if strand == "+" else "+"
            
            # 处理TSS区域
            tmpl_data = imotif_data[template_strand]
            nontmpl_data = imotif_data[strand]
            
            # 模板链计算 - TSS
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TSS"]["template"][i] += count
            
            # 非模板链计算 - TSS
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tss + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TSS"]["non_template"][i] += count
            
            # 模板链计算 - TES
            if len(tmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((tmpl_data['starts'] <= pos) & (tmpl_data['ends'] >= pos))
                    results["TES"]["template"][i] += count
            
            # 非模板链计算 - TES
            if len(nontmpl_data['starts']) > 0:
                for i, rel_pos in enumerate(relative_positions):
                    pos = tes + rel_pos
                    count = np.sum((nontmpl_data['starts'] <= pos) & (nontmpl_data['ends'] >= pos))
                    results["TES"]["non_template"][i] += count
            
            # 显示进度
            if (idx + 1) % 10000 == 0:
                print(f"  {species_name}: processed {idx+1}/{total_genes} genes")
        
        # 归一化处理
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
        
        # 保存该物种的结果
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
    """绘制所有物种的对比图"""
    
    # 设置绘图样式
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_palette("husl")
    
    # 1. TSS模板链对比
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    # TSS模板链
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
    
    # TSS非模板链
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
    
    # TES模板链
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
    
    # TES非模板链
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
    
    # 2. 分开绘制TSS和TES，每个物种单独显示
    for region in ["TSS", "TES"]:
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        axes = axes.flatten()
        
        for idx, result in enumerate(all_results):
            if result and idx < len(axes):
                ax = axes[idx]
                
                # 模板链
                ax.plot(relative_positions, result["data"][region]["template"],
                       label="Template Strand", color='blue', linewidth=2, alpha=0.8)
                # 非模板链
                ax.plot(relative_positions, result["data"][region]["non_template"],
                       label="Non-template Strand", color='red', linewidth=2, alpha=0.8)
                
                ax.set_title(f"{result['species']} - {region}", fontsize=12, fontweight='bold')
                ax.set_xlabel("Relative Position (bp)", fontsize=10)
                ax.set_ylabel("Enrichment", fontsize=10)
                ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
                ax.legend(fontsize=9)
                ax.grid(True, alpha=0.3)
                
                # 添加统计信息
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
    print("\n" + "="*60)
    print("PRIMATE SPECIES ANALYSIS SUMMARY")
    print("="*60)
    
    summary_data = []
    for result in all_results:
        if result:
            # 计算最大富集值位置
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
    
    # 显示表格
    df_summary = pd.DataFrame(summary_data)
    print("\nDetailed Statistics:")
    print(df_summary.to_string(index=False))
    
    # 保存为CSV
    df_summary.to_csv(os.path.join(output_dir, "summary_statistics.csv"), index=False)
    
    # 计算平均富集曲线
    print("\n" + "="*60)
    print("AVERAGE ENRICHMENT PROFILES")
    print("="*60)
    
    for region in ["TSS", "TES"]:
        for strand_type in ["template", "non_template"]:
            # 收集所有物种的数据
            all_curves = []
            for result in all_results:
                if result:
                    all_curves.append(result["data"][region][strand_type])
            
            if all_curves:
                avg_curve = np.mean(all_curves, axis=0)
                std_curve = np.std(all_curves, axis=0)
                
                # 保存平均曲线
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
                print(f"  Average enrichment: {avg_curve.mean():.4f} ± {std_curve.mean():.4f}")

def main():
    """主函数"""
    print("="*60)
    print("PRIMATE i-MOTIF ENRICHMENT ANALYSIS")
    print("="*60)
    print(f"Species to analyze: {list(species_info.keys())}")
    print(f"Window size: ±{window_size} bp")
    print(f"Output directory: {output_dir}")
    print("="*60)
    
    # 并行处理所有物种
    print("\nProcessing species in parallel...")
    with Pool(processes=min(n_processes, len(species_info))) as pool:
        args = [(name, data) for name, data in species_info.items()]
        all_results = pool.starmap(calculate_species_enrichment, args)
    
    # 过滤掉失败的结果
    valid_results = [r for r in all_results if r is not None]
    
    if valid_results:
        print(f"\nSuccessfully processed {len(valid_results)} out of {len(species_info)} species")
        
        # 绘制比较图
        print("\nGenerating plots...")
        plot_species_comparison(valid_results)
        
        # 计算统计信息
        calculate_combined_statistics(valid_results)
        
        # 生成最终报告
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
