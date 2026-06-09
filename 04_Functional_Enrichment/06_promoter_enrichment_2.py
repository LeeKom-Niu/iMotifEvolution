"""
计算并绘制pG4在启动子区域的富集图（图5A的人类启动子部分）
修正版本：纵轴从上到下为古老到年轻（Great ape → Homininae → Hominini → Human-specific）
符合 OUP 插图指南：字体 ≥7pt，线条粗细 0.25-1pt，色盲友好配色，PDF 文字可编辑。
已移除标题和 y 轴标签，基线仅保留虚线无文字，纵轴名称简写（去掉"IMs"）。
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
import os
import gzip
def set_oup_style():
    """设置符合 OUP 指南的绘图风格（字体统一加大至16pt）"""
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 16,
        'axes.linewidth': 0.5,
        'axes.edgecolor': 'black',
        'axes.labelpad': 10,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        'legend.fontsize': 12,
        'legend.frameon': False,
        'lines.linewidth': 1.0,
        'figure.dpi': 300,
        'savefig.dpi': 600,
        'figure.figsize': (11, 6.5),
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
set_oup_style()
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
INPUT_DIR = os.path.join(BASE_DIR, "functionalOutputs/Homo_sapiens")
OUTPUT_DIR = os.path.join(BASE_DIR, "enrichment_plots_oup")
os.makedirs(OUTPUT_DIR, exist_ok=True)
COLOR_MAP = {
    "Great ape": "
    "Homininae": "
    "Hominini": "
    "Human-specific": "
}
GROUP_MAP = {
    "hominid": "Great ape",
    "homininae": "Homininae",
    "hominini": "Hominini",
    "humanSpecific": "Human-specific"
}
DISPLAY_GROUP = GROUP_MAP
GENOME_LENGTH = 3117275501
def calculate_promoter_length():
    """计算启动子区域总长度"""
    promoter_file = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/functionalOutputs/Homo_sapiens/promoter_regions.bed.gz"
    
    if not os.path.exists(promoter_file):
        return 18842577
    
    total_length = 0
    with gzip.open(promoter_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                start = int(parts[1])
                end = int(parts[2])
                total_length += (end - start)
    
    return total_length
PROMOTER_LENGTH = calculate_promoter_length()
GROUP_TOTALS = {
    "hominid": 314942,
    "homininae": 124431,
    "hominini": 34964,
    "humanSpecific": 104483
}
def count_pg4_in_promoters():
    """统计每个群组在启动子中的pG4数量，输出使用简写名称"""
    
    counts = {}
    
    for group_short, group_display in GROUP_MAP.items():
        input_file = os.path.join(INPUT_DIR, f"allhsaG.intersected.betn.human_promoter.{group_short}G4s.bed.gz")
        
        if not os.path.exists(input_file):
            counts[group_display] = 0
            continue
        
        pg4_ids = set()
        with gzip.open(input_file, 'rt') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pg4_id = f"{parts[0]}:{parts[1]}-{parts[2]}:{parts[5]}"
                    pg4_ids.add(pg4_id)
        
        counts[group_display] = len(pg4_ids)
    
    return counts
def calculate_fold_enrichment(promoter_counts):
    """
    计算富集倍数
    """
    
    results = []
    
    for group_display, promoter_count in promoter_counts.items():
        group_short = [k for k, v in GROUP_MAP.items() if v == group_display][0]
        total_pg4 = GROUP_TOTALS[group_short]
        
        pG4_ratio = promoter_count / total_pg4
        promoter_ratio = PROMOTER_LENGTH / GENOME_LENGTH
        fold_enrichment = pG4_ratio / promoter_ratio
        
        percentage_in_promoter = (promoter_count / total_pg4) * 100
        percentage_genome = (PROMOTER_LENGTH / GENOME_LENGTH) * 100
        
        results.append({
            'Group': group_display,
            'Short_Name': group_short,
            'Promoter_Count': promoter_count,
            'Total_pG4': total_pg4,
            '%_in_Promoter': percentage_in_promoter,
            '%_Genome_Promoter': percentage_genome,
            'Fold_Enrichment': fold_enrichment,
            'Color': COLOR_MAP[group_display]
        })
        
    
    return pd.DataFrame(results)
def plot_enrichment(results_df):
    """
    绘制水平棒棒糖图
    纵轴顺序：从上到下为 Great ape → Homininae → Hominini → Human-specific
    """
    
    order = ["Great ape", "Homininae", "Hominini", "Human-specific"]
    results_df['Order'] = results_df['Group'].map({g: i for i, g in enumerate(order)})
    results_df = results_df.sort_values('Order', ascending=True)
    
    fig, ax = plt.subplots(figsize=(11, 6.5))
    
    y_pos = np.arange(len(results_df))
    
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax.hlines(
            y=i, 
            xmin=0, 
            xmax=row['Fold_Enrichment'],
            color=row['Color'],
            linewidth=1.0,
            alpha=0.7
        )
    
    ax.scatter(
        results_df['Fold_Enrichment'], 
        y_pos,
        c=results_df['Color'],
        s=200,
        edgecolor='black',
        linewidth=0.5,
        zorder=5
    )
    
    for i, (_, row) in enumerate(results_df.iterrows()):
        label_text = f"{row['Promoter_Count']:,} ({row['Fold_Enrichment']:.2f}x)"
        ax.text(
            row['Fold_Enrichment'] + 0.2, 
            i,
            label_text,
            va='center',
            ha='left',
            fontsize=12,
            fontweight='normal',
            color=row['Color']
        )
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Group'], fontsize=16)
    ax.set_ylabel('')
    ax.invert_yaxis()
    
    ax.set_xlabel('Fold Enrichment', fontsize=16)
    ax.set_title('')
    
    ax.axvline(x=1, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    
    ax.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.3)
    ax.set_axisbelow(True)
    
    x_max = results_df['Fold_Enrichment'].max() * 1.45
    ax.set_xlim([0, x_max])
    
    legend_elements = []
    for _, row in results_df.iterrows():
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', 
                   markerfacecolor=row['Color'], markersize=12,
                   label=f"{row['Group']}\n{row['Promoter_Count']:,} in promoters")
        )
    ax.legend(handles=legend_elements, loc='lower right', fontsize=12, frameon=False)
    
    plt.tight_layout()
    
    pdf_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.pdf")
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    
    tiff_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.tiff")
    plt.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    
    png_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    svg_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    
    plt.close()
    
    return fig, ax
def compare_with_paper(results_df):
    """与原文图5A比较（保持不变）"""
    
    paper_values = {
        "Great ape": {"fold": 6.71, "count": "~?"},
        "Homininae": {"fold": 5.77, "count": "~?"},
        "Hominini": {"fold": 4.67, "count": "~?"},
        "Human-specific": {"fold": 4.11, "count": 2753}
    }
    
    
    for _, row in results_df.iterrows():
        group = row['Group']
        your_fold = row['Fold_Enrichment']
        your_count = row['Promoter_Count']
        
        if group in paper_values:
            paper_fold = paper_values[group]['fold']
            paper_count = paper_values[group]['count']
            diff_fold = your_fold - paper_fold
            diff_percent = (diff_fold / paper_fold) * 100
        else:
def save_results(results_df):
    """保存计算结果（保持不变）"""
    
    output_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_results.csv")
    results_df.to_csv(output_csv, index=False)
    
    summary_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_summary.csv")
    summary = results_df[['Group', 'Promoter_Count', 'Total_pG4', '%_in_Promoter', 'Fold_Enrichment']].copy()
    summary.to_csv(summary_csv, index=False)
    
    
    for _, row in results_df.iterrows():
def main():
    """主函数"""
    
    promoter_counts = count_pg4_in_promoters()
    results_df = calculate_fold_enrichment(promoter_counts)
    plot_enrichment(results_df)
    compare_with_paper(results_df)
    save_results(results_df)
    
if __name__ == "__main__":
    main()
