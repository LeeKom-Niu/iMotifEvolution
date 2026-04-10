#!/usr/bin/env python3
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

# ===== OUP 插图指南风格设置（字体统一加大） =====
def set_oup_style():
    """设置符合 OUP 指南的绘图风格（字体统一加大至16pt）"""
    rcParams.update({
        # 字体设置
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 14,                     # 基础字体（稍大）
        # 坐标轴
        'axes.labelsize': 16,                 # 轴标签（Fold Enrichment）
        'axes.titlesize': 16,                 # 标题（已移除）
        'axes.linewidth': 0.5,                 # 轴线粗细
        'axes.edgecolor': 'black',
        'axes.labelpad': 10,
        # 刻度
        'xtick.labelsize': 16,                 # x 轴刻度标签
        'ytick.labelsize': 16,                 # y 轴刻度标签（组名）
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        # 图例
        'legend.fontsize': 12,                  # 图例文字
        'legend.frameon': False,
        # 线条
        'lines.linewidth': 1.0,                # 水平线宽
        # 图形尺寸
        'figure.dpi': 300,
        'savefig.dpi': 600,
        'figure.figsize': (11, 6.5),            # 稍增大以容纳大字体
        # PDF 文字可编辑
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

set_oup_style()

# ===== 配置参数 =====
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
INPUT_DIR = os.path.join(BASE_DIR, "functionalOutputs/Homo_sapiens")
OUTPUT_DIR = os.path.join(BASE_DIR, "enrichment_plots_oup")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 颜色映射（保留）
COLOR_MAP = {
    "Great ape": "#004d71",          # 深蓝
    "Homininae": "#8759a1",          # 紫色
    "Hominini": "#f75a78",            # 粉色
    "Human-specific": "#ffa600"       # 橙色
}

# 群组映射（注意去掉 IMs）
GROUP_MAP = {
    "hominid": "Great ape",
    "homininae": "Homininae",
    "hominini": "Hominini",
    "humanSpecific": "Human-specific"
}

# 为方便内部处理，保留完整名称与颜色的映射，但显示时用简写
DISPLAY_GROUP = GROUP_MAP

# ===== 基因组基本信息 =====
GENOME_LENGTH = 3117275501

def calculate_promoter_length():
    """计算启动子区域总长度"""
    promoter_file = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset/GreatApeT2T-G4s-main/datasets/functionalOutputs/Homo_sapiens/promoter_regions.bed.gz"
    
    if not os.path.exists(promoter_file):
        print(f"警告: 启动子文件不存在 {promoter_file}")
        print("使用作者代码中的值: 18,842,577 bp")
        return 18842577
    
    total_length = 0
    with gzip.open(promoter_file, 'rt') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                start = int(parts[1])
                end = int(parts[2])
                total_length += (end - start)
    
    print(f"启动子总长度: {total_length:,} bp")
    return total_length

PROMOTER_LENGTH = calculate_promoter_length()

# ===== 各群组pG4总数 =====
GROUP_TOTALS = {
    "hominid": 314942,
    "homininae": 124431,
    "hominini": 34964,
    "humanSpecific": 104483
}

# ===== 主要函数 =====
def count_pg4_in_promoters():
    """统计每个群组在启动子中的pG4数量，输出使用简写名称"""
    print("=== 统计启动子中的pG4数量 ===")
    
    counts = {}
    
    for group_short, group_display in GROUP_MAP.items():
        input_file = os.path.join(INPUT_DIR, f"allhsaG.intersected.betn.human_promoter.{group_short}G4s.bed.gz")
        
        if not os.path.exists(input_file):
            print(f"警告: 文件不存在 {input_file}")
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
        print(f"  {group_display}: {len(pg4_ids):,} 个唯一pG4")
    
    return counts

def calculate_fold_enrichment(promoter_counts):
    """
    计算富集倍数
    """
    print("\n=== 计算富集倍数 ===")
    
    results = []
    
    for group_display, promoter_count in promoter_counts.items():
        group_short = [k for k, v in GROUP_MAP.items() if v == group_display][0]
        total_pg4 = GROUP_TOTALS[group_short]
        
        pG4_ratio = promoter_count / total_pg4
        promoter_ratio = PROMOTER_LENGTH / GENOME_LENGTH
        fold_enrichment = pG4_ratio / promoter_ratio
        
        percentage_in_promoter = (promoter_count / total_pg4) * 100
        percentage_genome = (PROMOTER_LENGTH / GENOME_LENGTH) * 100
        
        # 使用简写名称作为组名
        results.append({
            'Group': group_display,
            'Short_Name': group_short,
            'Promoter_Count': promoter_count,
            'Total_pG4': total_pg4,
            '%_in_Promoter': percentage_in_promoter,
            '%_Genome_Promoter': percentage_genome,
            'Fold_Enrichment': fold_enrichment,
            'Color': COLOR_MAP[group_display]  # 颜色直接用简写对应的颜色
        })
        
        print(f"  {group_display}:")
        print(f"    启动子中pG4: {promoter_count:,} / {total_pg4:,} = {percentage_in_promoter:.2f}%")
        print(f"    启动子占基因组: {PROMOTER_LENGTH:,} / {GENOME_LENGTH:,} = {percentage_genome:.4f}%")
        print(f"    富集倍数: {fold_enrichment:.2f}")
    
    return pd.DataFrame(results)

def plot_enrichment(results_df):
    """
    绘制水平棒棒糖图
    纵轴顺序：从上到下为 Great ape → Homininae → Hominini → Human-specific
    """
    print("\n=== 绘制富集图（水平棒棒糖图） ===")
    
    # 按演化顺序排序
    order = ["Great ape", "Homininae", "Hominini", "Human-specific"]
    results_df['Order'] = results_df['Group'].map({g: i for i, g in enumerate(order)})
    results_df = results_df.sort_values('Order', ascending=True)
    
    fig, ax = plt.subplots(figsize=(11, 6.5))
    
    y_pos = np.arange(len(results_df))
    
    # 绘制水平线
    for i, (_, row) in enumerate(results_df.iterrows()):
        ax.hlines(
            y=i, 
            xmin=0, 
            xmax=row['Fold_Enrichment'],
            color=row['Color'],
            linewidth=1.0,
            alpha=0.7
        )
    
    # 绘制点
    ax.scatter(
        results_df['Fold_Enrichment'], 
        y_pos,
        c=results_df['Color'],
        s=200,
        edgecolor='black',
        linewidth=0.5,
        zorder=5
    )
    
    # 添加数量标签（字体稍大）
    for i, (_, row) in enumerate(results_df.iterrows()):
        label_text = f"{row['Promoter_Count']:,} ({row['Fold_Enrichment']:.2f}x)"
        ax.text(
            row['Fold_Enrichment'] + 0.2, 
            i,
            label_text,
            va='center',
            ha='left',
            fontsize=12,          # 数值标签12pt
            fontweight='normal',
            color=row['Color']
        )
    
    # 设置Y轴刻度标签（简写名称）
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results_df['Group'], fontsize=16)  # 与刻度标签大小一致
    ax.set_ylabel('')                 # 无 y 轴标签
    ax.invert_yaxis()
    
    # 设置X轴（字体由 axes.labelsize 控制，这里显式指定确保一致）
    ax.set_xlabel('Fold Enrichment', fontsize=16)
    ax.set_title('')
    
    # 基线虚线
    ax.axvline(x=1, color='black', linestyle='--', linewidth=0.5, alpha=0.5)
    
    # 网格线
    ax.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.3)
    ax.set_axisbelow(True)
    
    # X轴范围
    x_max = results_df['Fold_Enrichment'].max() * 1.45
    ax.set_xlim([0, x_max])
    
    # 图例（使用简写名称）
    legend_elements = []
    for _, row in results_df.iterrows():
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', 
                   markerfacecolor=row['Color'], markersize=12,
                   label=f"{row['Group']}\n{row['Promoter_Count']:,} in promoters")
        )
    ax.legend(handles=legend_elements, loc='lower right', fontsize=12, frameon=False)
    
    plt.tight_layout()
    
    # 保存文件
    pdf_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.pdf")
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    print(f"  PDF 已保存: {pdf_path}")
    
    tiff_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.tiff")
    plt.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight',
                pil_kwargs={'compression': 'tiff_lzw'})
    print(f"  TIFF 已保存: {tiff_path}")
    
    png_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"  PNG 已保存: {png_path}")
    
    svg_path = os.path.join(OUTPUT_DIR, "promoter_enrichment_human_lollipop_oup.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    print(f"  SVG 已保存: {svg_path}")
    
    plt.close()
    
    return fig, ax

def compare_with_paper(results_df):
    """与原文图5A比较（保持不变）"""
    print("\n=== 与原文图5A比较 ===")
    
    paper_values = {
        "Great ape": {"fold": 6.71, "count": "~?"},
        "Homininae": {"fold": 5.77, "count": "~?"},
        "Hominini": {"fold": 4.67, "count": "~?"},
        "Human-specific": {"fold": 4.11, "count": 2753}
    }
    
    print(f"{'Group':<20} {'Your Fold':<12} {'Paper Fold':<12} {'Difference':<12} {'Your Count':<12} {'Paper Count':<12}")
    print("-" * 85)
    
    for _, row in results_df.iterrows():
        group = row['Group']
        your_fold = row['Fold_Enrichment']
        your_count = row['Promoter_Count']
        
        if group in paper_values:
            paper_fold = paper_values[group]['fold']
            paper_count = paper_values[group]['count']
            diff_fold = your_fold - paper_fold
            diff_percent = (diff_fold / paper_fold) * 100
            print(f"{group:<20} {your_fold:<12.2f} {paper_fold:<12.2f} {diff_percent:<12.1f}% {your_count:<12} {paper_count:<12}")
        else:
            print(f"{group:<20} {your_fold:<12.2f} {'N/A':<12} {'N/A':<12} {your_count:<12} {'N/A':<12}")

def save_results(results_df):
    """保存计算结果（保持不变）"""
    print("\n=== 保存计算结果 ===")
    
    output_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_results.csv")
    results_df.to_csv(output_csv, index=False)
    print(f"  详细结果: {output_csv}")
    
    summary_csv = os.path.join(OUTPUT_DIR, "promoter_enrichment_summary.csv")
    summary = results_df[['Group', 'Promoter_Count', 'Total_pG4', '%_in_Promoter', 'Fold_Enrichment']].copy()
    summary.to_csv(summary_csv, index=False)
    print(f"  摘要结果: {summary_csv}")
    
    print("\n富集分析摘要:")
    print("-" * 80)
    print(f"{'Group':<20} {'Promoter':<12} {'Total':<12} {'% in Promoter':<15} {'Fold Enrichment':<15}")
    print("-" * 80)
    
    for _, row in results_df.iterrows():
        print(f"{row['Group']:<20} {row['Promoter_Count']:<12,} {row['Total_pG4']:<12,} {row['%_in_Promoter']:<15.2f} {row['Fold_Enrichment']:<15.2f}")

def main():
    """主函数"""
    print("=" * 80)
    print("pG4启动子富集分析 (人类) - 水平棒棒糖图版 (字体统一16pt，纵轴简写)")
    print("=" * 80)
    
    promoter_counts = count_pg4_in_promoters()
    results_df = calculate_fold_enrichment(promoter_counts)
    plot_enrichment(results_df)
    compare_with_paper(results_df)
    save_results(results_df)
    
    print("\n" + "=" * 80)
    print(f"分析完成! 结果保存在: {OUTPUT_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
