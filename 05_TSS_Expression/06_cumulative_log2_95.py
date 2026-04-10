#!/usr/bin/env python3
"""
针对 log₂(TPM+1) 数据，尝试剔除高表达端极端值（仅剔除大于某百分位数的基因），
并绘制 ECDF 图，观察不同阈值下“有 i-motif”与“无 i-motif”两组的差异。
物种显示顺序：human, bonobo, chimp, sumatran, gorilla。
依赖：已生成的 binary 版本 *_imotif_analysis_data.csv 文件。
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import rcParams

# ==================== 配置 ====================
BASE_DIR = os.getcwd()
DATA_DIR = os.path.join(BASE_DIR, "imotif_analysis_results_binary")
OUTPUT_DIR = os.path.join(BASE_DIR, "imotif_analysis_results_binary", "ecdf_trim_high")

ALL_SPECIES = ["sumatran", "gorilla", "human", "chimp", "bonobo"]
DISPLAY_ORDER = ["human", "bonobo", "chimp", "sumatran", "gorilla"]  # 绘图顺序

GROUPS = ['With i-Motif', 'Without i-Motif']
COLORS = {'With i-Motif': '#E64B35', 'Without i-Motif': '#4DBBD5'}
LINEWIDTH = 1.5

# 要尝试的百分位数阈值（剔除大于该百分位数的值）
THRESHOLDS = [95]  # 可自行增删

# 图形风格
def set_style():
    plt.style.use('default')
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 7,
        'axes.linewidth': 0.5,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def load_species_data(species):
    """读取指定物种的数据文件，返回 DataFrame（包含 log2_TPM 列）"""
    data_file = os.path.join(DATA_DIR, species, f"{species}_imotif_analysis_data.csv")
    if not os.path.exists(data_file):
        print(f"  ✗ 数据文件不存在: {data_file}")
        return None
    df = pd.read_csv(data_file)
    required = ['GeneID', 'group', 'log2_TPM']
    if not all(col in df.columns for col in required):
        print(f"  ✗ 数据文件缺少必要列: {data_file}")
        return None
    # 过滤表达量为0的基因（log2_TPM > 0 即 TPM>0）
    df = df[df['log2_TPM'] > 0].copy()
    if df.empty:
        print(f"  ✗ {species}: 无表达量>0的基因")
        return None
    return df

def trim_high_outliers(df, percentile):
    """
    剔除大于指定百分位数的值（基于 log2_TPM 列）
    percentile: 百分位数，如 95 表示剔除 > 95% 分位数的值
    返回新的 DataFrame
    """
    all_vals = df['log2_TPM']
    threshold = np.percentile(all_vals, percentile)
    df_trim = df[all_vals <= threshold].copy()
    return df_trim

def plot_ecdf(data_dict, percentile, x_label="log₂(TPM+1)"):
    """
    绘制一行五列的 ECDF 图
    data_dict: 物种 -> DataFrame（必须包含 'log2_TPM' 列）
    percentile: 用于文件名的阈值
    """
    n_species = len(DISPLAY_ORDER)
    fig, axes = plt.subplots(1, n_species, figsize=(15, 4))

    for idx, species in enumerate(DISPLAY_ORDER):
        ax = axes[idx]
        df = data_dict.get(species)
        if df is None or df.empty:
            ax.text(0.5, 0.5, f"No data for {species}", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(species.capitalize())
            continue

        df_with = df[df['group'] == 'With i-Motif']['log2_TPM']
        df_without = df[df['group'] == 'Without i-Motif']['log2_TPM']
        if len(df_with)==0 or len(df_without)==0:
            ax.text(0.5, 0.5, "Insufficient data", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(species.capitalize())
            continue

        # 计算 Mann-Whitney p 值
        _, p_val = stats.mannwhitneyu(df_with, df_without)

        all_data = df['log2_TPM']
        x_min = all_data.min()
        x_max = all_data.max()  # 此时 x_max 就是剔除后的最大值

        # 绘制 ECDF 曲线
        for grp in GROUPS:
            data = df[df['group'] == grp]['log2_TPM']
            if len(data) > 0:
                x = np.sort(data)
                y = np.arange(1, len(x)+1) / len(x)
                ax.step(x, y, where='post', label=f"{grp} (n={len(x)})",
                        color=COLORS[grp], linewidth=LINEWIDTH)

        # 中位数垂直线
        median_with = np.median(df_with) if len(df_with) else np.nan
        median_without = np.median(df_without) if len(df_without) else np.nan
        if not np.isnan(median_with):
            ax.axvline(median_with, color=COLORS['With i-Motif'], linestyle='--', linewidth=1, alpha=0.7)
        if not np.isnan(median_without):
            ax.axvline(median_without, color=COLORS['Without i-Motif'], linestyle='--', linewidth=1, alpha=0.7)

        ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.set_xlim(left=x_min, right=x_max)
        ax.set_ylim(0, 1)
        ax.set_xlabel(x_label)
        ax.set_ylabel('Cumulative frequency')
        ax.set_title(species.capitalize())

        ax.legend(loc='upper left', frameon=False, fontsize=6)
        ax.text(0.95, 0.05, f'p = {p_val:.2e}', transform=ax.transAxes, fontsize=7,
                ha='right', va='bottom', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.tight_layout()
    fname = f"ecdf_log2_trim_high_{percentile}.pdf"
    out_path = os.path.join(OUTPUT_DIR, fname)
    plt.savefig(out_path, format='pdf')
    plt.savefig(out_path.replace('.pdf', '.tiff'), format='tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(out_path.replace('.pdf', '.png'), format='png', dpi=300)
    plt.close()
    print(f"  已保存: {out_path}")

def main():
    set_style()
    ensure_dir(OUTPUT_DIR)

    # 1. 加载原始数据
    print("加载原始数据...")
    raw_data = {}
    for species in ALL_SPECIES:
        df = load_species_data(species)
        if df is not None:
            raw_data[species] = df
    if not raw_data:
        print("错误：未加载任何数据。")
        return

    # 2. 对每个阈值进行剔除并绘图
    for percentile in THRESHOLDS:
        print(f"\n处理阈值: 剔除 > {percentile}% 分位数")
        trimmed_data = {}
        for species, df in raw_data.items():
            df_trim = trim_high_outliers(df, percentile)
            trimmed_data[species] = df_trim
            print(f"  {species}: 原始 {len(df)} 个基因, 剔除后 {len(df_trim)} 个基因")
        plot_ecdf(trimmed_data, percentile, x_label="log₂(TPM+1)")

    print(f"\n所有结果已保存至: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
