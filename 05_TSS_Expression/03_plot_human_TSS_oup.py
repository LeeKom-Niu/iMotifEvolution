#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_human_TSS_oup.py - 绘制人类 TSS 区域 i-Motif 富集曲线（符合 OUP 插图指南）
修复：增加边距、字体嵌入设置，确保 PDF 在 AI 中完整显示。
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# ========== OUP 插图指南风格设置 ==========
plt.rcParams.update({
    # 字体（基础大小 8pt，所有文本 ≥7pt）
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,

    # 线条粗细（严格控制在 0.25-1 pt）
    'lines.linewidth': 0.8,
    'axes.linewidth': 0.5,
    'grid.linewidth': 0.3,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,

    # PDF 字体可编辑
    'pdf.fonttype': 42,
    'ps.fonttype': 42,

    # 其他
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.linestyle': ':',
    'grid.alpha': 0.3,
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})

# ========== 色盲友好配色 ==========
COLOR_TEMPLATE = '#1F77B4'      # 蓝色
COLOR_NONTEMPLATE = '#D62728'   # 红色

# ========== 参数设置 ==========
DATA_DIR = "primate_TSS_TES_enrichment_results"
SPECIES = "Human"
REGION = "TSS"
SMOOTH_SIGMA = 60
WINDOW_SIZE = 1000
OUTPUT_DIR = "human_TSS_figure"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_human_tss_data():
    """加载人类 TSS 数据"""
    file_path = os.path.join(DATA_DIR, SPECIES, f"{REGION}_results.tsv")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"数据文件不存在：{file_path}")
    df = pd.read_csv(file_path, sep='\t')
    if 'position' not in df.columns:
        raise ValueError("数据文件中缺少 'position' 列")
    df = df[(df['position'] >= -WINDOW_SIZE) & (df['position'] <= WINDOW_SIZE)].copy()
    return df

def apply_smoothing(df, sigma):
    """对富集度列应用高斯平滑"""
    if sigma > 0:
        df['template_enrich_smooth'] = gaussian_filter1d(df['template_enrich'], sigma=sigma)
        df['non_template_enrich_smooth'] = gaussian_filter1d(df['non_template_enrich'], sigma=sigma)
    else:
        df['template_enrich_smooth'] = df['template_enrich']
        df['non_template_enrich_smooth'] = df['non_template_enrich']
    return df

def plot_human_tss(df):
    """绘制人类 TSS 富集曲线"""
    fig, ax = plt.subplots(figsize=(4.5, 3.5))

    # 绘制两条曲线
    ax.plot(df['position'], df['template_enrich_smooth'],
            label='Template strand', color=COLOR_TEMPLATE, linewidth=0.8)
    ax.plot(df['position'], df['non_template_enrich_smooth'],
            label='Non-template strand', color=COLOR_NONTEMPLATE, linewidth=0.8)

    # 标记 TSS 位置
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.5, alpha=0.7, label='TSS')

    # 坐标轴范围
    ax.set_xlim(-WINDOW_SIZE, WINDOW_SIZE)
    ax.set_xticks(np.arange(-1000, 1001, 200))

    # 计算 y 轴范围，避免无效数据
    y_min = min(df['template_enrich_smooth'].min(), df['non_template_enrich_smooth'].min())
    y_max = max(df['template_enrich_smooth'].max(), df['non_template_enrich_smooth'].max())
    if np.isnan(y_min) or np.isnan(y_max):
        raise ValueError("平滑后数据包含 NaN，请检查原始数据")
    y_range = y_max - y_min
    ax.set_ylim(max(0, y_min - 0.1*y_range), y_max + 0.1*y_range)

    # 标签
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Normalized enrichment")
    ax.set_title("Human i-Motif enrichment around TSS", fontweight='normal')

    # 图例
    ax.legend(loc='upper right', frameon=False)

    # 网格线
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)

    plt.tight_layout()
    return fig

def save_figure(fig, base_name):
    """保存为 PDF（AI可编辑）和 TIFF（印刷）"""
    pdf_path = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
    tiff_path = os.path.join(OUTPUT_DIR, f"{base_name}.tiff")

    # PDF：增加边距，确保内容不被裁剪
    fig.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
    print(f"✓ PDF 已保存：{pdf_path}")

    # TIFF：600dpi，LZW 压缩
    fig.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight', pad_inches=0.1,
                pil_kwargs={'compression': 'tiff_lzw'})
    print(f"✓ TIFF 已保存：{tiff_path}")

def main():
    print("="*60)
    print("绘制人类 TSS 区域 i-Motif 富集曲线（符合 OUP 指南）")
    print("="*60)

    try:
        # 1. 加载数据
        print(f"正在加载 {SPECIES} {REGION} 数据...")
        df_raw = load_human_tss_data()
        print(f"  数据点数量：{len(df_raw)}")
        print(f"  位置范围：{df_raw['position'].min()} ～ {df_raw['position'].max()} bp")

        # 2. 应用平滑
        if SMOOTH_SIGMA > 0:
            print(f"应用高斯平滑，sigma = {SMOOTH_SIGMA}")
            df = apply_smoothing(df_raw, SMOOTH_SIGMA)
        else:
            df = apply_smoothing(df_raw, 0)

        # 检查平滑后数据是否有效
        if df['template_enrich_smooth'].isnull().all() or df['non_template_enrich_smooth'].isnull().all():
            raise ValueError("平滑后数据全部为 NaN")

        # 3. 绘图
        print("生成图形...")
        fig = plot_human_tss(df)

        # 4. 保存
        save_figure(fig, "human_TSS_enrichment_oup")
        plt.close(fig)

        print("\n✅ 完成！")
        print(f"输出目录：{os.path.abspath(OUTPUT_DIR)}")

    except Exception as e:
        print(f"\n❌ 错误：{e}")
        return 1

    return 0

if __name__ == "__main__":
    exit(main())
