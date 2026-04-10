#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
热图绘制脚本（OUP适配版，纵向窄版，颜色条紧贴底部）
从CSV文件读取i-Motif密度数据，生成可编辑的PDF矢量图。
颜色条水平置于底部，紧贴热图；图形宽度缩小以突出纵向。
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rcParams
import os
import sys

# ==================== OUP风格设置 ====================
def set_oup_style():
    """
    设置符合OUP插图指南的绘图风格：
    - 字体：Arial/Helvetica，TrueType嵌入（PDF可编辑）
    - 字号：基础12pt，其他元素适当放大
    - 线条：粗细适中
    - 颜色映射：viridis（科学通用）
    """
    # 字体嵌入设置（确保PDF文本可编辑）
    rcParams['pdf.fonttype'] = 42        # TrueType字体
    rcParams['ps.fonttype'] = 42         # PostScript TrueType
    
    # 字体族
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['font.size'] = 12           # 基础字号
    
    # 图形尺寸：宽度减小至5.5英寸，突出纵向布局；高度保持11英寸以容纳所有染色体
    rcParams['figure.figsize'] = (5.5, 11)
    rcParams['figure.dpi'] = 300
    
    # 线条和标记
    rcParams['lines.linewidth'] = 1.0
    rcParams['lines.markersize'] = 4
    
    # 坐标轴
    rcParams['axes.linewidth'] = 0.8
    rcParams['axes.labelpad'] = 6
    rcParams['axes.titlepad'] = 15
    rcParams['axes.labelsize'] = 14      # 轴标签字号
    
    # 刻度
    rcParams['xtick.major.width'] = 0.8
    rcParams['ytick.major.width'] = 0.8
    rcParams['xtick.minor.width'] = 0.6
    rcParams['ytick.minor.width'] = 0.6
    rcParams['xtick.labelsize'] = 10     # X轴刻度标签字号（略减，适应窄宽度）
    rcParams['ytick.labelsize'] = 11     # Y轴刻度标签字号
    
    # 图例
    rcParams['legend.fontsize'] = 11
    rcParams['legend.frameon'] = False

# ==================== 数据加载 ====================
def load_heatmap_data(file_path):
    """从CSV文件加载热图数据，返回DataFrame"""
    try:
        df = pd.read_csv(file_path, index_col=0)
        print(f"成功加载数据文件: {file_path}")
        print(f"数据形状: {df.shape}")
        print(f"行索引: {list(df.index)}")
        print(f"列名: {list(df.columns)}")
        print("\n数据预览:")
        print(df.head())
        df = df.astype(float)
        return df
    except Exception as e:
        print(f"加载数据文件时出错: {e}")
        return None

# ==================== 数据统计分析 ====================
def analyze_data(data_df):
    """打印数据统计信息"""
    print("\n" + "=" * 50)
    print("数据统计分析:")
    print("=" * 50)
    print(f"数据形状: {data_df.shape}")
    print(f"行数（染色体）: {len(data_df.index)}")
    print(f"列数（物种）: {len(data_df.columns)}")
    
    # 缺失值
    missing_counts = data_df.isnull().sum()
    print("\n缺失值统计:")
    for species, count in missing_counts.items():
        if count > 0:
            print(f"{species}: {count} 个缺失值")
    
    # 全局统计
    print("\n密度值统计 (/Mb):")
    print(f"全局最小值: {data_df.min().min():.6f}")
    print(f"全局最大值: {data_df.max().max():.6f}")
    print(f"全局平均值: {data_df.mean().mean():.6f}")
    print(f"全局中位数: {data_df.stack().median():.6f}")
    print(f"全局标准差: {data_df.stack().std():.6f}")
    
    # 物种统计
    print("\n各物种统计 (/Mb):")
    for species in data_df.columns:
        species_data = data_df[species].dropna()
        if len(species_data) > 0:
            print(f"{species}: 平均值={species_data.mean():.6f}, 最大值={species_data.max():.6f}")
    
    # 染色体排名
    print("\n染色体平均密度排名 (前5):")
    mean_by_chromosome = data_df.mean(axis=1)
    top_chromosomes = mean_by_chromosome.sort_values(ascending=False).head(5)
    for idx, (chr_name, density) in enumerate(top_chromosomes.items(), 1):
        print(f"  {idx}. Chr {chr_name}: {density:.6f}")

# ==================== 热图绘制 ====================
def plot_heatmap_from_data(data_df, output_dir, filename_prefix="iMotif_density_heatmap", annotate=True):
    """
    绘制热图并保存为可编辑PDF
    annotate: 是否在单元格中显示数值
    颜色条水平放置于底部，紧贴热图（pad=0.03）
    图形宽度缩小以突出纵向布局
    """
    set_oup_style()
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备数据
    data = data_df.copy()
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(5.5, 11))  # 窄版纵向尺寸
    
    # 颜色映射
    cmap = 'viridis'
    
    # 注释格式
    if annotate:
        annot = True
        annot_kws = {
            'size': 9,                 # 单元格内数字字号（稍减小以适应窄格）
            'color': 'white',
            'fontweight': 'bold'
        }
        fmt = '.3f'
    else:
        annot = False
        annot_kws = None
        fmt = None
    
    # 绘制热图
    heatmap = sns.heatmap(
        data,
        cmap=cmap,
        linewidths=0.5,
        linecolor='white',
        square=False,                   # 不强制正方形，纵向自然拉伸
        annot=annot,
        fmt=fmt,
        annot_kws=annot_kws,
        cbar_kws={
            'orientation': 'horizontal',   # 水平颜色条
            'location': 'bottom',          # 放在底部
            'label': 'i-Motif density (/Mb)',
            'shrink': 0.7,                 # 颜色条长度比例（略窄于热图宽度）
            'pad': 0.03,                    # ★ 关键修改：间距大幅减小，紧贴热图
            'ticks': plt.MaxNLocator(5)
        },
        ax=ax
    )
    
    # 颜色条标签字体
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('i-Motif density (/Mb)', fontsize=14)
    cbar.ax.tick_params(labelsize=11)
    
    # 标题已注释
    # ax.set_title(...)
    
    # Y轴标签（染色体）
    y_labels = []
    for label in data_df.index:
        if label in ['2a', '2b']:
            y_labels.append(f'Chr {label}')
        else:
            y_labels.append(f'Chr {label}')
    ax.set_yticklabels(y_labels, rotation=0, fontsize=11)
    
    # X轴标签（物种）映射
    species_names = {
        'bonobo': 'Bonobo',
        'chimp': 'Chimpanzee',
        'human': 'Human',
        'gorilla': 'Gorilla',
        'sumatran': 'S. orangutan',
        'bornean': 'B. orangutan'
    }
    x_labels = []
    for col in data_df.columns:
        if col in species_names:
            x_labels.append(species_names[col])
        else:
            x_labels.append(col)
    # 保持水平，字体大小10（已在rcParams中设置，此处显式指定以覆盖可能的变化）
    ax.set_xticklabels(x_labels, rotation=0, ha='center', fontsize=10)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存为可编辑PDF
    suffix = "_with_values" if annotate else "_color_only"
    pdf_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', format='pdf')
    
    # 保存PNG预览
    png_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight', format='png')
    
    print(f"可编辑PDF已保存至: {pdf_path}")
    print(f"PNG预览已保存至: {png_path}")
    
    plt.close(fig)
    return heatmap

# ==================== 主程序 ====================
def main():
    print("=" * 60)
    print("i-Motif密度热图生成脚本 (纵向窄版，颜色条紧贴底部)")
    print("=" * 60)
    
    # 文件路径
    data_file = "heatmap_data_fixed.csv"
    output_dir = "./heatmap_output_narrow"
    
    if not os.path.exists(data_file):
        print(f"错误: 数据文件 '{data_file}' 不存在！")
        sys.exit(1)
    
    heatmap_df = load_heatmap_data(data_file)
    if heatmap_df is None:
        sys.exit(1)
    
    analyze_data(heatmap_df)
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n正在绘制带数值的热图...")
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=True)
    
    print("\n正在绘制不带数值的热图...")
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=False)
    
    print("\n" + "=" * 60)
    print("热图生成完成！")
    print(f"输出目录: {output_dir}")
    print("=" * 60)

if __name__ == "__main__":
    main()
