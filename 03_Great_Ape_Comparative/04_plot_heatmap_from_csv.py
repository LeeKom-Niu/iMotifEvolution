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
def set_oup_style():
    """
    设置符合OUP插图指南的绘图风格：
    - 字体：Arial/Helvetica，TrueType嵌入（PDF可编辑）
    - 字号：基础12pt，其他元素适当放大
    - 线条：粗细适中
    - 颜色映射：viridis（科学通用）
    """
    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    rcParams['font.size'] = 12
    
    rcParams['figure.figsize'] = (5.5, 11)
    rcParams['figure.dpi'] = 300
    
    rcParams['lines.linewidth'] = 1.0
    rcParams['lines.markersize'] = 4
    
    rcParams['axes.linewidth'] = 0.8
    rcParams['axes.labelpad'] = 6
    rcParams['axes.titlepad'] = 15
    rcParams['axes.labelsize'] = 14
    
    rcParams['xtick.major.width'] = 0.8
    rcParams['ytick.major.width'] = 0.8
    rcParams['xtick.minor.width'] = 0.6
    rcParams['ytick.minor.width'] = 0.6
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 11
    
    rcParams['legend.fontsize'] = 11
    rcParams['legend.frameon'] = False
def load_heatmap_data(file_path):
    """从CSV文件加载热图数据，返回DataFrame"""
    try:
        df = pd.read_csv(file_path, index_col=0)
        df = df.astype(float)
        return df
    except Exception as e:
        return None
def analyze_data(data_df):
    """打印数据统计信息"""
    
    missing_counts = data_df.isnull().sum()
    for species, count in missing_counts.items():
        if count > 0:
    
    
    for species in data_df.columns:
        species_data = data_df[species].dropna()
        if len(species_data) > 0:
    
    mean_by_chromosome = data_df.mean(axis=1)
    top_chromosomes = mean_by_chromosome.sort_values(ascending=False).head(5)
    for idx, (chr_name, density) in enumerate(top_chromosomes.items(), 1):
def plot_heatmap_from_data(data_df, output_dir, filename_prefix="iMotif_density_heatmap", annotate=True):
    """
    绘制热图并保存为可编辑PDF
    annotate: 是否在单元格中显示数值
    颜色条水平放置于底部，紧贴热图（pad=0.03）
    图形宽度缩小以突出纵向布局
    """
    set_oup_style()
    os.makedirs(output_dir, exist_ok=True)
    
    data = data_df.copy()
    
    fig, ax = plt.subplots(figsize=(5.5, 11))
    
    cmap = 'viridis'
    
    if annotate:
        annot = True
        annot_kws = {
            'size': 9,
            'color': 'white',
            'fontweight': 'bold'
        }
        fmt = '.3f'
    else:
        annot = False
        annot_kws = None
        fmt = None
    
    heatmap = sns.heatmap(
        data,
        cmap=cmap,
        linewidths=0.5,
        linecolor='white',
        square=False,
        annot=annot,
        fmt=fmt,
        annot_kws=annot_kws,
        cbar_kws={
            'orientation': 'horizontal',
            'location': 'bottom',
            'label': 'i-Motif density (/Mb)',
            'shrink': 0.7,
            'pad': 0.03,
            'ticks': plt.MaxNLocator(5)
        },
        ax=ax
    )
    
    cbar = heatmap.collections[0].colorbar
    cbar.set_label('i-Motif density (/Mb)', fontsize=14)
    cbar.ax.tick_params(labelsize=11)
    
    
    y_labels = []
    for label in data_df.index:
        if label in ['2a', '2b']:
            y_labels.append(f'Chr {label}')
        else:
            y_labels.append(f'Chr {label}')
    ax.set_yticklabels(y_labels, rotation=0, fontsize=11)
    
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
    ax.set_xticklabels(x_labels, rotation=0, ha='center', fontsize=10)
    
    plt.tight_layout()
    
    suffix = "_with_values" if annotate else "_color_only"
    pdf_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.pdf")
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', format='pdf')
    
    png_path = os.path.join(output_dir, f"{filename_prefix}{suffix}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight', format='png')
    
    
    plt.close(fig)
    return heatmap
def main():
    
    data_file = "heatmap_data_fixed.csv"
    output_dir = "./heatmap_output_narrow"
    
    if not os.path.exists(data_file):
        sys.exit(1)
    
    heatmap_df = load_heatmap_data(data_file)
    if heatmap_df is None:
        sys.exit(1)
    
    analyze_data(heatmap_df)
    
    os.makedirs(output_dir, exist_ok=True)
    
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=True)
    
    plot_heatmap_from_data(heatmap_df, output_dir,
                           filename_prefix="iMotif_density_heatmap",
                           annotate=False)
    
if __name__ == "__main__":
    main()
