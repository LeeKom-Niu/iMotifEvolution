#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_two_bar_grouped.py - 双轴分组柱状图（数量左侧一组，密度右侧一组）
功能：将 i-Motif 数量和密度分为左右两组柱状图，共用 x 轴，密度轴在右侧。
      x 轴显示所有物种名（左右两组均显示），右侧组无柱顶标记。
输出：AI可编辑PDF + 高分辨率TIFF
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.ticker import FuncFormatter, StrMethodFormatter
from io import StringIO
import os

# ================== Nature 配色方案 ==================
NATURE_COLORS = {
    'blue': '#0072B2',
    'red': '#D55E00',
}

# ================== 绘图风格设置 ==================
def set_nature_style():
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'axes.linewidth': 0.5,
        'axes.edgecolor': 'black',
        'axes.labelpad': 6,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.3,
        'ytick.minor.width': 0.3,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
        'xtick.minor.size': 2,
        'ytick.minor.size': 2,
        'legend.fontsize': 10,
        'legend.frameon': False,
        'legend.loc': 'upper left',
        'lines.linewidth': 0.8,
        'figure.dpi': 300,
        'savefig.dpi': 600,
        'figure.figsize': (10, 6),
        'grid.color': '#cccccc',
        'grid.linestyle': ':',
        'grid.linewidth': 0.3,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

set_nature_style()

# ================== 数据读取与处理 ==================
data = """Species	TotalBedCount	TotalHaploidLength(bp)	TotalDiploidLength(bp)	TotalDensity(per_kb)
chimp	812713	3177739762	6355479524	0.127876
sumatran	810775	3259084148	6518168296	0.124387
human	791011	3117275501	6234551002	0.126875
bornean	804894	3220935163	6441870326	0.124947
gorilla	977152	3545834224	7091668448	0.137789
bonobo	824831	3244508021	6489016042	0.127112"""

df = pd.read_csv(StringIO(data), sep='\t')

# 物种名称映射
species_names = {
    'bonobo': 'Bonobo',
    'chimp': 'Chimpanzee',
    'human': 'Human',
    'gorilla': 'Gorilla',
    'bornean': 'B. orangutan',
    'sumatran': 'S. orangutan'
}
df['Species_Name'] = df['Species'].map(species_names)

# 设置顺序
species_order = ['Bonobo', 'Chimpanzee', 'Human', 'Gorilla', 'S. orangutan', 'B. orangutan']
df['Species_Name'] = pd.Categorical(df['Species_Name'], categories=species_order, ordered=True)
df = df.sort_values('Species_Name')

# 计算每 Mb 密度
df['Density_per_Mb'] = df['TotalDensity(per_kb)'] * 1000

# ================== 绘图 ==================
fig, ax1 = plt.subplots(figsize=(10, 6))

N = len(df)                     # 物种数
gap = 1                         # 左右组之间的空隙（空一个柱子宽度）
width = 0.6                     # 每个柱子的宽度

# 左侧组：数量柱子，x 坐标 0..N-1
x_left = np.arange(N)

# 右侧组：密度柱子，x 坐标 N+gap .. N+gap+N-1
x_right = np.arange(N) + N + gap

# 绘制数量柱（左侧组）
bars_left = ax1.bar(x_left, df['TotalBedCount'],
                    width=width, color=NATURE_COLORS['blue'],
                    edgecolor='black', linewidth=0.5,
                    label='i-Motif Count', zorder=3)

# 创建右侧 y 轴（密度轴）
ax2 = ax1.twinx()

# 绘制密度柱（右侧组）
bars_right = ax2.bar(x_right, df['Density_per_Mb'],
                     width=width, color=NATURE_COLORS['red'],
                     edgecolor='black', linewidth=0.5,
                     label='Density (per Mb)', zorder=3)

# ========== 坐标轴修饰 ==========
# 隐藏右侧轴脊（因为我们有右侧 y 轴）
ax2.spines['right'].set_visible(True)   # 保留右侧轴
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)

# 设置 x 轴刻度位置为左右两组中心，并显示两组标签（左右均显示物种名）
all_centers = np.concatenate([x_left + width/2, x_right + width/2])
all_labels = list(df['Species_Name']) + list(df['Species_Name'])   # 右侧也显示物种名称
ax1.set_xticks(all_centers)
ax1.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=15)

# 设置 y 轴
ax1.set_ylim(0, df['TotalBedCount'].max() * 1.15)
ax2.set_ylim(0, df['Density_per_Mb'].max() * 1.15)

ax1.tick_params(axis='y', labelsize=15, labelcolor=NATURE_COLORS['blue'])
ax2.tick_params(axis='y', labelsize=15, labelcolor=NATURE_COLORS['red'])

ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{int(x):,}'))
ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))

ax1.set_ylabel('i-Motif Count', color=NATURE_COLORS['blue'], fontsize=15)
ax2.set_ylabel('Density (per Mb)', color=NATURE_COLORS['red'], fontsize=15)

# ========== 基线（右侧密度轴 y=125 虚线） ==========
ax2.axhline(y=126, linestyle='--', color='gray', linewidth=0.8, alpha=0.7, zorder=0)

# ========== 图例（上移，避免与柱子重叠） ==========
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
legend = ax1.legend(lines1 + lines2, labels1 + labels2,
                    loc='upper left', frameon=False, fontsize=15,
                    bbox_to_anchor=(0.02, 1.05))   # 调整为合适位置

# ========== 输出与保存 ==========
# 自动调整布局，无需额外留白（已删除底部文字）
plt.tight_layout()

output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)

pdf_path = os.path.join(output_dir, "iMotif_distribution_grouped.pdf")
plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
print(f"✓ PDF 已保存: {pdf_path}")

tiff_path = os.path.join(output_dir, "iMotif_distribution_grouped.tiff")
plt.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight',
            pil_kwargs={'compression': 'tiff_lzw'})
print(f"✓ TIFF 已保存: {tiff_path}")

png_path = os.path.join(output_dir, "iMotif_distribution_grouped.png")
plt.savefig(png_path, format='png', dpi=150, bbox_inches='tight')
print(f"✓ PNG 预览: {png_path}")

# ========== 统计信息 ==========
print("\n" + "="*50)
print("统计数据（左右分组柱状图）:")
print("="*50)
print(df[['Species_Name', 'TotalBedCount', 'Density_per_Mb']].to_string(index=False))
print(f"\n总计 i-Motif 数量: {df['TotalBedCount'].sum():,}")
print(f"密度范围: {df['Density_per_Mb'].min():.2f} - {df['Density_per_Mb'].max():.2f} per Mb")
print("="*50)
