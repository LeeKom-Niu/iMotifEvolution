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
NATURE_COLORS = {
    'blue': '
    'red': '
}
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
        'grid.color': '
        'grid.linestyle': ':',
        'grid.linewidth': 0.3,
        'grid.alpha': 0.3,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
set_nature_style()
data = """Species	TotalBedCount	TotalHaploidLength(bp)	TotalDiploidLength(bp)	TotalDensity(per_kb)
chimp	812713	3177739762	6355479524	0.127876
sumatran	810775	3259084148	6518168296	0.124387
human	791011	3117275501	6234551002	0.126875
bornean	804894	3220935163	6441870326	0.124947
gorilla	977152	3545834224	7091668448	0.137789
bonobo	824831	3244508021	6489016042	0.127112"""
df = pd.read_csv(StringIO(data), sep='\t')
species_names = {
    'bonobo': 'Bonobo',
    'chimp': 'Chimpanzee',
    'human': 'Human',
    'gorilla': 'Gorilla',
    'bornean': 'B. orangutan',
    'sumatran': 'S. orangutan'
}
df['Species_Name'] = df['Species'].map(species_names)
species_order = ['Bonobo', 'Chimpanzee', 'Human', 'Gorilla', 'S. orangutan', 'B. orangutan']
df['Species_Name'] = pd.Categorical(df['Species_Name'], categories=species_order, ordered=True)
df = df.sort_values('Species_Name')
df['Density_per_Mb'] = df['TotalDensity(per_kb)'] * 1000
fig, ax1 = plt.subplots(figsize=(10, 6))
N = len(df)
gap = 1
width = 0.6
x_left = np.arange(N)
x_right = np.arange(N) + N + gap
bars_left = ax1.bar(x_left, df['TotalBedCount'],
                    width=width, color=NATURE_COLORS['blue'],
                    edgecolor='black', linewidth=0.5,
                    label='i-Motif Count', zorder=3)
ax2 = ax1.twinx()
bars_right = ax2.bar(x_right, df['Density_per_Mb'],
                     width=width, color=NATURE_COLORS['red'],
                     edgecolor='black', linewidth=0.5,
                     label='Density (per Mb)', zorder=3)
ax2.spines['right'].set_visible(True)
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
all_centers = np.concatenate([x_left + width/2, x_right + width/2])
all_labels = list(df['Species_Name']) + list(df['Species_Name'])
ax1.set_xticks(all_centers)
ax1.set_xticklabels(all_labels, rotation=45, ha='right', fontsize=15)
ax1.set_ylim(0, df['TotalBedCount'].max() * 1.15)
ax2.set_ylim(0, df['Density_per_Mb'].max() * 1.15)
ax1.tick_params(axis='y', labelsize=15, labelcolor=NATURE_COLORS['blue'])
ax2.tick_params(axis='y', labelsize=15, labelcolor=NATURE_COLORS['red'])
ax1.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{int(x):,}'))
ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))
ax1.set_ylabel('i-Motif Count', color=NATURE_COLORS['blue'], fontsize=15)
ax2.set_ylabel('Density (per Mb)', color=NATURE_COLORS['red'], fontsize=15)
ax2.axhline(y=126, linestyle='--', color='gray', linewidth=0.8, alpha=0.7, zorder=0)
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
legend = ax1.legend(lines1 + lines2, labels1 + labels2,
                    loc='upper left', frameon=False, fontsize=15,
                    bbox_to_anchor=(0.02, 1.05))
plt.tight_layout()
output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, "iMotif_distribution_grouped.pdf")
plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
tiff_path = os.path.join(output_dir, "iMotif_distribution_grouped.tiff")
plt.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight',
            pil_kwargs={'compression': 'tiff_lzw'})
png_path = os.path.join(output_dir, "iMotif_distribution_grouped.png")
plt.savefig(png_path, format='png', dpi=150, bbox_inches='tight')

