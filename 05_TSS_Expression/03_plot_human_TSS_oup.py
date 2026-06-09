"""
plot_human_TSS_oup.py - 绘制人类 TSS 区域 i-Motif 富集曲线（符合 OUP 插图指南）
修复：增加边距、字体嵌入设置，确保 PDF 在 AI 中完整显示。
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'lines.linewidth': 0.8,
    'axes.linewidth': 0.5,
    'grid.linewidth': 0.3,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'grid.linestyle': ':',
    'grid.alpha': 0.3,
    'legend.frameon': False,
    'figure.dpi': 300,
    'savefig.dpi': 600,
})
COLOR_TEMPLATE = '
COLOR_NONTEMPLATE = '
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
    ax.plot(df['position'], df['template_enrich_smooth'],
            label='Template strand', color=COLOR_TEMPLATE, linewidth=0.8)
    ax.plot(df['position'], df['non_template_enrich_smooth'],
            label='Non-template strand', color=COLOR_NONTEMPLATE, linewidth=0.8)
    ax.axvline(x=0, color='black', linestyle='--', linewidth=0.5, alpha=0.7, label='TSS')
    ax.set_xlim(-WINDOW_SIZE, WINDOW_SIZE)
    ax.set_xticks(np.arange(-1000, 1001, 200))
    y_min = min(df['template_enrich_smooth'].min(), df['non_template_enrich_smooth'].min())
    y_max = max(df['template_enrich_smooth'].max(), df['non_template_enrich_smooth'].max())
    if np.isnan(y_min) or np.isnan(y_max):
        raise ValueError("平滑后数据包含 NaN，请检查原始数据")
    y_range = y_max - y_min
    ax.set_ylim(max(0, y_min - 0.1*y_range), y_max + 0.1*y_range)
    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Normalized enrichment")
    ax.set_title("Human i-Motif enrichment around TSS", fontweight='normal')
    ax.legend(loc='upper right', frameon=False)
    ax.grid(True, linestyle=':', linewidth=0.3, alpha=0.3)
    plt.tight_layout()
    return fig
def save_figure(fig, base_name):
    """保存为 PDF（AI可编辑）和 TIFF（印刷）"""
    pdf_path = os.path.join(OUTPUT_DIR, f"{base_name}.pdf")
    tiff_path = os.path.join(OUTPUT_DIR, f"{base_name}.tiff")
    fig.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
    fig.savefig(tiff_path, format='tiff', dpi=600, bbox_inches='tight', pad_inches=0.1,
                pil_kwargs={'compression': 'tiff_lzw'})
def main():
    try:
        df_raw = load_human_tss_data()
        if SMOOTH_SIGMA > 0:
            df = apply_smoothing(df_raw, SMOOTH_SIGMA)
        else:
            df = apply_smoothing(df_raw, 0)
        if df['template_enrich_smooth'].isnull().all() or df['non_template_enrich_smooth'].isnull().all():
            raise ValueError("平滑后数据全部为 NaN")
        fig = plot_human_tss(df)
        save_figure(fig, "human_TSS_enrichment_oup")
        plt.close(fig)
    except Exception as e:
        return 1
    return 0
if __name__ == "__main__":
    exit(main())
