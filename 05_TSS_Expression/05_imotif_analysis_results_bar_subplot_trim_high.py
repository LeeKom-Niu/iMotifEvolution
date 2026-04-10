#!/usr/bin/env python3
"""
全物种启动子区 i‑Motif 存在与基因表达关联分析（原始TPM柱状图版）
分组规则：
  With i‑Motif:  有 i‑motif 重叠（max_imotif_score 非空）
  Without i‑Motif: 无 i‑motif 重叠
输出：
  - 每个阈值一个总目录，内含各物种子目录及对应的单物种柱状图、总图
"""

import os
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
from matplotlib import rcParams

# ==================== 配置 ====================
BASE_DIR = os.getcwd()
GENE_BED_DIR = os.path.join(BASE_DIR, "gene_bed_files")
IMOTIF_BED_DIR = os.path.join(BASE_DIR, "imotif_bed")
TPM_DIR = os.path.join(BASE_DIR, "TPM_matrices")
BASE_OUTPUT_DIR = os.path.join(BASE_DIR, "imotif_analysis_results_barplot_raw_trim_high")

BEDTOOLS = "/datapool/home/2023200496/envs/bedtools/bin/bedtools"

SPECIES = ["sumatran", "gorilla", "human", "chimp", "bonobo"]
# 绘图顺序（可按需修改）
DISPLAY_ORDER = ["human", "bonobo", "chimp", "sumatran", "gorilla"]

# 物种 -> i‑motif BED 文件映射
SPECIES_TO_IMOTIF_BED = {
    "sumatran": "Pongo_abelii_all.bed",
    "gorilla": "Gorilla_gorilla_all.bed",
    "human": "Homo_sapiens_all.bed",
    "chimp": "Pan_troglodytes_all.bed",
    "bonobo": "Pan_paniscus_all.bed",
}

# 误差棒类型：'sem' 或 'std'
ERRORBAR_TYPE = 'sem'   # 标准误

# ========== 新增配置：要尝试的剔除高表达阈值列表 ==========
TRIM_PERCENTILES = [95]   # 剔除大于该分位数的值
# =====================================================

def set_oup_style():
    plt.style.use('default')
    rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 11,
        'axes.linewidth': 0.5,
        'grid.linewidth': 0.3,
        'lines.linewidth': 1.0,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def ensure_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def print_file_head(file_path, n=5):
    if not os.path.exists(file_path):
        print(f"    文件不存在: {file_path}")
        return
    with open(file_path) as f:
        lines = [next(f) for _ in range(n)]
    for i, line in enumerate(lines):
        print(f"      Line {i+1}: {line.strip()}")

def create_promoter_bed(gene_bed_path, output_bed):
    df = pd.read_csv(gene_bed_path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
    promoters = []
    for _, row in df.iterrows():
        chrom = row['chr']
        gene_id = row['gene_id']
        strand = row['strand']
        if strand == '+':
            prom_start = max(0, row['start'] - 1000)
            prom_end = row['start']
        else:
            prom_start = row['end']
            prom_end = row['end'] + 1000
        promoters.append([chrom, prom_start, prom_end, gene_id, 0, strand])
    promoter_df = pd.DataFrame(promoters)
    promoter_df.to_csv(output_bed, sep='\t', header=False, index=False)
    print(f"  → 生成启动子区域: {len(promoter_df)} 个区域 -> {output_bed}")

def run_bedtools_intersect(promoter_bed, imotif_bed, output_intersect):
    cmd = [BEDTOOLS, "intersect", "-a", promoter_bed, "-b", imotif_bed,
           "-wa", "-wb", "-s"]
    try:
        with open(output_intersect, 'w') as fout:
            subprocess.run(cmd, stdout=fout, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"  ✗ bedtools 命令失败: {e.stderr.decode()}")
        return False
    if os.path.getsize(output_intersect) == 0:
        print(f"  ⚠ 重叠结果为空，没有 i‑motif 落在启动子区")
        return False
    else:
        print(f"  → 重叠结果非空: {output_intersect}")
        return True

def get_gene_max_imotif_score(intersect_file):
    df_intersect = pd.read_csv(intersect_file, sep='\t', header=None)
    if df_intersect.shape[1] != 12:
        print(f"   警告: 重叠文件列数 {df_intersect.shape[1]}，预期 12，尝试继续...")
    gene_col = 3
    score_col = 10
    if df_intersect.shape[1] <= score_col:
        print(f"   错误: 重叠文件列数不足，无法提取评分")
        return pd.DataFrame(columns=['gene_id', 'max_imotif_score'])
    df_max = df_intersect.groupby(gene_col)[score_col].max().reset_index()
    df_max.columns = ['gene_id', 'max_imotif_score']
    return df_max

def load_tpm_matrix(species):
    tpm_file = os.path.join(TPM_DIR, f"{species}_TPM_matrix.csv")
    if not os.path.exists(tpm_file):
        raise FileNotFoundError(f"TPM 矩阵不存在: {tpm_file}")
    tpm_df = pd.read_csv(tpm_file, index_col=0)
    tpm_df['mean_TPM'] = tpm_df.mean(axis=1)
    tpm_df['log2_TPM'] = np.log2(tpm_df['mean_TPM'] + 1)
    tpm_df = tpm_df.reset_index().rename(columns={'index': 'GeneID'})
    return tpm_df

def analyze_species(species, output_dir, trim_percentile, all_species_data):
    """
    分析单个物种，保存结果到 all_species_data 列表，并生成该物种的柱状图（原始TPM）
    trim_percentile: 剔除高表达极端值的分位数（例如 95），若为 None 则不剔除
    """
    print(f"\n{'='*60}")
    print(f"分析物种: {species} (剔除 >{trim_percentile}% 分位数)")
    print(f"{'='*60}")

    species_out = os.path.join(output_dir, species)
    ensure_dir(species_out)

    # 1. 基因 BED
    gene_bed = os.path.join(GENE_BED_DIR, f"{species}.genes.bed")
    if not os.path.exists(gene_bed):
        print(f"  ✗ 基因 BED 文件不存在: {gene_bed}，跳过")
        return

    # 2. i‑motif BED
    imotif_bed_name = SPECIES_TO_IMOTIF_BED.get(species)
    if not imotif_bed_name:
        print(f"  ✗ 未找到 {species} 的 i‑motif BED 映射，跳过")
        return
    imotif_bed = os.path.join(IMOTIF_BED_DIR, imotif_bed_name)
    if not os.path.exists(imotif_bed):
        print(f"  ✗ i‑motif BED 文件不存在: {imotif_bed}，跳过")
        return

    print("\n  基因 BED 前 3 行:")
    print_file_head(gene_bed, 3)
    print("\n  i‑motif BED 前 3 行:")
    print_file_head(imotif_bed, 3)

    # 3. 生成启动子 BED
    promoter_bed = os.path.join(species_out, f"{species}_promoter.bed")
    create_promoter_bed(gene_bed, promoter_bed)

    # 4. bedtools 重叠
    intersect_file = os.path.join(species_out, f"{species}_promoter_imotif_intersect.bed")
    has_overlap = run_bedtools_intersect(promoter_bed, imotif_bed, intersect_file)
    if not has_overlap:
        print(f"  → 由于没有重叠，跳过 {species} 的后续分析")
        return

    # 5. 提取每个基因的最大 i‑motif 评分
    max_score_df = get_gene_max_imotif_score(intersect_file)
    print(f"  → 找到 {len(max_score_df)} 个基因有 i‑motif 重叠")

    # 6. 加载 TPM 矩阵
    tpm_df = load_tpm_matrix(species)
    print(f"  → TPM 矩阵基因数: {len(tpm_df)}")

    # 7. 合并
    merged = tpm_df.merge(max_score_df, left_on='GeneID', right_on='gene_id', how='left')
    merged = merged.drop(columns=['gene_id'])

    # 分组：有 i‑motif vs 无 i‑motif
    merged['group'] = merged['max_imotif_score'].apply(
        lambda x: 'With i-Motif' if pd.notna(x) else 'Without i-Motif'
    )

    # 过滤掉表达量为 0 的基因
    merged_filtered = merged[merged['mean_TPM'] > 0].copy()
    print(f"  → 过滤后（表达量>0）基因数: {len(merged_filtered)}")

    # 剔除高表达极端值（如果 trim_percentile 不是 None）
    if trim_percentile is not None:
        all_tpm = merged_filtered['mean_TPM']
        high_thresh = np.percentile(all_tpm, trim_percentile)
        n_before = len(merged_filtered)
        merged_filtered = merged_filtered[merged_filtered['mean_TPM'] <= high_thresh].copy()
        n_after = len(merged_filtered)
        print(f"  → 剔除 >{trim_percentile}% 分位数 (阈值 {high_thresh:.2f}) 后，基因数从 {n_before} 降至 {n_after}")

    cat_counts = merged_filtered['group'].value_counts()
    print("\n  分组统计:")
    for cat in ['With i-Motif', 'Without i-Motif']:
        cnt = cat_counts.get(cat, 0)
        print(f"    {cat}: {cnt} 个基因")

    # 提取各组原始 TPM
    with_data = merged_filtered[merged_filtered['group'] == 'With i-Motif']
    without_data = merged_filtered[merged_filtered['group'] == 'Without i-Motif']

    if len(with_data) == 0 or len(without_data) == 0:
        print("  → 某一分组为空，无法进行统计检验和绘图，跳过")
        # 仍保存数据文件
        data_out = os.path.join(species_out, f"{species}_imotif_analysis_data.csv")
        merged_filtered[['GeneID', 'group', 'max_imotif_score', 'mean_TPM', 'log2_TPM']].to_csv(data_out, index=False)
        print(f"  → 数据已保存至: {data_out}")
        return

    # 统计检验（Mann-Whitney U）和原始 TPM 倍数计算
    tpm_with = with_data['mean_TPM']
    tpm_without = without_data['mean_TPM']
    p_tpm = stats.mannwhitneyu(tpm_with, tpm_without).pvalue
    median_with = tpm_with.median()
    median_without = tpm_without.median()
    fold_change = median_with / median_without if median_without > 0 else np.inf

    print("\n Mann-Whitney U 检验 p 值 (原始 TPM):")
    print(f"  Raw TPM: With vs Without: p = {p_tpm:.2e}")
    print(f"\n 原始 TPM 中位数:")
    print(f"   With i-Motif:  {median_with:.2f}")
    print(f"   Without i-Motif: {median_without:.2f}")
    print(f"  倍数 (With/Without): {fold_change:.2f}")

    # 计算原始 TPM 的均值和标准误（或标准差）
    mean_with = tpm_with.mean()
    mean_without = tpm_without.mean()
    if ERRORBAR_TYPE == 'sem':
        error_with = tpm_with.sem()
        error_without = tpm_without.sem()
    else:  # std
        error_with = tpm_with.std()
        error_without = tpm_without.std()

    # 保存数据供总图使用（包括均值和误差）
    all_species_data.append({
        'species': species,
        'mean_with': mean_with,
        'mean_without': mean_without,
        'error_with': error_with,
        'error_without': error_without,
        'n_with': len(tpm_with),
        'n_without': len(tpm_without),
        'p_tpm': p_tpm,
        'fold_change': fold_change
    })

    # ========== 绘制单物种柱状图（原始 TPM） ==========
    colors = {'With i-Motif': '#E64B35', 'Without i-Motif': '#4DBBD5'}
    categories = ['With i-Motif', 'Without i-Motif']
    means = [mean_with, mean_without]
    errors = [error_with, error_without]
    ns = [len(tpm_with), len(tpm_without)]

    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)

    # 绘制柱状图
    bars = ax.bar(categories, means, yerr=errors, capsize=5,
                  color=[colors[cat] for cat in categories],
                  edgecolor='black', linewidth=1, alpha=0.8,
                  error_kw={'linewidth': 1.5, 'ecolor': 'black'})

    # 添加显著性标记和倍数标记
    y_max = max(means + errors) * 1.1  # 为标注留出空间
    if p_tpm < 0.001: sig = '***'
    elif p_tpm < 0.01: sig = '**'
    elif p_tpm < 0.05: sig = '*'
    else: sig = 'ns'

    if p_tpm < 0.05:
        # 画一条连接两个柱子的横线
        x1, x2 = 0, 1
        y = y_max
        ax.plot([x1, x1, x2, x2], [y, y+0.05*y_max, y+0.05*y_max, y], 'k-', linewidth=0.8)
        ax.text((x1+x2)/2, y+0.08*y_max, sig, ha='center', va='bottom', fontsize=12)
        # 显示倍数
        ax.text((x1+x2)/2, y_max*0.9, f'Fold = {fold_change:.2f}', ha='center', va='bottom', fontsize=10, style='italic')
    else:
        # 仅显示倍数和 p 值
        ax.text(0.5, y_max*0.9, f'Fold = {fold_change:.2f}\np = {p_tpm:.2e}',
                ha='center', va='bottom', fontsize=9, style='italic')

    # 设置标签和标题
    ax.set_ylabel('Mean TPM (raw)', fontsize=14)
    ax.set_title(f'{species.capitalize()}', loc='left', fontweight='bold', fontsize=16)
    ax.set_xticks(range(len(categories)))
    ax.set_xticklabels([f'{cat}\n(n={n})' for cat, n in zip(categories, ns)], fontsize=12)
    ax.set_ylim(0, y_max * 1.05)  # 留一点顶部空间

    # 可选：添加网格线
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    # 保存图形
    pdf_path = os.path.join(species_out, f"{species}_barplot_raw.pdf")
    plt.savefig(pdf_path, format='pdf')
    plt.savefig(os.path.join(species_out, f"{species}_barplot_raw.tiff"), format='tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.savefig(os.path.join(species_out, f"{species}_barplot_raw.png"), format='png', dpi=150)
    plt.close()
    print(f"  → 柱状图保存至: {species_out}/")

    # 保存数据文件
    data_out = os.path.join(species_out, f"{species}_imotif_analysis_data.csv")
    merged_filtered[['GeneID', 'group', 'max_imotif_score', 'mean_TPM', 'log2_TPM']].to_csv(data_out, index=False)
    print(f"  → 数据已保存至: {data_out}")

    print(f"\n✓ {species} 分析完成！")

def plot_combined_barplot(all_species_data, output_dir, trim_percentile):
    """
    绘制所有物种的分组柱状图总图（单一图，x轴为物种，横向长条形状）
    按 DISPLAY_ORDER 排序，优化标注位置
    trim_percentile: 用于标题说明
    """
    # 按 DISPLAY_ORDER 重新排序数据
    data_dict = {data['species']: data for data in all_species_data}
    sorted_data = [data_dict[sp] for sp in DISPLAY_ORDER if sp in data_dict]
    n = len(sorted_data)
    if n == 0:
        print("没有有效数据，无法绘制总图。")
        return

    # 横向长条形状：宽度与物种数量成正比，但固定高度较小，形成横向长条
    fig_width = max(8, n * 1.2)   # 每个物种约1.2英寸宽
    fig_height = 5                # 固定高度，使图形扁长
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)

    # 颜色配置
    colors = {'With i-Motif': '#E64B35', 'Without i-Motif': '#4DBBD5'}

    # 准备数据
    species_names = [data['species'].capitalize() for data in sorted_data]
    means_with = [data['mean_with'] for data in sorted_data]
    means_without = [data['mean_without'] for data in sorted_data]
    errors_with = [data['error_with'] for data in sorted_data]
    errors_without = [data['error_without'] for data in sorted_data]
    p_values = [data['p_tpm'] for data in sorted_data]
    folds = [data['fold_change'] for data in sorted_data]

    # 设置分组柱状图位置
    x = np.arange(n)  # 物种位置
    width = 0.35      # 柱子宽度
    # 绘制柱子（With i-Motif 在左侧，Without i-Motif 在右侧）
    bars_with = ax.bar(x - width/2, means_with, width, yerr=errors_with, capsize=5,
                       color=colors['With i-Motif'], edgecolor='black', linewidth=1, alpha=0.8,
                       label='With i-Motif', error_kw={'linewidth': 1.5, 'ecolor': 'black'})
    bars_without = ax.bar(x + width/2, means_without, width, yerr=errors_without, capsize=5,
                          color=colors['Without i-Motif'], edgecolor='black', linewidth=1, alpha=0.8,
                          label='Without i-Motif', error_kw={'linewidth': 1.5, 'ecolor': 'black'})

    # 计算每个柱组的最大高度（用于放置标记）
    max_heights = [max(means_with[i] + errors_with[i], means_without[i] + errors_without[i]) for i in range(n)]
    # 为标记留出空间：取最大高度的1.3倍，但至少比最大值高20%
    y_max_plot = max(max_heights) * 1.3
    if y_max_plot < 1e-5:
        y_max_plot = 1.0

    # 添加基线（y=0线）
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, zorder=0)

    # 添加显著性标记和倍数标记
    for i in range(n):
        # 当前柱组最高点（考虑误差）
        top_with = means_with[i] + errors_with[i]
        top_without = means_without[i] + errors_without[i]
        max_top = max(top_with, top_without)
        # 标记起始y坐标：柱组最高点 + 间距（动态，取最大高度的5%）
        y_start = max_top + 0.05 * y_max_plot
        # 如果y_start超出y轴上限，适当调整y轴上限
        if y_start > y_max_plot:
            y_max_plot = y_start + 0.05 * y_max_plot

        # 画横线连接两个柱子（从 With 柱右侧到 Without 柱左侧）
        x_left = x[i] - width/2 + width
        x_right = x[i] + width/2
        ax.plot([x_left, x_right], [y_start, y_start], 'k-', linewidth=0.8)

        # 根据p值确定星号
        p = p_values[i]
        if p < 0.001:
            sig = '***'
        elif p < 0.01:
            sig = '**'
        elif p < 0.05:
            sig = '*'
        else:
            sig = 'ns'

        # 标注星号（位于横线上方）
        ax.text(x[i], y_start + 0.02 * y_max_plot, sig, ha='center', va='bottom', fontsize=10)

        # 标注倍数（位于横线下方，柱组内部）
        fold = folds[i]
        # 选择两个柱子中较高的一个，在柱子内部上方放置倍数文本
        higher_bar_top = max(means_with[i], means_without[i])
        # 如果柱子高度过小，则放在横线下方
        if higher_bar_top > 0:
            y_fold = higher_bar_top + 0.02 * y_max_plot
            # 确保不超出横线
            if y_fold >= y_start:
                y_fold = y_start - 0.03 * y_max_plot
        else:
            y_fold = y_start - 0.03 * y_max_plot
        ax.text(x[i], y_fold, f'Fold={fold:.2f}', ha='center', va='bottom', fontsize=8, style='italic')

    # 设置坐标轴
    ax.set_xticks(x)
    ax.set_xticklabels(species_names, rotation=0, ha='center', fontsize=11)
    ax.set_ylabel('Mean TPM', fontsize=14)
    if trim_percentile is not None:
        title = f'Expression comparison (genes with/without i-Motif in promoters)\n(high expression trimmed >{trim_percentile}% quantile)'
    else:
        title = 'Expression comparison (genes with/without i-Motif in promoters)\n(no trimming)'
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(loc='upper left', frameon=True, fontsize=10)
    ax.set_ylim(0, y_max_plot)
    ax.yaxis.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    # 保存图形
    outfile = os.path.join(output_dir, "combined_barplot_raw.pdf")
    plt.savefig(outfile, format='pdf')
    plt.savefig(os.path.join(output_dir, "combined_barplot_raw.png"), format='png', dpi=150)
    plt.close()
    print(f"\n总图（横向长条）已保存: {outfile}")

def main():
    set_oup_style()
    ensure_dir(BASE_OUTPUT_DIR)

    # 对每个剔除阈值进行循环
    for trim_percentile in TRIM_PERCENTILES:
        print(f"\n{'='*60}")
        print(f"开始处理剔除阈值: >{trim_percentile}% 分位数")
        print(f"{'='*60}")

        # 创建该阈值对应的输出目录
        output_dir = os.path.join(BASE_OUTPUT_DIR, f"trim_high_{trim_percentile}")
        ensure_dir(output_dir)

        all_species_data = []  # 收集该阈值下的所有物种数据

        for sp in SPECIES:
            try:
                analyze_species(sp, output_dir, trim_percentile, all_species_data)
            except Exception as e:
                print(f"  ✗ 分析 {sp} 时出错: {e}")
                import traceback
                traceback.print_exc()
                continue

        # 绘制该阈值下的总图
        if all_species_data:
            plot_combined_barplot(all_species_data, output_dir, trim_percentile)
        else:
            print("没有成功分析任何物种，无法生成总图。")

        print(f"\n阈值 {trim_percentile}% 处理完成，结果位于: {output_dir}")

    print("\n" + "="*60)
    print("所有阈值处理完成！总结果目录: " + BASE_OUTPUT_DIR)
    print("="*60)

if __name__ == "__main__":
    main()

