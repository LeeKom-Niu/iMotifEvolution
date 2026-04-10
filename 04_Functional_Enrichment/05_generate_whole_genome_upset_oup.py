#!/usr/bin/env python3
"""
绘制全基因组G4共享模式的Upset图
专注于全基因组分析，保持与原脚本完全相同的样式
修改：
  1. 输出可编辑PDF（pdf.fonttype=42）
  2. 同时输出高分辨率PNG（300 dpi）
  3. 新增高分辨率TIFF（600 dpi，LZW压缩）用于出版印刷
  4. 提供图例位置调整的示例注释
"""

import os
import sys
from itertools import combinations
import matplotlib
matplotlib.use('Agg')  # 非交互式环境使用
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import upsetplot as upsplt
import warnings
from collections import defaultdict

# 设置PDF字体为可编辑（Type 42 TrueType）
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42  # 同时设置PostScript

warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=UserWarning)

# 设置项目根目录
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"

# 设置matplotlib参数
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'

# 色盲友好调色板
CBPalette = {
    "Black": "#000000",
    "Orange": "#F4A637",
    "Light blue": "#B6DBFF",
    "Vermilion": "#DB5829",
    "Mid blue": "#7BB0DF",
    "Maroon": "#894B45",
    "Dark blue": "#1964B0",
    "Light purple": "#D2BBD7",
    "Light teal": "#00C992",
    "Purple": "#AE75A2",
    "Teal": "#008A69",
    "Dark purple": "#882D71",
    "Dark teal": "#386350",
    "Grey": "#DEDEDE",
    "Yellow": "#E9DC6D"
}

# 物种名称映射字典
speciesnos = {
    1: ['human', 'Homo_sapiens', 'hs1', 'H. sapiens'],
    2: ['bonobo', 'Pan_paniscus', 'pan', 'P. paniscus'],
    3: ['chimp', 'Pan_troglodytes', 'pan', 'P. troglodytes'],
    4: ['gorilla', 'Gorilla_gorilla', 'gor', 'G. gorilla'],
    5: ['sorang', 'Pongo_abelii', 'pon', 'P. abelii'],
    6: ['borang', 'Pongo_pygmaeus', 'pon', 'P. pygmaeus']
}

def comma_formatter(x, pos):
    '''千位分隔符格式化'''
    return '{:,.0f}'.format(x)

def splitIDsandarrange(df):
    '''分割ID并整理列顺序'''
    dfOut = df[0].str.split('|', expand=True)
    dfOut.columns = ['species', 'chrom', 'start', 'strand', 'length', 'score']
    dfOut['species'] = dfOut['species'].astype(int)
    dfOut['start'] = dfOut['start'].astype(int)
    dfOut['length'] = dfOut['length'].astype(int)
    dfOut['end'] = dfOut['start'] + dfOut['length']
    dfOut.sort_values(by=['species', 'chrom', 'start'], inplace=True)
    dfOut['chrom'] = dfOut['chrom'].apply(lambda x: f'chr{x}')
    dfOut.drop(columns=['length'], inplace=True)
    dfOut.reset_index(drop=True, inplace=True)
    dfOut["dummy"] = '.'
    dfOut = dfOut[['species', 'chrom', 'start', 'end', 'dummy', 'score', 'strand']]
    return dfOut

def stackedBarUpset(upsetDatadf, species, alignedUnique):
    '''为堆叠条形Upset图调整数据'''
    upsetDatadf.reset_index(inplace=True)
    
    for nos, target in enumerate(species):
        # 找到该物种特有的行
        condition = ((upsetDatadf[target] == True) & 
                     (upsetDatadf[[s for s in species if s != target]].eq(False).all(axis=1)))
        condition_indices = upsetDatadf[condition].index
        
        # 分配对齐和未对齐的物种特异性pG4s
        if len(condition_indices) > 0:
            # 前alignedUnique[nos]个为对齐的
            upsetDatadf.loc[condition_indices[:alignedUnique[nos]], "forUnique01"] = "Aligned species-specific IMs"
            # 其余的为未对齐的
            upsetDatadf.loc[condition_indices[alignedUnique[nos]:], "forUnique01"] = "Unaligned species-specific IMs"
    
    # 将共享的pG4s标记出来
    upsetDatadf["forUnique01"] = upsetDatadf["forUnique01"].fillna("Shared IMs")
    upsetDatadf.set_index(species, inplace=True)
    return upsetDatadf

def analyze_intersection_patterns(presAbs_df, species_list, title="Whole Genome"):
    """分析并输出交集模式的数量统计"""
    
    print(f"\n{'='*80}")
    print(f"Whole Genome Intersection Pattern Analysis")
    print('='*80)
    
    # 统计每种模式的数量
    pattern_counts = {}
    total_ims = len(presAbs_df)
    
    # 将DataFrame转换为列表以便处理
    data_matrix = presAbs_df[species_list].values
    
    # 统计所有可能的交集模式
    for i in range(len(data_matrix)):
        # 创建模式字符串，例如："111000" 表示前三个物种存在，后三个不存在
        pattern_str = ''.join(['1' if x == 1 else '0' for x in data_matrix[i]])
        
        if pattern_str not in pattern_counts:
            pattern_counts[pattern_str] = 0
        pattern_counts[pattern_str] += 1
    
    # 按数量排序
    sorted_patterns = sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)
    
    # 输出控制台信息
    print(f"\nTotal IMs: {total_ims:,}")
    print(f"\nDifferent intersection patterns: {len(pattern_counts)}")
    print("\nDetailed statistics:")
    print("-" * 100)
    print(f"{'Pattern':<15} {'Binary':<25} {'Species Combination':<40} {'Count':>10} {'Percentage':>10}")
    print("-" * 100)
    
    pattern_data = []
    for pattern_str, count in sorted_patterns:
        # 解析模式字符串
        species_present = []
        pattern_binary = ' '.join(pattern_str[i:i+1] for i in range(0, len(pattern_str), 1))
        
        for i, species in enumerate(species_list):
            if pattern_str[i] == '1':
                species_present.append(species)
        
        # 生成易读的物种组合描述
        if len(species_present) == 0:
            species_combo = "No species"
        elif len(species_present) == 1:
            species_combo = f"{species_present[0]} (specific)"
        elif len(species_present) == len(species_list):
            species_combo = "All species shared"
        else:
            species_combo = ' + '.join(species_present)
        
        percentage = (count / total_ims) * 100
        print(f"{pattern_str:<15} {pattern_binary:<25} {species_combo:<40} {count:>10,} {percentage:>9.2f}%")
        
        pattern_data.append({
            'pattern': pattern_str,
            'pattern_binary': pattern_binary,
            'species_count': len(species_present),
            'species_combination': species_combo,
            'count': count,
            'percentage': percentage
        })
    
    # 按物种数量分组统计
    print("\nStatistics by number of shared species:")
    print("-" * 60)
    
    species_count_stats = defaultdict(int)
    for pattern_str, count in pattern_counts.items():
        species_count = pattern_str.count('1')
        species_count_stats[species_count] += count
    
    for count in sorted(species_count_stats.keys()):
        total = species_count_stats[count]
        percentage = (total / total_ims) * 100
        if count == 1:
            print(f"{count} species (species-specific): {total:>15,} ({percentage:>6.2f}%)")
        elif count == len(species_list):
            print(f"{count} species (all species shared): {total:>12,} ({percentage:>6.2f}%)")
        else:
            print(f"{count} species: {total:>25,} ({percentage:>6.2f}%)")
    
    # 输出摘要统计
    print("\nSummary statistics:")
    print("-" * 60)
    print(f"Total IMs: {total_ims:,}")
    print(f"Total species: {len(species_list)}")
    print(f"Total patterns: {len(pattern_counts)}")
    print(f"Average IMs per pattern: {total_ims/len(pattern_counts):,.1f}")
    
    return pd.DataFrame(pattern_data), total_ims

def generate_whole_genome_upset(df, alignedUniqueGQs, output_dir="plots/whole_genome", stats_dir="stats/whole_genome"):
    '''生成全基因组Upset图，输出可编辑PDF、高分辨率PNG和TIFF'''
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)
    
    # 物种顺序 - 与原脚本完全一致
    order = ["P. paniscus", "P. troglodytes", "H. sapiens", "G. gorilla", "P. pygmaeus", "P. abelii"]
    
    print(f"\n{'='*100}")
    print(f"Processing: Whole Genome")
    print('='*100)
    
    # 生成存在-缺失矩阵
    print("Generating presence-absence matrix...")
    presabs = []
    grouped_df = df.groupby("ID")
    
    total_groups = len(df.groupby("ID"))
    print(f"Total IMs groups: {total_groups:,}")
    
    for idx, (name, group) in enumerate(grouped_df):
        if (idx + 1) % 100000 == 0:
            print(f"  Processing progress: {idx + 1:,}/{total_groups:,} ({((idx+1)/total_groups*100):.1f}%)")
        
        specieslist = [0] * 6
        presspecies = group["SPECIES"].unique()
        for index in presspecies:
            specieslist[index-1] = 1
        presabs.append(specieslist)
    
    presAbs_df = pd.DataFrame(presabs, columns=[i[3] for i in speciesnos.values()])
    presAbs_dfUpset = presAbs_df[order]
    speciesUpset = order
    
    # 分析交集模式
    print("\nAnalyzing intersection patterns...")
    pattern_df, total_ims = analyze_intersection_patterns(presAbs_dfUpset, speciesUpset, "Whole Genome")
    
    # 保存模式统计到CSV
    stats_filename = "whole_genome_intersection_patterns.csv"
    stats_path = os.path.join(stats_dir, stats_filename)
    pattern_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
    print(f"\nDetailed statistics saved to: {stats_path}")
    
    # 获取对齐的物种特异性pG4s数量
    print("\nGetting aligned species-specific pG4s counts...")
    alignedUniqueUpset = []
    
    for s in range(1, 7):
        # 全基因组：选择所有该物种的记录
        specAlignedUniqueGQs = alignedUniqueGQs[alignedUniqueGQs['species'] == s]
        count = specAlignedUniqueGQs.shape[0]
        alignedUniqueUpset.append(count)
        print(f"  {speciesnos[s][3]}: {count:,} aligned species-specific pG4s")
    
    # 根据order列表重新排列 - 与原脚本完全一致
    human_element = alignedUniqueUpset.pop(0)
    alignedUniqueUpset.insert(2, human_element)
    sorang_element = alignedUniqueUpset.pop(4)
    alignedUniqueUpset.insert(5, sorang_element)
    alignedUniqueUpset = np.array(alignedUniqueUpset)
    
    # 生成Upset图数据
    print("\nGenerating Upset plot data...")
    presAbsMatrixNormUpsetDict = {}
    for column in presAbs_dfUpset.columns:
        indices = [i for i, value in enumerate(presAbs_dfUpset[column]) if value == 1]
        presAbsMatrixNormUpsetDict[column] = indices
    
    # 创建Upset图数据
    upsetData = upsplt.from_contents(presAbsMatrixNormUpsetDict)
    upsetData["forUnique01"] = "Shared IMs"
    upsetData = stackedBarUpset(upsetData, speciesUpset, alignedUniqueUpset)
    
    # 创建Upset图 - 与原脚本完全相同的参数
    print("Creating Upset plot...")
    upset = upsplt.UpSet(
        upsetData, 
        sort_by="cardinality", 
        sort_categories_by="-input", 
        facecolor=CBPalette["Black"],
        show_counts=False, 
        totals_plot_elements=6, 
        intersection_plot_elements=0
    )
    
    # 添加堆叠条形图 - 与原脚本完全相同的参数
    upset.add_stacked_bars(
        by="forUnique01", 
        colors=[CBPalette["Dark blue"], CBPalette["Vermilion"], CBPalette["Light blue"]], 
        elements=10
    )
    
    # 绘制图形 - 与原脚本相同的尺寸
    fig = plt.figure(figsize=(12, 8))
    plot_result = upset.plot(fig=fig)
    
    # 设置堆叠条形图区域的格式 - 与原脚本完全一致
    plot_result["extra0"].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["extra0"].set_yticklabels(plot_result["extra0"].get_yticklabels(), fontsize=11)
    
    # 设置矩阵区域的格式 - 与原脚本完全一致
    plot_result["matrix"].set_yticklabels(
        ["S. orangutan", "B. orangutan", "Gorilla", "Human", "Chimpanzee", "Bonobo"], 
        fontsize=12
    )
    
    # 设置总数区域的格式 - 与原脚本完全一致
    plot_result["totals"].xaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["totals"].set_xlabel("\nTotal IMs in species", fontsize=16)
    plot_result["totals"].set_xticklabels(plot_result["totals"].get_xticklabels(), fontsize=7)
    
    # 设置整体标签 - 与原脚本完全一致
    plt.ylabel("Number of IMs\n", fontsize=18)
    plt.grid(alpha=0.5, linestyle="--")
    
    # ========== 图例位置调整说明 ==========
    # upsetplot的图例通常位于堆叠条形图所在的子图（"extra0"）中。
    # 您可以通过以下方式获取图例句柄并调整位置：
    #   legend = plot_result["extra0"].get_legend()
    #   if legend:
    #       legend.set_bbox_to_anchor((1.05, 1))  # 调整到子图右侧
    # 或者使用更精细的参数：legend.set_bbox_to_anchor((x, y), loc='upper left')
    # 如果您想完全控制，也可以直接创建新图例：
    #   handles, labels = plot_result["extra0"].get_legend_handles_labels()
    #   plot_result["extra0"].legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1))
    # 注意：需要在plt.savefig之前修改。
    # 默认情况下，upsetplot会自动放置图例，您可以根据需要取消下面代码的注释进行调整。
    #
    # 示例：将图例移到图形右上角外部
    legend = plot_result["extra0"].get_legend()
    if legend:
        legend.set_bbox_to_anchor((0.25, 1))
    # =====================================
    
    # 保存可编辑PDF
    pdf_path = os.path.join(output_dir, "whole_genome.filtered_upset.pdf")
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    print(f"\nPDF saved (editable text): {pdf_path}")
    
    # 保存SVG（矢量格式）
    svg_path = os.path.join(output_dir, "whole_genome.filtered_upset.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    print(f"SVG format saved: {svg_path}")
    
    # 保存高分辨率PNG（300 dpi）
    png_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_300dpi.png")
    plt.savefig(png_highres_path, format='png', dpi=300, bbox_inches='tight')
    print(f"High-resolution PNG saved (300 dpi): {png_highres_path}")
    
    # 保存高分辨率TIFF（600 dpi，LZW压缩）- 出版首选
    # 注意：需要安装pillow库以支持TIFF格式和LZW压缩
    try:
        tiff_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_600dpi.tiff")
        plt.savefig(tiff_highres_path, format='tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
        print(f"High-resolution TIFF saved (600 dpi, LZW): {tiff_highres_path}")
    except Exception as e:
        print(f"Warning: Could not save TIFF. Please ensure pillow is installed: pip install pillow")
        print(f"Error: {e}")
    
    plt.close()
    
    return pattern_df, total_ims

def main():
    """主函数 - 专门处理全基因组数据，保持与原脚本完全相同的逻辑"""
    
    print("=== Whole Genome G4 Sharing Pattern Upset Plot Generation Script ===")
    print(f"Project directory: {BASE_DIR}")
    
    # 创建输出目录 - 使用与原脚本相同的目录结构
    stats_dir = os.path.join(BASE_DIR, "output/stats")
    plots_dir = os.path.join(BASE_DIR, "output/plots/upsetPlots")
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    # 加载对齐的物种特异性pG4s数据集
    print("\n1. Loading aligned species-specific pG4s data...")
    aligned_unique_path = os.path.join(BASE_DIR, "output/datasets/alignedUniquehsaG.egs")
    
    if not os.path.exists(aligned_unique_path):
        print(f"Error: File does not exist - {aligned_unique_path}")
        sys.exit(1)
    
    try:
        alignedUniquedf = pd.read_csv(aligned_unique_path, sep='\t', header=None)
        alignedUniqueGQs = splitIDsandarrange(alignedUniquedf)
        print(f"  Successfully loaded: {len(alignedUniqueGQs)} records")
    except Exception as e:
        print(f"  Failed to load: {e}")
        sys.exit(1)
    
    # 加载全基因组数据
    print("\n2. Loading whole genome data...")
    whole_genome_path = os.path.join(BASE_DIR, "output/datasets/allhsaG.graph.df")
    
    if not os.path.exists(whole_genome_path):
        print(f"Error: Whole genome data file does not exist - {whole_genome_path}")
        sys.exit(1)
    
    try:
        # 读取完整数据，与原脚本保持一致
        df_whole = pd.read_csv(whole_genome_path, header=0, sep="\t", low_memory=False)
        print(f"  Successfully loaded: {len(df_whole)} rows")
    except Exception as e:
        print(f"  Failed to load: {e}")
        sys.exit(1)
    
    # 检查数据质量
    print("\n3. Checking data quality...")
    print(f"  Unique ID count: {df_whole['ID'].nunique():,}")
    print(f"  Species distribution:")
    species_counts = df_whole['SPECIES'].value_counts().sort_index()
    for species_id, count in species_counts.items():
        species_name = speciesnos.get(species_id, ['Unknown'])[3]
        print(f"    {species_name}: {count:,} records")
    
    # 生成全基因组Upset图
    print("\n4. Generating whole genome Upset plot...")
    try:
        pattern_df, total_ims = generate_whole_genome_upset(
            df_whole, alignedUniqueGQs, plots_dir, stats_dir
        )
        
        # 输出最终摘要
        print(f"\n{'='*100}")
        print("Whole genome analysis completed!")
        print('='*100)
        print(f"Total IMs: {total_ims:,}")
        print(f"Unique sharing patterns: {len(pattern_df)}")
        print(f"Average IMs per pattern: {pattern_df['count'].mean():,.1f}")
        
        print(f"\nOutput directories:")
        print(f"  Plots: {plots_dir}")
        print(f"  Statistics: {stats_dir}")
        print('='*100)
        
    except Exception as e:
        print(f"\nError generating Upset plot: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
