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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd
import upsetplot as upsplt
import warnings
from collections import defaultdict
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
warnings.simplefilter("ignore", category=FutureWarning)
warnings.simplefilter("ignore", category=UserWarning)
BASE_DIR = "/datapool/life-zhanghk/niulk/my_project/18_mammalian/03_upset"
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 9
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
CBPalette = {
    "Black": "
    "Orange": "
    "Light blue": "
    "Vermilion": "
    "Mid blue": "
    "Maroon": "
    "Dark blue": "
    "Light purple": "
    "Light teal": "
    "Purple": "
    "Teal": "
    "Dark purple": "
    "Dark teal": "
    "Grey": "
    "Yellow": "
}
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
        condition = ((upsetDatadf[target] == True) & 
                     (upsetDatadf[[s for s in species if s != target]].eq(False).all(axis=1)))
        condition_indices = upsetDatadf[condition].index
        
        if len(condition_indices) > 0:
            upsetDatadf.loc[condition_indices[:alignedUnique[nos]], "forUnique01"] = "Aligned species-specific IMs"
            upsetDatadf.loc[condition_indices[alignedUnique[nos]:], "forUnique01"] = "Unaligned species-specific IMs"
    
    upsetDatadf["forUnique01"] = upsetDatadf["forUnique01"].fillna("Shared IMs")
    upsetDatadf.set_index(species, inplace=True)
    return upsetDatadf
def analyze_intersection_patterns(presAbs_df, species_list, title="Whole Genome"):
    """分析并输出交集模式的数量统计"""
    
    
    pattern_counts = {}
    total_ims = len(presAbs_df)
    
    data_matrix = presAbs_df[species_list].values
    
    for i in range(len(data_matrix)):
        pattern_str = ''.join(['1' if x == 1 else '0' for x in data_matrix[i]])
        
        if pattern_str not in pattern_counts:
            pattern_counts[pattern_str] = 0
        pattern_counts[pattern_str] += 1
    
    sorted_patterns = sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)
    
    
    pattern_data = []
    for pattern_str, count in sorted_patterns:
        species_present = []
        pattern_binary = ' '.join(pattern_str[i:i+1] for i in range(0, len(pattern_str), 1))
        
        for i, species in enumerate(species_list):
            if pattern_str[i] == '1':
                species_present.append(species)
        
        if len(species_present) == 0:
            species_combo = "No species"
        elif len(species_present) == 1:
            species_combo = f"{species_present[0]} (specific)"
        elif len(species_present) == len(species_list):
            species_combo = "All species shared"
        else:
            species_combo = ' + '.join(species_present)
        
        percentage = (count / total_ims) * 100
        
        pattern_data.append({
            'pattern': pattern_str,
            'pattern_binary': pattern_binary,
            'species_count': len(species_present),
            'species_combination': species_combo,
            'count': count,
            'percentage': percentage
        })
    
    
    species_count_stats = defaultdict(int)
    for pattern_str, count in pattern_counts.items():
        species_count = pattern_str.count('1')
        species_count_stats[species_count] += count
    
    for count in sorted(species_count_stats.keys()):
        total = species_count_stats[count]
        percentage = (total / total_ims) * 100
        if count == 1:
        elif count == len(species_list):
        else:
    
    
    return pd.DataFrame(pattern_data), total_ims
def generate_whole_genome_upset(df, alignedUniqueGQs, output_dir="plots/whole_genome", stats_dir="stats/whole_genome"):
    '''生成全基因组Upset图，输出可编辑PDF、高分辨率PNG和TIFF'''
    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)
    
    order = ["P. paniscus", "P. troglodytes", "H. sapiens", "G. gorilla", "P. pygmaeus", "P. abelii"]
    
    
    presabs = []
    grouped_df = df.groupby("ID")
    
    total_groups = len(df.groupby("ID"))
    
    for idx, (name, group) in enumerate(grouped_df):
        if (idx + 1) % 100000 == 0:
        
        specieslist = [0] * 6
        presspecies = group["SPECIES"].unique()
        for index in presspecies:
            specieslist[index-1] = 1
        presabs.append(specieslist)
    
    presAbs_df = pd.DataFrame(presabs, columns=[i[3] for i in speciesnos.values()])
    presAbs_dfUpset = presAbs_df[order]
    speciesUpset = order
    
    pattern_df, total_ims = analyze_intersection_patterns(presAbs_dfUpset, speciesUpset, "Whole Genome")
    
    stats_filename = "whole_genome_intersection_patterns.csv"
    stats_path = os.path.join(stats_dir, stats_filename)
    pattern_df.to_csv(stats_path, index=False, encoding='utf-8-sig')
    
    alignedUniqueUpset = []
    
    for s in range(1, 7):
        specAlignedUniqueGQs = alignedUniqueGQs[alignedUniqueGQs['species'] == s]
        count = specAlignedUniqueGQs.shape[0]
        alignedUniqueUpset.append(count)
    
    human_element = alignedUniqueUpset.pop(0)
    alignedUniqueUpset.insert(2, human_element)
    sorang_element = alignedUniqueUpset.pop(4)
    alignedUniqueUpset.insert(5, sorang_element)
    alignedUniqueUpset = np.array(alignedUniqueUpset)
    
    presAbsMatrixNormUpsetDict = {}
    for column in presAbs_dfUpset.columns:
        indices = [i for i, value in enumerate(presAbs_dfUpset[column]) if value == 1]
        presAbsMatrixNormUpsetDict[column] = indices
    
    upsetData = upsplt.from_contents(presAbsMatrixNormUpsetDict)
    upsetData["forUnique01"] = "Shared IMs"
    upsetData = stackedBarUpset(upsetData, speciesUpset, alignedUniqueUpset)
    
    upset = upsplt.UpSet(
        upsetData, 
        sort_by="cardinality", 
        sort_categories_by="-input", 
        facecolor=CBPalette["Black"],
        show_counts=False, 
        totals_plot_elements=6, 
        intersection_plot_elements=0
    )
    
    upset.add_stacked_bars(
        by="forUnique01", 
        colors=[CBPalette["Dark blue"], CBPalette["Vermilion"], CBPalette["Light blue"]], 
        elements=10
    )
    
    fig = plt.figure(figsize=(12, 8))
    plot_result = upset.plot(fig=fig)
    
    plot_result["extra0"].yaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["extra0"].set_yticklabels(plot_result["extra0"].get_yticklabels(), fontsize=11)
    
    plot_result["matrix"].set_yticklabels(
        ["S. orangutan", "B. orangutan", "Gorilla", "Human", "Chimpanzee", "Bonobo"], 
        fontsize=12
    )
    
    plot_result["totals"].xaxis.set_major_formatter(FuncFormatter(comma_formatter))
    plot_result["totals"].set_xlabel("\nTotal IMs in species", fontsize=16)
    plot_result["totals"].set_xticklabels(plot_result["totals"].get_xticklabels(), fontsize=7)
    
    plt.ylabel("Number of IMs\n", fontsize=18)
    plt.grid(alpha=0.5, linestyle="--")
    
    legend = plot_result["extra0"].get_legend()
    if legend:
        legend.set_bbox_to_anchor((0.25, 1))
    
    pdf_path = os.path.join(output_dir, "whole_genome.filtered_upset.pdf")
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    
    svg_path = os.path.join(output_dir, "whole_genome.filtered_upset.svg")
    plt.savefig(svg_path, format='svg', transparent=True, bbox_inches='tight')
    
    png_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_300dpi.png")
    plt.savefig(png_highres_path, format='png', dpi=300, bbox_inches='tight')
    
    try:
        tiff_highres_path = os.path.join(output_dir, "whole_genome.filtered_upset_600dpi.tiff")
        plt.savefig(tiff_highres_path, format='tiff', dpi=600, bbox_inches='tight',
                    pil_kwargs={'compression': 'tiff_lzw'})
    except Exception as e:
    
    plt.close()
    
    return pattern_df, total_ims
def main():
    """主函数 - 专门处理全基因组数据，保持与原脚本完全相同的逻辑"""
    
    
    stats_dir = os.path.join(BASE_DIR, "output/stats")
    plots_dir = os.path.join(BASE_DIR, "output/plots/upsetPlots")
    os.makedirs(stats_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    
    aligned_unique_path = os.path.join(BASE_DIR, "output/datasets/alignedUniquehsaG.egs")
    
    if not os.path.exists(aligned_unique_path):
        sys.exit(1)
    
    try:
        alignedUniquedf = pd.read_csv(aligned_unique_path, sep='\t', header=None)
        alignedUniqueGQs = splitIDsandarrange(alignedUniquedf)
    except Exception as e:
        sys.exit(1)
    
    whole_genome_path = os.path.join(BASE_DIR, "output/datasets/allhsaG.graph.df")
    
    if not os.path.exists(whole_genome_path):
        sys.exit(1)
    
    try:
        df_whole = pd.read_csv(whole_genome_path, header=0, sep="\t", low_memory=False)
    except Exception as e:
        sys.exit(1)
    
    species_counts = df_whole['SPECIES'].value_counts().sort_index()
    for species_id, count in species_counts.items():
        species_name = speciesnos.get(species_id, ['Unknown'])[3]
    
    try:
        pattern_df, total_ims = generate_whole_genome_upset(
            df_whole, alignedUniqueGQs, plots_dir, stats_dir
        )
        
        
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        sys.exit(1)
if __name__ == "__main__":
    main()
