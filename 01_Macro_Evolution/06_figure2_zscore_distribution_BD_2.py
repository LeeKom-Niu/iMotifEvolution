#!/usr/bin/env python3
"""
figure2_zscore_distribution_BD.py - 生成稳健Z-score核密度估计分布图（BD版本）
修改：从第一张图（四分类图）中去掉病毒分类
包含：
B. 三分类Z-score分布（Bacteria, Archaea, Eukaryota）
D. 真核分类Z-score分布（原生动物在最前面）
保持原始风格，所有类别都有阴影
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import spearmanr, norm, gaussian_kde
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.signal import savgol_filter
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
import json

# 设置简洁样式
def set_clean_style():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'font.size': 9,
        'pdf.fonttype': 42,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })

# 颜色方案（保持与原脚本一致，但去掉病毒颜色）
three_category_colors = {
    'Bacteria': '#2CA02C',     # 绿色
    'Archaea': '#1F77B4',      # 蓝色
    'Eukaryota': '#D62728',    # 红色
}

# 病毒颜色（备用，但不在第一张图中使用）
viral_color = '#17BECF'        # 青色

nine_category_colors = {
    'Bacteria': '#2CA02C',     # 绿色
    'Archaea': '#1F77B4',      # 蓝色
    'Viral': '#17BECF',        # 青色
    'Fungi': '#9467BD',        # 紫色
    'Plant': '#8C564B',        # 棕色
    'Invertebrate': '#E377C2', # 粉色
    'Protozoa': '#BCBD22',     # 黄绿色
    'Vertebrate Other': '#7F7F7F',  # 灰色
    'Mammalian': '#D62728',    # 红色（哺乳动物保持红色）
}

class ZScoreDistributionBD:
    def __init__(self, simulation_file, real_genome_file, z_threshold=1.96):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.z_threshold = z_threshold
        
        # 三分类（去掉病毒）
        self.three_categories = ['Bacteria', 'Archaea', 'Eukaryota']
        
        # 九分类（保持原样）
        self.nine_categories = [
            'Bacteria', 'Archaea', 'Viral', 
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        # 真核细分类群（原生动物在最前面，保持原始风格）
        self.eukaryotic_categories = [
            'Protozoa',     # 原生动物放在最前面
            'Fungi', 
            'Plant', 
            'Invertebrate',
            'Vertebrate Other', 
            'Mammalian'     # 哺乳动物保持红色
        ]
        
        # 数据
        self.simulation_data = None
        self.real_data = None
        self.real_data_3cat = None  # 三分类数据（修改：3cat代替4cat）
        self.real_data_9cat = None  # 九分类数据
        self.background_model = None
    
    def load_and_preprocess_data(self):
        """加载并预处理数据"""
        print("=== 加载数据 ===")
        
        # 加载模拟数据（理论背景）
        self.simulation_data = pd.read_csv(self.simulation_file)
        print(f"理论背景数据: {len(self.simulation_data)} 行")
        
        # 加载真实基因组数据
        print(f"加载真实基因组数据: {self.real_genome_file}")
        self.real_data = pd.read_csv(self.real_genome_file, sep='\t')
        
        # 标准化列名
        column_mapping = {}
        for col in self.real_data.columns:
            col_lower = col.lower()
            if 'genome' in col_lower or 'species' in col_lower:
                column_mapping['Genome'] = col
            elif 'class' in col_lower or 'domain' in col_lower:
                column_mapping['Classification'] = col
            elif 'gc' in col_lower:
                column_mapping['GC_content'] = col
            elif 'density' in col_lower:
                column_mapping['Genomic_density'] = col
        
        for new_name, old_name in column_mapping.items():
            if old_name in self.real_data.columns:
                self.real_data = self.real_data.rename(columns={old_name: new_name})
        
        # 清理数据
        if 'Classification' in self.real_data.columns:
            self.real_data = self.real_data[self.real_data['Classification'] != 'Classification']
            self.real_data['Classification'] = self.real_data['Classification'].str.strip()
            
            # 将Vertebrate Mammalian替换为Mammalian
            self.real_data['Classification'] = self.real_data['Classification'].replace(
                'Vertebrate Mammalian', 'Mammalian'
            )
        
        # 转换数据类型
        self.real_data['GC_content'] = pd.to_numeric(self.real_data['GC_content'], errors='coerce')
        self.real_data['Genomic_density'] = pd.to_numeric(self.real_data['Genomic_density'], errors='coerce')
        
        # GC含量百分比转小数
        if self.real_data['GC_content'].max() > 1:
            self.real_data['GC_content'] = self.real_data['GC_content'] / 100.0
        
        # 创建三分类数据（修改：3cat代替4cat，去掉病毒）
        self.real_data_3cat = self.real_data.copy()
        
        # 定义真核细分类群
        eukaryotic_subcategories = [
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        # 将真核细分类群映射到Eukaryota
        eukaryotic_mask = self.real_data_3cat['Classification'].isin(eukaryotic_subcategories)
        self.real_data_3cat.loc[eukaryotic_mask, 'Classification_3cat'] = 'Eukaryota'
        
        # 其他类别保持不变（只包括Bacteria和Archaea，去掉Viral）
        other_categories = ['Bacteria', 'Archaea']  # 修改：去掉Viral
        for cat in other_categories:
            mask = self.real_data_3cat['Classification'] == cat
            self.real_data_3cat.loc[mask, 'Classification_3cat'] = cat
        
        # 注意：病毒数据被排除在三分类之外
        
        # 创建九分类数据（保持原样，包含病毒）
        self.real_data_9cat = self.real_data[self.real_data['Classification'].isin(self.nine_categories)].copy()
        self.real_data_9cat['Classification_9cat'] = self.real_data_9cat['Classification']
        
        print(f"成功加载 {len(self.real_data):,} 个真实基因组")
        print(f"三分类数据（无病毒）: {len(self.real_data_3cat.dropna(subset=['Classification_3cat'])):,} 个基因组")
        print(f"九分类数据（包含病毒）: {len(self.real_data_9cat):,} 个基因组")
        
        # 统计各分类数量
        if 'Classification_3cat' in self.real_data_3cat.columns:
            three_cat_counts = self.real_data_3cat['Classification_3cat'].value_counts()
            print(f"\n三分类统计:")
            for cat, count in three_cat_counts.items():
                print(f"  {cat}: {count:,}")
        
        if 'Classification_9cat' in self.real_data_9cat.columns:
            nine_cat_counts = self.real_data_9cat['Classification_9cat'].value_counts()
            print(f"\n九分类统计:")
            for cat, count in nine_cat_counts.items():
                if cat in self.eukaryotic_categories:  # 只打印真核分类
                    print(f"  {cat}: {count:,}")
        
        return True
    
    def build_background_model(self):
        """从模拟数据构建背景模型"""
        print("\n=== 构建背景模型 ===")
        
        # 查找列名
        possible_gc_columns = ['GC_content', 'GC', 'gc_content', 'gc']
        possible_density_columns = ['Density_IM_per_Mb', 'density_IM_per_Mb', 'density', 'Density', 'median_density']
        
        gc_col = None
        density_col = None
        
        for col in possible_gc_columns:
            if col in self.simulation_data.columns:
                gc_col = col
                break
        
        for col in possible_density_columns:
            if col in self.simulation_data.columns:
                density_col = col
                break
        
        if not gc_col or not density_col:
            # 尝试通过模式查找列
            for col in self.simulation_data.columns:
                if 'gc' in col.lower():
                    gc_col = col
                if 'density' in col.lower() or 'median' in col.lower():
                    density_col = col
        
        if not gc_col or not density_col:
            raise ValueError(f"在模拟数据中找不到GC和密度列。可用列: {list(self.simulation_data.columns)}")
        
        print(f"使用GC列: '{gc_col}', 密度列: '{density_col}'")
        
        # 提取数据
        background_data = self.simulation_data[[gc_col, density_col]].copy()
        background_data = background_data.rename(columns={gc_col: 'GC_content', density_col: 'density'})
        
        # 转换GC含量
        if background_data['GC_content'].max() > 1:
            background_data['GC_content'] = background_data['GC_content'] / 100.0
        
        # 移除NaN值
        background_data = background_data.dropna()
        
        # 按GC含量排序
        background_data = background_data.sort_values('GC_content')
        
        # 提取值
        gc_values = background_data['GC_content'].values
        density_values = background_data['density'].values
        
        # 按GC含量排序
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]
        
        # 创建背景模型
        self.background_model = interp1d(
            gc_values, density_values,
            kind='linear', 
            bounds_error=False,
            fill_value=(density_values[0], density_values[-1])
        )
        
        return background_data
    
    def calculate_deviations(self):
        """计算稳健Z-score"""
        print(f"\n=== 计算稳健Z-score ===")
        
        # 为三分类数据计算偏离
        self.real_data_3cat['expected_density'] = self.background_model(self.real_data_3cat['GC_content'])
        deviations_3cat = self.real_data_3cat['Genomic_density'] - self.real_data_3cat['expected_density']
        
        # 稳健Z-score（使用三分类数据的整体分布）
        median_deviation = deviations_3cat.median()
        mad = (deviations_3cat - median_deviation).abs().median()
        
        if mad > 0:
            self.real_data_3cat['robust_z_score'] = (deviations_3cat - median_deviation) / (1.4826 * mad)
            print(f"稳健Z-score计算: median_deviation={median_deviation:.2f}, MAD={mad:.2f}")
            
            # 对九分类数据使用相同的Z-score计算方法
            self.real_data_9cat['expected_density'] = self.background_model(self.real_data_9cat['GC_content'])
            deviations_9cat = self.real_data_9cat['Genomic_density'] - self.real_data_9cat['expected_density']
            self.real_data_9cat['robust_z_score'] = (deviations_9cat - median_deviation) / (1.4826 * mad)
        else:
            self.real_data_3cat['robust_z_score'] = 0
            self.real_data_9cat['robust_z_score'] = 0
        
        return True
    
    def create_figure_2_BD(self):
        """创建图形2：稳健Z-score核密度估计分布（BD版本） - 一行两图"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        
        # 2B: 三分类Z-score分布（修改：从四分类改为三分类）
        self._plot_zscore_by_3categories(axes[0])
        
        # 2D: 真核分类Z-score分布（原生动物在最前面）
        self._plot_zscore_by_eukaryotic_categories_protozoa_first(axes[1])
        
        #plt.suptitle('Figure : Robust Z-score Kernel Density Estimation Distribution', fontsize=14, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
    
    def _plot_zscore_by_3categories(self, ax):
        """绘制三分类Z-score分布（去掉病毒）"""
        # 收集数据
        categories = []
        zscore_data = []
        
        for category in self.three_categories:
            cat_data = self.real_data_3cat[self.real_data_3cat['Classification_3cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 10:
                categories.append(category)
                zscore_data.append(z_scores)
        
        if len(categories) >= 2:
            # 创建KDE曲线
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)
            
            # 首先绘制填充区域（在背景）
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)
            
            # 然后绘制曲线（在前景）
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2.5, label=category)
            
            # 添加阈值线（富集和缺失的图注）
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')
            
            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            #ax.set_title('A. Z-score Distribution by Three Categories')
            
            # 图例：只显示分类和阈值线
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            #ax.set_title('A. Z-score Distribution by Three Categories')
    
    def _plot_zscore_by_eukaryotic_categories_protozoa_first(self, ax):
        """绘制真核分类Z-score分布（原生动物在最前面，保持原始风格）"""
        # 使用自定义顺序：原生动物在最前面
        categories_ordered = self.eukaryotic_categories.copy()
        
        # 收集数据
        categories = []
        zscore_data = []
        
        for category in categories_ordered:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) >= 3:  # 降低阈值以适应真核类别
                categories.append(category)
                zscore_data.append(z_scores)
        
        if len(categories) >= 2:
            # 创建KDE曲线
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)
            
            # 首先绘制填充区域（在背景）
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)
            
            # 然后绘制曲线（在前景）
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2, label=category)
            
            # 添加阈值线（富集和缺失的图注）
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')
            
            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            #ax.set_title('B. Z-score Distribution by Eukaryotic Categories')
            
            # 图例：只显示分类和阈值线（使用两列）
            ax.legend(fontsize=8, loc='upper right', ncol=2)
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            #ax.set_title('B. Z-score Distribution by Eukaryotic Categories')
    
    def perform_statistical_analysis(self):
        """执行统计分析"""
        results = {}
        
        print("\n=== 统计分析 ===")
        
        # 1. 三分类Z-score统计（修改：3cat代替4cat）
        print("\n1. 三分类Z-score统计（无病毒）:")
        for category in self.three_categories:
            cat_data = self.real_data_3cat[self.real_data_3cat['Classification_3cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 0:
                mean_val = z_scores.mean()
                std_val = z_scores.std()
                median_val = z_scores.median()
                n_val = len(z_scores)
                results[f'{category}_zscore'] = {
                    'mean': float(mean_val),
                    'std': float(std_val),
                    'median': float(median_val),
                    'n': n_val
                }
                print(f"   {category}: mean={mean_val:.2f}, std={std_val:.2f}, median={median_val:.2f}, n={n_val:,}")
        
        # 2. 真核分类Z-score统计（按显示顺序）
        print("\n2. 真核分类Z-score统计（原生动物在最前面）:")
        for category in self.eukaryotic_categories:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 0:
                mean_val = z_scores.mean()
                median_val = z_scores.median()
                n_val = len(z_scores)
                results[f'{category}_zscore'] = {
                    'mean': float(mean_val),
                    'median': float(median_val),
                    'n': n_val
                }
                print(f"   {category}: mean={mean_val:.2f}, median={median_val:.2f}, n={n_val:,}")
        
        return results
    
    def export_results(self, output_dir):
        """导出结果"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n=== 导出结果到 {output_dir} ===")
        
        # 保存三分类Z-score数据（修改：3cat代替4cat）
        results_file_3cat = output_dir / "zscore_data_3categories.csv"
        self.real_data_3cat[['Genome', 'Classification_3cat', 'GC_content', 
                           'Genomic_density', 'expected_density', 'robust_z_score']].to_csv(
                               results_file_3cat, index=False)
        print(f"✓ 三分类Z-score数据（无病毒）: {results_file_3cat}")
        
        # 保存真核分类Z-score数据（按显示顺序）
        eukaryotic_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'].isin(self.eukaryotic_categories)]
        # 按显示顺序排序
        eukaryotic_data['Classification_9cat'] = pd.Categorical(
            eukaryotic_data['Classification_9cat'], 
            categories=self.eukaryotic_categories, 
            ordered=True
        )
        eukaryotic_data = eukaryotic_data.sort_values('Classification_9cat')
        
        results_file_eukaryotic = output_dir / "zscore_data_eukaryotic_categories.csv"
        eukaryotic_data[['Genome', 'Classification_9cat', 'GC_content', 
                       'Genomic_density', 'expected_density', 'robust_z_score']].to_csv(
                           results_file_eukaryotic, index=False)
        print(f"✓ 真核分类Z-score数据（原生动物在前）: {results_file_eukaryotic}")
        
        # 保存统计结果
        stats_results = self.perform_statistical_analysis()
        if stats_results:
            stats_file = output_dir / "zscore_statistics.json"
            with open(stats_file, 'w') as f:
                json.dump(stats_results, f, indent=2, ensure_ascii=False)
            print(f"✓ Z-score统计结果: {stats_file}")
        
        return output_dir

def main():
    """主函数"""
    set_clean_style()
    
    # 文件路径
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"
    
    print("="*80)
    print("稳健Z-score核密度估计分布分析 (BD版本)")
    print(f"Z阈值: {1.96} (95%置信区间)")
    print("主要修改：从第一张图（四分类图）中去掉病毒分类")
    print("现在显示：")
    print("  图A: 三分类Z-score分布 (Bacteria, Archaea, Eukaryota)")
    print("  图B: 真核分类Z-score分布（原生动物在最前面）")
    print("="*80)
    
    # 初始化分析器
    analyzer = ZScoreDistributionBD(simulation_file, real_genome_file, z_threshold=1.96)
    
    try:
        # 执行分析
        analyzer.load_and_preprocess_data()
        analyzer.build_background_model()
        analyzer.calculate_deviations()
        
        # 创建BD版本图形
        print("\n=== 创建BD版本图形 ===")
        print("创建图2 (BD版本): 稳健Z-score核密度估计分布...")
        fig = analyzer.create_figure_2_BD()
        
        # 导出结果
        output_dir = Path(base_dir) / "figure2_BD_analysis"
        analyzer.export_results(output_dir)
        
        # 保存图形
        print("\n=== 保存图形 ===")
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.png", dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.pdf", bbox_inches='tight')
        print(f"✓ 图2 (BD版本): {output_dir}/figure2_zscore_distribution_BD.png/pdf")
        
        # 显示统计摘要
        print("\n" + "="*80)
        print("📊 统计摘要")
        print("="*80)
        
        stats = analyzer.perform_statistical_analysis()
        
        print("\n三分类Z-score摘要（无病毒）:")
        for category in analyzer.three_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
                print(f"  {category}:")
                print(f"    均值: {cat_stats['mean']:.2f}, 中位数: {cat_stats['median']:.2f}")
                print(f"    标准差: {cat_stats['std']:.2f}, 样本数: {cat_stats['n']:,}")
        
        print("\n真核分类Z-score摘要（按显示顺序）:")
        for category in analyzer.eukaryotic_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
                print(f"  {category}:")
                print(f"    均值: {cat_stats['mean']:.2f}, 中位数: {cat_stats['median']:.2f}")
                print(f"    样本数: {cat_stats['n']:,}")
        
        print("\n" + "="*80)
        print("✅ 分析完成! BD版本图形和结果已保存。")
        print("="*80)
        
    except Exception as e:
        print(f"\n❌ 分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
