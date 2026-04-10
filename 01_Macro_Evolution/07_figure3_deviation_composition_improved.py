#!/usr/bin/env python3
"""
figure3_deviation_composition.py - 创建Figure 3: 各生物类别偏离物种组成百分比
修改内容：
1. 对于细菌的富集比例太小的问题，用红色字体标记在上面
2. 柱状图直接顶格，不需要标记n的数量
3. 排序改为从低到高
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
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'legend.fontsize': 9,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'axes.linewidth': 0.8,
        'figure.dpi': 150,
        'savefig.dpi': 150,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
    })

# 颜色方案
four_category_colors = {
    'Bacteria': '#2CA02C',     # 绿色
    'Archaea': '#1F77B4',      # 蓝色
    'Eukaryota': '#D62728',    # 红色
}

nine_category_colors = {
    'Bacteria': '#2CA02C',     # 绿色
    'Archaea': '#1F77B4',      # 蓝色
    'Fungi': '#9467BD',        # 紫色
    'Plant': '#8C564B',        # 棕色
    'Invertebrate': '#E377C2', # 粉色
    'Protozoa': '#BCBD22',     # 黄绿色
    'Vertebrate Other': '#7F7F7F',  # 灰色
    'Mammalian': '#D62728',    # 红色
}

# 偏离类别颜色
deviation_colors = {
    'Enriched': '#D62728',     # 红色
    'Normal': '#7F7F7F',       # 灰色
    'Depleted': '#1F77B4'      # 蓝色
}

class GCBackgroundComparison:
    def __init__(self, simulation_file, real_genome_file, z_threshold=1.96):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.z_threshold = z_threshold
        self.confidence_level = (1 - 2 * (1 - norm.cdf(z_threshold))) * 100
        
        # 四分类（去掉病毒）
        self.four_categories = ['Bacteria', 'Archaea', 'Eukaryota']
        
        # 九分类（更新为Mammalian，去掉病毒）
        self.nine_categories = [
            'Bacteria', 'Archaea', 
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        # 真核细分类群（用于映射到Eukaryota）
        self.eukaryotic_subcategories = [
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        # 数据
        self.simulation_data = None
        self.real_data = None
        self.real_data_4cat = None  # 四分类数据
        self.real_data_9cat = None  # 九分类数据
        self.background_model = None
        self.background_stats = None
        
        # 用于可视化
        self.smooth_background_x = None
        self.smooth_background_y = None
    
    def load_and_preprocess_data(self):
        """加载并预处理数据"""
        print("=== 加载数据 ===")
        
        # 加载模拟数据（理论背景）
        self.simulation_data = pd.read_csv(self.simulation_file)
        print(f"理论背景数据: {len(self.simulation_data)} 行")
        print(f"模拟数据列: {list(self.simulation_data.columns)}")
        
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
            print("将GC含量从百分比转换为小数")
        
        # 创建四分类数据（去掉病毒）
        self.real_data_4cat = self.real_data.copy()
        
        # 将真核细分类群映射到Eukaryota
        eukaryotic_mask = self.real_data_4cat['Classification'].isin(self.eukaryotic_subcategories)
        self.real_data_4cat.loc[eukaryotic_mask, 'Classification_4cat'] = 'Eukaryota'
        
        # 其他类别保持不变，但排除病毒
        other_categories = ['Bacteria', 'Archaea']
        for cat in other_categories:
            mask = self.real_data_4cat['Classification'] == cat
            self.real_data_4cat.loc[mask, 'Classification_4cat'] = cat
        
        # 过滤掉病毒和其他未分类的
        self.real_data_4cat = self.real_data_4cat[self.real_data_4cat['Classification_4cat'].isin(self.four_categories)].copy()
        
        # 创建九分类数据（去掉病毒）
        self.real_data_9cat = self.real_data[self.real_data['Classification'].isin(self.nine_categories)].copy()
        self.real_data_9cat['Classification_9cat'] = self.real_data_9cat['Classification']
        
        print(f"成功加载 {len(self.real_data):,} 个真实基因组")
        print(f"四分类数据(无病毒): {len(self.real_data_4cat):,} 个基因组")
        print(f"九分类数据(无病毒): {len(self.real_data_9cat):,} 个基因组")
        
        # 打印各类别数量
        print("\n四分类数量统计(无病毒):")
        for category in self.four_categories:
            count = (self.real_data_4cat['Classification_4cat'] == category).sum()
            if count > 0:
                print(f"  {category}: {count:,} genomes")
        
        print("\n九分类数量统计(无病毒):")
        for category in self.nine_categories:
            count = (self.real_data_9cat['Classification_9cat'] == category).sum()
            if count > 0:
                print(f"  {category}: {count:,} genomes")
        
        return True
    
    def build_background_model(self):
        """从模拟数据构建平滑背景模型"""
        print("\n=== 构建背景模型 ===")
        
        # 检查模拟数据是否有必要的列
        print(f"模拟数据列: {list(self.simulation_data.columns)}")
        
        # 尝试不同的列名
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
            print("将模拟GC含量从百分比转换为小数")
        
        # 移除NaN值
        background_data = background_data.dropna()
        
        # 按GC含量排序
        background_data = background_data.sort_values('GC_content')
        
        print(f"背景数据点: {len(background_data)}")
        print(f"GC含量范围: {background_data['GC_content'].min():.3f} - {background_data['GC_content'].max():.3f}")
        print(f"密度范围: {background_data['density'].min():.3f} - {background_data['density'].max():.3f}")
        
        # 提取值
        gc_values = background_data['GC_content'].values
        density_values = background_data['density'].values
        
        # 按GC含量排序
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]
        
        # 创建更多点的平滑曲线
        self.smooth_background_x = np.linspace(gc_values.min(), gc_values.max(), 500)
        
        # 使用Savitzky-Golay滤波器平滑曲线
        if len(gc_values) > 10:
            # 应用Savitzky-Golay滤波器
            window_length = min(15, len(gc_values) - 1)
            if window_length % 2 == 0:
                window_length -= 1  # 确保为奇数
            
            smoothed_density = savgol_filter(density_values, window_length=window_length, polyorder=3)
            
            # 使用UnivariateSpline进行非常平滑的插值
            try:
                # 使用小的平滑因子
                spline = UnivariateSpline(gc_values, smoothed_density, s=0.001, k=3)
                self.smooth_background_y = spline(self.smooth_background_x)
                print("使用UnivariateSpline + Savitzky-Golay创建非常平滑的背景模型")
            except:
                # 回退到三次样条
                cs = interp1d(gc_values, smoothed_density, kind='cubic', fill_value='extrapolate')
                self.smooth_background_y = cs(self.smooth_background_x)
                print("使用三次样条创建平滑背景模型")
        else:
            # 直接使用三次样条
            cs = interp1d(gc_values, density_values, kind='cubic', fill_value='extrapolate')
            self.smooth_background_y = cs(self.smooth_background_x)
            print("使用三次样条创建平滑背景模型")
        
        # 存储用于预测（使用线性插值避免过拟合）
        self.background_model = interp1d(
            gc_values, density_values,
            kind='linear', 
            bounds_error=False,
            fill_value=(density_values[0], density_values[-1])
        )
        
        # 保存背景统计
        self.background_stats = background_data
        
        print(f"背景模型构建成功，覆盖GC范围: {gc_values.min():.3f} - {gc_values.max():.3f}")
        
        return background_data
    
    def calculate_deviations(self):
        """计算偏离指标"""
        print(f"\n=== 计算偏离指标 (Z阈值={self.z_threshold}, {self.confidence_level:.1f}%置信区间) ===")
        
        # 为四分类数据计算偏离
        self.real_data_4cat['expected_density'] = self.background_model(self.real_data_4cat['GC_content'])
        deviations_4cat = self.real_data_4cat['Genomic_density'] - self.real_data_4cat['expected_density']
        
        # 稳健Z-score（使用四分类数据的整体分布）
        median_deviation = deviations_4cat.median()
        mad = (deviations_4cat - median_deviation).abs().median()
        
        if mad > 0:
            self.real_data_4cat['robust_z_score'] = (deviations_4cat - median_deviation) / (1.4826 * mad)
            print(f"稳健Z-score计算: median_deviation={median_deviation:.2f}, MAD={mad:.2f}")
            
            # 对九分类数据使用相同的Z-score计算方法
            self.real_data_9cat['expected_density'] = self.background_model(self.real_data_9cat['GC_content'])
            deviations_9cat = self.real_data_9cat['Genomic_density'] - self.real_data_9cat['expected_density']
            self.real_data_9cat['robust_z_score'] = (deviations_9cat - median_deviation) / (1.4826 * mad)
        else:
            self.real_data_4cat['robust_z_score'] = 0
            self.real_data_9cat['robust_z_score'] = 0
        
        # 分类（四分类）
        conditions = [
            self.real_data_4cat['robust_z_score'] > self.z_threshold,
            self.real_data_4cat['robust_z_score'] < -self.z_threshold,
            (self.real_data_4cat['robust_z_score'] >= -self.z_threshold) & 
            (self.real_data_4cat['robust_z_score'] <= self.z_threshold)
        ]
        choices = ['Enriched', 'Depleted', 'Normal']
        self.real_data_4cat['deviation_category'] = np.select(conditions, choices, default='Normal')
        
        # 分类（九分类）
        conditions_9cat = [
            self.real_data_9cat['robust_z_score'] > self.z_threshold,
            self.real_data_9cat['robust_z_score'] < -self.z_threshold,
            (self.real_data_9cat['robust_z_score'] >= -self.z_threshold) & 
            (self.real_data_9cat['robust_z_score'] <= self.z_threshold)
        ]
        self.real_data_9cat['deviation_category'] = np.select(conditions_9cat, choices, default='Normal')
        
        # 统计
        print("\n四分类偏离统计(无病毒):")
        category_counts = self.real_data_4cat['deviation_category'].value_counts()
        for category, count in category_counts.items():
            percentage = count / len(self.real_data_4cat) * 100
            print(f"  {category}: {count:,} genomes ({percentage:.1f}%)")
        
        print("\n九分类偏离统计(无病毒):")
        category_counts_9 = self.real_data_9cat['deviation_category'].value_counts()
        for category, count in category_counts_9.items():
            percentage = count / len(self.real_data_9cat) * 100
            print(f"  {category}: {count:,} genomes ({percentage:.1f}%)")
        
        return True
    
    def create_figure_3_deviation_composition(self):
        """创建图形3：各生物类别偏离物种组成百分比（AC版本）"""
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # 3A: 四分类偏离组成（堆积柱状图）
        self._plot_deviation_composition_4cat_stacked(axes[0])
        
        # 3B: 九分类偏离组成（堆积柱状图）
        self._plot_deviation_composition_9cat_stacked(axes[1])
        
        plt.suptitle('Figure 3: Percentage Composition of Deviation Categories (without Viral)', fontsize=14, y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
    
    def _plot_deviation_composition_4cat_stacked(self, ax):
        """绘制四分类偏离组成（堆积柱状图，在富集部分标记百分比）"""
        # 计算各分类的偏离组成
        composition_data = {}
        
        for category in self.four_categories:
            cat_data = self.real_data_4cat[self.real_data_4cat['Classification_4cat'] == category]
            if len(cat_data) > 0:
                total = len(cat_data)
                composition = {}
                for dev_cat in ['Enriched', 'Normal', 'Depleted']:
                    count = (cat_data['deviation_category'] == dev_cat).sum()
                    composition[dev_cat] = count / total * 100
                composition_data[category] = composition
        
        if composition_data:
            categories = list(composition_data.keys())
            enriched_pcts = [composition_data[c]['Enriched'] for c in categories]
            normal_pcts = [composition_data[c]['Normal'] for c in categories]
            depleted_pcts = [composition_data[c]['Depleted'] for c in categories]
            
            x = np.arange(len(categories))
            bottom = np.zeros(len(categories))
            
            # 绘制堆积柱状图
            enriched_bars = ax.bar(x, enriched_pcts, label='Enriched', color=deviation_colors['Enriched'], 
                                  alpha=0.8, bottom=bottom)
            bottom += enriched_pcts
            
            normal_bars = ax.bar(x, normal_pcts, label='Normal', color=deviation_colors['Normal'], 
                                alpha=0.8, bottom=bottom)
            bottom += normal_pcts
            
            depleted_bars = ax.bar(x, depleted_pcts, label='Depleted', color=deviation_colors['Depleted'], 
                                  alpha=0.8, bottom=bottom)
            
            # 在富集部分标记富集百分比
            for i, (bar, pct) in enumerate(zip(enriched_bars, enriched_pcts)):
                height = bar.get_height()
                if height > 0:
                    # 对于富集比例很小的情况（如Bacteria: 0.9%），用红色字体标记在柱子上方
                    if pct < 5:  # 小于5%的情况
                        # 在柱子顶部上方标记，使用红色字体
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height + 1,
                               f'{pct:.0f}%', ha='center', va='bottom', fontsize=9, 
                               fontweight='bold', color='red')
                    else:
                        # 在柱子内部标记，使用白色字体
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height/2,
                               f'{pct:.0f}%', ha='center', va='center', fontsize=9, 
                               fontweight='bold', color='white')
            
            ax.set_xlabel('Biological Category')
            ax.set_ylabel('Percentage (%)')
            ax.set_title('A. Four-Category Deviation Composition (without Viral)')
            ax.set_xticks(x)
            ax.set_xticklabels(categories)
            ax.legend(fontsize=8, loc='upper left')
            
            # 设置y轴范围，让柱子直接顶格
            ax.set_ylim(0, 105)  # 留5%的空间给顶部的标签
            
            ax.grid(True, alpha=0.3, axis='y', linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('A. Four-Category Deviation Composition')
    
    def _plot_deviation_composition_9cat_stacked(self, ax):
        """绘制九分类偏离组成（堆积柱状图，在富集部分标记百分比，按富集比例从低到高排序）"""
        # 计算各分类的偏离组成
        composition_data = {}
        
        for category in self.nine_categories:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            if len(cat_data) > 0:
                total = len(cat_data)
                composition = {}
                for dev_cat in ['Enriched', 'Normal', 'Depleted']:
                    count = (cat_data['deviation_category'] == dev_cat).sum()
                    composition[dev_cat] = count / total * 100
                composition_data[category] = composition
        
        if composition_data:
            # 按富集比例从低到高排序
            categories_sorted = sorted(composition_data.keys(), 
                                      key=lambda x: composition_data[x]['Enriched'], 
                                      reverse=False)  # reverse=False表示从低到高
            
            categories = categories_sorted
            enriched_pcts = [composition_data[c]['Enriched'] for c in categories]
            normal_pcts = [composition_data[c]['Normal'] for c in categories]
            depleted_pcts = [composition_data[c]['Depleted'] for c in categories]
            
            x = np.arange(len(categories))
            bottom = np.zeros(len(categories))
            
            # 绘制堆积柱状图
            enriched_bars = ax.bar(x, enriched_pcts, label='Enriched', color=deviation_colors['Enriched'], 
                                  alpha=0.8, bottom=bottom)
            bottom += enriched_pcts
            
            normal_bars = ax.bar(x, normal_pcts, label='Normal', color=deviation_colors['Normal'], 
                                alpha=0.8, bottom=bottom)
            bottom += normal_pcts
            
            depleted_bars = ax.bar(x, depleted_pcts, label='Depleted', color=deviation_colors['Depleted'], 
                                  alpha=0.8, bottom=bottom)
            
            # 在富集部分标记富集百分比
            for i, (bar, pct) in enumerate(zip(enriched_bars, enriched_pcts)):
                height = bar.get_height()
                if height > 0:
                    # 对于富集比例很小的情况（如Bacteria: 0.9%），用红色字体标记在柱子上方
                    if pct < 5:  # 小于5%的情况
                        # 在柱子顶部上方标记，使用红色字体
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height + 1,
                               f'{pct:.0f}%', ha='center', va='bottom', fontsize=9, 
                               fontweight='bold', color='red')
                    else:
                        # 在柱子内部标记，使用白色字体
                        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + height/2,
                               f'{pct:.0f}%', ha='center', va='center', fontsize=9, 
                               fontweight='bold', color='white')
            
            ax.set_xlabel('Biological Category')
            ax.set_ylabel('Percentage (%)')
            ax.set_title('B. Nine-Category Deviation Composition (Sorted by Enrichment from Low to High, without Viral)')
            ax.set_xticks(x)
            ax.set_xticklabels(categories, rotation=45, ha='right')
            ax.legend(fontsize=8, loc='upper left')
            
            # 设置y轴范围，让柱子直接顶格
            ax.set_ylim(0, 105)  # 留5%的空间给顶部的标签
            
            ax.grid(True, alpha=0.3, axis='y', linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('B. Nine-Category Deviation Composition')

def main():
    """主函数"""
    set_clean_style()
    
    # 文件路径
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"
    
    print("="*80)
    print("创建Figure 3: 各生物类别偏离物种组成百分比")
    print(f"Z阈值: {1.96} (95%置信区间)")
    print("="*80)
    
    # 初始化分析器
    analyzer = GCBackgroundComparison(simulation_file, real_genome_file, z_threshold=1.96)
    
    try:
        # 执行分析
        analyzer.load_and_preprocess_data()
        analyzer.build_background_model()
        analyzer.calculate_deviations()
        
        # 创建Figure 3
        print("\n=== 创建Figure 3 ===")
        print("创建图3: 各生物类别偏离物种组成百分比（AC版本）...")
        fig3 = analyzer.create_figure_3_deviation_composition()
        
        # 导出结果
        output_dir = Path(base_dir) / "figure3_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存图形
        print("\n=== 保存图形 ===")
        fig3.savefig(output_dir / "figure3_deviation_composition_AC.png", dpi=300, bbox_inches='tight')
        fig3.savefig(output_dir / "figure3_deviation_composition_AC.pdf", bbox_inches='tight')
        print(f"✓ 图3 (AC版本): {output_dir}/figure3_deviation_composition_AC.png/pdf")
        
        # 显示图形
        plt.show()
        
        print("\n" + "="*80)
        print("✅ Figure 3 创建完成! 图形已保存。")
        print("="*80)
        
    except Exception as e:
        print(f"\n❌ 分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
