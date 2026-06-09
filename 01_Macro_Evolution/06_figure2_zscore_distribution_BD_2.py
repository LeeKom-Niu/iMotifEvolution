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
three_category_colors = {
    'Bacteria': '
    'Archaea': '
    'Eukaryota': '
}
viral_color = '
nine_category_colors = {
    'Bacteria': '
    'Archaea': '
    'Viral': '
    'Fungi': '
    'Plant': '
    'Invertebrate': '
    'Protozoa': '
    'Vertebrate Other': '
    'Mammalian': '
}
class ZScoreDistributionBD:
    def __init__(self, simulation_file, real_genome_file, z_threshold=1.96):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.z_threshold = z_threshold
        
        self.three_categories = ['Bacteria', 'Archaea', 'Eukaryota']
        
        self.nine_categories = [
            'Bacteria', 'Archaea', 'Viral', 
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        self.eukaryotic_categories = [
            'Protozoa',
            'Fungi', 
            'Plant', 
            'Invertebrate',
            'Vertebrate Other', 
            'Mammalian'
        ]
        
        self.simulation_data = None
        self.real_data = None
        self.real_data_3cat = None
        self.real_data_9cat = None
        self.background_model = None
    
    def load_and_preprocess_data(self):
        """加载并预处理数据"""
        
        self.simulation_data = pd.read_csv(self.simulation_file)
        
        self.real_data = pd.read_csv(self.real_genome_file, sep='\t')
        
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
        
        if 'Classification' in self.real_data.columns:
            self.real_data = self.real_data[self.real_data['Classification'] != 'Classification']
            self.real_data['Classification'] = self.real_data['Classification'].str.strip()
            
            self.real_data['Classification'] = self.real_data['Classification'].replace(
                'Vertebrate Mammalian', 'Mammalian'
            )
        
        self.real_data['GC_content'] = pd.to_numeric(self.real_data['GC_content'], errors='coerce')
        self.real_data['Genomic_density'] = pd.to_numeric(self.real_data['Genomic_density'], errors='coerce')
        
        if self.real_data['GC_content'].max() > 1:
            self.real_data['GC_content'] = self.real_data['GC_content'] / 100.0
        
        self.real_data_3cat = self.real_data.copy()
        
        eukaryotic_subcategories = [
            'Fungi', 'Plant', 'Invertebrate', 
            'Protozoa', 'Vertebrate Other', 'Mammalian'
        ]
        
        eukaryotic_mask = self.real_data_3cat['Classification'].isin(eukaryotic_subcategories)
        self.real_data_3cat.loc[eukaryotic_mask, 'Classification_3cat'] = 'Eukaryota'
        
        other_categories = ['Bacteria', 'Archaea']
        for cat in other_categories:
            mask = self.real_data_3cat['Classification'] == cat
            self.real_data_3cat.loc[mask, 'Classification_3cat'] = cat
        
        
        self.real_data_9cat = self.real_data[self.real_data['Classification'].isin(self.nine_categories)].copy()
        self.real_data_9cat['Classification_9cat'] = self.real_data_9cat['Classification']
        
        
        if 'Classification_3cat' in self.real_data_3cat.columns:
            three_cat_counts = self.real_data_3cat['Classification_3cat'].value_counts()
            for cat, count in three_cat_counts.items():
        
        if 'Classification_9cat' in self.real_data_9cat.columns:
            nine_cat_counts = self.real_data_9cat['Classification_9cat'].value_counts()
            for cat, count in nine_cat_counts.items():
                if cat in self.eukaryotic_categories:
        
        return True
    
    def build_background_model(self):
        """从模拟数据构建背景模型"""
        
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
            for col in self.simulation_data.columns:
                if 'gc' in col.lower():
                    gc_col = col
                if 'density' in col.lower() or 'median' in col.lower():
                    density_col = col
        
        if not gc_col or not density_col:
            raise ValueError(f"在模拟数据中找不到GC和密度列。可用列: {list(self.simulation_data.columns)}")
        
        
        background_data = self.simulation_data[[gc_col, density_col]].copy()
        background_data = background_data.rename(columns={gc_col: 'GC_content', density_col: 'density'})
        
        if background_data['GC_content'].max() > 1:
            background_data['GC_content'] = background_data['GC_content'] / 100.0
        
        background_data = background_data.dropna()
        
        background_data = background_data.sort_values('GC_content')
        
        gc_values = background_data['GC_content'].values
        density_values = background_data['density'].values
        
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]
        
        self.background_model = interp1d(
            gc_values, density_values,
            kind='linear', 
            bounds_error=False,
            fill_value=(density_values[0], density_values[-1])
        )
        
        return background_data
    
    def calculate_deviations(self):
        """计算稳健Z-score"""
        
        self.real_data_3cat['expected_density'] = self.background_model(self.real_data_3cat['GC_content'])
        deviations_3cat = self.real_data_3cat['Genomic_density'] - self.real_data_3cat['expected_density']
        
        median_deviation = deviations_3cat.median()
        mad = (deviations_3cat - median_deviation).abs().median()
        
        if mad > 0:
            self.real_data_3cat['robust_z_score'] = (deviations_3cat - median_deviation) / (1.4826 * mad)
            
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
        
        self._plot_zscore_by_3categories(axes[0])
        
        self._plot_zscore_by_eukaryotic_categories_protozoa_first(axes[1])
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        return fig
    
    def _plot_zscore_by_3categories(self, ax):
        """绘制三分类Z-score分布（去掉病毒）"""
        categories = []
        zscore_data = []
        
        for category in self.three_categories:
            cat_data = self.real_data_3cat[self.real_data_3cat['Classification_3cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) > 10:
                categories.append(category)
                zscore_data.append(z_scores)
        
        if len(categories) >= 2:
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)
            
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)
            
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) > 10:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = three_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2.5, label=category)
            
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')
            
            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            
            ax.legend(fontsize=8, loc='upper right')
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
    
    def _plot_zscore_by_eukaryotic_categories_protozoa_first(self, ax):
        """绘制真核分类Z-score分布（原生动物在最前面，保持原始风格）"""
        categories_ordered = self.eukaryotic_categories.copy()
        
        categories = []
        zscore_data = []
        
        for category in categories_ordered:
            cat_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'] == category]
            z_scores = cat_data['robust_z_score'].dropna()
            if len(z_scores) >= 3:
                categories.append(category)
                zscore_data.append(z_scores)
        
        if len(categories) >= 2:
            x_min = min([data.min() for data in zscore_data]) - 1
            x_max = max([data.max() for data in zscore_data]) + 1
            x_range = np.linspace(x_min, x_max, 200)
            
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.fill_between(x_range, y_values, alpha=0.2, color=color)
            
            for i, (category, data) in enumerate(zip(categories, zscore_data)):
                if len(data) >= 3:
                    kde = gaussian_kde(data)
                    y_values = kde(x_range)
                    color = nine_category_colors[category]
                    ax.plot(x_range, y_values, color=color, linewidth=2, label=category)
            
            ax.axvline(x=self.z_threshold, color='red', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Enriched (Z>{self.z_threshold})')
            ax.axvline(x=-self.z_threshold, color='blue', linestyle='--', alpha=0.6, 
                      linewidth=1.2, label=f'Depleted (Z<-{self.z_threshold})')
            
            ax.set_xlabel('Robust Z-score')
            ax.set_ylabel('Probability Density')
            
            ax.legend(fontsize=8, loc='upper right', ncol=2)
            ax.grid(True, alpha=0.3, linestyle=':')
        else:
            ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center', transform=ax.transAxes)
    
    def perform_statistical_analysis(self):
        """执行统计分析"""
        results = {}
        
        
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
        
        return results
    
    def export_results(self, output_dir):
        """导出结果"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        
        results_file_3cat = output_dir / "zscore_data_3categories.csv"
        self.real_data_3cat[['Genome', 'Classification_3cat', 'GC_content', 
                           'Genomic_density', 'expected_density', 'robust_z_score']].to_csv(
                               results_file_3cat, index=False)
        
        eukaryotic_data = self.real_data_9cat[self.real_data_9cat['Classification_9cat'].isin(self.eukaryotic_categories)]
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
        
        stats_results = self.perform_statistical_analysis()
        if stats_results:
            stats_file = output_dir / "zscore_statistics.json"
            with open(stats_file, 'w') as f:
                json.dump(stats_results, f, indent=2, ensure_ascii=False)
        
        return output_dir
def main():
    """主函数"""
    set_clean_style()
    
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"
    
    
    analyzer = ZScoreDistributionBD(simulation_file, real_genome_file, z_threshold=1.96)
    
    try:
        analyzer.load_and_preprocess_data()
        analyzer.build_background_model()
        analyzer.calculate_deviations()
        
        fig = analyzer.create_figure_2_BD()
        
        output_dir = Path(base_dir) / "figure2_BD_analysis"
        analyzer.export_results(output_dir)
        
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.png", dpi=300, bbox_inches='tight')
        fig.savefig(output_dir / "figure2_zscore_distribution_BD.pdf", bbox_inches='tight')
        
        
        stats = analyzer.perform_statistical_analysis()
        
        for category in analyzer.three_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
        
        for category in analyzer.eukaryotic_categories:
            if f'{category}_zscore' in stats:
                cat_stats = stats[f'{category}_zscore']
        
        
    except Exception as e:
        import traceback
        traceback.print_exc()
if __name__ == "__main__":
    main()
