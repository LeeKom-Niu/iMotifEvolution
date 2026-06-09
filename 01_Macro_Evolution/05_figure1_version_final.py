"""
figure1_two_versions_final.py - 最终版本：双图对比的两个版本
版本1: 完整信息版（含n值）
版本2: 发表纯净版（简洁，无n值）
理论背景曲线描述更准确
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
def set_font_settings():
    """设置字体，兼容Linux环境"""
    import matplotlib
    available_fonts = matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    
    preferred_fonts = ['DejaVu Sans', 'Liberation Sans', 'FreeSans', 'Nimbus Sans']
    
    selected_font = None
    for font_name in preferred_fonts:
        for font_path in available_fonts:
            if font_name.lower() in font_path.lower():
                selected_font = font_name
                break
        if selected_font:
            break
    
    if not selected_font:
        for font_path in available_fonts:
            font_prop = matplotlib.font_manager.FontProperties(fname=font_path)
            font_name = font_prop.get_name()
            if 'sans' in font_name.lower():
                selected_font = font_name
                break
    
    if not selected_font:
        selected_font = 'sans-serif'
    
    
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': [selected_font],
        'font.size': 8,
        'pdf.fonttype': 42,
        
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'axes.linewidth': 0.6,
        'axes.unicode_minus': False,
        'axes.labelweight': 'normal',
        
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        
        'legend.fontsize': 7,
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '
        'legend.fancybox': False,
        
        'lines.linewidth': 1.8,
        'lines.markersize': 4.0,
        'lines.markeredgewidth': 0.5,
        
        'scatter.marker': 'o',
        'scatter.edgecolors': 'white',
        
        'figure.dpi': 300,
        'figure.constrained_layout.use': True,
        'savefig.dpi': 600,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        
        'grid.alpha': 0.2,
        'grid.linestyle': ':',
        'grid.linewidth': 0.5,
    })
two_category_colors = {
    'Prokaryote': '
    'Eukaryota': '
}
eukaryotic_subgroup_colors = {
    'Protozoa': '
    'Fungi': '
    'Plant': '
    'Invertebrate': '
    'Vertebrate Other': '
    'Mammalian': '
}
BACKGROUND_COLOR = '
BACKGROUND_FILL_COLOR = '
class FinalDensityPlotter:
    def __init__(self, simulation_file, real_genome_file):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.background_data = None
        self.real_data = None
        self.smooth_background = None
        
        self.two_categories = ['Prokaryote', 'Eukaryota']
        self.eukaryotic_subgroups = [
            'Protozoa', 'Fungi', 'Plant', 
            'Invertebrate', 'Vertebrate Other', 'Mammalian'
        ]
        
        self.stats = {}
        
        self.bg_gc_min = None
        self.bg_gc_max = None
    
    def load_and_preprocess_data(self):
        """加载并预处理所有数据"""
        
        try:
            real_data = pd.read_csv(self.real_genome_file, sep='\t')
        except Exception as e:
            return False
        
        column_mapping = {}
        for col in real_data.columns:
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
            if old_name in real_data.columns:
                real_data = real_data.rename(columns={old_name: new_name})
        
        if 'Classification' in real_data.columns:
            real_data = real_data[real_data['Classification'] != 'Classification']
            real_data['Classification'] = real_data['Classification'].str.strip()
            real_data['Classification'] = real_data['Classification'].replace({
                'Vertebrate Mammalian': 'Mammalian',
                'Vertebrate_Other': 'Vertebrate Other'
            })
        
        real_data['GC_content'] = pd.to_numeric(real_data['GC_content'], errors='coerce')
        real_data['Genomic_density'] = pd.to_numeric(real_data['Genomic_density'], errors='coerce')
        
        if real_data['GC_content'].max() > 1:
            real_data['GC_content'] = real_data['GC_content'] / 100.0
        
        original_len = len(real_data)
        real_data = real_data.dropna(subset=['GC_content', 'Genomic_density', 'Classification'])
        
        prokaryotic_mask = real_data['Classification'].isin(['Bacteria', 'Archaea'])
        gc_zero_mask = (real_data['GC_content'] == 0) & prokaryotic_mask
        if gc_zero_mask.any():
            real_data = real_data[~gc_zero_mask]
        
        self._create_classification_columns(real_data)
        
        self.real_data = real_data
        
        
        self._calculate_adaptive_background_range()
        
        return True
    
    def _create_classification_columns(self, data):
        """创建分类列"""
        data['Classification_2cat'] = np.nan
        
        prokaryotic_mask = data['Classification'].isin(['Bacteria', 'Archaea'])
        data.loc[prokaryotic_mask, 'Classification_2cat'] = 'Prokaryote'
        
        eukaryotic_categories = ['Protozoa', 'Fungi', 'Plant', 
                               'Invertebrate', 'Vertebrate Other', 'Mammalian']
        eukaryotic_mask = data['Classification'].isin(eukaryotic_categories)
        data.loc[eukaryotic_mask, 'Classification_2cat'] = 'Eukaryota'
        
        data['Classification_eukaryote'] = np.nan
        data.loc[eukaryotic_mask, 'Classification_eukaryote'] = data.loc[eukaryotic_mask, 'Classification']
        
        for cat_type in ['2cat', 'eukaryote']:
            col = f'Classification_{cat_type}'
            if col in data.columns:
                counts = data[col].value_counts()
                self.stats[cat_type] = counts.to_dict()
                for cat, count in counts.items():
    
    def _calculate_adaptive_background_range(self):
        """根据真实数据计算自适应背景曲线范围"""
        if self.real_data is None:
            return
        
        all_gc = self.real_data['GC_content'].values
        gc_min = np.percentile(all_gc, 2)
        gc_max = np.percentile(all_gc, 98)
        
        gc_min = max(0.0, gc_min - 0.05)
        gc_max = min(0.9, gc_max + 0.05)
        
        self.bg_gc_min = gc_min
        self.bg_gc_max = gc_max
        
    
    def load_theoretical_background(self):
        """加载随机模拟背景曲线"""
        
        try:
            sim_data = pd.read_csv(self.simulation_file)
            
            gc_col = None
            density_col = None
            
            for col in sim_data.columns:
                col_lower = col.lower()
                if 'gc' in col_lower:
                    gc_col = col
                if 'density' in col_lower or ('im' in col_lower and 'mb' in col_lower):
                    density_col = col
            
            if not gc_col or not density_col:
                gc_col = sim_data.columns[0]
                density_col = sim_data.columns[1]
            
            bg_data = sim_data[[gc_col, density_col]].copy()
            bg_data.columns = ['GC_content', 'density']
        except Exception as e:
            return False
        
        if bg_data['GC_content'].max() > 1:
            bg_data['GC_content'] = bg_data['GC_content'] / 100.0
        
        bg_data = bg_data.groupby('GC_content', as_index=False)['density'].mean()
        
        bg_data = bg_data.sort_values('GC_content').dropna().reset_index(drop=True)
        
        
        self.background_data = bg_data
        return True
    
    def create_adaptive_smooth_curve(self):
        """创建自适应的平滑随机模拟背景曲线"""
        
        if self.background_data is None:
            return False
        
        gc_values = self.background_data['GC_content'].values
        density_values = self.background_data['density'].values
        
        unique_indices = np.unique(gc_values, return_index=True)[1]
        gc_values = gc_values[unique_indices]
        density_values = density_values[unique_indices]
        
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]
        
        
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            extended_min = max(0.0, self.bg_gc_min - 0.05)
            extended_max = min(0.9, self.bg_gc_max + 0.05)
            
            mask = (gc_values >= extended_min) & (gc_values <= extended_max)
            if mask.any():
                gc_values = gc_values[mask]
                density_values = density_values[mask]
            else:
        
        if len(gc_values) < 3:
            gc_values = self.background_data['GC_content'].values
            density_values = self.background_data['density'].values
            gc_values = gc_values[:]
            density_values = density_values[:]
        
        if self.bg_gc_min is not None and gc_values.min() > self.bg_gc_min:
            low_gc = np.array([self.bg_gc_min])
            low_density = np.interp(low_gc, gc_values, density_values)
            gc_values = np.concatenate([low_gc, gc_values])
            density_values = np.concatenate([low_density, density_values])
        
        if self.bg_gc_max is not None and gc_values.max() < self.bg_gc_max:
            high_gc = np.array([self.bg_gc_max])
            high_density = np.interp(high_gc, gc_values, density_values)
            gc_values = np.concatenate([gc_values, high_gc])
            density_values = np.concatenate([density_values, high_density])
        
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            smooth_gc = np.linspace(self.bg_gc_min, self.bg_gc_max, 300)
        else:
            smooth_gc = np.linspace(gc_values.min(), gc_values.max(), 300)
        
        
        try:
            interp_func = interp1d(gc_values, density_values, kind='cubic', 
                                 fill_value='extrapolate', bounds_error=False)
            smooth_density = interp_func(smooth_gc)
        except Exception as e:
            try:
                interp_func = interp1d(gc_values, density_values, kind='quadratic', 
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
            except:
                interp_func = interp1d(gc_values, density_values, kind='linear', 
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
        
        smooth_density = np.maximum(smooth_density, 0)
        
        smooth_density = gaussian_filter1d(smooth_density, sigma=1.5)
        
        
        self.smooth_background = pd.DataFrame({
            'GC_content': smooth_gc,
            'density': smooth_density
        })
        
        return True
    
    def create_complete_version_figure(self):
        """创建完整信息版图形（含n值）"""
        
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        
        self._plot_prok_euk_comparison(axes[0], show_n=True)
        
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=True)
        
        
        return fig
    
    def create_clean_version_figure(self):
        """创建发表纯净版图形（无n值）"""
        
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        
        self._plot_prok_euk_comparison(axes[0], show_n=False)
        
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=False)
        
        
        return fig
    
    def _plot_prok_euk_comparison(self, ax, show_n=True):
        """绘制原核与真核对比图"""
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            ax.fill_between(bg_gc, 0, bg_density, 
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')
            
            ax.plot(bg_gc, bg_density, 
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)
        
        if self.real_data is not None and 'Classification_2cat' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_2cat'])
            
            for category in self.two_categories:
                cat_data = plot_data[plot_data['Classification_2cat'] == category]
                if len(cat_data) > 0:
                    color = two_category_colors[category]
                    
                    label = f"{category}"
                    if show_n and '2cat' in self.stats and category in self.stats['2cat']:
                        label += f" (n={self.stats['2cat'][category]:,})"
                    
                    ax.scatter(
                        cat_data['GC_content'], 
                        cat_data['Genomic_density'],
                        color=color, 
                        alpha=0.75,
                        s=12,
                        label=label,
                        edgecolors='white',
                        linewidth=0.3,
                        zorder=10
                    )
        
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')
        
        self._set_adaptive_axis_limits(ax, '2cat')
        
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)
        
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)
        
        self._add_legend_to_axis(ax, show_n=show_n)
    
    def _plot_eukaryote_subgroup_comparison(self, ax, show_n=True):
        """绘制真核细分对比图"""
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            ax.fill_between(bg_gc, 0, bg_density, 
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')
            
            ax.plot(bg_gc, bg_density, 
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)
        
        if self.real_data is not None and 'Classification_eukaryote' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_eukaryote'])
            
            for category in self.eukaryotic_subgroups:
                cat_data = plot_data[plot_data['Classification_eukaryote'] == category]
                if len(cat_data) > 0:
                    color = eukaryotic_subgroup_colors[category]
                    
                    label = f"{category}"
                    if show_n and 'eukaryote' in self.stats and category in self.stats['eukaryote']:
                        label += f" (n={self.stats['eukaryote'][category]:,})"
                    
                    ax.scatter(
                        cat_data['GC_content'], 
                        cat_data['Genomic_density'],
                        color=color, 
                        alpha=0.8,
                        s=14,
                        label=label,
                        edgecolors='white',
                        linewidth=0.3,
                        zorder=10
                    )
        
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')
        
        self._set_adaptive_axis_limits(ax, 'eukaryote')
        
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)
        
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)
        
        self._add_eukaryotic_legend_to_axis(ax, show_n=show_n)
    
    def _set_adaptive_axis_limits(self, ax, cat_type):
        """自适应设置坐标轴范围"""
        if self.real_data is None:
            return
        
        col = f'Classification_{cat_type}'
        if col not in self.real_data.columns:
            return
        
        plot_data = self.real_data.dropna(subset=[col])
        if len(plot_data) == 0:
            return
        
        x_data = plot_data['GC_content'].values
        
        x_min = np.percentile(x_data, 2) - 0.05
        x_max = np.percentile(x_data, 98) + 0.05
        
        x_min = max(-0.02, x_min)
        x_max = min(0.92, x_max)
        
        ax.set_xlim(x_min, x_max)
        
        y_data = plot_data['Genomic_density'].values
        
        y_99 = np.percentile(y_data, 99)
        
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            mask = (bg_gc >= x_min) & (bg_gc <= x_max)
            if mask.any():
                bg_max_in_range = bg_density[mask].max()
            else:
                bg_max_in_range = bg_density.max()
            
            y_max = max(y_99, bg_max_in_range) * 1.15
        else:
            y_max = y_99 * 1.15
        
        ax.set_ylim(-0.02 * y_max, y_max)
        
    
    def _add_legend_to_axis(self, ax, show_n=True):
        """添加图例到坐标轴（用于原核与真核图）"""
        handles, labels = ax.get_legend_handles_labels()
        
        if len(handles) == 0:
            return
        
        bg_handles = []
        bg_labels = []
        data_handles = []
        data_labels = []
        
        for handle, label in zip(handles, labels):
            if 'background' in label.lower():
                bg_handles.append(handle)
                bg_labels.append(label)
            else:
                data_handles.append(handle)
                data_labels.append(label)
        
        if data_handles:
            leg1 = ax.legend(data_handles, data_labels, 
                           loc='upper left', fontsize=7.5,
                           handletextpad=0.5, 
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=1)
            
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR, 
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])
                
                ax.legend([bg_patch], [bg_labels[0]], 
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)
                
                ax.add_artist(leg1)
    
    def _add_eukaryotic_legend_to_axis(self, ax, show_n=True):
        """添加图例到坐标轴（用于真核细分图）- 修正位置：左上角"""
        handles, labels = ax.get_legend_handles_labels()
        
        if len(handles) == 0:
            return
        
        bg_handles = []
        bg_labels = []
        data_handles = []
        data_labels = []
        
        for handle, label in zip(handles, labels):
            if 'background' in label.lower():
                bg_handles.append(handle)
                bg_labels.append(label)
            else:
                data_handles.append(handle)
                data_labels.append(label)
        
        if data_handles:
            leg1 = ax.legend(data_handles, data_labels, 
                           loc='upper left', fontsize=7,
                           handletextpad=0.5, 
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=2)
            
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR, 
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])
                
                ax.legend([bg_patch], [bg_labels[0]], 
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)
                
                ax.add_artist(leg1)
    
    def save_both_versions(self, output_dir):
        """保存两个版本的图形"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        
        saved_files = []
        
        fig_complete = self.create_complete_version_figure()
        saved_complete = self._save_single_figure(fig_complete, output_dir, "figure1_complete_version")
        saved_files.extend(saved_complete)
        
        fig_clean = self.create_clean_version_figure()
        saved_clean = self._save_single_figure(fig_clean, output_dir, "figure1_clean_version")
        saved_files.extend(saved_clean)
        
        plt.close('all')
        
        return saved_files
    
    def _save_single_figure(self, fig, output_dir, filename):
        """保存单个图形"""
        formats = [
            ('tiff', 'tiff', {'pil_kwargs': {'compression': 'tiff_lzw'}}),
            ('pdf', 'pdf', {}),
            ('png', 'png', {}),
            ('svg', 'svg', {}),
        ]
        
        saved_files = []
        for fmt_name, fmt_ext, kwargs in formats:
            file_path = output_dir / f"{filename}.{fmt_ext}"
            try:
                fig.savefig(file_path, dpi=600, format=fmt_ext, **kwargs)
                saved_files.append(file_path)
            except Exception as e:
        
        return saved_files
    
    def export_statistics(self, output_dir):
        """导出统计数据"""
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        stats_data = []
        
        for cat_type in ['2cat', 'eukaryote']:
            col = f'Classification_{cat_type}'
            if self.real_data is not None and col in self.real_data.columns:
                cat_data = self.real_data.dropna(subset=[col])
                
                for category in cat_data[col].unique():
                    sub_data = cat_data[cat_data[col] == category]
                    
                    stats_row = {
                        'Category_Type': cat_type,
                        'Category': category,
                        'Count': len(sub_data),
                        'GC_Mean': sub_data['GC_content'].mean(),
                        'GC_Std': sub_data['GC_content'].std(),
                        'GC_Min': sub_data['GC_content'].min(),
                        'GC_Max': sub_data['GC_content'].max(),
                        'Density_Mean': sub_data['Genomic_density'].mean(),
                        'Density_Std': sub_data['Genomic_density'].std(),
                        'Density_Min': sub_data['Genomic_density'].min(),
                        'Density_Max': sub_data['Genomic_density'].max(),
                        'Density_Median': sub_data['Genomic_density'].median()
                    }
                    stats_data.append(stats_row)
        
        if stats_data:
            stats_df = pd.DataFrame(stats_data)
            stats_path = output_dir / "figure1_statistics.csv"
            stats_df.to_csv(stats_path, index=False)
            
            for cat_type in ['2cat', 'eukaryote']:
                cat_stats = stats_df[stats_df['Category_Type'] == cat_type]
                if not cat_stats.empty:
                    for _, row in cat_stats.iterrows():
                              f"GC={row['GC_Mean']:.3f}±{row['GC_Std']:.3f} "
                              f"({row['GC_Min']:.3f}-{row['GC_Max']:.3f}), "
                              f"Density={row['Density_Mean']:.2f}±{row['Density_Std']:.2f} IM/Mb")
def main():
    """主函数"""
    set_font_settings()
    
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"
    
    
    plotter = FinalDensityPlotter(simulation_file, real_genome_file)
    
    try:
        if not plotter.load_and_preprocess_data():
            return
        
        if not plotter.load_theoretical_background():
            return
        
        if not plotter.create_adaptive_smooth_curve():
            return
        
        output_dir = Path(base_dir) / "figures" / "figure1_two_versions"
        plotter.export_statistics(output_dir)
        
        saved_files = plotter.save_both_versions(output_dir)
        
        
        for file_path in saved_files:
            if "complete_version" in str(file_path):
        
        for file_path in saved_files:
            if "clean_version" in str(file_path):
        
        
        
    except Exception as e:
        import traceback
        traceback.print_exc()
if __name__ == "__main__":
    main()
