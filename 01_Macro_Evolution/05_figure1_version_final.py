#!/usr/bin/env python3
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

# 设置字体 - 在Linux上使用可用的字体
def set_font_settings():
    """设置字体，兼容Linux环境"""
    import matplotlib
    # 获取可用字体
    available_fonts = matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf')
    
    # 优先选择DejaVu Sans，它是Linux上通常可用的优秀字体
    preferred_fonts = ['DejaVu Sans', 'Liberation Sans', 'FreeSans', 'Nimbus Sans']
    
    selected_font = None
    for font_name in preferred_fonts:
        # 检查是否有这些字体
        for font_path in available_fonts:
            if font_name.lower() in font_path.lower():
                selected_font = font_name
                break
        if selected_font:
            break
    
    if not selected_font:
        # 如果没有找到首选字体，使用第一个可用的无衬线字体
        for font_path in available_fonts:
            font_prop = matplotlib.font_manager.FontProperties(fname=font_path)
            font_name = font_prop.get_name()
            if 'sans' in font_name.lower():
                selected_font = font_name
                break
    
    if not selected_font:
        selected_font = 'sans-serif'  # 使用系统默认
    
    print(f"使用字体: {selected_font}")
    
    # 设置Nature风格，使用找到的字体
    plt.rcParams.update({
        # 字体设置
        'font.family': 'sans-serif',
        'font.sans-serif': [selected_font],
        'font.size': 8,
        'pdf.fonttype': 42,
        
        # 坐标轴
        'axes.labelsize': 9,
        'axes.titlesize': 10,
        'axes.linewidth': 0.6,
        'axes.unicode_minus': False,
        'axes.labelweight': 'normal',
        
        # 刻度
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        
        # 图例
        'legend.fontsize': 7,
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '#CCCCCC',
        'legend.fancybox': False,
        
        # 线条
        'lines.linewidth': 1.8,
        'lines.markersize': 4.0,
        'lines.markeredgewidth': 0.5,
        
        # 散点
        'scatter.marker': 'o',
        'scatter.edgecolors': 'white',
        
        # 图形
        'figure.dpi': 300,
        'figure.constrained_layout.use': True,
        'savefig.dpi': 600,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        
        # 网格
        'grid.alpha': 0.2,
        'grid.linestyle': ':',
        'grid.linewidth': 0.5,
    })

# 优化颜色方案
two_category_colors = {
    'Prokaryote': '#1F77B4',     # 深蓝色 - 原核生物
    'Eukaryota': '#D62728',      # 红色 - 真核生物
}

# 真核细分颜色方案
eukaryotic_subgroup_colors = {
    'Protozoa': '#BCBD22',       # 橄榄绿 - 原生动物
    'Fungi': '#9467BD',          # 紫色 - 真菌
    'Plant': '#8C564B',          # 棕色 - 植物
    'Invertebrate': '#E377C2',   # 粉色 - 无脊椎动物
    'Vertebrate Other': '#7F7F7F',  # 灰色 - 脊椎动物其他
    'Mammalian': '#D62728',      # 红色 - 哺乳动物
}

# 理论背景曲线颜色 - 使用更中性的颜色
BACKGROUND_COLOR = '#2CA02C'      # 柔和的绿色，与数据点区分
BACKGROUND_FILL_COLOR = '#C7E9C0' # 淡绿色填充

class FinalDensityPlotter:
    def __init__(self, simulation_file, real_genome_file):
        self.simulation_file = simulation_file
        self.real_genome_file = real_genome_file
        self.background_data = None
        self.real_data = None
        self.smooth_background = None
        
        # 分类定义
        self.two_categories = ['Prokaryote', 'Eukaryota']
        # 真核类群
        self.eukaryotic_subgroups = [
            'Protozoa', 'Fungi', 'Plant', 
            'Invertebrate', 'Vertebrate Other', 'Mammalian'
        ]
        
        # 统计信息
        self.stats = {}
        
        # 自适应背景曲线范围
        self.bg_gc_min = None
        self.bg_gc_max = None
    
    def load_and_preprocess_data(self):
        """加载并预处理所有数据"""
        print("=== 加载数据 ===")
        
        # 1. 加载真实基因组数据
        print(f"1. 加载真实基因组数据: {self.real_genome_file}")
        try:
            real_data = pd.read_csv(self.real_genome_file, sep='\t')
            print(f"   ✓ 成功加载 {len(real_data)} 行数据")
        except Exception as e:
            print(f"   ✗ 加载失败: {e}")
            return False
        
        # 标准化列名
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
        
        # 清理数据
        if 'Classification' in real_data.columns:
            real_data = real_data[real_data['Classification'] != 'Classification']
            real_data['Classification'] = real_data['Classification'].str.strip()
            real_data['Classification'] = real_data['Classification'].replace({
                'Vertebrate Mammalian': 'Mammalian',
                'Vertebrate_Other': 'Vertebrate Other'
            })
        
        # 转换数据类型
        real_data['GC_content'] = pd.to_numeric(real_data['GC_content'], errors='coerce')
        real_data['Genomic_density'] = pd.to_numeric(real_data['Genomic_density'], errors='coerce')
        
        # GC含量百分比转小数
        if real_data['GC_content'].max() > 1:
            real_data['GC_content'] = real_data['GC_content'] / 100.0
            print("   ✓ 将GC含量从百分比转换为小数")
        
        # 移除无效数据
        original_len = len(real_data)
        real_data = real_data.dropna(subset=['GC_content', 'Genomic_density', 'Classification'])
        print(f"   ✓ 数据清洗: {original_len} → {len(real_data)} 行")
        
        # 过滤掉GC含量为0的原核数据点
        prokaryotic_mask = real_data['Classification'].isin(['Bacteria', 'Archaea'])
        gc_zero_mask = (real_data['GC_content'] == 0) & prokaryotic_mask
        if gc_zero_mask.any():
            print(f"   ⚠ 移除 {gc_zero_mask.sum()} 个GC含量为0的原核数据点")
            real_data = real_data[~gc_zero_mask]
        
        # 创建分类列
        self._create_classification_columns(real_data)
        
        self.real_data = real_data
        
        print(f"\n2. 数据统计:")
        print(f"   - 总基因组数: {len(real_data):,}")
        print(f"   - GC含量范围: {real_data['GC_content'].min():.3f} - {real_data['GC_content'].max():.3f}")
        print(f"   - 密度范围: {real_data['Genomic_density'].min():.3f} - {real_data['Genomic_density'].max():.3f}")
        
        # 计算自适应背景曲线范围
        self._calculate_adaptive_background_range()
        
        return True
    
    def _create_classification_columns(self, data):
        """创建分类列"""
        # 两分类（原核和真核）
        data['Classification_2cat'] = np.nan
        
        # 原核：细菌 + 古菌
        prokaryotic_mask = data['Classification'].isin(['Bacteria', 'Archaea'])
        data.loc[prokaryotic_mask, 'Classification_2cat'] = 'Prokaryote'
        
        # 真核
        eukaryotic_categories = ['Protozoa', 'Fungi', 'Plant', 
                               'Invertebrate', 'Vertebrate Other', 'Mammalian']
        eukaryotic_mask = data['Classification'].isin(eukaryotic_categories)
        data.loc[eukaryotic_mask, 'Classification_2cat'] = 'Eukaryota'
        
        # 真核细分
        data['Classification_eukaryote'] = np.nan
        data.loc[eukaryotic_mask, 'Classification_eukaryote'] = data.loc[eukaryotic_mask, 'Classification']
        
        # 统计各分类数量
        for cat_type in ['2cat', 'eukaryote']:
            col = f'Classification_{cat_type}'
            if col in data.columns:
                counts = data[col].value_counts()
                self.stats[cat_type] = counts.to_dict()
                print(f"   - {cat_type}分类: {len(counts)} 类")
                for cat, count in counts.items():
                    print(f"     {cat}: {count:,}")
    
    def _calculate_adaptive_background_range(self):
        """根据真实数据计算自适应背景曲线范围"""
        if self.real_data is None:
            return
        
        # 获取所有数据的GC范围
        all_gc = self.real_data['GC_content'].values
        gc_min = np.percentile(all_gc, 2)  # 2%分位数作为下限
        gc_max = np.percentile(all_gc, 98)  # 98%分位数作为上限
        
        # 确保合理的范围
        gc_min = max(0.0, gc_min - 0.05)  # 扩大一点范围
        gc_max = min(0.9, gc_max + 0.05)  # 限制在0.9以内
        
        self.bg_gc_min = gc_min
        self.bg_gc_max = gc_max
        
        print(f"\n3. 自适应背景曲线范围:")
        print(f"   - GC范围: {self.bg_gc_min:.3f} - {self.bg_gc_max:.3f}")
    
    def load_theoretical_background(self):
        """加载随机模拟背景曲线"""
        print("\n=== 加载随机模拟背景曲线 ===")
        
        print(f"1. 加载模拟数据: {self.simulation_file}")
        try:
            sim_data = pd.read_csv(self.simulation_file)
            
            # 查找GC和密度列
            gc_col = None
            density_col = None
            
            for col in sim_data.columns:
                col_lower = col.lower()
                if 'gc' in col_lower:
                    gc_col = col
                if 'density' in col_lower or ('im' in col_lower and 'mb' in col_lower):
                    density_col = col
            
            if not gc_col or not density_col:
                # 如果没找到，使用前两列
                gc_col = sim_data.columns[0]
                density_col = sim_data.columns[1]
                print(f"   ⚠ 使用默认列: GC={gc_col}, Density={density_col}")
            
            bg_data = sim_data[[gc_col, density_col]].copy()
            bg_data.columns = ['GC_content', 'density']
            print(f"   ✓ 从模拟数据创建背景曲线 ({len(bg_data)} 行)")
        except Exception as e:
            print(f"   ✗ 加载模拟数据失败: {e}")
            return False
        
        # 确保GC含量为小数
        if bg_data['GC_content'].max() > 1:
            bg_data['GC_content'] = bg_data['GC_content'] / 100.0
            print("   ✓ 将GC含量从百分比转换为小数")
        
        # 移除重复的GC值（取平均值）
        print(f"   - 处理前: {len(bg_data)} 行")
        bg_data = bg_data.groupby('GC_content', as_index=False)['density'].mean()
        print(f"   - 移除重复GC值后: {len(bg_data)} 行")
        
        # 确保数据排序
        bg_data = bg_data.sort_values('GC_content').dropna().reset_index(drop=True)
        
        print(f"\n2. 随机模拟背景数据:")
        print(f"   - 数据点数: {len(bg_data)}")
        print(f"   - GC范围: {bg_data['GC_content'].min():.3f} - {bg_data['GC_content'].max():.3f}")
        print(f"   - 密度范围: {bg_data['density'].min():.3f} - {bg_data['density'].max():.3f}")
        
        self.background_data = bg_data
        return True
    
    def create_adaptive_smooth_curve(self):
        """创建自适应的平滑随机模拟背景曲线"""
        print("\n=== 创建自适应平滑随机模拟背景曲线 ===")
        
        if self.background_data is None:
            print("   ✗ 没有背景数据")
            return False
        
        # 提取数据
        gc_values = self.background_data['GC_content'].values
        density_values = self.background_data['density'].values
        
        # 确保严格递增且无重复
        unique_indices = np.unique(gc_values, return_index=True)[1]
        gc_values = gc_values[unique_indices]
        density_values = density_values[unique_indices]
        
        # 确保排序
        sort_idx = np.argsort(gc_values)
        gc_values = gc_values[sort_idx]
        density_values = density_values[sort_idx]
        
        print(f"1. 原始数据:")
        print(f"   - 点数: {len(gc_values)}")
        print(f"   - GC范围: {gc_values.min():.3f} - {gc_values.max():.3f}")
        
        # 根据自适应范围筛选数据
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            # 扩展一点范围以获得更好的插值
            extended_min = max(0.0, self.bg_gc_min - 0.05)
            extended_max = min(0.9, self.bg_gc_max + 0.05)
            
            mask = (gc_values >= extended_min) & (gc_values <= extended_max)
            if mask.any():
                gc_values = gc_values[mask]
                density_values = density_values[mask]
                print(f"   - 自适应筛选后: {len(gc_values)} 点")
                print(f"   - 筛选范围: {extended_min:.3f} - {extended_max:.3f}")
            else:
                print(f"   ⚠ 无数据在自适应范围内，使用原始数据")
        
        # 确保有足够的数据点
        if len(gc_values) < 3:
            print(f"   ⚠ 数据点不足({len(gc_values)})，使用原始范围")
            gc_values = self.background_data['GC_content'].values
            density_values = self.background_data['density'].values
            gc_values = gc_values[:]
            density_values = density_values[:]
        
        # 扩展数据到自适应范围的边界
        if self.bg_gc_min is not None and gc_values.min() > self.bg_gc_min:
            # 在低端添加点
            low_gc = np.array([self.bg_gc_min])
            low_density = np.interp(low_gc, gc_values, density_values)
            gc_values = np.concatenate([low_gc, gc_values])
            density_values = np.concatenate([low_density, density_values])
            print(f"   ✓ 扩展低GC端到 {self.bg_gc_min:.3f}")
        
        if self.bg_gc_max is not None and gc_values.max() < self.bg_gc_max:
            # 在高端添加点
            high_gc = np.array([self.bg_gc_max])
            high_density = np.interp(high_gc, gc_values, density_values)
            gc_values = np.concatenate([gc_values, high_gc])
            density_values = np.concatenate([density_values, high_density])
            print(f"   ✓ 扩展高GC端到 {self.bg_gc_max:.3f}")
        
        # 创建自适应采样点
        if self.bg_gc_min is not None and self.bg_gc_max is not None:
            smooth_gc = np.linspace(self.bg_gc_min, self.bg_gc_max, 300)
        else:
            smooth_gc = np.linspace(gc_values.min(), gc_values.max(), 300)
        
        print(f"2. 创建自适应平滑曲线:")
        print(f"   - 采样点数: {len(smooth_gc)}")
        print(f"   - 最终GC范围: {smooth_gc.min():.3f} - {smooth_gc.max():.3f}")
        
        try:
            # 使用三次样条插值
            interp_func = interp1d(gc_values, density_values, kind='cubic', 
                                 fill_value='extrapolate', bounds_error=False)
            smooth_density = interp_func(smooth_gc)
            print("   ✓ 使用三次样条插值")
        except Exception as e:
            print(f"   ✗ 三次样条失败: {e}")
            try:
                # 尝试二次样条
                interp_func = interp1d(gc_values, density_values, kind='quadratic', 
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
                print("   ✓ 使用二次样条插值")
            except:
                # 使用线性插值
                interp_func = interp1d(gc_values, density_values, kind='linear', 
                                     fill_value='extrapolate', bounds_error=False)
                smooth_density = interp_func(smooth_gc)
                print("   ✓ 使用线性插值")
        
        # 确保非负值
        smooth_density = np.maximum(smooth_density, 0)
        
        # 轻微平滑处理
        smooth_density = gaussian_filter1d(smooth_density, sigma=1.5)
        
        print(f"3. 平滑曲线统计:")
        print(f"   - 密度范围: {smooth_density.min():.3f} - {smooth_density.max():.3f}")
        print(f"   - 平均密度: {smooth_density.mean():.3f}")
        
        # 保存平滑曲线
        self.smooth_background = pd.DataFrame({
            'GC_content': smooth_gc,
            'density': smooth_density
        })
        
        return True
    
    def create_complete_version_figure(self):
        """创建完整信息版图形（含n值）"""
        print("\n=== 创建完整信息版Figure ===")
        
        # 创建图形：两张图并排
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        
        # 1. 原核与真核对比图
        print("1. 创建原核与真核对比图...")
        self._plot_prok_euk_comparison(axes[0], show_n=True)
        
        # 2. 真核细分对比图
        print("2. 创建真核细分对比图...")
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=True)
        
        # 添加面板标签 (A, B)
        #for i, ax in enumerate(axes):
        #    label = chr(65 + i)  # A, B
        #    ax.text(-0.12, 1.02, label, transform=ax.transAxes,
        #           fontsize=11, fontweight='bold', va='bottom', ha='right',
        #           fontfamily=plt.rcParams['font.sans-serif'][0])
        
        print("✓ 完整信息版Figure创建完成")
        return fig
    
    def create_clean_version_figure(self):
        """创建发表纯净版图形（无n值）"""
        print("\n=== 创建发表纯净版Figure ===")
        
        # 创建图形：两张图并排
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
        
        # 1. 原核与真核对比图
        print("1. 创建原核与真核对比图...")
        self._plot_prok_euk_comparison(axes[0], show_n=False)
        
        # 2. 真核细分对比图
        print("2. 创建真核细分对比图...")
        self._plot_eukaryote_subgroup_comparison(axes[1], show_n=False)
        
        # 添加面板标签 (a, b)
        #for i, ax in enumerate(axes):
        #    label = chr(65 + i)  # a, b
        #    ax.text(-0.12, 1.02, label, transform=ax.transAxes,
        #           fontsize=11, fontweight='bold', va='bottom', ha='right',
        #           fontfamily=plt.rcParams['font.sans-serif'][0])
        
        print("✓ 发表纯净版Figure创建完成")
        return fig
    
    def _plot_prok_euk_comparison(self, ax, show_n=True):
        """绘制原核与真核对比图"""
        # 绘制随机模拟背景曲线（自适应范围）
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            # 使用半透明填充创建柔和背景
            ax.fill_between(bg_gc, 0, bg_density, 
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')
            
            # 绘制背景曲线
            ax.plot(bg_gc, bg_density, 
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)
        
        # 绘制散点数据
        if self.real_data is not None and 'Classification_2cat' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_2cat'])
            
            for category in self.two_categories:
                cat_data = plot_data[plot_data['Classification_2cat'] == category]
                if len(cat_data) > 0:
                    color = two_category_colors[category]
                    
                    # 为图例创建标签
                    label = f"{category}"
                    if show_n and '2cat' in self.stats and category in self.stats['2cat']:
                        label += f" (n={self.stats['2cat'][category]:,})"
                    
                    # 使用稍大的点以便更好展示
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
        
        # 设置坐标轴
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')
       #ax.set_title('Prokaryotes vs eukaryotes', fontsize=11, fontweight='bold', pad=12)
        
        # 自适应设置坐标轴范围
        self._set_adaptive_axis_limits(ax, '2cat')
        
        # 设置刻度
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)
        
        # 添加网格
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)
        
        # 添加图例
        self._add_legend_to_axis(ax, show_n=show_n)
    
    def _plot_eukaryote_subgroup_comparison(self, ax, show_n=True):
        """绘制真核细分对比图"""
        # 绘制随机模拟背景曲线（自适应范围）
        if self.smooth_background is not None:
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            ax.fill_between(bg_gc, 0, bg_density, 
                          color=BACKGROUND_FILL_COLOR, alpha=0.2, zorder=1,
                          label='Random simulation background')
            
            ax.plot(bg_gc, bg_density, 
                   color=BACKGROUND_COLOR, linewidth=2.0, zorder=2,
                   linestyle='-', alpha=0.8)
        
        # 绘制散点数据
        if self.real_data is not None and 'Classification_eukaryote' in self.real_data.columns:
            plot_data = self.real_data.dropna(subset=['Classification_eukaryote'])
            
            for category in self.eukaryotic_subgroups:
                cat_data = plot_data[plot_data['Classification_eukaryote'] == category]
                if len(cat_data) > 0:
                    color = eukaryotic_subgroup_colors[category]
                    
                    # 为图例创建标签
                    label = f"{category}"
                    if show_n and 'eukaryote' in self.stats and category in self.stats['eukaryote']:
                        label += f" (n={self.stats['eukaryote'][category]:,})"
                    
                    ax.scatter(
                        cat_data['GC_content'], 
                        cat_data['Genomic_density'],
                        color=color, 
                        alpha=0.8,
                        s=14,  # 真核数据点更大以便区分
                        label=label,
                        edgecolors='white',
                        linewidth=0.3,
                        zorder=10
                    )
        
        # 设置坐标轴
        ax.set_xlabel('GC content', fontsize=10, fontweight='medium')
        ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10, fontweight='medium')
       #ax.set_title('Eukaryotic subgroups', fontsize=11, fontweight='bold', pad=12)
        
        # 自适应设置坐标轴范围
        self._set_adaptive_axis_limits(ax, 'eukaryote')
        
        # 设置刻度
        ax.set_xticks(np.arange(0, 1.0, 0.2))
        ax.tick_params(axis='both', which='major', labelsize=9)
        
        # 添加网格
        ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5, zorder=0)
        
        # 添加图例 - 修正位置：现在在左上角
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
        
        # X轴范围：基于数据分布
        x_data = plot_data['GC_content'].values
        
        # 使用数据的2%和98%分位数，并扩大一点
        x_min = np.percentile(x_data, 2) - 0.05
        x_max = np.percentile(x_data, 98) + 0.05
        
        # 确保合理范围
        x_min = max(-0.02, x_min)
        x_max = min(0.92, x_max)
        
        ax.set_xlim(x_min, x_max)
        
        # Y轴范围：基于数据和背景曲线
        y_data = plot_data['Genomic_density'].values
        
        # 使用数据的99%分位数作为参考
        y_99 = np.percentile(y_data, 99)
        
        # 如果有背景曲线，考虑其最大值
        if self.smooth_background is not None:
            # 只考虑在当前X范围内的背景曲线
            bg_gc = self.smooth_background['GC_content'].values
            bg_density = self.smooth_background['density'].values
            
            # 筛选在X范围内的背景数据
            mask = (bg_gc >= x_min) & (bg_gc <= x_max)
            if mask.any():
                bg_max_in_range = bg_density[mask].max()
            else:
                bg_max_in_range = bg_density.max()
            
            y_max = max(y_99, bg_max_in_range) * 1.15
        else:
            y_max = y_99 * 1.15
        
        # 确保最小值为0
        ax.set_ylim(-0.02 * y_max, y_max)
        
        print(f"   - {cat_type} 坐标轴范围: X={x_min:.3f}-{x_max:.3f}, Y=0-{y_max:.1f}")
    
    def _add_legend_to_axis(self, ax, show_n=True):
        """添加图例到坐标轴（用于原核与真核图）"""
        handles, labels = ax.get_legend_handles_labels()
        
        if len(handles) == 0:
            return
        
        # 分离背景和数据的图例
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
        
        # 添加数据图例
        if data_handles:
            # 原核与真核使用单列，放在左上角
            leg1 = ax.legend(data_handles, data_labels, 
                           loc='upper left', fontsize=7.5,
                           handletextpad=0.5, 
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=1)
            
            # 添加背景图例
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR, 
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])
                
                # 在右上角添加背景图例
                ax.legend([bg_patch], [bg_labels[0]], 
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)
                
                # 恢复第一个图例
                ax.add_artist(leg1)
    
    def _add_eukaryotic_legend_to_axis(self, ax, show_n=True):
        """添加图例到坐标轴（用于真核细分图）- 修正位置：左上角"""
        handles, labels = ax.get_legend_handles_labels()
        
        if len(handles) == 0:
            return
        
        # 分离背景和数据的图例
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
        
        # 添加数据图例
        if data_handles:
            # 真核细分使用两列，放在左上角
            leg1 = ax.legend(data_handles, data_labels, 
                           loc='upper left', fontsize=7,
                           handletextpad=0.5, 
                           borderaxespad=0.3,
                           framealpha=0.95,
                           ncol=2)
            
            # 添加背景图例
            if bg_handles:
                from matplotlib.patches import Patch
                bg_patch = Patch(facecolor=BACKGROUND_FILL_COLOR, 
                               edgecolor=BACKGROUND_COLOR,
                               linewidth=1,
                               alpha=0.7,
                               label=bg_labels[0])
                
                # 在右上角添加背景图例
                ax.legend([bg_patch], [bg_labels[0]], 
                        loc='upper right', fontsize=7.5,
                        handletextpad=0.5,
                        borderaxespad=0.3,
                        framealpha=0.95)
                
                # 恢复第一个图例
                ax.add_artist(leg1)
    
    def save_both_versions(self, output_dir):
        """保存两个版本的图形"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n=== 保存两个版本的图形 ===")
        
        saved_files = []
        
        # 1. 保存完整信息版
        print("\n1. 保存完整信息版...")
        fig_complete = self.create_complete_version_figure()
        saved_complete = self._save_single_figure(fig_complete, output_dir, "figure1_complete_version")
        saved_files.extend(saved_complete)
        
        # 2. 保存发表纯净版
        print("\n2. 保存发表纯净版...")
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
                print(f"   ✓ {fmt_name.upper()}格式: {file_path}")
                saved_files.append(file_path)
            except Exception as e:
                print(f"   ✗ 保存{fmt_name.upper()}失败: {e}")
        
        return saved_files
    
    def export_statistics(self, output_dir):
        """导出统计数据"""
        print("\n=== 导出统计数据 ===")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        stats_data = []
        
        # 收集所有统计数据
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
        
        # 保存为CSV
        if stats_data:
            stats_df = pd.DataFrame(stats_data)
            stats_path = output_dir / "figure1_statistics.csv"
            stats_df.to_csv(stats_path, index=False)
            print(f"✓ 统计数据: {stats_path}")
            
            # 打印摘要
            print("\n📊 数据摘要:")
            for cat_type in ['2cat', 'eukaryote']:
                cat_stats = stats_df[stats_df['Category_Type'] == cat_type]
                if not cat_stats.empty:
                    print(f"\n{cat_type}分类:")
                    for _, row in cat_stats.iterrows():
                        print(f"  {row['Category']}: {row['Count']:,} genomes, "
                              f"GC={row['GC_Mean']:.3f}±{row['GC_Std']:.3f} "
                              f"({row['GC_Min']:.3f}-{row['GC_Max']:.3f}), "
                              f"Density={row['Density_Mean']:.2f}±{row['Density_Std']:.2f} IM/Mb")

def main():
    """主函数"""
    # 先设置字体
    set_font_settings()
    
    # 文件路径
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    simulation_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    real_genome_file = Path(base_dir) / "all_species_with_GCcontent_and_density.txt"
    
    print("=" * 70)
    print("Figure 1: 双版本双图对比分析 - 最终版")
    print("=" * 70)
    print("版本1: 完整信息版（含n值）")
    print("版本2: 发表纯净版（无n值）")
    print("=" * 70)
    print("图形内容:")
    print("1. 原核与真核比较")
    print("2. 真核细分: 原生动物、真菌、植物、无脊椎动物、脊椎动物其他、哺乳动物")
    print("=" * 70)
    print("改进点:")
    print("1. 真核细分图例位置：左上角")
    print("2. 背景曲线描述：Random simulation background")
    print("3. 背景曲线颜色：绿色（与数据点更好区分）")
    print("=" * 70)
    
    # 初始化绘图器
    plotter = FinalDensityPlotter(simulation_file, real_genome_file)
    
    try:
        # 1. 加载真实数据
        if not plotter.load_and_preprocess_data():
            print("❌ 真实数据加载失败")
            return
        
        # 2. 加载/创建随机模拟背景曲线
        if not plotter.load_theoretical_background():
            print("❌ 随机模拟背景加载失败")
            return
        
        # 3. 创建自适应的平滑随机模拟背景曲线
        if not plotter.create_adaptive_smooth_curve():
            print("❌ 自适应平滑随机模拟背景曲线创建失败")
            return
        
        # 4. 导出统计数据
        output_dir = Path(base_dir) / "figures" / "figure1_two_versions"
        plotter.export_statistics(output_dir)
        
        # 5. 保存两个版本的图形
        saved_files = plotter.save_both_versions(output_dir)
        
        print("\n" + "=" * 70)
        print("✅ 两个版本Figure生成完成!")
        print("=" * 70)
        print(f"主要输出文件:")
        
        # 按版本分类显示
        print("\n📁 完整信息版:")
        for file_path in saved_files:
            if "complete_version" in str(file_path):
                print(f"  {file_path.suffix.upper()}: {file_path.name}")
        
        print("\n📁 发表纯净版:")
        for file_path in saved_files:
            if "clean_version" in str(file_path):
                print(f"  {file_path.suffix.upper()}: {file_path.name}")
        
        print("=" * 70)
        
        # 显示关键信息
        print("\n🎨 版本对比:")
        print("1. 完整信息版:")
        print("   - 包含各分类的样本数量 (n=)")
        print("   - 适合内部报告和审稿人查看")
        print("2. 发表纯净版:")
        print("   - 简洁，只显示分类名称")
        print("   - 适合最终发表图形")
        print("3. 统一改进:")
        print("   - 图例位置：都在左上角（真核细分使用两列）")
        print("   - 背景描述：Random simulation background")
        print("   - 颜色方案：绿色背景曲线，与数据点形成良好对比")
        print("=" * 70)
        
    except Exception as e:
        print(f"\n❌ 分析过程中出现错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
