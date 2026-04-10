#!/usr/bin/env python3
"""
plot_pure_gc_density.py - 纯净版GC含量与i-motif密度关系图
生成Nature期刊风格的核心图表
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns

# 设置Nature期刊纯净版样式
def set_nature_pure_style():
    """设置Nature期刊的纯净版绘图样式"""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 9,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.6,
        'ytick.minor.width': 0.6,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'figure.autolayout': False,
        'axes.grid': True,
        'grid.alpha': 0.2,
        'grid.linewidth': 0.5
    })

def load_data():
    """加载数据"""
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    combined_file = Path(base_dir) / "results" / "simulation_results_GC_combined.csv"
    
    if not combined_file.exists():
        raise FileNotFoundError(f"未找到数据文件: {combined_file}")
    
    df = pd.read_csv(combined_file)
    
    # 计算汇总统计
    summary = df.groupby('GC_content').agg({
        'Density_IM_per_Mb': ['mean', 'std', 'count', 'sem'],
        'Actual_C_content': 'mean',
        'Actual_G_content': 'mean'
    }).round(4)
    
    summary.columns = ['Density_mean', 'Density_std', 'Replicate_count', 'Density_sem',
                      'Actual_C_mean', 'Actual_G_mean']
    summary = summary.reset_index()
    
    # 计算置信区间 (95% CI)
    summary['ci_lower'] = summary['Density_mean'] - 1.96 * summary['Density_sem']
    summary['ci_upper'] = summary['Density_mean'] + 1.96 * summary['Density_sem']
    
    return df, summary

def create_confidence_interval_plot(df, summary):
    """创建带置信区间的散点图"""
    fig, ax = plt.subplots(figsize=(4.0, 3.0))  # Nature标准尺寸
    
    # 使用viridis颜色映射
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary)))
    
    # 绘制所有重复的数据点（带轻微抖动）
    for idx, (_, row) in enumerate(summary.iterrows()):
        gc_content = row['GC_content']
        subset = df[df['GC_content'] == gc_content]
        
        # 添加轻微水平抖动避免重叠
        jitter = np.random.normal(0, 0.0015, len(subset))
        ax.scatter(gc_content + jitter, subset['Density_IM_per_Mb'],
                  alpha=0.5, s=25, color=colors[idx], edgecolors='white', linewidth=0.3,
                  zorder=2)
    
    # 绘制均值曲线
    ax.plot(summary['GC_content'], summary['Density_mean'],
           color='#2E4057', linewidth=2, marker='o', markersize=6,
           markerfacecolor='white', markeredgecolor='#2E4057', markeredgewidth=1.5,
           zorder=3)
    
    # 绘制95%置信区间
    ax.fill_between(summary['GC_content'],
                   summary['ci_lower'], summary['ci_upper'],
                   alpha=0.15, color='#2E4057', zorder=1)
    
    # 设置坐标轴
    ax.set_xlabel('GC content', fontsize=10)
    ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10)
    
    # 设置网格
    ax.grid(True, alpha=0.2, linewidth=0.5, linestyle='-', zorder=0)
    
    # 设置刻度
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xticklabels(['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])
    
    # 设置坐标轴范围
    ax.set_xlim(0.08, 0.92)
    
    # 添加置信区间说明
    ax.text(0.02, 0.98, 'Mean ± 95% CI', transform=ax.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, pad=0.3))
    
    plt.tight_layout()
    return fig

def create_fitted_curve_plot(df, summary):
    """创建带拟合曲线的图"""
    fig, ax = plt.subplots(figsize=(4.0, 3.0))  # Nature标准尺寸
    
    # 使用viridis颜色映射
    colors = plt.cm.viridis(np.linspace(0, 1, len(summary)))
    
    # 绘制所有重复的数据点
    for idx, (_, row) in enumerate(summary.iterrows()):
        gc_content = row['GC_content']
        subset = df[df['GC_content'] == gc_content]
        
        jitter = np.random.normal(0, 0.0015, len(subset))
        ax.scatter(gc_content + jitter, subset['Density_IM_per_Mb'],
                  alpha=0.4, s=25, color=colors[idx], edgecolors='white', linewidth=0.3,
                  zorder=2)
    
    # 绘制均值曲线
    ax.plot(summary['GC_content'], summary['Density_mean'],
           color='#2E4057', linewidth=2, marker='o', markersize=6,
           markerfacecolor='white', markeredgecolor='#2E4057', markeredgewidth=1.5,
           zorder=3, label='Mean')
    
    # 尝试多种拟合，选择最佳
    x = summary['GC_content'].values
    y = summary['Density_mean'].values
    
    # 1. 二次多项式拟合
    try:
        coeff_quad = np.polyfit(x, y, 2)
        poly_quad = np.poly1d(coeff_quad)
        x_fit = np.linspace(x.min(), x.max(), 200)
        y_fit_quad = poly_quad(x_fit)
        r2_quad = 1 - np.sum((y - poly_quad(x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_quad = -np.inf
    
    # 2. 三次多项式拟合
    try:
        coeff_cubic = np.polyfit(x, y, 3)
        poly_cubic = np.poly1d(coeff_cubic)
        y_fit_cubic = poly_cubic(x_fit)
        r2_cubic = 1 - np.sum((y - poly_cubic(x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_cubic = -np.inf
    
    # 3. 指数拟合
    try:
        # y = a * exp(b*x)
        log_y = np.log(y)
        coeff_exp = np.polyfit(x, log_y, 1)
        a = np.exp(coeff_exp[1])
        b = coeff_exp[0]
        y_fit_exp = a * np.exp(b * x_fit)
        r2_exp = 1 - np.sum((y - a*np.exp(b*x))**2) / np.sum((y - np.mean(y))**2)
    except:
        r2_exp = -np.inf
    
    # 选择R²最高的拟合
    fits = {
        'Quadratic': (r2_quad, y_fit_quad, '#D62728'),
        'Cubic': (r2_cubic, y_fit_cubic, '#FF7F0E'),
        'Exponential': (r2_exp, y_fit_exp, '#2CA02C')
    }
    
    best_fit_name = max(fits, key=lambda k: fits[k][0])
    best_r2, best_y_fit, best_color = fits[best_fit_name]
    
    # 绘制最佳拟合曲线
    ax.plot(x_fit, best_y_fit, '--', color=best_color, linewidth=2, 
           alpha=0.8, zorder=4, label=f'{best_fit_name} fit')
    
    # 添加R²值
    if best_r2 > 0:
        ax.text(0.98, 0.98, f'R² = {best_r2:.3f}\n({best_fit_name})', 
                transform=ax.transAxes, fontsize=8, verticalalignment='top',
                horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, pad=0.3))
    
    # 设置坐标轴
    ax.set_xlabel('GC content', fontsize=10)
    ax.set_ylabel('i-Motif density (IM/Mb)', fontsize=10)
    
    # 设置网格
    ax.grid(True, alpha=0.2, linewidth=0.5, linestyle='-', zorder=0)
    
    # 设置刻度
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    ax.set_xticklabels(['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])
    
    # 设置坐标轴范围
    ax.set_xlim(0.08, 0.92)
    
    # 添加图例
    ax.legend(loc='upper left', frameon=True, fancybox=True, 
              framealpha=0.8, edgecolor='none', fontsize=8)
    
    plt.tight_layout()
    return fig, best_fit_name, best_r2

def create_correlation_heatmap(df):
    """创建相关性热图（可选）"""
    # 计算Spearman和Pearson相关性
    corr_cols = ['Actual_GC_content', 'Actual_C_content', 'Actual_G_content', 'Density_IM_per_Mb']
    
    # Pearson相关性
    pearson_corr = df[corr_cols].corr(method='pearson')
    
    # Spearman相关性
    spearman_corr = df[corr_cols].corr(method='spearman')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.8))
    
    # 重命名列名
    display_names = ['GC', 'C', 'G', 'Density']
    pearson_corr.index = display_names
    pearson_corr.columns = display_names
    spearman_corr.index = display_names
    spearman_corr.columns = display_names
    
    # Pearson热图
    im1 = ax1.imshow(pearson_corr.values, cmap='coolwarm', vmin=-1, vmax=1)
    for i in range(len(pearson_corr)):
        for j in range(len(pearson_corr)):
            ax1.text(j, i, f'{pearson_corr.iloc[i, j]:.2f}', 
                    ha='center', va='center', 
                    color='white' if abs(pearson_corr.iloc[i, j]) > 0.5 else 'black',
                    fontsize=8)
    
    ax1.set_xticks(range(len(display_names)))
    ax1.set_yticks(range(len(display_names)))
    ax1.set_xticklabels(display_names)
    ax1.set_yticklabels(display_names)
    ax1.set_title('Pearson Correlation', fontsize=10, pad=10)
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    # Spearman热图
    im2 = ax2.imshow(spearman_corr.values, cmap='coolwarm', vmin=-1, vmax=1)
    for i in range(len(spearman_corr)):
        for j in range(len(spearman_corr)):
            ax2.text(j, i, f'{spearman_corr.iloc[i, j]:.2f}', 
                    ha='center', va='center', 
                    color='white' if abs(spearman_corr.iloc[i, j]) > 0.5 else 'black',
                    fontsize=8)
    
    ax2.set_xticks(range(len(display_names)))
    ax2.set_yticks(range(len(display_names)))
    ax2.set_xticklabels(display_names)
    ax2.set_yticklabels(display_names)
    ax2.set_title('Spearman Correlation', fontsize=10, pad=10)
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    return fig

def save_figures(figures_dict, base_dir):
    """保存所有图形"""
    output_dir = Path(base_dir) / "figures_pure"
    output_dir.mkdir(exist_ok=True)
    
    for name, fig in figures_dict.items():
        fig.savefig(output_dir / f"{name}.png", dpi=600, bbox_inches='tight')
        fig.savefig(output_dir / f"{name}.pdf", bbox_inches='tight')
        fig.savefig(output_dir / f"{name}.svg", bbox_inches='tight')
        print(f"✓ 保存: {name}.png/pdf/svg")
    
    print(f"\n所有纯净版图形已保存到: {output_dir}")

def main():
    """主函数"""
    # 设置纯净版Nature样式
    set_nature_pure_style()
    
    print("=== 纯净版GC含量与i-motif密度关系图 ===")
    
    # 加载数据
    try:
        df, summary = load_data()
        print(f"✓ 成功加载数据，共{len(df)}次模拟")
        print(f"✓ GC含量梯度: {len(summary)}个 (10%-90%)")
        
        # 显示基本统计
        print(f"\n密度统计摘要:")
        print(f"  最高密度: {summary['Density_mean'].max():.1f} IM/Mb @ GC{summary.loc[summary['Density_mean'].idxmax(), 'GC_content']*100:.0f}%")
        print(f"  最低密度: {summary['Density_mean'].min():.1f} IM/Mb @ GC{summary.loc[summary['Density_mean'].idxmin(), 'GC_content']*100:.0f}%")
        print(f"  平均变异系数: {((summary['Density_std'] / summary['Density_mean']) * 100).mean():.1f}%")
        
    except FileNotFoundError as e:
        print(f"✗ 错误: {e}")
        print("请先运行 combine_results_GC.py 合并数据")
        return
    
    # 创建图形
    print("\n创建图形...")
    
    # 1. 带置信区间的散点图
    fig_ci = create_confidence_interval_plot(df, summary)
    
    # 2. 带拟合曲线的图
    fig_fit, fit_name, fit_r2 = create_fitted_curve_plot(df, summary)
    print(f"✓ 最佳拟合: {fit_name}, R² = {fit_r2:.3f}")
    
    # 3. 可选：相关性热图
    try:
        fig_corr = create_correlation_heatmap(df)
        print("✓ 创建相关性热图")
        include_corr = True
    except:
        print("⚠ 跳过相关性热图（数据可能存在问题）")
        include_corr = False
    
    # 保存图形
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    figures = {
        'gc_density_confidence_interval': fig_ci,
        'gc_density_fitted_curve': fig_fit,
    }
    
    if include_corr:
        figures['correlation_heatmap'] = fig_corr
    
    save_figures(figures, base_dir)
    
    # 计算并显示相关性
    print("\n=== 相关性分析 ===")
    
    # Pearson相关性（线性关系）
    pearson_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'], method='pearson')
    pearson_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'], method='pearson')
    
    # Spearman相关性（单调关系）
    spearman_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'], method='spearman')
    spearman_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'], method='spearman')
    
    print(f"Pearson相关系数:")
    print(f"  GC含量 vs 密度: {pearson_gc:.3f}")
    print(f"  C含量 vs 密度: {pearson_c:.3f}")
    
    print(f"\nSpearman相关系数:")
    print(f"  GC含量 vs 密度: {spearman_gc:.3f}")
    print(f"  C含量 vs 密度: {spearman_c:.3f}")
    
    # 判断使用哪种相关性
    if abs(spearman_gc - pearson_gc) > 0.1:
        print("\n⚠ 建议使用Spearman相关性（数据可能非正态分布）")
    else:
        print("\n✓ Pearson和Spearman相关性接近，均可使用")
    
    # 显示图形
    plt.show()
    
    print("\n=== 完成 ===")
    print("纯净版图形已生成，适用于论文发表。")

if __name__ == "__main__":
    main()
