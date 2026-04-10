#!/usr/bin/env python3
"""
combine_results_GC.py - 合并GC含量模拟实验结果
"""

import pandas as pd
from pathlib import Path
import glob
import numpy as np
from datetime import datetime

def combine_all_results():
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    results_dir = Path(base_dir) / "results"
    
    # 硬编码配置
    GC_LEVELS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
    REPLICATES = 30
    GENOME_SIZE = 10000000
    
    print(f"=== 合并所有GC含量模拟结果 ===")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"GC含量梯度: {len(GC_LEVELS)}个 (10% - 90%)")
    print(f"每个梯度重复次数: {REPLICATES}")
    print(f"预期总模拟次数: {len(GC_LEVELS) * REPLICATES}")
    
    # 查找所有单个GC含量的结果文件
    result_files = glob.glob(str(results_dir / "simulation_GC*_detailed.csv"))
    
    print(f"\n找到 {len(result_files)} 个结果文件")
    
    all_data = []
    missing_levels = []
    
    for gc_level in GC_LEVELS:
        gc_percent = int(gc_level * 100)
        expected_file = results_dir / f"simulation_GC{gc_percent}_detailed.csv"
        
        if expected_file.exists():
            try:
                df = pd.read_csv(expected_file)
                if len(df) > 0:
                    all_data.append(df)
                    print(f"✓ GC{gc_percent}%: {len(df)}/{REPLICATES} 个重复")
                else:
                    print(f"✗ GC{gc_percent}%: 文件为空")
                    missing_levels.append(gc_level)
            except Exception as e:
                print(f"✗ GC{gc_percent}%: 读取错误 - {e}")
                missing_levels.append(gc_level)
        else:
            print(f"✗ GC{gc_percent}%: 文件不存在")
            missing_levels.append(gc_level)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_file = results_dir / "simulation_results_GC_combined.csv"
        combined_df.to_csv(combined_file, index=False)
        
        total_simulations = len(combined_df)
        expected_total = len(GC_LEVELS) * REPLICATES
        completion_rate = total_simulations / expected_total * 100
        
        print(f"\n合并完成!")
        print(f"总数据行: {total_simulations}")
        print(f"完成率: {completion_rate:.1f}% ({total_simulations}/{expected_total})")
        print(f"保存到: {combined_file}")
        
        # 计算汇总统计
        summary = combined_df.groupby('GC_content').agg({
            'Density_IM_per_Mb': ['mean', 'std', 'count', 'min', 'max', 'sem'],
            'iMotif_count': ['mean', 'sum'],
            'Actual_GC_content': 'mean',
            'Actual_C_content': 'mean',
            'Actual_G_content': 'mean',
            'Pos_strand_C': 'mean',
            'Neg_strand_C': 'mean',
            'Run_time_seconds': 'mean'
        }).round(4)
        
        # 展平多级列索引
        summary.columns = [
            'Density_mean', 'Density_std', 'Replicate_count', 
            'Density_min', 'Density_max', 'Density_sem',
            'iMotif_mean', 'iMotif_total',
            'Actual_GC_mean', 'Actual_C_mean', 'Actual_G_mean',
            'Pos_C_mean', 'Neg_C_mean', 'Run_time_mean'
        ]
        summary = summary.reset_index()
        
        # 计算置信区间
        summary['ci_lower'] = summary['Density_mean'] - 1.96 * summary['Density_sem']
        summary['ci_upper'] = summary['Density_mean'] + 1.96 * summary['Density_sem']
        
        # 计算变异系数
        summary['cv'] = (summary['Density_std'] / summary['Density_mean']) * 100
        
        # 计算C含量与GC含量的比例
        summary['C_to_GC_ratio'] = summary['Actual_C_mean'] / summary['GC_content']
        
        print(f"\n=== 汇总统计 ===")
        print(summary.to_string())
        
        # 计算置信区间和相关性
        print("\n=== 95% 置信区间 ===")
        summary_stats = []
        for _, row in summary.iterrows():
            gc_content = row['GC_content']
            mean_density = row['Density_mean']
            std_density = row['Density_std']
            n = row['Replicate_count']
            gc_percent = int(gc_content * 100)
            
            # 计算标准误和置信区间
            se = std_density / np.sqrt(n)
            ci_lower = mean_density - 1.96 * se
            ci_upper = mean_density + 1.96 * se
            
            print(f"GC{gc_percent}%: {mean_density:.2f} ± {1.96*se:.2f} IM/Mb "
                  f"[{ci_lower:.2f}, {ci_upper:.2f}] (n={n})")
            
            # 提取对应分组的数据计算相关性
            group_data = combined_df[combined_df['GC_content'] == gc_content]
            corr_c = group_data['Actual_C_content'].corr(group_data['Density_IM_per_Mb']) if len(group_data) > 1 else np.nan
            
            summary_stats.append({
                'GC_content': gc_content,
                'GC_percent': gc_percent,
                'Density_mean': mean_density,
                'Density_std': std_density,
                'CI_lower': ci_lower,
                'CI_upper': ci_upper,
                'Replicates': n,
                'SE': se,
                'CV': row['cv'],
                'Actual_C_mean': row['Actual_C_mean'],
                'Actual_G_mean': row['Actual_G_mean'],
                'C_to_GC_ratio': row['C_to_GC_ratio'],
                'Corr_C_Density': corr_c
            })
        
        # 保存汇总统计
        summary_df = pd.DataFrame(summary_stats)
        summary_file = results_dir / "simulation_summary_statistics_GC.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"\n汇总统计保存到: {summary_file}")
        
        # 计算总体统计和相关性
        print("\n=== 总体统计 ===")
        print(f"总i-motif数量: {combined_df['iMotif_count'].sum():,}")
        print(f"平均密度: {combined_df['Density_IM_per_Mb'].mean():.2f} IM/Mb")
        print(f"密度标准差: {combined_df['Density_IM_per_Mb'].std():.2f} IM/Mb")
        print(f"总运行时间: {combined_df['Run_time_seconds'].sum()/3600:.1f} 小时")
        
        # 计算整体相关性
        print("\n=== 整体相关性分析 ===")
        print(f"C含量与密度相关性: {combined_df['Actual_C_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"G含量与密度相关性: {combined_df['Actual_G_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"GC含量与密度相关性: {combined_df['Actual_GC_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"正链C与密度相关性: {combined_df['Pos_strand_C'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"负链C与密度相关性: {combined_df['Neg_strand_C'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        
        # 检查缺失的GC含量
        if missing_levels:
            print(f"\n警告: 以下GC含量数据缺失或不全:")
            for gc in missing_levels:
                print(f"  GC{int(gc*100)}%")
        
        # 保存关键发现
        key_findings = results_dir / "key_findings.txt"
        with open(key_findings, 'w') as f:
            f.write("=== GC含量模拟实验关键发现 ===\n\n")
            f.write(f"实验完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总模拟次数: {total_simulations}/{expected_total}\n")
            f.write(f"完成率: {completion_rate:.1f}%\n\n")
            
            # 找到最高和最低密度
            max_idx = summary_df['Density_mean'].idxmax()
            min_idx = summary_df['Density_mean'].idxmin()
            
            f.write(f"最高i-motif密度: {summary_df.loc[max_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[max_idx, 'GC_percent']}%\n")
            f.write(f"最低i-motif密度: {summary_df.loc[min_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[min_idx, 'GC_percent']}%\n\n")
            
            # 计算密度变化范围
            density_range = summary_df['Density_mean'].max() - summary_df['Density_mean'].min()
            f.write(f"密度变化范围: {density_range:.1f} IM/Mb\n")
            
            # 计算C含量与密度的整体相关性
            overall_corr = combined_df['Actual_C_content'].corr(combined_df['Density_IM_per_Mb'])
            f.write(f"C含量与密度整体相关性: {overall_corr:.3f}\n")
        
        print(f"\n关键发现保存到: {key_findings}")
        
        return combined_df, summary_df
    else:
        print("错误: 没有找到任何有效的结果文件")
        return None, None

if __name__ == "__main__":
    combine_all_results()
