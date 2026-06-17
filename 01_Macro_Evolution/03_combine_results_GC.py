#!/usr/bin/env python3
"""
combine_results_GC.py - Combine GC content simulation results
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
    
    print(f"=== Combining all GC content simulation results ===")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"GC content gradients: {len(GC_LEVELS)} levels (10% - 90%)")
    print(f"Replicates per gradient: {REPLICATES}")
    print(f"Expected total simulations: {len(GC_LEVELS) * REPLICATES}")
    
    # 查找所有单个GC含量的结果文件
    result_files = glob.glob(str(results_dir / "simulation_GC*_detailed.csv"))
    
    print(f"\nFound {len(result_files)} result files")
    
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
                    print(f"✓ GC{gc_percent}%: {len(df)}/{REPLICATES} replicates")
                else:
                    print(f"✗ GC{gc_percent}%: File is empty")
                    missing_levels.append(gc_level)
            except Exception as e:
                print(f"✗ GC{gc_percent}%: Read error - {e}")
                missing_levels.append(gc_level)
        else:
            print(f"✗ GC{gc_percent}%: File does not exist")
            missing_levels.append(gc_level)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_file = results_dir / "simulation_results_GC_combined.csv"
        combined_df.to_csv(combined_file, index=False)
        
        total_simulations = len(combined_df)
        expected_total = len(GC_LEVELS) * REPLICATES
        completion_rate = total_simulations / expected_total * 100
        
        print(f"\nMerge completed!")
        print(f"Total data rows: {total_simulations}")
        print(f"Completion rate: {completion_rate:.1f}% ({total_simulations}/{expected_total})")
        print(f"Saved to: {combined_file}")
        
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
        
        # Flatten multi-level column index
        summary.columns = [
            'Density_mean', 'Density_std', 'Replicate_count', 
            'Density_min', 'Density_max', 'Density_sem',
            'iMotif_mean', 'iMotif_total',
            'Actual_GC_mean', 'Actual_C_mean', 'Actual_G_mean',
            'Pos_C_mean', 'Neg_C_mean', 'Run_time_mean'
        ]
        summary = summary.reset_index()
        
        # Calculate confidence interval
        summary['ci_lower'] = summary['Density_mean'] - 1.96 * summary['Density_sem']
        summary['ci_upper'] = summary['Density_mean'] + 1.96 * summary['Density_sem']
        
        # Calculate coefficient of variation
        summary['cv'] = (summary['Density_std'] / summary['Density_mean']) * 100
        
        # Calculate C content to GC content ratio
        summary['C_to_GC_ratio'] = summary['Actual_C_mean'] / summary['GC_content']
        
        print(f"\n=== Summary Statistics ===")
        print(summary.to_string())

        # Calculate confidence intervals and correlations
        print("\n=== 95% Confidence Intervals ===")
        summary_stats = []
        for _, row in summary.iterrows():
            gc_content = row['GC_content']
            mean_density = row['Density_mean']
            std_density = row['Density_std']
            n = row['Replicate_count']
            gc_percent = int(gc_content * 100)
            
            # Calculate standard error and confidence interval
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
        
        # Save summary statistics
        summary_df = pd.DataFrame(summary_stats)
        summary_file = results_dir / "simulation_summary_statistics_GC.csv"
        summary_df.to_csv(summary_file, index=False)
        print(f"\nSummary statistics saved to: {summary_file}")
        
        # Calculate overall statistics and correlations
        print("\n=== Overall Statistics ===")
        print(f"Total i-motif count: {combined_df['iMotif_count'].sum():,}")
        print(f"Mean density: {combined_df['Density_IM_per_Mb'].mean():.2f} IM/Mb")
        print(f"Density SD: {combined_df['Density_IM_per_Mb'].std():.2f} IM/Mb")
        print(f"Total runtime: {combined_df['Run_time_seconds'].sum()/3600:.1f} hours")
        
        # Calculate overall correlations
        print("\n=== Overall Correlation Analysis ===")
        print(f"C content vs density: {combined_df['Actual_C_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"G content vs density: {combined_df['Actual_G_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"GC content vs density: {combined_df['Actual_GC_content'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"Positive strand C vs density: {combined_df['Pos_strand_C'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        print(f"Negative strand C vs density: {combined_df['Neg_strand_C'].corr(combined_df['Density_IM_per_Mb']):.3f}")
        
        # Check for missing GC levels
        if missing_levels:
            print(f"\nWarning: Missing or incomplete data for the following GC levels:")
            for gc in missing_levels:
                print(f"  GC{int(gc*100)}%")
        
        # Save key findings
        key_findings = results_dir / "key_findings.txt"
        with open(key_findings, 'w') as f:
            f.write("=== GC Content Simulation Key Findings ===\n\n")
            f.write(f"Experiment completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总模拟次数: {total_simulations}/{expected_total}\n")
            f.write(f"Completion rate: {completion_rate:.1f}%\n\n")

            # Find highest and lowest density
            max_idx = summary_df['Density_mean'].idxmax()
            min_idx = summary_df['Density_mean'].idxmin()
            
            f.write(f"Highest i-motif density: {summary_df.loc[max_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[max_idx, 'GC_percent']}%\n")
            f.write(f"Lowest i-motif density: {summary_df.loc[min_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[min_idx, 'GC_percent']}%\n\n")
            
            # Calculate density range
            density_range = summary_df['Density_mean'].max() - summary_df['Density_mean'].min()
            f.write(f"Density range: {density_range:.1f} IM/Mb\n")

            # Calculate overall C content vs density correlation
            overall_corr = combined_df['Actual_C_content'].corr(combined_df['Density_IM_per_Mb'])
            f.write(f"C content vs density overall correlation: {overall_corr:.3f}\n")
        
        print(f"\nKey findings saved to: {key_findings}")
        
        return combined_df, summary_df
    else:
        print("Error: No valid result files found")
        return None, None

if __name__ == "__main__":
    combine_all_results()
