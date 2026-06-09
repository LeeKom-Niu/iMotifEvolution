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
    
    GC_LEVELS = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
    REPLICATES = 30
    GENOME_SIZE = 10000000
    
    
    result_files = glob.glob(str(results_dir / "simulation_GC*_detailed.csv"))
    
    
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
                else:
                    missing_levels.append(gc_level)
            except Exception as e:
                missing_levels.append(gc_level)
        else:
            missing_levels.append(gc_level)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        combined_file = results_dir / "simulation_results_GC_combined.csv"
        combined_df.to_csv(combined_file, index=False)
        
        total_simulations = len(combined_df)
        expected_total = len(GC_LEVELS) * REPLICATES
        completion_rate = total_simulations / expected_total * 100
        
        
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
        
        summary.columns = [
            'Density_mean', 'Density_std', 'Replicate_count', 
            'Density_min', 'Density_max', 'Density_sem',
            'iMotif_mean', 'iMotif_total',
            'Actual_GC_mean', 'Actual_C_mean', 'Actual_G_mean',
            'Pos_C_mean', 'Neg_C_mean', 'Run_time_mean'
        ]
        summary = summary.reset_index()
        
        summary['ci_lower'] = summary['Density_mean'] - 1.96 * summary['Density_sem']
        summary['ci_upper'] = summary['Density_mean'] + 1.96 * summary['Density_sem']
        
        summary['cv'] = (summary['Density_std'] / summary['Density_mean']) * 100
        
        summary['C_to_GC_ratio'] = summary['Actual_C_mean'] / summary['GC_content']
        
        
        summary_stats = []
        for _, row in summary.iterrows():
            gc_content = row['GC_content']
            mean_density = row['Density_mean']
            std_density = row['Density_std']
            n = row['Replicate_count']
            gc_percent = int(gc_content * 100)
            
            se = std_density / np.sqrt(n)
            ci_lower = mean_density - 1.96 * se
            ci_upper = mean_density + 1.96 * se
            
                  f"[{ci_lower:.2f}, {ci_upper:.2f}] (n={n})")
            
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
        
        summary_df = pd.DataFrame(summary_stats)
        summary_file = results_dir / "simulation_summary_statistics_GC.csv"
        summary_df.to_csv(summary_file, index=False)
        
        
        
        if missing_levels:
            for gc in missing_levels:
        
        key_findings = results_dir / "key_findings.txt"
        with open(key_findings, 'w') as f:
            f.write("=== GC含量模拟实验关键发现 ===\n\n")
            f.write(f"实验完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总模拟次数: {total_simulations}/{expected_total}\n")
            f.write(f"完成率: {completion_rate:.1f}%\n\n")
            
            max_idx = summary_df['Density_mean'].idxmax()
            min_idx = summary_df['Density_mean'].idxmin()
            
            f.write(f"最高i-motif密度: {summary_df.loc[max_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[max_idx, 'GC_percent']}%\n")
            f.write(f"最低i-motif密度: {summary_df.loc[min_idx, 'Density_mean']:.1f} IM/Mb "
                   f"@ GC{summary_df.loc[min_idx, 'GC_percent']}%\n\n")
            
            density_range = summary_df['Density_mean'].max() - summary_df['Density_mean'].min()
            f.write(f"密度变化范围: {density_range:.1f} IM/Mb\n")
            
            overall_corr = combined_df['Actual_C_content'].corr(combined_df['Density_IM_per_Mb'])
            f.write(f"C含量与密度整体相关性: {overall_corr:.3f}\n")
        
        
        return combined_df, summary_df
    else:
        return None, None
if __name__ == "__main__":
    combine_all_results()
