#!/usr/bin/env python3
"""
parallel_simulation_GC.py - GC含量背景模拟实验并行脚本
生成不同GC含量的随机序列，运行iM-Seeker检测i-motif密度
"""

import sys
import os
import numpy as np
import multiprocessing as mp
import subprocess
import pandas as pd
from pathlib import Path
import random
import json
from datetime import datetime

# 硬编码配置参数
CONFIG = {
    'simulation': {
        'at_ratio': 1.0,  # A和T的比例 (A:T = 1:1)
        'gc_equal': True,  # G和C是否均等分配
        'im_seeker': {
            'model_dir': "/datapool/home/2023200496/niulk/my_project/genome/iM_Seeker_model"
        }
    }
}

def load_config():
    """返回硬编码的配置"""
    return CONFIG

def generate_sequence_gc(args):
    """生成单条序列的任务函数 - GC含量版本"""
    gc_level, replicate, genome_size, base_dir, config = args
    
    # 为每个进程设置不同的随机种子
    random_seed = hash((gc_level, replicate, datetime.now().timestamp())) % (2**32)
    np.random.seed(random_seed)
    random.seed(random_seed)
    
    gc_percent = int(gc_level * 100)
    work_dir = Path(base_dir) / "sequences" / f"GC{gc_percent}_rep{replicate:03d}"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    sequence_file = work_dir / "random_seq.fa"
    output_dir = work_dir / "output"
    log_file = work_dir / "generation.log"
    
    log_content = []
    log_content.append(f"=== GC{gc_percent}_rep{replicate:03d} ===")
    log_content.append(f"随机种子: {random_seed}")
    log_content.append(f"GC含量目标: {gc_level}")
    log_content.append(f"基因组大小: {genome_size} bp")
    
    print(f"生成 GC{gc_percent}_rep{replicate:03d} (种子: {random_seed})")
    
    # 1. 生成随机序列 - 根据GC含量和配置
    total_gc = int(genome_size * gc_level)
    total_at = genome_size - total_gc
    
    # 从配置中获取参数
    at_ratio = config['simulation']['at_ratio']
    gc_equal = config['simulation']['gc_equal']
    
    # 计算G和C的数量
    if gc_equal:
        # G和C均等分配
        num_g = total_gc // 2
        num_c = total_gc - num_g  # 处理奇数情况
    else:
        # 可以修改为G和C不同比例，默认均等
        num_g = total_gc // 2
        num_c = total_gc - num_g
    
    # 计算A和T的数量
    total_ratio = at_ratio + 1.0  # A:T = at_ratio:1
    num_a = int(total_at * at_ratio / total_ratio)
    num_t = total_at - num_a
    
    # 创建包含所有碱基的列表
    bases = (['G'] * num_g + 
             ['C'] * num_c + 
             ['A'] * num_a + 
             ['T'] * num_t)
    
    # 充分打乱碱基顺序
    np.random.shuffle(bases)
    random_sequence = ''.join(bases)
    
    # 验证碱基含量
    total_len = len(random_sequence)
    actual_g = random_sequence.count('G') / total_len
    actual_c = random_sequence.count('C') / total_len
    actual_a = random_sequence.count('A') / total_len
    actual_t = random_sequence.count('T') / total_len
    actual_gc = actual_g + actual_c
    
    # 计算正链C含量（i-motif形成位点）
    pos_strand_c = actual_c
    # 计算负链C含量（正链的G含量）
    neg_strand_c = actual_g
    # 计算正链G含量
    pos_strand_g = actual_g
    
    log_content.append(f"目标GC含量: {gc_level:.4f}")
    log_content.append(f"实际GC含量: {actual_gc:.4f}")
    log_content.append(f"实际G含量: {actual_g:.4f}")
    log_content.append(f"实际C含量: {actual_c:.4f}")
    log_content.append(f"实际A含量: {actual_a:.4f}")
    log_content.append(f"实际T含量: {actual_t:.4f}")
    log_content.append(f"正链C含量: {pos_strand_c:.4f}")
    log_content.append(f"正链G含量: {pos_strand_g:.4f}")
    log_content.append(f"负链C含量: {neg_strand_c:.4f}")
    log_content.append(f"序列长度: {total_len} bp")
    log_content.append(f"G/C比例: {actual_g/actual_c:.3f}" if actual_c > 0 else "G/C比例: N/A (C=0)")
    log_content.append(f"A/T比例: {actual_a/actual_t:.3f}" if actual_t > 0 else "A/T比例: N/A (T=0)")
    
    print(f"GC{gc_percent}_rep{replicate:03d}: 目标GC={gc_level:.3f}, 实际GC={actual_gc:.3f}")
    print(f"  正链C={pos_strand_c:.3f}, G={pos_strand_g:.3f}, 负链C={neg_strand_c:.3f}")
    
    # 保存序列
    with open(sequence_file, 'w') as f:
        f.write(f'>random_GC_{gc_level:.3f}_rep_{replicate:03d}_seed_{random_seed}\n')
        for i in range(0, len(random_sequence), 80):
            f.write(random_sequence[i:i+80] + '\n')
    
    # 2. 运行im-seeker
    output_dir.mkdir(exist_ok=True)
    
    model_dir = config['simulation']['im_seeker']['model_dir']
    
    cmd = [
        'iM-Seeker.py',
        '--sequence', str(sequence_file),
        '--classification_model', f"{model_dir}/pickle_model_classification.pkl",
        '--regression_model', f"{model_dir}/pickle_model_regression.pkl", 
        '--output_folder', str(output_dir)
    ]
    
    log_content.append(f"运行命令: {' '.join(cmd)}")
    
    print(f"运行im-seeker: GC{gc_percent}_rep{replicate:03d}")
    
    start_time = datetime.now()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        run_time = (datetime.now() - start_time).total_seconds()
        
        log_content.append(f"运行时间: {run_time:.1f}秒")
        log_content.append(f"退出码: {result.returncode}")
        
        if result.returncode != 0:
            log_content.append(f"stderr: {result.stderr[:500]}")
            print(f"im-seeker失败 GC{gc_percent}_rep{replicate:03d}")
            # 保存日志
            with open(log_file, 'w') as f:
                f.write('\n'.join(log_content))
            return None
        else:
            log_content.append("im-seeker运行成功")
            print(f"im-seeker成功 GC{gc_percent}_rep{replicate:03d}")
            
    except Exception as e:
        log_content.append(f"运行异常: {str(e)}")
        print(f"运行错误 GC{gc_percent}_rep{replicate:03d}: {e}")
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_content))
        return None
    
    # 3. 解析结果
    imotif_count, motif_details = parse_results(output_dir, work_dir)
    
    log_content.append(f"找到i-motif数量: {imotif_count}")
    
    # 密度计算：基于20Mb（正负链）
    effective_length = genome_size * 2  # 10Mb单链 → 20Mb双链处理
    density = (imotif_count / effective_length) * 1000000
    
    # 计算连续C的分布（重要指标）
    c_distribution = analyze_c_distribution(random_sequence)
    log_content.append(f"连续C分布: {json.dumps(c_distribution)}")
    
    # 计算连续G的分布（研究G的潜在干扰）
    g_distribution = analyze_g_distribution(random_sequence)
    log_content.append(f"连续G分布: {json.dumps(g_distribution)}")
    
    log_content.append(f"i-motif密度: {density:.2f} IM/Mb")
    
    # 保存日志
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_content))
    
    print(f"GC{gc_percent}_rep{replicate:03d}: {imotif_count}个i-motif, 密度{density:.2f} IM/Mb")
    
    return {
        'GC_content': gc_level,
        'Replicate': replicate,
        'Random_seed': random_seed,
        'Genome_size': genome_size,
        'Effective_length_Mb': effective_length / 1000000,
        'iMotif_count': imotif_count,
        'Density_IM_per_Mb': density,
        'Actual_GC_content': actual_gc,
        'Actual_G_content': actual_g,
        'Actual_C_content': actual_c,
        'Actual_A_content': actual_a,
        'Actual_T_content': actual_t,
        'Pos_strand_C': pos_strand_c,
        'Pos_strand_G': pos_strand_g,
        'Neg_strand_C': neg_strand_c,
        'G_to_C_ratio': actual_g/actual_c if actual_c > 0 else None,
        'A_to_T_ratio': actual_a/actual_t if actual_t > 0 else None,
        'Run_time_seconds': run_time if 'run_time' in locals() else None,
        'C_distribution': json.dumps(c_distribution),
        'G_distribution': json.dumps(g_distribution)
    }

def parse_results(output_dir, work_dir):
    """解析im-seeker输出结果"""
    output_dir = Path(output_dir)
    
    # 查找最终的预测文件
    result_file = output_dir / "iM-seeker_final_prediction.txt"
    
    imotif_count = 0
    motif_details = []
    
    if result_file.exists():
        try:
            with open(result_file, 'r') as f:
                lines = f.readlines()
            
            # 解析i-motif详细信息
            if len(lines) > 1:  # 有标题行和至少一个结果
                # 标题行
                header = lines[0].strip().split('\t')
                
                # 解析每个motif
                for line in lines[1:]:
                    if line.strip():
                        imotif_count += 1
                        parts = line.strip().split('\t')
                        if len(parts) >= 6:
                            motif_details.append({
                                'chromosome': parts[0],
                                'start': int(parts[1]),
                                'end': int(parts[2]),
                                'strand': parts[3],
                                'score': float(parts[4]) if parts[4] != 'NA' else None,
                                'sequence': parts[5] if len(parts) > 5 else ''
                            })
            
            print(f"找到 {imotif_count} 个i-motif预测")
            
            # 保存详细结果
            if motif_details:
                details_file = work_dir / "motif_details.json"
                with open(details_file, 'w') as f:
                    json.dump(motif_details, f, indent=2)
                    
        except Exception as e:
            print(f"解析文件错误 {result_file}: {e}")
            # 尝试简单计数
            with open(result_file, 'r') as f:
                lines = f.readlines()
            if len(lines) > 0:
                imotif_count = len(lines) - 1  # 减去标题行
    else:
        print(f"结果文件不存在: {result_file}")
    
    return imotif_count, motif_details

def analyze_c_distribution(sequence):
    """分析连续C的分布"""
    c_distribution = {1: 0, 2: 0, 3: 0, '4+': 0}
    
    current_c_streak = 0
    
    for base in sequence:
        if base == 'C':
            current_c_streak += 1
        else:
            if current_c_streak > 0:
                if current_c_streak >= 4:
                    c_distribution['4+'] += 1
                elif current_c_streak in c_distribution:
                    c_distribution[current_c_streak] += 1
                current_c_streak = 0
    
    # 检查末尾的C连续
    if current_c_streak > 0:
        if current_c_streak >= 4:
            c_distribution['4+'] += 1
        elif current_c_streak in c_distribution:
            c_distribution[current_c_streak] += 1
    
    return c_distribution

def analyze_g_distribution(sequence):
    """分析连续G的分布"""
    g_distribution = {1: 0, 2: 0, 3: 0, '4+': 0}
    
    current_g_streak = 0
    
    for base in sequence:
        if base == 'G':
            current_g_streak += 1
        else:
            if current_g_streak > 0:
                if current_g_streak >= 4:
                    g_distribution['4+'] += 1
                elif current_g_streak in g_distribution:
                    g_distribution[current_g_streak] += 1
                current_g_streak = 0
    
    # 检查末尾的G连续
    if current_g_streak > 0:
        if current_g_streak >= 4:
            g_distribution['4+'] += 1
        elif current_g_streak in g_distribution:
            g_distribution[current_g_streak] += 1
    
    return g_distribution

def main():
    if len(sys.argv) != 4:
        print("用法: python parallel_simulation_GC.py <GC_level> <replicates> <genome_size>")
        print("示例: python parallel_simulation_GC.py 0.40 30 10000000")
        sys.exit(1)
    
    gc_level = float(sys.argv[1])
    replicates = int(sys.argv[2])
    genome_size = int(sys.argv[3])
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"
    
    # 加载配置
    config = load_config()
    
    # 准备任务参数
    tasks = [(gc_level, rep, genome_size, base_dir, config) 
             for rep in range(1, replicates + 1)]
    
    gc_percent = int(gc_level * 100)
    print(f"\n{'='*60}")
    print(f"开始并行处理 GC{gc_percent}%")
    print(f"重复次数: {replicates}")
    print(f"基因组大小: {genome_size} bp")
    print(f"有效长度（正负链）: {genome_size * 2} bp = {genome_size * 2 / 1000000:.1f} Mb")
    print(f"碱基分配策略: G和C{'均等' if config['simulation']['gc_equal'] else '不等'}分配")
    print(f"A/T比例: {config['simulation']['at_ratio']}:1")
    print(f"可用CPU数量: {mp.cpu_count()}")
    print(f"开始时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print('='*60)
    
    # 使用进程池并行执行
    process_count = min(mp.cpu_count(), replicates, 60)  # 不超过SLURM分配的CPU数
    print(f"使用 {process_count} 个进程并行执行")
    
    with mp.Pool(processes=process_count) as pool:
        results = pool.map(generate_sequence_gc, tasks)
    
    # 保存结果
    valid_results = [r for r in results if r is not None]
    if valid_results:
        df = pd.DataFrame(valid_results)
        
        # 确保结果目录存在
        results_dir = Path(base_dir) / "results"
        results_dir.mkdir(exist_ok=True)
        
        # 保存详细结果
        results_file = results_dir / f"simulation_GC{gc_percent}_detailed.csv"
        df.to_csv(results_file, index=False)
        
        # 保存汇总结果
        summary_file = results_dir / f"simulation_GC{gc_percent}_summary.txt"
        
        # 计算统计信息
        success_rate = len(valid_results) / replicates * 100
        avg_density = df['Density_IM_per_Mb'].mean()
        std_density = df['Density_IM_per_Mb'].std()
        avg_gc = df['Actual_GC_content'].mean()
        avg_c = df['Actual_C_content'].mean()
        avg_g = df['Actual_G_content'].mean()
        
        summary_content = [
            f"GC含量: {gc_level} ({gc_percent}%)",
            f"重复次数: {replicates}",
            f"成功次数: {len(valid_results)} ({success_rate:.1f}%)",
            f"平均实际GC含量: {avg_gc:.4f}",
            f"平均实际C含量: {avg_c:.4f}",
            f"平均实际G含量: {avg_g:.4f}",
            f"平均正链C含量: {df['Pos_strand_C'].mean():.4f}",
            f"平均负链C含量: {df['Neg_strand_C'].mean():.4f}",
            f"平均i-motif密度: {avg_density:.2f} ± {std_density:.2f} IM/Mb",
            f"密度范围: {df['Density_IM_per_Mb'].min():.2f} - {df['Density_IM_per_Mb'].max():.2f} IM/Mb",
            f"平均运行时间: {df['Run_time_seconds'].mean():.1f} 秒",
            f"总i-motif数量: {df['iMotif_count'].sum():,}",
        ]
        
        if avg_density > 0:
            cv = (std_density / avg_density) * 100
            summary_content.append(f"变异系数: {cv:.1f}%")
        
        # 计算C含量与密度的相关性
        if len(df) > 1:
            corr_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'])
            corr_g = df['Actual_G_content'].corr(df['Density_IM_per_Mb'])
            corr_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'])
            summary_content.append(f"C含量与密度相关性: {corr_c:.3f}")
            summary_content.append(f"G含量与密度相关性: {corr_g:.3f}")
            summary_content.append(f"GC含量与密度相关性: {corr_gc:.3f}")
        
        # 写入汇总文件
        with open(summary_file, 'w') as f:
            f.write('\n'.join(summary_content))
        
        # 打印结果
        print(f"\n{'='*60}")
        print(f"GC{gc_percent}% 完成")
        print(f"成功: {len(valid_results)}/{replicates} ({success_rate:.1f}%)")
        print(f"平均GC含量: {avg_gc:.4f} (目标: {gc_level:.4f})")
        print(f"平均C含量: {avg_c:.4f}, 平均G含量: {avg_g:.4f}")
        print(f"平均密度: {avg_density:.2f} ± {std_density:.2f} IM/Mb")
        print(f"密度范围: {df['Density_IM_per_Mb'].min():.2f} - {df['Density_IM_per_Mb'].max():.2f} IM/Mb")
        print(f"详细结果保存到: {results_file}")
        print(f"汇总结果保存到: {summary_file}")
        
        if avg_density > 0:
            cv = (std_density / avg_density) * 100
            print(f"变异系数: {cv:.1f}%")
        
        if len(df) > 1:
            print(f"C含量与密度相关性: {corr_c:.3f}")
            print(f"GC含量与密度相关性: {corr_gc:.3f}")
        
        print(f"结束时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print('='*60)
        
        return df
    else:
        print("\n警告：所有任务都失败了！")
        return None

if __name__ == "__main__":
    main()
