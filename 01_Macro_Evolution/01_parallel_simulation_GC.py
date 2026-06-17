#!/usr/bin/env python3
"""
parallel_simulation_GC.py - Parallel script for GC content background simulation
Generate random sequences with different GC levels, run iM-Seeker to detect i-motif density
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

# Hardcoded configuration parameters
CONFIG = {
    'simulation': {
        'at_ratio': 1.0,  # A and T ratio (A:T = 1:1)
        'gc_equal': True,  # Whether G and C are equally distributed
        'im_seeker': {
            'model_dir': "/datapool/home/2023200496/niulk/my_project/genome/iM_Seeker_model"
        }
    }
}

def load_config():
    """Return hardcoded config"""
    return CONFIG

def generate_sequence_gc(args):
    """Task function for generating a single sequence - GC content version"""
    gc_level, replicate, genome_size, base_dir, config = args

    # Set different random seeds for each process
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
    log_content.append(f"Random seed: {random_seed}")
    log_content.append(f"Target GC content: {gc_level}")
    log_content.append(f"Genome size: {genome_size} bp")

    # 1. Generate random sequence - based on GC content and config
    total_gc = int(genome_size * gc_level)
    total_at = genome_size - total_gc

    # Get parameters from config
    at_ratio = config['simulation']['at_ratio']
    gc_equal = config['simulation']['gc_equal']

    # Calculate G and C counts
    if gc_equal:
        # G and C equally distributed
        num_g = total_gc // 2
        num_c = total_gc - num_g  # Handle odd case
    else:
        # Can modify for unequal G/C ratio, default to equal
        num_g = total_gc // 2
        num_c = total_gc - num_g

    # Calculate A and T counts
    total_ratio = at_ratio + 1.0  # A:T = at_ratio:1
    num_a = int(total_at * at_ratio / total_ratio)
    num_t = total_at - num_a

    # Create list containing all bases
    bases = (['G'] * num_g +
             ['C'] * num_c +
             ['A'] * num_a +
             ['T'] * num_t)

    # Fully shuffle the base order
    np.random.shuffle(bases)
    random_sequence = ''.join(bases)

    # Validate base content
    total_len = len(random_sequence)
    actual_g = random_sequence.count('G') / total_len
    actual_c = random_sequence.count('C') / total_len
    actual_a = random_sequence.count('A') / total_len
    actual_t = random_sequence.count('T') / total_len
    actual_gc = actual_g + actual_c

    # Calculate positive strand C content (i-motif formation site)
    pos_strand_c = actual_c
    # Calculate negative strand C content (positive strand G content)
    neg_strand_c = actual_g
    # Calculate positive strand G content
    pos_strand_g = actual_g

    log_content.append(f"Target GC content: {gc_level:.4f}")
    log_content.append(f"Actual GC content: {actual_gc:.4f}")
    log_content.append(f"Actual G content: {actual_g:.4f}")
    log_content.append(f"Actual C content: {actual_c:.4f}")
    log_content.append(f"Actual A content: {actual_a:.4f}")
    log_content.append(f"Actual T content: {actual_t:.4f}")
    log_content.append(f"Positive strand C content: {pos_strand_c:.4f}")
    log_content.append(f"Positive strand G content: {pos_strand_g:.4f}")
    log_content.append(f"Negative strand C content: {neg_strand_c:.4f}")
    log_content.append(f"Sequence length: {total_len} bp")
    log_content.append(f"G/C ratio: {actual_g/actual_c:.3f}" if actual_c > 0 else "G/C ratio: N/A (C=0)")
    log_content.append(f"A/T ratio: {actual_a/actual_t:.3f}" if actual_t > 0 else "A/T ratio: N/A (T=0)")

    # Save sequence
    with open(sequence_file, 'w') as f:
        f.write(f'>random_GC_{gc_level:.3f}_rep_{replicate:03d}_seed_{random_seed}\n')
        for i in range(0, len(random_sequence), 80):
            f.write(random_sequence[i:i+80] + '\n')

    # 2. Run im-seeker
    output_dir.mkdir(exist_ok=True)

    model_dir = config['simulation']['im_seeker']['model_dir']

    cmd = [
        'iM-Seeker.py',
        '--sequence', str(sequence_file),
        '--classification_model', f"{model_dir}/pickle_model_classification.pkl",
        '--regression_model', f"{model_dir}/pickle_model_regression.pkl",
        '--output_folder', str(output_dir)
    ]

    log_content.append(f"Running command: {' '.join(cmd)}")

    print(f"Running im-seeker: GC{gc_percent}_rep{replicate:03d}")

    start_time = datetime.now()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        run_time = (datetime.now() - start_time).total_seconds()

        log_content.append(f"Runtime: {run_time:.1f} seconds")
        log_content.append(f"Exit code: {result.returncode}")

        if result.returncode != 0:
            log_content.append(f"stderr: {result.stderr[:500]}")
            print(f"im-seeker failed GC{gc_percent}_rep{replicate:03d}")
            # Save log
            with open(log_file, 'w') as f:
                f.write('\n'.join(log_content))
            return None
        else:
            log_content.append("im-seeker ran successfully")
            print(f"im-seeker succeeded GC{gc_percent}_rep{replicate:03d}")

    except Exception as e:
        log_content.append(f"Runtime error: {str(e)}")
        print(f"Runtime error GC{gc_percent}_rep{replicate:03d}: {e}")
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_content))
        return None

    # 3. Parse results
    imotif_count, motif_details = parse_results(output_dir, work_dir)

    log_content.append(f"i-motif count: {imotif_count}")

    # Density calculation: based on 20Mb (both strands)
    effective_length = genome_size * 2  # 10Mb single strand -> 20Mb double strand
    density = (imotif_count / effective_length) * 1000000

    # Calculate consecutive C distribution (important metric)
    c_distribution = analyze_c_distribution(random_sequence)
    log_content.append(f"Consecutive C distribution: {json.dumps(c_distribution)}")

    # Calculate consecutive G distribution (study potential G interference)
    g_distribution = analyze_g_distribution(random_sequence)
    log_content.append(f"Consecutive G distribution: {json.dumps(g_distribution)}")

    log_content.append(f"i-motif density: {density:.2f} IM/Mb")

    # Save log
    with open(log_file, 'w') as f:
        f.write('\n'.join(log_content))

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
    """Parse im-seeker output results"""
    output_dir = Path(output_dir)

    # Find the final prediction file
    result_file = output_dir / "iM-seeker_final_prediction.txt"

    imotif_count = 0
    motif_details = []

    if result_file.exists():
        try:
            with open(result_file, 'r') as f:
                lines = f.readlines()

            # Parse i-motif details
            if len(lines) > 1:  # Has header and at least one result
                # Header line
                header = lines[0].strip().split('\t')

                # Parse each motif
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

            # Save detailed results
            if motif_details:
                details_file = work_dir / "motif_details.json"
                with open(details_file, 'w') as f:
                    json.dump(motif_details, f, indent=2)

        except Exception as e:
            print(f"Parse file error {result_file}: {e}")
            # Try simple count
            with open(result_file, 'r') as f:
                lines = f.readlines()
            if len(lines) > 0:
                imotif_count = len(lines) - 1  # Subtract header
    else:
        print(f"Result file not found: {result_file}")

    return imotif_count, motif_details

def analyze_c_distribution(sequence):
    """Analyze consecutive C distribution"""
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

    # Check trailing C streak
    if current_c_streak > 0:
        if current_c_streak >= 4:
            c_distribution['4+'] += 1
        elif current_c_streak in c_distribution:
            c_distribution[current_c_streak] += 1

    return c_distribution

def analyze_g_distribution(sequence):
    """Analyze consecutive G distribution"""
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

    # Check trailing G streak
    if current_g_streak > 0:
        if current_g_streak >= 4:
            g_distribution['4+'] += 1
        elif current_g_streak in g_distribution:
            g_distribution[current_g_streak] += 1

    return g_distribution

def main():
    if len(sys.argv) != 4:
        print("Usage: python parallel_simulation_GC.py <GC_level> <replicates> <genome_size>")
        print("Example: python parallel_simulation_GC.py 0.40 30 10000000")
        sys.exit(1)

    gc_level = float(sys.argv[1])
    replicates = int(sys.argv[2])
    genome_size = int(sys.argv[3])
    base_dir = "/datapool/home/2023200496/niulk/my_project/1555_simulation_study_extended_GC"

    # Load config
    config = load_config()

    # Prepare task parameters
    tasks = [(gc_level, rep, genome_size, base_dir, config)
             for rep in range(1, replicates + 1)]

    gc_percent = int(gc_level * 100)
    print(f"\n{'='*60}")
    print(f"Starting parallel processing GC{gc_percent}%")
    print(f"Replicates: {replicates}")
    print(f"Genome size: {genome_size} bp")
    print(f"Effective length (both strands): {genome_size * 2} bp = {genome_size * 2 / 1000000:.1f} Mb")
    print(f"Base allocation strategy: G and C {'equal' if config['simulation']['gc_equal'] else 'unequal'} allocation")
    print(f"A/T ratio: {config['simulation']['at_ratio']}:1")
    print(f"Available CPUs: {mp.cpu_count()}")
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print('='*60)

    # Execute in parallel using process pool
    process_count = min(mp.cpu_count(), replicates, 60)  # Do not exceed SLURM allocated CPUs
    print(f"Using {process_count} processes for parallel execution")

    with mp.Pool(processes=process_count) as pool:
        results = pool.map(generate_sequence_gc, tasks)

    # Save results
    valid_results = [r for r in results if r is not None]
    if valid_results:
        df = pd.DataFrame(valid_results)

        # Ensure results directory exists
        results_dir = Path(base_dir) / "results"
        results_dir.mkdir(exist_ok=True)

        # Save detailed results
        results_file = results_dir / f"simulation_GC{gc_percent}_detailed.csv"
        df.to_csv(results_file, index=False)

        # Save summary results
        summary_file = results_dir / f"simulation_GC{gc_percent}_summary.txt"

        # Calculate statistics
        success_rate = len(valid_results) / replicates * 100
        avg_density = df['Density_IM_per_Mb'].mean()
        std_density = df['Density_IM_per_Mb'].std()
        avg_gc = df['Actual_GC_content'].mean()
        avg_c = df['Actual_C_content'].mean()
        avg_g = df['Actual_G_content'].mean()

        summary_content = [
            f"GC content: {gc_level} ({gc_percent}%)",
            f"Replicates: {replicates}",
            f"Successful: {len(valid_results)} ({success_rate:.1f}%)",
            f"Mean actual GC content: {avg_gc:.4f}",
            f"Mean actual C content: {avg_c:.4f}",
            f"Mean actual G content: {avg_g:.4f}",
            f"Mean positive strand C content: {df['Pos_strand_C'].mean():.4f}",
            f"Mean negative strand C content: {df['Neg_strand_C'].mean():.4f}",
            f"Mean i-motif density: {avg_density:.2f} ± {std_density:.2f} IM/Mb",
            f"Density range: {df['Density_IM_per_Mb'].min():.2f} - {df['Density_IM_per_Mb'].max():.2f} IM/Mb",
            f"Mean runtime: {df['Run_time_seconds'].mean():.1f} sec",
            f"Total i-motif count: {df['iMotif_count'].sum():,}",
        ]

        if avg_density > 0:
            cv = (std_density / avg_density) * 100
            summary_content.append(f"Coefficient of variation: {cv:.1f}%")

        # Calculate C content vs density correlation
        if len(df) > 1:
            corr_c = df['Actual_C_content'].corr(df['Density_IM_per_Mb'])
            corr_g = df['Actual_G_content'].corr(df['Density_IM_per_Mb'])
            corr_gc = df['Actual_GC_content'].corr(df['Density_IM_per_Mb'])
            summary_content.append(f"C content vs density correlation: {corr_c:.3f}")
            summary_content.append(f"G content vs density correlation: {corr_g:.3f}")
            summary_content.append(f"GC content vs density correlation: {corr_gc:.3f}")

        # Write summary file
        with open(summary_file, 'w') as f:
            f.write('\n'.join(summary_content))

        # Print results
        print(f"\n{'='*60}")
        print(f"GC{gc_percent}% Complete")
        print(f"Successful: {len(valid_results)}/{replicates} ({success_rate:.1f}%)")
        print(f"Mean GC content: {avg_gc:.4f} (target: {gc_level:.4f})")
        print(f"Mean C content: {avg_c:.4f}, Mean G content: {avg_g:.4f}")
        print(f"Mean density: {avg_density:.2f} ± {std_density:.2f} IM/Mb")
        print(f"Density range: {df['Density_IM_per_Mb'].min():.2f} - {df['Density_IM_per_Mb'].max():.2f} IM/Mb")
        print(f"Detailed results saved to: {results_file}")
        print(f"Summary results saved to: {summary_file}")

        if avg_density > 0:
            cv = (std_density / avg_density) * 100
            print(f"Coefficient of variation: {cv:.1f}%")

        if len(df) > 1:
            print(f"C content vs density correlation: {corr_c:.3f}")
            print(f"GC content vs density correlation: {corr_gc:.3f}")

        print(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print('='*60)

        return df
    else:
        print("\nWarning: All tasks failed!")
        return None

if __name__ == "__main__":
    main()
