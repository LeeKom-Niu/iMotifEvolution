#!/bin/bash

#SBATCH --job-name=plot_mam
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --time=30-00:00:00
#SBATCH --partition=life-zhanghk
#SBATCH --output=plot_mam_%j.out
#SBATCH --error=plot_mam_%j.err

echo "开始时间: $(date)"

python plot_primate_TSS_TES.py > plot_primate_TSS_TES.log 2>&1

echo "任务完成时间: $(date)"

