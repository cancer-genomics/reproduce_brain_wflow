#!/bin/bash

#SBATCH --mem=1GB
#SBATCH --time=7-1:00:00
#SBATCH --job-name=iu
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm-%x-%J-%a.output.txt
#SBATCH --error=logs/slurm-%x-%J-%a.error.txt


module load conda_R/4.4

path=intended_use/scripts
cd $path

echo ${SLURM_ARRAY_TASK_ID}
time Rscript run_simulation_thread.r ${SLURM_ARRAY_TASK_ID} ../results/runs/results_${SLURM_ARRAY_TASK_ID}.txt
echo 'Done'
