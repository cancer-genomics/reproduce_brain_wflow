#!/bin/bash
#SBATCH --partition=cancergen
#SBATCH --time=96:00:00
#SBATCH --job-name=LOO
#SBATCH --output=./logs/LOO-%A-%a.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-484

module load conda_R/4.3.x

#fold="fold"$SGE_TASK_ID

export R_LIBS_USER=/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17-bioc-release

Rscript A_Ensemble_ARTEMIS1_delfi_LOO.r $SLURM_ARRAY_TASK_ID
