#!/bin/bash
#SBATCH --partition=cancergen
#SBATCH --time=96:00:00
#SBATCH --job-name=Lock
#SBATCH --output=./logs/Lock
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

module load conda_R/4.3.x

export R_LIBS_USER=/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17-bioc-release

Rscript B_Train_Lock.r
