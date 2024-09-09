#!/bin/bash

#SBATCH 
#SBATCH 
#SBATCH 
#SBATCH --mem=1GB 
#SBATCH           
#SBATCH --time=7-1:00:00
#SBATCH --partition=cancergen 
#SBATCH --job-name=ctcf       
#SBATCH                       
#SBATCH --cpus-per-task=2     
#SBATCH                       
#SBATCH --output=logs/slurm-%x-%J-%a.output.txt      
#SBATCH --error=logs/slurm-%x-%J-%a.error.txt        

# module load conda_R/4.3
module load R/4.3
# export R_LIBS_USER=/dcl01/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.10-bioc-release/
export R_LIBS_USER=/dcs04/scharpf/data/pipeline-hub/pipeline-lib/R-libs/3.17.bioc-release

path=/dcs07/scharpf/data/nniknafs/delfi-brain/data/manuscript/analysis/tfbs-expr/submission-23-July-2024
cd $path

echo ${SLURM_ARRAY_TASK_ID}
Rscript scripts/05_brain_null.r -s  ${SLURM_ARRAY_TASK_ID}

echo 'Done'

