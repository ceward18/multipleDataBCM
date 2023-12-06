#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=SIR_montreal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=16:00:00
#SBATCH --mem=248G
#SBATCH --array=5-10
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
module load openblas/0.3.5_gcc8.2.0_multiarch
module load R

####### Run your script #########################
Rscript run_models.R $SLURM_ARRAY_TASK_ID
