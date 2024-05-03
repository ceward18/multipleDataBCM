#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=SIHRD_nyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --array=1-12
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
module load openblas/0.3.5_gcc8.2.0_multiarch
module load R/4.1.0

####### Run your script #########################
Rscript run_models.R $SLURM_ARRAY_TASK_ID
