#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=BCM_SIR
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=8:00:00
#SBATCH --mem=64000M
#SBATCH --array=71-80,101-110,161-170
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
source /etc/profile
module load R/4.2.0-openblas-rocky8

####### Run your script #########################
Rscript run_models.R $SLURM_ARRAY_TASK_ID
