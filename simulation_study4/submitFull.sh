#!/bin/bash
####### Reserve computing resources #############
#SBATCH --job-name=BCM_full
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time=24:00:00
#SBATCH --mem 64000M
#SBATCH --array=16-30
#SBATCH --output=./out/Array.%A_%a.out
#SBATCH --error=./err/Array.%A_%a.error

####### Set environment variables ###############
module load gcc/8.2.0
module load openblas/0.3.5_gcc8.2.0_multiarch
module load java/9.0.4
module load lib/hdf5/1.10.0.1
export PATH=/home/ward-c/R/bin:$PATH
export R_LIBS=/home/ward-c/R/lib64:$R_LIBS

####### Run your script #########################
Rscript run_models.R $SLURM_ARRAY_TASK_ID
