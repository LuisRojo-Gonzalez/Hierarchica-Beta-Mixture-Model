#!/bin/bash
##---------------SLURM Parameters - NLHPC ----------------
#SBATCH -J Mixture
#SBATCH -p slims
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=1000
#SBATCH --array=[1-36]%18
#SBATCH -o Output/Mixture_%a.out
#SBATCH -e Fail/Mixture_%a.err
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=luis.rojo.g@usach.cl
#SBATCH --mail-type=NONE

#-----------------Toolchain---------------------------
# ----------------Modules----------------------------
ml  R/4.0.0

# ----------------Command--------------------------
Rscript --vanilla MixModel.R $SLURM_ARRAY_TASK_ID
