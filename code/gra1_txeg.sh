#!/bin/bash
#SBATCH --job-name="txeg_st"
#SBATCH --time=10:00:00   # walltime 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32   # 16 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # partition type
#SBATCH --mail-user=briana.whitaker@usda.gov
#SBATCH --mail-type=BEGIN,END,FAIL

date
module load r/3.6.1
Rscript $HOME/txeg/txeg_sourcetracker_gradient_noMonoc_1.R
date
