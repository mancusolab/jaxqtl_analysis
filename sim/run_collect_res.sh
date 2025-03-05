#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem=15Gb
#SBATCH --array=1
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

Rscript combine_jaxqtl_saigeqtl.R --nsim 500 --params params_add_Va0

