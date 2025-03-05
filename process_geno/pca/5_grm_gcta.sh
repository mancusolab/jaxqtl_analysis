#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=40Gb
#SBATCH --array=1
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

gcta=/project/nmancuso_8/elezhang/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1

$gcta --bfile ./ldprune/allchr_pruned --make-grm --out ./grm/allchr

