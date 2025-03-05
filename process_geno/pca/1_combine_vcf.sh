#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=16Gb
#SBATCH --array=1
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

indir="/project/gazal_569/DATA/yazar2022/genotype/filter_vcf_r09"

bcftools concat ${indir}/chr{1..22}.dose.filtered.R2_0.9.vcf.gz -Oz -o ${indir}/allchr.dose.filtered.R2_0.9

