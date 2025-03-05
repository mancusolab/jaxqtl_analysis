#!/bin/bash

plink2='/project/nmancuso_8/elezhang/software/plink2'

${plink2} --keep IID.txt --make-bed  --out ./allchr --vcf /project/gazal_569/DATA/yazar2022/genotype/filter_vcf_r09/allchr.dose.filtered.R2_0.9
