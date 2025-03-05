plink2='/project/nmancuso_8/elezhang/software/plink2'

${plink2} --bfile ./allchr --hwe 1e-6 --make-bed --maf 0.01 --out ./hwe_0.01_snp/allchr --max-alleles 2 --snps-only

