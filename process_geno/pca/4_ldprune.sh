plink2='/project/nmancuso_8/elezhang/software/plink2'

${plink2} --bfile ./hwe_0.01_snp/allchr --indep-pairwise 250 100 0.3 --out ./ldprune/allchr --rm-dup exclude-mismatch

${plink2} --bfile ./hwe_0.01_snp/allchr --extract ./ldprune/allchr.prune.in --make-bed --out ./ldprune/allchr_pruned
