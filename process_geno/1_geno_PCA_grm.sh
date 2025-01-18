# use the same dataset that check relatedness to calculate PC

# steps overview:
# 1. Filter variants: INFO > 0.9, maf > 0.01, biallelic SNPs, HWE p > 1e-6, LD pruned r2=0.3
# 2. Use GCTA to generate GRM matrix


# 1. concat all chromosome files and out put bed file
bcftools concat chr{1..22}.dose.filtered.R2_0.9.vcf.gz -o  allchr.dose.filtered.R2_0.9.vcf.gz
plink2 --vcf allchr.dose.filtered.R2_0.9.vcf.gz --make-bed --out /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n981_info_09/allchr

# 2. convert vcf to bed file (allchr)
# Filter variants for cis-eQTL analysis
# INFO > 0.9 (done), MAF > 0.01, HWE p < 1e-6
cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n981_info_09

# calculate allele frequency
MAF=0.01
plink2 --bfile allchr --maf ${MAF} --snps-only --max-alleles 2 --make-bed  --out ./maf_0.01_snp/allchr

# calculate HWE
plink2 --bfile ./maf_0.01_snp/allchr --hardy --out ./hwe_0.01_snp/allchr

# create snplist with HWE p < 1e-6 in R to make hwe_nopass.snplist to exclude; 
# use fread() to read *.hardy if no ID found (data type is not right)

# filter out HWE < 1e-6 and make bed (use this to calculate PC)
plink2 --bfile ./maf_0.01_snp/allchr --exclude ./hwe_0.01_snp/hwe_nopass.snplist --make-bed --out ./hwe_0.01_snp/allchr_pass_hwe

# LD pruning r2=0.3
plink2 --bfile ./hwe_0.01_snp/allchr_pass_hwe --indep-pairwise 250 100 0.3 --rm-dup exclude-mismatch --out ./ldprune/allchr

plink2 --bfile ./hwe_0.01_snp/allchr_pass_hwe --extract ./ldprune/allchr.prune.in --make-bed --out ./ldprune/allchr_pruned

# calculate PC
plink2 --bfile ./ldprune/allchr_pruned --pca 10 --out ./PCA/allchr_pruned

# calculate grm
gcta=/project/nmancuso_8/elezhang/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
$gcta --bfile ./hwe_0.01_snp/allchr_pass_hwe --make-grm --out ./allchr_pass_hwe

gcta=/project/nmancuso_8/elezhang/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
$gcta --bfile ./ldprune/allchr_pruned --make-grm --out ./grm/allchr

