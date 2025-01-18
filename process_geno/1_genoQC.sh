# Filter SNP variants for cis-eQTL analysis
# INFO > 0.8 (done), keep SNPs, MAF > 0.05, HWE p < 1e-6

# genotype directory: genotype with INFO>0.8 by chr
cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n981
outdir="hwe_maf0.05_snp"

# 1. Remove variants that not pass HWE test < 1e-6
for i in {1..22}
do
plink --bfile ./chr${i} --hwe 1e-6 --out ./${outdir}/hwe/chr${i} --make-bed
done


# 2. Filter SNP variants with MAF > 0.05
MAF=0.05
for i in {1..22}
do
plink2 --bfile ./${outdir}/hwe/chr${i} --maf ${MAF} --snps-only --max-alleles 2 --make-bed  --out ./${outdir}/chr${i}
done
