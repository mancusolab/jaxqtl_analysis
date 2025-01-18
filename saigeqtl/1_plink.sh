# start with SNPs after filtering by HWE and MAF > 0.05
# filter for SNPs with MAC>20 and do LD pruning

cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp
mkdir -p ldprune_mac20
cd ldprune_mac20

plink=/project/nmancuso_8/elezhang/software/plink2_20240105/plink2

# get prune.in list
# here we can feed step1 with all pruned snps
for chr in {1..22};do
${plink} --bfile ../chr${chr} --mac 20 --indep-pairwise 250 100 0.3 --rm-dup exclude-mismatch --out ./chr${chr}
${plink} --bfile ../chr${chr} --extract ./chr${chr}.prune.in --make-bed --out chr${chr}
done
