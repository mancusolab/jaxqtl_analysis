#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=30Gb
#SBATCH --array=1
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu
#SBATCH --output=CD4_NC.out

source /home1/zzhang39/renv_saigeqtl.sh
cd /project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl

step1_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step1_fitNULLGLMM_qtl.R
step2_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step2_tests_qtl.R
step3_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step3_gene_pvalue_qtl.R

geno_dir="../data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldprune_mac20"

celltype="CD4_NC"

for IDX in `seq 1 50`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./input/params_50`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

gene=${1}
chr=${2}
tss=${3}

# 1. run null model
# SNPs are filtered by maf > 0.05, so guranteed MAC>=20
# ignore GRM, use qr transformation on covariates,
Rscript ${step1_script} \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=./input/pheno/${celltype}.pheno.tsv.gz 	\
        --phenoCol=${gene}      \
        --covarColList=sex,age,PC1,PC2,PC3,PC4,PC5,PC6,ct_pc1,ct_pc2    \
        --sampleCovarColList=sex,age,PC1,PC2,PC3,PC4,PC5,PC6      \
        --sampleIDColinphenoFile=individual \
        --offsetCol=log_offset  \
        --traitType=count \
        --outputPrefix=./output/${celltype}_${gene} \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=TRUE  \
        --skipModelFitting=FALSE  \
        --tol=0.00001  \
        --maxiter=1000   \
        --plinkFile=${geno_dir}/chr${chr}      \
        --IsOverwriteVarianceRatioFile=TRUE
        
# 2. score test

# cis region file, window size = 500000 (total 1Mb)
regionFile=./input/region/${gene}_cis_region.txt

tss_start=`python -c "print(int(${tss}-500000))"`
tss_start=`python -c "print(${tss_start} if ${tss_start} > 0 else 1)"`
tss_end=`python -c "print(int(${tss}+500000))"`

echo -e "${chr}\t${tss_start}\t${tss_end}" > ${regionFile}

step1prefix=./output/${celltype}_${gene}
step2prefix=./output/${celltype}_${gene}_cis

Rscript ${step2_script}    \
        --bedFile=${geno_dir}/chr${chr}.bed      \
        --bimFile=${geno_dir}/chr${chr}.bim      \
        --famFile=${geno_dir}/chr${chr}.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=${chr}       \
        --minMAF=0 \
        --minMAC=20 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=2 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=10000

Rscript ${step3_script} \
        --assocFile=${step2prefix}        \
        --geneName=${gene}       \
        --genePval_outputFile=./output/${celltype}_${gene}_cis_genePval

done
