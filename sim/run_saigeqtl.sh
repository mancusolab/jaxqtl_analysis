#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=20Gb
#SBATCH --array=1-375
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl/code/sim_sc

params_file="params_add_Va0"

step1_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step1_fitNULLGLMM_qtl.R
step2_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step2_tests_qtl.R
step3_script=/project/nmancuso_8/elezhang/software/qtl/extdata/step3_gene_pvalue_qtl.R

# eg. start = 1, stop = 10
start=`python -c "print(1 + 300*int(int($idx-1)))"`
stop=$((start + 299))

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./params/${params_file}_sqtl`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

maf=${2}
nobs=${6}
sim_idx=${7}
pheno_idx=${8}

geno_file="./input/geno/chr1_${maf}_n${nobs}"

step1prefix="./output/${params_file}/sim${sim_idx}.pheno${pheno_idx}_sqtl"
step2prefix="${step1prefix}_cis"

chr=1
regionFile="./input/geno/chr${chr}_${maf}_n${nobs}_cis_region.txt"

pheno_file="./tmp/sim${sim_idx}.pheno${pheno_idx}.tsv.gz"
gene="gene"

# 1. fit null model
timeout 500s Rscript ${step1_script} \
        --useSparseGRMtoFitNULL=FALSE  \
        --useGRMtoFitNULL=FALSE \
        --phenoFile=${pheno_file} \
        --phenoCol=${gene}     \
        --covarColList=sex,age   \
        --sampleCovarColList=sex,age   \
        --sampleIDColinphenoFile=iid \
        --offsetCol=log_offset  \
        --traitType=count \
        --outputPrefix=${step1prefix} \
        --skipVarianceRatioEstimation=FALSE  \
        --isRemoveZerosinPheno=FALSE \
        --isCovariateOffset=FALSE  \
        --isCovariateTransform=FALSE  \
        --skipModelFitting=FALSE  \
        --tol=0.001   \
        --maxiter=200   \
        --plinkFile=${geno_file}      \
        --SPAcutoff=10000  \
        --IsOverwriteVarianceRatioFile=TRUE

# 2. cis test
Rscript ${step2_script}    \
        --bedFile=${geno_file}.bed      \
        --bimFile=${geno_file}.bim      \
        --famFile=${geno_file}.fam      \
        --SAIGEOutputFile=${step2prefix}     \
        --chrom=${chr}       \
        --minMAF=0 \
        --minMAC=1 \
        --LOCO=FALSE    \
        --GMMATmodelFile=${step1prefix}.rda     \
        --SPAcutoff=10000 \
        --varianceRatioFile=${step1prefix}.varianceRatio.txt    \
        --rangestoIncludeFile=${regionFile}     \
        --markers_per_chunk=500

# 3. get causal SNP pval
Rscript 2_get_pval_sqtl.R \
  --idx ${pheno_idx} \
  --gene gene${pheno_idx} \
  --prefix ${step1prefix} \
  -o ${step1prefix}

# clean up
done
