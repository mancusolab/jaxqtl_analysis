#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=20Gb
#SBATCH --array=1-54
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl

# pull the idx'th line from file params and store in variable named params
params=`sed "${idx}q;d" ./result/pval_beta_acat/code/params`

# split strings in params variable by space and store in bash arg variables
set -- junk $params
shift

celltype=$1
chr=$2
genelist=${3}
test_method="score"
mode="cis_acat"

geno="./data/geno_n982_info_08_updateid/hwe_maf0.05_snp/chr${chr}"
covar="./data/features_new/donor_features.all.6PC.tsv"
pheno="./data/pheno/celltype16_new/${celltype}.bed.gz"
chunk="./result/pval_beta_acat/input/${genelist}"

model="NB"
nperm=1000

out="./result/pval_beta_acat/output/${genelist}_jaxqtl_nb_pacat"

jaxqtl \
 --geno ${geno} \
 --covar ${covar} \
 --pheno ${pheno} \
 --model ${model} \
 --mode ${mode} \
 --genelist ${chunk} \
 --test-method ${test_method} \
 --nperm ${nperm} \
 --addpc 2 \
 --standardize \
 --out ${out}

