#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=20Gb
#SBATCH --array=1-22
#SBATCH --partition=main
#SBATCH --account=nmancuso_8
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl

celltype="NK"

# pull the idx'th line from file params and store in variable named params
# params=`sed "${idx}q;d" ./data/pheno/celltype16/metadata/genelist/${celltype}/params`

# split strings in params variable by space and store in bash arg variables
# set -- junk $params
# shift

chr=${idx}

geno="./data/geno_n982_info_08/hwe_maf0.05_snp/chr${chr}"
covar="./data/features/donor_features.all.6PC.NK.PC.bed"
pheno="./data/pheno/celltype16/${celltype}.tmm.bed.gz"

nperm=1000
window=500000
seed=${idx}
# robust="True"

out="./result/cis/celltype16_tensorqtl/${celltype}/chr${chr}"

python3 -m tensorqtl ${geno} ${pheno} ${out} \
    --covariates ${covar} \
    --window ${window} \
    --seed ${seed} \
    --mode cis_nominal

