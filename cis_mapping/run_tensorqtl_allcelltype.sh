#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=15Gb
#SBATCH --array=1-66
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl

# eg. start = 1, stop = 10
start=`python -c "print(1 + 5*int(int($idx-1)))"`
stop=$((start + 4))

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./data/pheno/celltype16_new/metadata/genelist/params_all_tqtl`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

celltype=${1}
chr=${2}
mode="cis_nominal" #"cis"

geno="./data/geno_n982_info_08_updateid/hwe_maf0.05_snp/chr${chr}"
covar="./data/features_new/donor_features.all.6PC.${celltype}.PC.bed"
pheno="./data/pheno/celltype16_new/${celltype}.tmm.bed.gz"

nperm=1000
window=500000
seed=${idx}

out="./result/cis/celltype16_tensorqtl_new/${celltype}/chr${chr}"

python3 -m tensorqtl ${geno} ${pheno} ${out} \
    --covariates ${covar} \
    --window ${window} \
    --seed ${seed} \
    --mode ${mode} \
    --permutations ${nperm}
done

