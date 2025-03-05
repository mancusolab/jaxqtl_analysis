#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=15Gb
#SBATCH --array=1-1334
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
params=`sed "${IDX}q;d" ./data/pheno/celltype16_new/metadata/genelist/params_all_celltype`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

celltype=${1}
chr=${2}
chunk_file=${3}
test_method="score"
mode="cis"

geno="./data/geno_n982_info_08_updateid/hwe_maf0.05_snp/chr${chr}"
covar="./data/features_new/donor_features.all.6PC.tsv"
pheno="./data/pheno/celltype16_new/${celltype}.bed.gz"
genelist="./data/pheno/celltype16_new/metadata/genelist/${celltype}/chr${chr}/${chunk_file}"

model="poisson"
nperm=1000
# robust="True"

out="./result/cis/celltype16_new/${celltype}/chr${chr}/${chunk_file}.newbeta.pois"

jaxqtl \
 -geno ${geno} \
 -covar ${covar} \
 -pheno ${pheno} \
 -model ${model} \
 -mode ${mode} \
 -genelist ${genelist} \
 -test-method ${test_method} \
 -nperm ${nperm} \
 -addpc 2 \
 --standardize \
 -out ${out}

done

