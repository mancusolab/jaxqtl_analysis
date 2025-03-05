#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=10Gb
#SBATCH --array=1-65
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich

# eg. start = 1, stop = 10
start=`python -c "print(1 + 5*int(int($idx-1)))"`
stop=$((start + 4))

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
#params=`sed "${IDX}q;d" params_sm_jaxqtl_tqtl`
#params=`sed "${IDX}q;d" params_sm_jaxqtl_nb_lm`
params=`sed "${IDX}q;d" params_sm_lm`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

chunk=$1
chr=$2
celltype=$3
software=$4
model=$5

geno_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp"

pip_prefix="/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/result_wald_label/${software}"
pip_suffix="${model}.L10.estvarFALSE.wald.tsv.gz"
annot_dir="${geno_dir}/ldsc"
genelist="./params_sm/${software}/${chunk}"
gene_meta="./gtf.82.bed"
sm_model="logit"
pip_cs="pip"

out_path=/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/enrich_sm/${software}

python enrich_annot_sm.py \
  --geno ${geno_dir}/chr${chr} \
  --chr ${chr} \
  --celltype ${celltype} \
  --model ${sm_model} \
  --pip-cs ${pip_cs} \
  --pip-prefix ${pip_prefix} \
  --pip-suffix ${pip_suffix} \
  --annot ${annot_dir}/chr${chr}.annot.gz \
  --genelist ${genelist} \
  --gene-meta ${gene_meta} \
  --out ${out_path}
done

