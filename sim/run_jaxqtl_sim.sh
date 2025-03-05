#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=20Gb
#SBATCH --array=1-75
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

params_file="params_add_Va0"

# eg. start = 1, stop = 10
start=`python -c "print(1 + 3*int(int($idx-1)))"`
stop=$((start + 2))

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./params/${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

sim_script="sim_sc.py"

CT=${1}
maf=${2}

num_sim=500

beta0=${3}
Va=${4}
Vre=${5}
nobs=${6}

geno="./input/geno/chr1_${maf}_n${nobs}"
libsize="./input/pheno/${CT}.features.tsv.gz"

model="poisson" #"NB"

out="./output/${params_file}/sim${IDX}"
out_sc="./tmp/sim${IDX}"

python ${sim_script} \
  --geno ${geno} \
  --CT ${CT} \
  --libsize-path ${libsize} \
  --nobs ${nobs} \
  --maf ${maf} \
  --model ${model} \
  --beta0 ${beta0} \
  --Va ${Va} \
  --Vre ${Vre} \
  --num-sim ${num_sim} \
  --seed ${IDX} \
  --write-sc \
  --out-sc ${out_sc} \
  --out ${out} 
done

