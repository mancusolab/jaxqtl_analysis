#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=10Gb
#SBATCH --array=1-20
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

# create annotation file using bed file

cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldsc

# eg. start = 1, stop = 10
start=`python -c "print(1 + 6*int(int($idx-1)))"`
stop=$((start + 5))

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" params`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

annot_dir=${1}
# no bed suffix
# annotation file need format chr as "chr1", wheareas bim file doesn't
annot_name=${2}

for chr in {1..22};do

annot_out=${annot_name}.${chr}.annot.gz

# copy over original bed file to tmp
# add "chr" to annotation file in place
if [ -e ${annot_dir}/${annot_name}.bed ]; then
cp ${annot_dir}/${annot_name}.bed ./tmp
else
continue
fi

find_chr=`sed -n '1{/^chr/p};q' ./tmp/${annot_name}.bed | wc --lines`

if [ $find_chr -lt 1 ];then
sed -i 's/^/chr/' ./tmp/${annot_name}.bed
fi

python make_annot.py \
 --windowsize 0 \
 --bed-file ./tmp/${annot_name}.bed \
 --bimfile ../chr${chr}.bim \
 --annot-file ${annot_out}

done

# clean up
rm -rf tmp/${annot_name}.bed
done
