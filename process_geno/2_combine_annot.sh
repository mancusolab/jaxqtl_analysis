#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=10Gb
#SBATCH --array=1-21
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

# for each chromosome, combine annotation

cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldsc

chr=${idx}
output="chr${chr}.annot"

# pull first file
params=`sed "1q;d" annot_list`
echo "Running instance 1 with params: ${params}"
set -- junk $params
shift
first_annot=$1

zcat ${first_annot}.${chr}.annot.gz > ${output} 

awk -v var1="$first_annot" 'BEGIN{FS=OFS="\t"} NR == 1 {$1 = var1; print; next} {print}' ${output} > ${output}.tmp && mv ${output}.tmp ${output} 


start=2
stop=`wc --lines < annot_list`

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" annot_list`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

# no bed suffix
# annotation file need format chr as "chr1", wheareas bim file doesn't
annot_name=${1}

awk 'BEGIN{FS=OFS="\t"}
     NR==FNR {
         if (NR==1) {names_header="First Name"; next}
         name[NR]=$1
     }
     NR!=FNR {
         if (FNR==1) {print names_header, "Years Old"; next}
         print name[FNR+1], $1
     }' names.txt ages.txt > combined.txt


done
