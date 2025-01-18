#!/bin/bash

# find eQTLs in peaks region of annotation

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0

ref_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed"

params_file=params_scent
start=1
stop=`wc --lines < ${params_file}`

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

query_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/${4}/all_celltype"

celltype=$1
ref=${2}.subtract.PLS.ENCODE
query=$3

bedtools intersect -wb -wa \
    -a ${ref_dir}/${ref} \
    -b ${query_dir}/${query}.bed > \
    ${query_dir}/${query}.intersect.${ref}
done
