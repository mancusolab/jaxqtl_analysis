#!/bin/bash

### remove PLS from annotation bed file

# PLS from ENCODE (hg38)
cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0

# find overlap
ref="PLS.ENCODE"
bed1="Bcell"
bed2="Myeloid"
bed3="Tnk"

# subtract any overlap of B from A
for bed in $bed1 $bed2 $bed3;do
bedtools subtract \
  -a ${bed}.gene_peaks.merged.hg19.bed \
  -b ${ref}.hg19.bed \
  -A > \
  ${bed}.subtract.${ref}
done

# format for qtltools (remove "chr" in the first column)
sed -i 's/chr//' PLS.ENCODE.hg19.bed 

done

# use R script to format


# merge all immune cells
# intersect multiple bed files
bedtools intersect -wa -wb \
    -a ${bed1}.gene_peaks.merged.hg19.bed \
    -b ${bed2}.gene_peaks.merged.hg19.bed  > \
    ${bed1}.intersect.${bed2}
    
# use R to filter genes are the same

# sort bed file before bedtools merge
for bed in $bed1 $bed2 $bed3;do
sort -k1,1 -k2,2n ${bed}.subtract.PLS.ENCODE > ${bed}.subtract.PLS.ENCODE
done

# merge 3 bed files for each gene

stop=`wc --lines < ./tmp/scent_immune_genelist`
start=1

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./tmp/scent_immune_genelist`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

gene=$1

mkdir ./tmp/$gene
awk -v var1="$gene" '$4 == var1' Bcell.subtract.PLS.ENCODE > ./tmp/$gene/Bcell_${gene}
awk -v var1="$gene" '$4 == var1' Myeloid.subtract.PLS.ENCODE > ./tmp/$gene/Myeloid_${gene}
awk -v var1="$gene" '$4 == var1' Tnk.subtract.PLS.ENCODE > ./tmp/$gene/Tnk_${gene}

cat tmp/${gene}/*_${gene} > tmp/$gene/${gene}.combine
sort -k1,1 -k2,2n ./tmp/$gene/${gene}.combine > ./tmp/$gene/${gene}.combine.sorted

bedtools merge -i ./tmp/$gene/${gene}.combine.sorted > ./tmp/$gene/${gene}.merge

awk -v var1="$gene" '{OFS = "\t" ; print $1,$2,$3,var1}' ./tmp/$gene/${gene}.merge > ./tmp/${gene}.merge

rm -rf ./tmp/$gene
rmdir ./tmp/$gene

done

# remove regions in all immune with any overlap with B
bedtools subtract \
    -a allimmune.intersect \
    -b ${ref}.hg19.bed  \
    -A > \
    allimmune.subtract.PLS.ENCODE
