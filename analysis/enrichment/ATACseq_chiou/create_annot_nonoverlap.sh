# create ATACseq annotation file without any overlap between cell types

cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/chiou

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0


# do B cell
celltype="B_all"
find . -maxdepth 1 -type f -name "*.bed" ! -name "ATACseq_B_*.bed" -exec cat {} + > combined.bed
sort -k1,1 -k2,2n combined.bed > sorted_combined.bed
bedtools subtract -a "ATACseq_${celltype}.bed" -b sorted_combined.bed -A > "./specific/ATACseq_${celltype}_sp.bed"
rm -rf combined.bed sorted_combined.bed

for celltype in B_IN B_Mem;do
find . -maxdepth 1 -type f -name "*.bed" ! -name "ATACseq_B_all.bed" ! -name "ATACseq_${celltype}.bed" \
  -exec cat {} + > combined.bed
sort -k1,1 -k2,2n combined.bed > sorted_combined.bed
bedtools subtract -a "ATACseq_${celltype}.bed" -b sorted_combined.bed -A > "./specific/ATACseq_${celltype}_sp.bed"
rm -rf combined.bed sorted_combined.bed
done

# for other cell types
for celltype in CD4_all CD8_all DC Mono_C Mono_NC NK NK_R; do
find . -maxdepth 1 -type f -name "*.bed" ! -name "ATACseq_${celltype}.bed" \
  -exec cat {} + > combined.bed
sort -k1,1 -k2,2n combined.bed > sorted_combined.bed
bedtools subtract -a "ATACseq_${celltype}.bed" -b sorted_combined.bed -A > "./specific/ATACseq_${celltype}_sp.bed"
rm -rf combined.bed sorted_combined.bed
done

# intersect with specific eQTLs SNP positions
> all_intersect_atac
for FILE in *_mash.bed;do
  for ANNOT in ATACseq_*;do
     #celltype=${ANNOT%.bed}
     #snp=$(echo "$filename" | cut -d'.' -f1)
     #celltype=$(echo "$filename" | cut -d'.' -f2)
     celltype=${ANNOT%.bed}
     snplist=${FILE%_mash.bed}
     
     bedtools intersect -wa -wb -a $FILE -b $ANNOT > tmp
     awk -v var1="$snplist" -v var2="$celltype" 'BEGIN {OFS="\t"}  {print $0, var2, var1}' tmp >> all_intersect_atac
     rm -rf tmp
done
done

# # concatenate results together 
# find . -maxdepth 1 -type f -name "*_mash*" ! -name "*.bed" -exec cat {} + > all_intersect.tsv
