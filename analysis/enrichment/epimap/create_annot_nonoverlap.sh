# create Epimap annotation file without any overlap between cell types

cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/epimap/specific

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0

# for other cell types
for reg in promoter enhancer;do
for celltype in Bcells CD4_Tcells CD8_Tcells DC Mono_C NK; do
find . -maxdepth 1 -type f -name "*.${reg}.merge.bed" ! -name "${celltype}.${reg}.merge.bed" \
  -exec cat {} + > combined.bed
sort -k1,1 -k2,2n combined.bed > sorted_combined.bed
bedtools subtract -a "${celltype}.${reg}.merge.bed" -b sorted_combined.bed -A > "${celltype}.${reg}.sp.bed"
rm -rf combined.bed sorted_combined.bed
done
done

# intersect with specific eQTLs SNP positions
> all_intersect_epimap
for FILE in *_mash.bed;do
  for ANNOT in *.merge.bed;do
     filename=${ANNOT%.merge.bed}
     celltype=$(echo "$filename" | cut -d'.' -f1)
     reg=$(echo "$filename" | cut -d'.' -f2)
     snplist=${FILE%_mash.bed}
     
     bedtools intersect -wa -wb -a $FILE -b $ANNOT > tmp
     awk -v var1="$celltype" -v var2="$reg" -v var3="$snplist" 'BEGIN {OFS="\t"}  {print $0, var1, var2, var3}' tmp >> all_intersect_epimap
     rm -rf tmp
done
done
