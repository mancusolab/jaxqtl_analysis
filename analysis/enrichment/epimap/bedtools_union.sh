# take union over multiple samples

# set working dir
# cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/epimap
cd /project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/epimap_onesample

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0

# params_file=params_bed
params_file=params_onesample_bed

start=1
stop=`wc --lines < ${params_file}`

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

bed=$1
celltype_out=$2
indir=$3

# cat files together (assume super group is in order)
if [ ! -e ./tmp/${celltype_out}.enhancer.bed ]; then
touch ./tmp/${celltype_out}.enhancer.bed
fi
cat ${indir}/${bed}.enhancer.bed >> ./tmp/${celltype_out}.enhancer.bed

if [ ! -e ./tmp/${celltype_out}.promoter.bed ]; then
touch ./tmp/${celltype_out}.promoter.bed
fi
cat ${indir}/${bed}.promoter.bed >> ./tmp/${celltype_out}.promoter.bed

done


# merge each type bed file
params_file=supergroup

start=1
stop=`wc --lines < ${params_file}`

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

bed=$1

# sort the bed files
sort -k1,1 -k2,2n ./tmp/${bed}.enhancer.bed > ./tmp/${bed}.enhancer.bed.sorted
sort -k1,1 -k2,2n ./tmp/${bed}.promoter.bed > ./tmp/${bed}.promoter.bed.sorted

bedtools merge -i ./tmp/${bed}.enhancer.bed.sorted > ${bed}.enhancer.merge.bed

bedtools merge -i ./tmp/${bed}.promoter.bed.sorted > ${bed}.promoter.merge.bed

sed -i 's/chr//' ${bed}.enhancer.merge.bed
sed -i 's/chr//' ${bed}.promoter.merge.bed

# # subtract promoter region from enhancer region
# bedtools subtract \
#     -a ${bed}.enhancer.merge.bed \
#     -b ${bed}.promoter.merge.bed  \
#     -A > \
#     ${bed}.enhancer.subtract.promoter.merge.bed
done