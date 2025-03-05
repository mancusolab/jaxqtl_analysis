#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=10Gb
#SBATCH --array=1-85
#SBATCH --partition=conti
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd /project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich

# eg. start = 1, stop = 10
start=`python -c "print(1 + 1000*int(int($idx-1)))"`
stop=$((start + 999))

module purge
module load gcc/11.3.0
module load bedtools2/2.30.0

plink2='/project/nmancuso_8/elezhang/software/plink2_20240105/plink2'

gen_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp"
annot_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed"
out_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/scent_0616"
# out_dir="/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/scent_new_fixalpha/bulk/"
tmp="${out_dir}/tmp"

# params_file="egene_lookup_new_fixalpha.tsv"
params_file="egene_lookup_0616.tsv"

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ${params_file}`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

# read in tss of gene
gene_id=$1
gene_symbol=$2
chr=$3
tss=$5
window=500000
ref_cell=$7

tss_start=`python -c "print(int(${tss}-${window}))"`
tss_start=`python -c "print(${tss_start} if ${tss_start} > 0 else 1)"`
tss_end=`python -c "print(int(${tss}+${window}))"`

bim=${tmp}/${gene_id}

if [ -e ${bim}.bim ]; then
echo "bim exist"
else
${plink2} --bfile ${gen_dir}/chr${chr} \
  --chr ${chr} \
  --from-bp ${tss_start} \
  --to-bp ${tss_end} \
  --make-just-bim \
  --out ${bim}
fi

if [ -e ${bim}.bim ]; then
# count how many rows in bim file
num_cis=`wc --lines < ${bim}.bim`

# convert .bim to .bed
awk '{OFS = "\t" ; print $1,$4-1,$4}' ${bim}.bim > ${bim}.bed

# do intersection

# cut out line that contains enhancer
ref_cell=${ref_cell}.subtract.PLS.ENCODE
annot_file=${annot_dir}/${ref_cell}

annot_gene=${tmp}/${ref_cell}.${gene_id}

awk -v var1="$gene_symbol" '$4 == var1' $annot_file > ${annot_gene}

out_name=${out_dir}/${gene_id}_${ref_cell}

# if not empty
if [ -s "${annot_gene}" ]; then
intersect=${bim}.intersect.${ref_cell}.${gene_id}
bedtools intersect -wb -wa \
    -a ${annot_gene} \
    -b ${bim}.bed > \
    ${intersect}
else
num_match=0
echo -e "${chr}\t${gene_id}\t${gene_symbol}\t${tss_start}\t${tss_end}\t${num_cis}\t${num_match}" > ${out_name}
continue
fi

# count lines
num_match=`wc --lines < ${intersect}`

echo -e "${chr}\t${gene_id}\t${gene_symbol}\t${tss_start}\t${tss_end}\t${num_cis}\t${num_match}" > ${out_name}
else
continue
fi
done



