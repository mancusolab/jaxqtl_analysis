### liftover bed files between hg19 and hg38
# arg1: input file
# arg2: 38 or 19 (old)
# arg3: 38 or 19 (new)
# arg4: out dir

liftover="/project/nmancuso_8/elezhang/software/liftOver"
 
infile=$1
old=$2
new=$3
outdir=$4

echo $infile
echo $old
echo $new
echo $outdir

prefix=${infile%".gz"}
prefix=${infile%".bed"}
prefix=${infile%".bed.gz"}

echo $prefix

chain="/project/nmancuso_8/data/liftover_chains/hg${old}ToHg${new}.over.chain.gz"

${liftover} \
    ${infile} \
    ${chain} \
    ${prefix}.hg${new}.bed \
    ${prefix}.unlifted.bed


