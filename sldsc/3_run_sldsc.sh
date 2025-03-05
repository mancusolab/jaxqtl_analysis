#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=15Gb
#SBATCH --array=1-91
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

# eg. start = 1, stop = 10
start=`python -c "print(1 + 20*int(int($idx-1)))"`
stop=$((start + 19))

outdir="sldsc_unioncells_allgtex_onek1k"

for IDX in `seq ${start} ${stop}`; do
# run w/e you need to run
params=`sed "${IDX}q;d" ./params`
echo "Running instance ${IDX} with params: ${params}"
set -- junk $params
shift

CT=${1}
ANNOT=${2}
TRAIT=${3}
sumstat_dir=${4}

ldsc="/project/nmancuso_8/elezhang/software/ldsc/ldsc.py"

OPTIONS="--overlap-annot --print-coefficients";
REF_DIR="/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3"
BASELINE="$REF_DIR/baseline_v1.2/baseline.";
BASELINELD="$REF_DIR/baselineLD_v2.2/baselineLD.";
FREQ="--frqfile-chr $REF_DIR/plink_files/1000G.EUR.QC.";
WEIGHTS="--w-ld-chr $REF_DIR/weights/weights.hm3_noMHC.";
    
SUMSTAT="${sumstat_dir}/$TRAIT.sumstats.gz"

echo $TRAIT - $ANNOT - $CT - baseline
python ${ldsc} \
  --h2 $SUMSTAT \
  --ref-ld-chr $BASELINE,../gtex/unioncells_onek1k/$ANNOT.,$CT/$ANNOT. \
  $FREQ $WEIGHTS $OPTIONS \
  --out $outdir/$TRAIT.baseline.$ANNOT.$CT

echo $TRAIT - $ANNOT - $CT - baselineLD
python ${ldsc} \
 --h2 $SUMSTAT \
 --ref-ld-chr $BASELINELD,../gtex/unioncells_onek1k/$ANNOT.,$CT/$ANNOT. \
 $FREQ $WEIGHTS $OPTIONS \
 --out ${outdir}/$TRAIT.baselineLD.$ANNOT.$CT

done

