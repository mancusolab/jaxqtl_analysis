#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=5Gb
#SBATCH --array=1
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

ldsc="/project/nmancuso_8/elezhang/software/ldsc/ldsc.py"

#for CT in unioncells allcells B_IN B_Mem CD4_NC CD8_ET CD8_NC Mono_C Mono_NC NK CD4_ET CD4_SOX4 CD8_S100B DC NK_R Plasma; do
for CT in Tcell_group Bcell_group Mono_group NK_group;do
    #for ANNOT in pip cs pip_cs; do
    for ANNOT in cs; do
        for CHR in {1..22}; do
            if [ ! -f "$CT/$ANNOT.$CHR.l2.ldscore.gz" ]; then
                echo "$CT $ANNOT $CHR"
                sbatch --mem=20000 -t 0-02:00 --wrap="python ${ldsc} \
                    --l2 \
                    --bfile /project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.$CHR \
                    --ld-wind-cm 1 \
                    --annot $CT/$ANNOT.$CHR.annot.gz \
                    --thin-annot \
                    --out $CT/$ANNOT.$CHR \
                    --print-snps /project/gazal_569/DATA/ldsc/reference_files/w_hm3.txt"
            fi
        done
    done
done

# check number of files: 352
