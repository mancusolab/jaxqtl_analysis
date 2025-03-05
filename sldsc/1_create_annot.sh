#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=10Gb
#SBATCH --array=1
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

finemap_dir=../../result_wald_label/jaxqtl
ldsc="/project/nmancuso_8/elezhang/software/ldsc/ldsc.py"

#create annotation files 
for CT in allcells B_IN B_Mem CD4_NC CD8_ET CD8_NC Mono_C Mono_NC NK CD4_ET CD4_SOX4 CD8_S100B DC NK_R Plasma; do
    echo $CT
    mkdir -p $CT
    for FILE in ${finemap_dir}/*.$CT.nb.L10.estvarFALSE.*.tsv.gz; do
        zcat $FILE | awk '{if (NR>1 && $4>0.0001) {print $0}}' >> $CT/$CT.txt
    done
    sbatch --mem=20000 -t 0-02:00 --wrap="perl create_annot.pl $CT"
done

