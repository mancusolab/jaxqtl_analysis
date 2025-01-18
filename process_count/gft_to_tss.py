# create tss file for collapse model
import qtl.io
import pandas as pd

# use collapose model will affect 409 genes (among those pass 0.01 threshold), 
# compared to v87 version "gene" feature start site
# summary of difference:
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0     8.0    25.0   996.6   146.0 83617.0

annotation_gtf="Homo_sapiens.GRCh37.82.genes.gtf"
gtf=qtl.io.gtf_to_tss_bed(annotation_gtf, feature='transcript')
gtf.to_csv("Homo_sapiens.GRCh37.82.bed.gz",index=False, sep="\t")
