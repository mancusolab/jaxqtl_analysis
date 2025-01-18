# map between hg38 and hg19 for filtered SNPs

library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_rsid")

dbsnp_dir <- "/project/nmancuso_8/data/dbSNP/dbSNP_v155/split/"

version <- "25" # "39": hg38; "25": hg19

hg19 <- "GCF_000001405.chr1.renamed.25.gz"

for (chr in 1:22){
  query <- fread(paste0("chr", chr, ".bim"), col.names = c("CHR", "ID", "cM", "POS", "ALT", "REF"))
  hg19 <- fread(paste0(dbsnp_dir, "GCF_000001405.chr", chr, ".renamed.25.gz"), skip="#CHROM")
  
  hg19 %>% select(CHR=`#CHROM`, POS_hg19=POS, ID, REF_hg19=REF, ALT_hg19=ALT) %>% 
    inner_join(query, by = c("CHR", "POS", "ID", "REF_hg19"="REF"))
  
  hg38 <- fread(paste0(dbsnp_dir, "GCF_000001405.chr", chr, ".renamed.39.gz"), skip="#CHROM")
  
  hg38 %>% select(CHR=`#CHROM`, POS_hg38=POS, ID, REF_hg38=REF, ALT_hg38=ALT) %>% 
    inner_join(hg19 %>% select(CHR=`#CHROM`, POS_hg19=POS, ID, REF_hg19=REF, ALT_hg19=ALT), 
               by = c("CHR", "ID"))
  
}
