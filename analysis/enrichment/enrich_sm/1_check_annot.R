# check pip over annotations
# calculate proportion of variants that fall into the annotation
library(tidyverse)
library(data.table)

library(optparse)


setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldsc")
out_dir <- "./summary/"

df <- fread(paste0("chr1.annot.gz"), header=T)
cols <- colnames(df)

num_var <- nrow(df)
allsum <- colSums(df)
for (chr in 2:22){
  # one annotation per column
  df <- fread(paste0("chr",chr, ".annot.gz"), header=T)
  if (all.equal(colnames(df), cols)){
    allsum <- allsum + colSums(df)
    num_var <- num_var + nrow(df)
    print(chr)
  }
}


df %>% fwrite(paste0(out_dir, "/", method, "_pip_annot.", chr_idx, ".tsv.gz"), sep="\t")

