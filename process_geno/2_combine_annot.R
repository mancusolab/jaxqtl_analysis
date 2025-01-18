# combine annotation in R

library(tidyverse)
library(data.table)
library(optparse)

# option_list <- list(
#   make_option(c("--chr"), type="character", default=NULL, 
#               help="which chr", metavar="character"),
#   make_option(c("--annotlist"), type="character", default=NULL, 
#               help="list of annotation, no header", metavar="character"),
#   make_option(c("--indir"), type="character", default=NULL, 
#               help="annotation dir", metavar="character"),
#   make_option(c("--outdir"), type="character", default=NULL, 
#               help="output dir", metavar="character"),
# )
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)
# 
# chr <- opt$chr

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldsc")

# read all annotation to combine
annot_list <- read_tsv("./annot_list", F) %>% pull(X1)

for (chr in 1:22){
  
  out <- fread(paste0(annot_list[1], ".", chr, ".annot.gz"), header=T)
  colnames(out) <- annot_list[1]
  
  for (annot in annot_list[-1]){
    df <- fread(paste0(annot, ".", chr, ".annot.gz"), header=T)
    colnames(df) <- annot
    out <- bind_cols(out, df)
    print(annot)
  }
  
  if (ncol(out) == length(annot_list)){
    out %>% fwrite(paste0("chr", chr, ".annot.gz"), sep = "\t") 
  }
  
  rm(out)
}
