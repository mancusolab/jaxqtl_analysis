## get p values 
library(optparse)
library(tidyverse)
library(data.table)

option_list <- list(
  make_option(c("--idx"), type="integer", default=NULL, 
              help="index of causal SNP", metavar="number"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene name", metavar="character"),
  make_option(c("--prefix"), type="character", default=NULL, 
              help="Cis file prefix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file prefix name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# load glmm
load(paste0(opt$prefix, ".rda"))

# ratio file
ratio <- fread(paste0(opt$prefix, ".varianceRatio.txt"), header = FALSE)

# cis file
data.table::fread(paste0(opt$prefix, "_cis"), header=TRUE) %>% 
  janitor::clean_names() %>% 
  filter(marker_id == paste0("rs", opt$idx)) %>% 
  mutate(phenotype_id = opt$gene,
         dispersion =  modglmm$theta[1],
         re_var = modglmm$theta[2],
         ratio_sparse = ratio$V1[1],
         ratio_null = ratio$V1[2],
         ratio_null_noXadj = ratio$V1[3],
         ratio_null_sample = ratio$V1[4],
         null_converged = modglmm$converged) %>% 
  write_tsv(paste0(opt$out, "_causal_snp"))

