# get permutation p values files
library(optparse)
library(tidyverse)
library(glue)

option_list <- list(
  make_option(c("--i"), type="integer", default=NULL, 
              help="row index of metadata")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

datasets <- read_tsv("CT_metadata")
i <- as.integer(opt$i)

out_name <- datasets$out_name[i]
if (!file.exists(glue("{out_name}.tsv.gz"))){
  qts <- datasets$study_id[i]
  qtd <- datasets$dataset_id[i]
  out_name <- datasets$out_name[i]
  print(i)
  system(glue("curl --output {out_name}.tsv.gz ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{qts}/{qtd}/{qtd}.permuted.tsv.gz"))
}