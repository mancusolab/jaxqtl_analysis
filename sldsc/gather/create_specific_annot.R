library(tidyverse)
library(glue)

library("optparse")

option_list = list(
  make_option(c("--params"), type="character", default=NULL, 
              help="specific annotation list", metavar="character"),
  make_option(c("--indir"), type="character", default=NULL, 
              help="old annotation dir", metavar="character"),
  make_option(c("--outdir"), type="character", default=NULL, 
              help="outdir name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

# read params for CT or tissues 
params <- read_tsv(opt$params, F) %>% pull(X1)
indir <- opt$indir
outdir <- opt$outdir

for (i in 1:22){
  # construct matrix
  annot <- list()
  for (CT in params){
    tmp <- read_tsv(glue("{indir}/{CT}/cs.{i}.annot.gz"))
    annot[CT] <- list(tmp$cs)
  }
  annot <- bind_rows(annot)
  
  # find snps specific to one annotation
  sp_loci <- rowSums(annot) == 1
  for (CT in params){
    annot_sp <- as.integer(annot[[CT]] * sp_loci > 0)
    if (!dir.exists(glue("{outdir}/{CT}"))){
      system(glue("mkdir -p {outdir}/{CT}"))
    }
    tibble(cs = annot_sp) %>% write_tsv(glue("{outdir}/{CT}/cs.{i}.annot"))
  }
  print(i)
}

system(glue("gzip {outdir}/*/*.annot"))
