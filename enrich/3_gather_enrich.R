# gather enrichment results
library(tidyverse)
library(data.table)

params <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich/params", F)
colnames(params) <- c("method", "celltype", "snplist", "annot", "prefix")

indir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/"

params <- params %>% 
  mutate(out_name = paste0(method, "_", snplist, "_", annot, ".txt"))

headers <- c("obs_eqtl", "total_eqtl", "expected_eqtl", "sd", 
             "pval", "OR_L", "OR", "OR_U")

allres <- data.frame()
for (i in params$out_name){
  if (file.exists(paste0(indir, i))){
    oneres <- fread(paste0(indir, i), header = FALSE) 
  }else{
    oneres <- as.data.frame(matrix(rep(NA, length(headers)), 1, length(headers)))
  }
  colnames(oneres) <- headers
  allres <- bind_rows(allres, oneres)
}

allres <- bind_cols(params, allres)

allres %>% fwrite(paste0(indir, "allenrich_res.tsv.gz"), sep="\t")


