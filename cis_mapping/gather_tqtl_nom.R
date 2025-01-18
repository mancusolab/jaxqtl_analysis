# gather wald p values from tensorqtl
library(tidyverse)
library(data.table)
library(arrow)

ct_files <- list.dirs("../")
ct_files <- ct_files[!grepl("all_celltype", ct_files)]
ct_files <- ct_files[-1]

query <- tibble(chr = c(5, 20), 
                phenotype_id = c("ENSG00000134352", "ENSG00000101017"))

chr <- 20
genelist <- c("ENSG00000101017")

df <- data.frame()
for (CT in ct_files){
  ct <- gsub("..//","", CT)
  onedf <- read_parquet(paste0(CT, "/chr", chr, ".cis_qtl_pairs.", chr, ".parquet"))
  onedf <- onedf %>% filter(phenotype_id %in% genelist) %>% mutate(celltype = ct)
  
  df <- bind_rows(df, onedf)
}

df %>% fwrite(paste0(genelist, ".cis_qtl_pairs.tsv.gz"), sep="\t")
