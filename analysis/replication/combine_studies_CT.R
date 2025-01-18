# combine 
library("tidyverse")
library("glue")
library("qvalue")

setwd("/project/nmancuso_8/data/eQTL-Catalogue/allsumstats")

# run local
datasets <- read_tsv("CT_metadata")
for (CT in unique(datasets$group)){
  allres <- data.frame()
  tmp <- datasets %>% filter(group == CT)
  
  for (file_name in tmp$out_name){
    df <- read_tsv(glue("{file_name}.allpairs.tsv.gz"), col_types="nncciiccccciinnnnn") %>%
      filter(chromosome %in% 1:22 & type == "SNP") %>% 
      mutate(qval = qvalue(p_beta, fdr.level = 0.05)$qvalue) %>%
      filter(qval < 0.05)
    allres <- bind_rows(allres, df)
  }
  allres %>% distinct() %>% 
    write_tsv(glue(CT, ".tsv.gz"))
}
