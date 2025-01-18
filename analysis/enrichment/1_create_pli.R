### create pLI score 

library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation")

ref <- fread("./gnomad.v2.1.1.lof_metrics.by_gene.txt") %>% 
  select(gene, transcript, pLI)

query <- fread("../pheno_meta/gene_lookup.tsv.gz"); nrow(query)

query <- query %>% left_join(ref, by = c("GeneSymbol" = "gene")) %>% 
  filter(!is.na(pLI)) %>% 
  group_by(Geneid) %>% 
  filter(pLI == max(pLI)) %>% 
  ungroup() %>% 
  distinct(Geneid, .keep_all = T)

query %>% fwrite("./genes_pLI.tsv.gz", sep="\t")
