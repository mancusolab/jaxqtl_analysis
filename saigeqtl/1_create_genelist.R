## find genes that pass expression threshold > 10% in all three representative CTs
## randomly select 50 genes from chr1 for speed benchmark
## write genelist file (no headers) with three cols: phenotype_id, chr, tss

library(tidyverse)
library(glue)
outdir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl/input/"
metadir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/"
gtf <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/gene_lookup.tsv.gz")

CTs <- c("CD4_NC", "B_IN", "Plasma")
thresh <- 0.1 # 10%
which_chr <- 1
n_genes <- 50

alldf <- c()
for (CT in CTs){
  meta <- read_tsv(glue(metadir, CT, "_gene_summary.tsv.gz")) %>%
    filter(express_percent > thresh & chr == which_chr) %>% 
    mutate(celltype = CT)
  alldf <- bind_rows(alldf, meta)
}

set.seed(2024)
alldf %>% count(Geneid) %>% filter(n == length(CTs)) %>% 
  sample_n(n_genes, replace = F) %>% 
  inner_join(gtf %>% select(Geneid, chr, end), by="Geneid") %>% 
  select(-n) %>%
  arrange(end) %>%
  write_tsv(glue(outdir, "params_50"), col_names = F)
  
