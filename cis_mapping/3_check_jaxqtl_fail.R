# check if failed cis-egenes in jaxqtl found by tqtl

library(tidyverse)
library(data.table)

jqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16/"
tqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl/"

gene_sum_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16/metadata/"
threshold <- 0.01 # threshold for jaxqtl (and tqtl)
fdr_val <- 0.05
gtex_1per <- "1per"

allcelltypes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv", F) %>%
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_")) %>% 
  pull(celltype)

for (cell_type in allcelltypes){
  print(cell_type)
  jaxqtl_genes_fail <- fread(paste0(jqtl_dir, cell_type, "/jaxqtl_allres_cis_score.", cell_type, ".genes.fail.tsv.gz"), header = TRUE)
  tqtl_genes <- fread(paste0(tqtl_dir, cell_type, "/tqtl_allres.", cell_type, ".", gtex_1per, ".cisgenes.tsv"), header = F) %>% 
    rename(phenotype_id = V1)
  
  found <- tqtl_genes %>% 
    filter(phenotype_id %in% jaxqtl_genes_fail$Geneid)
  
  print(paste0("Number of genes found: ", nrow(found)))
}
