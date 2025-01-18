# format Genesymbols in bed file

library(tidyverse)

gene_lookup <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/gene_lookup.tsv.gz")

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed")

celltypes <- c("Bcell", "Myeloid", "Tnk")

# for (cell in celltypes){
#   bed <- read_tsv(paste0(cell, ".subtract.PLS.ENCODE"), F)
#   colnames(bed) <- c("chr", "peak_left", "peak_right", "GeneSymbol", "region")
#   bed %>% left_join(gene_lookup, by = c("chr", "GeneSymbol"))
# }

scent_map %>% 
  rename(celltype=yazar_celltype, annot=annot_celltype) %>% 
  select(celltype, annot)


# sc findings
bind_rows(NB_egenes_df %>% distinct(phenotype_id, celltype),
          tqtl_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
  distinct() %>% 
  left_join(gene_lookup, by = "phenotype_id") %>% 
  select(phenotype_id, GeneSymbol, chr, start, end, strand, celltype) %>% 
  left_join(scent_map, by=c("celltype"="yazar_celltype")) %>% 
  select(-celltype) %>% 
  write_tsv("./egene_lookup_new_fixalpha.tsv",col_names = F)

# bulk findings
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

expand.grid.df(bind_rows(jaxqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype),
                         tqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
                 distinct() %>% 
                 left_join(gene_lookup, by = "phenotype_id") %>% 
                 select(phenotype_id, GeneSymbol, chr, start, end, strand, celltype), 
               data.frame(celltypes)) %>% 
  select(-celltype) %>% 
  write_tsv("./egene_bulk_lookup_fixalpha.tsv", col_names = F)

bind_rows(jaxqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype),
          tqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
  distinct() %>% 
  left_join(gene_lookup, by = "phenotype_id") %>% 
  select(phenotype_id, GeneSymbol, chr, start, end, strand, celltype) %>% 
  left_join(scent_map, by=c("celltype"="yazar_celltype")) %>% 
  select(-celltype) %>% 
  write_tsv("./egene_lookup_new_fixalpha.tsv",col_names = F)

# update 0616: combine all eGenes together

expand.grid.df(bind_rows(NB_egenes_df %>% distinct(phenotype_id, celltype),
                         jqtl_lm_score_egenes_df %>% distinct(phenotype_id, celltype),
                         jaxqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype),
                         jaxqtl_linear_allcell_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
                 distinct(phenotype_id, celltype) %>% 
                 left_join(gene_lookup, by = "phenotype_id") %>% 
                 select(phenotype_id, GeneSymbol, chr, start, end, strand, celltype), 
               data.frame(celltypes)) %>% 
  select(-celltype) %>% 
  write_tsv("./egene_lookup_0616.tsv",col_names = F)
