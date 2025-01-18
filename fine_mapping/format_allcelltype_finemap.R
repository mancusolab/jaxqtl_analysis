# format all cell type finemap
# add bins to tss

library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/result_wald_label")
cols_names <- c("chr", "pos_1", "pos", "snp", "pip", "phenotype_id", "celltype")

both_egenes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/enrich_sm/both_egenes")

gtf <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/Homo_sapiens.GRCh37.82.bed.gz")

# TODO: add file path for cis-region of egene

for (method in c("jaxqtl", "tqtl")){
  cs <- read_tsv(paste0(method, "/bed/allcelltype_finemap.cs.tsv.gz"), F) %>% 
    mutate(X7=gsub(".cs.bed", "", X7))
  colnames(cs) <- cols_names
  
  cs <- cs %>% left_join(gtf %>% select(phenotype_id=gene_id,tss=end),by=c("phenotype_id"))
  
  cs <- cs %>% inner_join(both_egenes, by=c("phenotype_id", "celltype"))
  
  window <- 500000
  
  tmp <- cs %>% distinct(celltype, phenotype_id, tss)
  
  cs_new <- data.frame()
  cissnp_df <- data.frame()
  for (idx in 1:nrow(tmp)){
    gene <- tmp$phenotype_id[idx]
    cell <- tmp$celltype[idx]
    tss <- tmp$tss[idx]
    break_points <- seq(tss - window, tss + window, 500)
    # cut into (break_point1, break_point2], 
    labels <- 1:(length(break_points)-1)
    
    cs_new <- bind_rows(cs_new, 
                        cs %>% filter(celltype == cell & phenotype_id == gene) %>% 
                          mutate(bins = cut(pos, breaks=break_points, labels=labels)))
    print(idx)
    
    # TODO: count number of cis-snps in this window
  }
  
  if (nrow(cs) == nrow(cs_new)){
    cs_new %>% write_tsv(paste0(method, "/bed/allcelltype_finemap_bothegenes_addtssbins.cs.tsv.gz"))
  }

}

