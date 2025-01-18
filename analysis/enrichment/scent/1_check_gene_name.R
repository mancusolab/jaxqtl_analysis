# check SCENT gene name
library(tidyverse)
annot_dir <- "../data/OneK1K/annotation/SCENT/"
annot_cells <- c("Bcell", "Tnk", "Myeloid")

annot_df <- data.frame()
for (cell in annot_cells){
  df <- read_tsv(paste0(annot_dir, cell, ".gene_peaks.merged.bed.gz"), F)
  annot_df <- bind_rows(annot_df, df)
}

annot_df

# 7852 unique Genes in SCENT result
# 7133 found in OneK1K
# this number doesn't change if replace special character and enforce lower case on gene name
annot_df %>% distinct(chr=X1, GeneSymbol=X4) %>% # nrow()
  mutate(GeneSymbol=tolower(GeneSymbol),
         GeneSymbol=gsub("[^[:alnum:]]", "_", GeneSymbol)) %>% 
  inner_join(gene_lookup %>% mutate(chr = paste0("chr", chr),
                                   GeneSymbol=tolower(GeneSymbol),
                                   GeneSymbol=gsub("[^[:alnum:]]", "_", GeneSymbol)) %>% 
               select(phenotype_id, chr, GeneSymbol)) %>% nrow()

gene_lookup %>% filter(grepl("CXorf38", GeneSymbol, ignore.case=TRUE))
