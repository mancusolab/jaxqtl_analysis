# create permutation of pseudobulk data

library(tidyverse)
library(data.table)

wkdir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/"
setwd(wkdir)

celltype_path="./pheno/celltype16_new/metadata/celltype_14.tsv"
indir = "./pheno/celltype16_new/"
out_dir = "./pheno/celltype16_new/perm/"

all_celltype <- read_tsv(celltype_path, F) %>% 
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_"))

set.seed(2024)
for (celltype in all_celltype$celltype){
  for (suffix in c(".bed", ".tmm.bed")){
    df <- fread(paste0(indir, celltype, suffix, ".gz"), header=TRUE) 
    new_id <- sample(colnames(df)[5:ncol(df)], replace=FALSE)
    colnames(df)[5:ncol(df)] <- new_id
    df %>% fwrite(paste0(out_dir, celltype, suffix, ".gz"), sep="\t") 
  }
 print(celltype)
}

# permute the labels in tensorqtl tmm.bed same as in above


