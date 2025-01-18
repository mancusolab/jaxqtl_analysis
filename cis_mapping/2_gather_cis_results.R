### gather results across chromosomes

library(tidyverse)
library(data.table)
library(qvalue)

# set working directory
setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha_perm")

model <- "lm" # pois, nb, lm
test <- "score" # wald
gene_sum_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new_fixalpha/metadata/"

allcelltypes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv", F) %>%
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_")) %>% 
  pull(celltype)

allcelltypes <- c(allcelltypes, "allcells")

for (cell_type in allcelltypes){
  print(cell_type)
  allres <- data.frame()
  for (chr_idx in 1:22){
    filelist <- list.files(paste0(cell_type, "/chr", chr_idx))
    res_files <- filelist[grepl(paste0("newbeta.", model,".cis_", test), filelist, fixed = TRUE)]
    print(length(res_files))
    
    for (filename in res_files){
      one_file <- fread(paste0(cell_type, "/chr", chr_idx, "/", filename), header = TRUE)
      allres <- bind_rows(allres, one_file)
    }
    
    print(chr_idx)
  }
  
  allres %>% fwrite(paste0("all_celltype/jaxqtl_allres_cis_", test, ".newbeta.", model, ".", cell_type,".tsv.gz"), sep="\t")
  
  print(mean(is.na(allres$pval_beta)))
  print(sum(allres$alpha_cov < 1.01e-8, na.rm = T))

}

## gather cis results tensorqtl
setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl_new")

for (cell_type in allcelltypes){
  print(cell_type)
  allres <- data.frame()
  for (chr_idx in 1:22){
    one_file <- fread(paste0(cell_type, "/chr", chr_idx, ".cis_qtl.txt.gz"), header = TRUE)
    allres <- bind_rows(allres, one_file)
  }
  
  allres %>% fwrite(paste0("all_celltype/tqtl_allres.cis.", cell_type,".tsv.gz"), sep="\t")
}

## gather add_covar results

model <- "nb" # pois, nb, lm

for (cell_type in allcelltypes){
  print(cell_type)
  allres <- data.frame()
  for (chr_idx in 1:22){
    filelist <- list.files(paste0(cell_type, "/chr", chr_idx))
    res_files <- filelist[grepl(paste0(model,".prs.robust.cis_wald.tsv.gz"), filelist, fixed = TRUE)]
    print(length(res_files))
    
    for (filename in res_files){
      one_file <- fread(paste0(cell_type, "/chr", chr_idx, "/", filename), header = TRUE)
      allres <- bind_rows(allres, one_file)
    }
    
    print(chr_idx)
  }
  
  allres %>% 
    filter(model_converged == TRUE) %>% 
    fwrite(paste0("all_celltype/jaxqtl_allres_", model, ".robust.prs.", cell_type,".tsv.gz"), sep="\t")
  
  print(mean(is.na(allres$pval_beta)))
  
}
