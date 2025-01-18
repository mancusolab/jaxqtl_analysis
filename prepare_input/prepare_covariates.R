## prepare covariates files

library(tidyverse)
setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/features")

# recode sex: female as 1
covars <- read_tsv("./new_features.tsv") %>% 
  mutate(sex=sex-1)

PCs <- read_tsv("../geno_n982_info_09/PCA/allchr_pruned.eigenvec") %>% 
  select(IID, PC1:PC6)

covars <- covars %>% left_join(PCs, by = c("individual" = "IID")) %>% 
  rename(iid = individual)

rownames_to_column(df, "VALUE")

# make covariate files for jaxqtl
covars %>% 
  write_tsv("./donor_features.all.6PC.tsv")

# make bed file for tqtl
t(covars) %>% 
  as.data.frame() %>% 
  rownames_to_column("X") %>% 
  write_tsv("./donor_features.all.6PC.bed", col_names = F)

