# prepare pheno input for saige qtl

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
celltype <- args[1] # NK

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl")

########## format phenotype ##########

# read sparse data
phe <- Matrix::readMM(paste0("./input/pheno/", celltype, ".sc.mtx")) # cell x gene

# individual covariates (cell has same order as phe)
cell_meta <- fread(paste0("./input/pheno/", celltype, ".sc.cellmeta.tsv.gz")) %>% 
  mutate(sex = sex - 1) %>% 
  select(individual, age, sex, ind_cell_offset)

# center and standardize
cell_meta <- cell_meta %>% 
  distinct(individual, .keep_all=T) %>% 
  mutate(age = age/sd(age),
         sex = sex/sd(sex)) %>% 
  select(individual, age, sex) %>% 
  right_join(cell_meta %>% select(-c(age, sex)), by = "individual")

# individual genotype pcs
geno_pcs <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/features/donor_features.all.6PC.tsv") %>% 
  select(-c(sex, age)) %>% 
  mutate(across(PC1:PC6, ~ .x/sd(.x)))

cell_meta <- cell_meta %>% 
  left_join(geno_pcs, by = c("individual" = "iid"))

# e-genes found by either jaxqtl and/or lm (union set)
gene_meta <- fread(paste0("./input/pheno/", celltype, ".sc.genelist.gz"))
# gene_chr <- fread(paste0("./input/params_", celltype), header = FALSE)
# names(gene_chr) <- c("Geneid", "chr", "tss")
# gene_meta <- gene_meta %>% left_join(gene_chr, by="Geneid")

# individual-cell expression pcs
pcs <- fread(paste0("./input/pheno/", celltype, ".sc.2PC.gz"), header=FALSE)
names(pcs) <- c("ct_pc1", "ct_pc2")
pcs <- pcs %>% mutate(across(ct_pc1:ct_pc2, ~ .x/sd(.x)))

phe_df <- as.data.frame(as.matrix(phe))
colnames(phe_df) <- gene_meta$Geneid

# bind covariates
dat <- bind_cols(cell_meta, pcs, phe_df) %>% 
  mutate(log_offset = log(ind_cell_offset)) %>% 
  select(-ind_cell_offset) %>% 
  select(individual:ct_pc2, log_offset, everything())

dat %>%
  select(individual:ct_pc2, log_offset, all_of(gene_meta$Geneid)) %>% 
  fwrite(paste0("./input/pheno/", celltype, ".pheno.tsv.gz"), sep="\t")


