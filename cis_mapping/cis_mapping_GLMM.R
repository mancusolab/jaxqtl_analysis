# cis mapping using poisson mixed effect model
library(tidyverse)
library(data.table)
# library(lme4)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/GLMM")

# read sparse data
phe <- Matrix::readMM('./data/pheno/NK.shared.sc.mtx') # cell x gene

# individual covariates (cell has same order as phe)
cell_meta <- fread("./data/pheno/NK.shared.cellmeta.tsv.gz") %>% 
  mutate(sex = sex - 1) %>% 
  group_by(individual) %>% 
  mutate(ind_offset = sum(ind_cell_offset)) %>% 
  ungroup()

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
gene_meta <- fread("./data/NK.shared.leadsnp.tsv.gz")

# read cis SNPs from "ENSG00000147813" (chr8)
gene <- "ENSG00000147813" # 8:144676862
geno <- fread(paste0("./data/geno/", gene, "_cis.G.tsv.gz"))  # NxM
geno_meta <- fread(paste0("./data/geno/", gene, "_cis.G.meta.tsv.gz")) %>% 
  mutate(SNP = paste0("V", 1:ncol(geno)))

geno_iid_order <- read_tsv("./data/geno/iid_order.tsv.gz") %>% pull(iid)
geno <- geno %>% mutate(individual = geno_iid_order) %>% select(individual, everything())

geno <- geno %>% gather(key="SNP",value="genotype",-c(individual))

# individual-cell expression pcs
ct_pcs <- fread("./data/pheno/NK.sc.2PC.gz", header=FALSE)
names(ct_pcs) <- c("ct_pc1", "ct_pc2")
ct_pcs <- ct_pcs %>% mutate(across(ct_pc1:ct_pc2, ~ .x/sd(.x)))

# set.seed(2023)
# # permute rows once and create pseudo bulk data
# sample_id <- sample(1:nrow(cell_meta), replace = FALSE)
# 
# phe <- phe[sample_id, ]
phe_df <- as.data.frame(as.matrix(phe)) # invert to individual x cell
colnames(phe_df) <- gene_meta$phenotype_id

# ct_pcs <- ct_pcs[sample_id, ]

# bind covariates
dat <- bind_cols(cell_meta, ct_pcs, phe_df) %>% 
  mutate(log_offset = log(ind_cell_offset))

# glmer
res <- geno_meta %>% mutate(glmm_pval = NA, glmm_slope = NA, glmm_se = NA)

for (snp in geno_meta$SNP){
  print(snp)
  dat <- dat %>% left_join(geno %>% 
                             filter(SNP == snp) %>% 
                             select("individual", "genotype"), 
                           by = "individual")
  
  # fit <- glmer(ENSG00000167077 ~ offset(log_offset) + (1|individual) + genotype + age + sex + ct_pc1 + ct_pc2,
  #              data = dat, family = poisson,
  #              nAGQ = 25,
  #              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  
  
  try({
    fit <- MASS::glmmPQL(as.formula(paste0(gene, " ~ offset(log_offset) + genotype + age + sex + ct_pc1 + ct_pc2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6")),
                         random = ~ 1 | individual,
                         data = dat, family = poisson, niter=100)
    
    res$glmm_pval[res$SNP == snp] <- summary(fit)$tTable[,"p-value"]["genotype"]
    res$glmm_slope[res$SNP == snp] <- summary(fit)$tTable[,"Value"]["genotype"]
    res$glmm_se[res$SNP == snp] <- summary(fit)$tTable[,"Std.Error"]["genotype"]
    
  })
  
  dat <- dat %>% select(-genotype)
}

# problem: ENSG00000128604
res %>% write_tsv(paste0("./result/glmm_pois_NK.", gene, ".cis.tsv.gz"))
