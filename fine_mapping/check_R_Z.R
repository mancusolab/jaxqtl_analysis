# check R and Z score convergence

library(tidyverse)
library(data.table)
library(arrow)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/code")
params <- read_tsv("./nb_allegenes.tsv", F)
colnames(params) <- c("phenotype_id", "chr", "celltype", "N", "chunk")

blacklist <- read_tsv("../result_wald/blacklist.tsv",F) %>% pull(X1)
method <- "wald"
model <- "nb"

params <- params %>% filter(!phenotype_id %in% blacklist) %>% 
  mutate(num_cis = NA,
         nocvg_Z = NA,
         na_Z = NA,
         na_R = NA)

for (idx in 1:nrow(params)){
  print(idx)
  
  chr <- params$chr[[idx]]
  celltype <- params$celltype[[idx]]
  gene <- params$phenotype_id[[idx]]
  chunk <- params$chunk[[idx]]
  
  res_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/"
  eqtl_dir <- paste0(res_dir, celltype, "/chr", chr, "/")
  
  R_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/", celltype, "/")
  
  zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")
  Z <- read_parquet(zscore_file) %>% select(-c("__index_level_0__"))
  print(paste("read parquet: ", zscore_file))
  
  Z <- Z %>% filter(phenotype_id == gene)
  
  # pull vector of sumstats # check if there is missing value and convergence
  z_scores <- Z$slope/Z$slope_se
  snp_to_rm <- which(Z$converged == FALSE)
  
  # read in correlation (r), cell-type specific because projecting out diff expression PCs
  R <- data.table::fread(paste0(R_dir, gene, ".ld_wt.tsv.gz"), header = F)
  R <- as.matrix(R)
  
  params$num_cis[[idx]] <- length(z_scores)
  params$na_R[[idx]] <- sum(is.na(R))
  params$nocvg_Z[[idx]] <- length(snp_to_rm)
  params$na_Z[[idx]] <- sum(is.na(z_scores))
}

params %>% write_tsv(paste0("../result_summary/check_Z_R.", model, ".", method, ".tsv.gz"))


# warning: now fixed
# ENSG00000162627 from chr1 in B_Mem has all missing value in correlation R (likely caused by GLM weight all NAs)
