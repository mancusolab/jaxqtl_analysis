## create up sampling 
library(tidyverse)
library(data.table)
library(glue)
library(genio)

library(optparse)

option_list <- list(
  make_option(c("-c", "--celltype"), type="character", default=NULL, 
              help="which cell type", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

CT <- opt$celltype
outdir_saigeqtl <- "/project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl/input/pheno"
outdir_jqtl <- "/project/nmancuso_8/elezhang/projects/jaxqtl/compare_speed/input"
indir_geno <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/ldprune_mac20"

# downsample individuals
CTs <- c("CD4_NC", "B_IN", "Plasma")
N_list <- c(10000)
max_cells <- 5000
CT_cts <- floor(c(0.37, 0.06, 0.003) * max_cells)
n_pcs <- 2

set.seed(2024)
num_cells <- CT_cts[CTs == CT]

for (N in N_list){
  # first read single cell data matrix
  file_name <- glue("/project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl/input/pheno/{CT}.pheno.tsv.gz")
  df <- read_tsv(file_name)
  
  # first filter individuals ids
  N_ids <- sample(c(unique(df$individual)), size=N, replace=TRUE)
  
  new_dat <- data.frame()
  for (i in 1:N){
    old_id <- N_ids[i]
    tmp <- df %>% filter(individual == old_id) %>%
      sample_n(size = num_cells, replace = TRUE) %>% 
      mutate(individual = glue("{i}_{i}"))
    new_dat <- bind_rows(new_dat, tmp)
    print(i)
  }
  
  new_dat %>% fwrite(glue("{outdir_saigeqtl}/{CT}.N{N}.pheno.tsv.gz"), sep="\t")
  new_ids <- unique(new_dat$individual)
  
  # create pseudo-bulk data for jaxqtl
  file_name <- glue("{outdir_jqtl}/{CT}.N100.bed.gz") # to get geneid info
  old_bulk <- read_tsv(file_name) %>% select(`#Chr`:Geneid)
  
  jqtl_dat <- new_dat %>% 
    rename(iid=individual) %>% 
    gather(key = Geneid, value = ct, starts_with("ENSG")) %>% 
    group_by(iid, Geneid) %>% 
    summarize(ct = sum(ct)) %>% ungroup() %>% 
    spread(key = iid, value = ct) %>% 
    left_join(old_bulk, by="Geneid") %>% 
    select(`#Chr`, start, end, Geneid, everything())
  
  jqtl_dat %>% fwrite(glue("{outdir_jqtl}/{CT}.N{N}.bed.gz"), sep="\t")
  print("write jaxqtl bulk data")
  
  # create covariate file for jaxqtl pseudo-bulk data
  new_cov <- new_dat %>% rename(iid=individual) %>% 
    distinct(iid, age, sex, PC1, PC2, PC3, PC4, PC5, PC6) 
  
  new_cov %>% 
    write_tsv(glue("{outdir_jqtl}/donor_features.N{N}.{CT}.6PC.tsv"))
  
  # create log(1+y) data for linear model (for simplicity)
  tqtl_mat <- log((jqtl_dat[,-c(1:4)])+1)
  
  jqtl_dat %>% select(`#Chr`, start, end, Geneid) %>% 
    bind_cols(tqtl_mat) %>% 
    fwrite(glue("{outdir_jqtl}/{CT}.N{N}.tmm.bed.gz"), sep="\t")
  
  print("write linear bulk data")
  
  # calculate PC and create covariate for tensorqtl linear model
  tqtl_mat_scale <- scale(t(tqtl_mat))
  svd_fit <- svd(tqtl_mat_scale)
  pcs <- svd_fit$u[,1:n_pcs] * svd_fit$d[1:n_pcs]
  pcs <- as.data.frame(pcs)
  colnames(pcs) <- paste0("EPC", 1:n_pcs)
  
  tqtl_cov <- bind_cols(new_cov, pcs)
  
  as.data.frame(t(tqtl_cov)) %>% 
    rownames_to_column() %>% 
    write_tsv(glue("{outdir_jqtl}/donor_features.all.6PC.{CT}.N{N}.PC.bed"),
              col_names = F)
}


# # sample for genotype data
# N <- 10000
# geno <- read_plink(glue("{indir_geno}/chr1"))
# bim <- geno$bim
# m_loci <- nrow(bim)
# 
# af <- rowSums(geno$X) / (2 * ncol(geno$X))
# X <- sapply(af, function(p) rbinom(n = N, size = 2, prob = p))
# new_ids <- paste0(1:N, "_", 1:N)
# rownames(X) <- new_ids
# X <- t(X)
# fam <- tibble(fam = 0, id = new_ids, pat = 0, mat = 0, sex = 0, pheno = -9)
# 
# # write new geno
# write_plink(glue("/project/nmancuso_8/elezhang/projects/jaxqtl/saigeqtl/input/geno/chr1_10k"), 
#             X, bim = bim, fam = fam)
