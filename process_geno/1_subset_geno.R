## pull out SNP data from plink file
library(tidyverse)
library(data.table)
library(genio)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/test_GLMM")

geno_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08/hwe_maf0.05_snp/"

snplist <- fread("/project/nmancuso_8/elezhang/projects/jaxqtl/test_GLMM/NK.both.genes.leadsnp.tsv.gz",
                 header=TRUE)

G <- data.frame()
for (chr_idx in 1:22){
  bim <- read_bim(paste0(geno_dir, "chr", chr_idx, ".bim"))
  fam <- read_fam(paste0(geno_dir, "chr", chr_idx, ".fam"))
  g <- read_bed(paste0(geno_dir, "chr", chr_idx), bim$id, fam$id)
  varid <- rownames(g)
  g <- as.data.frame(g) %>% mutate(snp = varid) %>% 
    filter(snp %in% snplist$variant_id) %>% 
    distinct(snp, .keep_all = TRUE) %>% 
    select(snp, everything())
  G <- bind_rows(G, g)
}

out <- snplist %>% left_join(G, by = c("variant_id" = "snp"))
