# gather genes with non-converged glm

library(tidyverse)
library(data.table)
library(arrow)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/code")
params <- read_tsv("./nb_allegenes.tsv", F)
colnames(params) <- c("phenotype_id", "chr", "celltype", "N", "chunk")

params <- read_tsv("./nb_wald_has_ncvg.tsv", F)
colnames(params) <- c("phenotype_id", "chr", "celltype", "N", "chunk")
method <- "wald"
model <- "nb"

outdf <- data.frame()
for (idx in 1:nrow(params)){
  print(idx)
  
  chr <- params$chr[[idx]]
  celltype <- params$celltype[[idx]]
  gene <- params$phenotype_id[[idx]]
  chunk <- params$chunk[[idx]]
  
  res_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/"
  eqtl_dir <- paste0(res_dir, celltype, "/chr", chr, "/")
  
  zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")
  Z <- read_parquet(zscore_file) %>% select(-c("__index_level_0__"))
  print(paste("read parquet: ", zscore_file))
  
  Z <- Z %>% filter(phenotype_id == gene)

  outdf <- bind_rows(outdf, Z %>% mutate(celltype = celltype))
}

outdf %>% write_tsv(paste0("../result_summary/jaxqtl_nom_noncvg_moreiter.", model, ".", method, ".tsv.gz"))

# local check for nocvg wald results (and their corresponding score test)
outdf_wald <- fread("../result/finemap/result_summary/jaxqtl_nom_noncvg_moreiter.nb.wald.tsv.gz")
outdf_score <- fread("../result/finemap/result_summary/jaxqtl_nom_noncvg.nb.score.tsv.gz")

joined <- outdf_wald %>% mutate(slope=as.numeric(slope), slope_se = as.numeric(slope_se)) %>% 
  mutate(Z_wald = slope /slope_se) %>% 
  left_join(outdf_score %>% mutate(Z_score = slope / slope_se) %>% 
              select(snp, phenotype_id, celltype, Z_score),
            by=c("snp", "phenotype_id", "celltype")) 

joined %>% 
  filter(converged == FALSE) %>%
  ggplot(aes(x = Z_score, y = Z_wald)) + geom_point(aes(color = converged)) +
  facet_grid(.~converged) +
  geom_abline(slope = 1, intercept = 0, color = "red") + axis_theme

joined %>% 
  mutate(cvg_Z = ifelse(converged == FALSE, "No", "Yes"),
         maf = ifelse(af < 0.5, af, 1 - 0.5)) %>% 
  ggplot(aes(x = cvg_Z, y = maf)) + geom_boxplot() + axis_theme

outdf_wald %>% filter(alpha < 1. & converged == FALSE)
summary(outdf_wald$alpha)

outdf_wald %>% distinct(phenotype_id, celltype)

outdf_wald %>% 
  filter(phenotype_id == "ENSG00000163563" & celltype == "DC") %>% 
  # filter(converged == FALSE) %>%
  mutate(maf = ifelse(af < 0.5, af, 1-af)) %>% 
  ggplot(aes(x = converged, y = alpha)) + geom_boxplot()
