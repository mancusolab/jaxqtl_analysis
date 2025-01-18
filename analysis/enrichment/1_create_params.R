# create parameter file for enrichment analysis

library(tidyverse)
library(data.table)

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich")

indir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/"

bedfiles <- read_tsv(paste0(indir, "bedfile_list"), col_names = F) %>% pull(X1)
bedfiles <- gsub(".bed*$","",bedfiles)

allcelltypes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv", F) %>%
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_")) %>% 
  pull(celltype)

# tqtl: NK.1per.leadSNP,

threshold <- 0.01
model <- "nb"

snplist_tqtl <- c(paste0(allcelltypes, ".threshold", threshold, ".leadSNP"), 
                  paste0(allcelltypes, ".", model, ".threshold", threshold, ".leadSNP.tqtl.only"))

tqtl_params <- expand.grid("celltype16_tensorqtl", snplist_tqtl, bedfiles)
colnames(tqtl_params) <- c("method", "snplist", "annot")

tqtl_params <- tqtl_params %>% 
  mutate(tss = gsub(".leadSNP|.leadSNP.tqtl.only", ".genespass.tss", snplist),
         tss = gsub(paste0(".", model), "", tss)) %>% 
  select(method, snplist, annot, tss)

# jqtl: NK.leadSNP

snplist_jqtl <- c(paste0(allcelltypes, ".", model, ".threshold", threshold, ".leadSNP"), 
                  paste0(allcelltypes, ".", model, ".threshold", threshold, ".leadSNP.jqtl.only"))

jqtl_params <- expand.grid("celltype16", snplist_jqtl, bedfiles)
colnames(jqtl_params) <- c("method", "snplist", "annot")

jqtl_params <- jqtl_params %>% 
  mutate(tss = gsub(".leadSNP|.leadSNP.jqtl.only", ".genespass.tss", snplist)) %>% 
  select(method, snplist, annot, tss)

bind_rows(tqtl_params, jqtl_params) %>% 
  write_tsv("./params", col_names = F)

