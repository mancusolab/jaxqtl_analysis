# find cis-SNPs from Onek1k
library(tidyverse)
library(arrow)
library(glue)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/sldsc_wald_label/jaxqtl/allcis")

indir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha"

params <- read_tsv("allgenes_pass_1per.tsv")

df <- data.frame()
for (i in 1:nrow(params)){
  chunk <- params$chunk[i]
  chr <- params$chr[i]
  CT <- params$celltype[i]
  
  file <- glue("{indir}/{CT}/chr{chr}/{chunk}.nb.cis_qtl_pairs.{chr}.wald.parquet")
  if (file.exists(file)){
    tmp <- read_parquet(file, col_select=c("snp", "chrom", "pos"))
    tmp <- tmp %>% 
      distinct() %>% 
      mutate(chrom = as.character(chrom), pos = as.integer(pos),
             pip = 0, cs = 1)
    df <- bind_rows(df, tmp)
  }else{
    print(glue("{file} not exist!"))
  }
  print(i)
}

df %>% 
  distinct() %>% 
  data.table::fwrite("allcis.txt", sep="\t", col.names = F)
