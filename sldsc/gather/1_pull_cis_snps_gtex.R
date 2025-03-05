# find cis-SNPs from GTEx sumstats (hg38)
library(tidyverse)
library(arrow)
library(glue)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/sldsc_wald_label/gtex_selection/allcis")

indir <- "/project/nmancuso_8/data/GTEx/GTEXv8/from_web/eQTL_results/allsumstats"

params <- read_tsv("../tissue.list",F) %>% 
  pull(X1)

df <- data.frame()
for (CT in params){
  for (i in 1:22){
    file <- glue("{indir}/{CT}/{CT}.v8.EUR.allpairs.chr{i}.parquet")
    if (file.exists(file)){
      tmp <- read_parquet(file, col_select = "variant_id") %>% distinct()
      
      tmp <- tmp %>% 
        separate(variant_id, into=c("chrom", "pos", "ref", "alt", "rm"), sep="_", remove = F) %>% 
        distinct(variant_id, chrom, pos) %>% 
        filter(chrom %in% paste0("chr", 1:22)) %>% 
        mutate(chrom = as.integer(gsub("chr", "", chrom)), 
               pos = as.integer(pos),
               pip = 0, cs = 1)
      df <- bind_rows(df, tmp)
    }else{
      print(glue("{file} not exist!"))
    } 
  }
  df <- df %>% distinct()
  print(i)
}

df %>% 
  distinct() %>% 
  data.table::fwrite("allcis_hg38.txt", sep="\t", col.names = F)

# bed file for liftover
df %>% 
  distinct() %>% 
  mutate(chrom = paste0("chr", chrom), pos_1 = pos-1) %>% 
  select(chrom, pos_1, pos) %>% 
  data.table::fwrite("allcis_hg38.bed", sep="\t", col.names = F)
