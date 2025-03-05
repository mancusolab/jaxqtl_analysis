# create annotation for gtex brain tissues
library(tidyverse)
library(data.table)
library(glue)

indir <- "/project/nmancuso_8/data/GTEx/GTEXv8/from_web/eQTL_results/susie_finemap"
outdir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/sldsc_wald_label/gtex"

df <- fread(glue("{indir}/GTEx_49tissues_release1.tsv.bgz"))

# cs=-1 means no in 95% CS
# susie results has cs_id from 1-10
# total 13 brain tissues
df <- df %>% 
  filter(method == "SUSIE") %>% 
  # filter(grepl("Brain",tissue)) %>% 
  filter(cs_id > 0 & chromosome %in% paste0("chr", 1:22)) %>% 
  filter(str_length(allele1) == 1 & str_length(allele2) == 1) %>% 
  mutate(chr = gsub("chr", "", chromosome)) 

tissues <- unique(df$tissue)

# write out per tissue 
for (CT in tissues){
  system(glue("mkdir -p {outdir}/{CT}"))
  df %>% filter(tissue == CT) %>%
    mutate(variant = paste0(variant, "_b37"),
           cs = 1) %>% 
    distinct(variant, chr, end, pip, cs) %>%
    write_tsv(glue("{outdir}/{CT}/{CT}.txt"), col_names = F)
  print(CT)
}

tibble(x=tissues) %>% write_tsv(glue("{outdir}/tissue.list"), col_names = F)

