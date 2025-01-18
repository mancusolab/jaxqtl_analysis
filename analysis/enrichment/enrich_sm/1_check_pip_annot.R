# check pip over annotations
# calculate proportion of variants that fall into the annotation
library(tidyverse)
library(data.table)

library(optparse)

option_list <- list(
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="jaxqtl or tqtl [default= %default]", metavar="character"),
  make_option(c("-c", "--chr"), type="character", default=NULL, 
              help="which chromosome [default= %default]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="out directory [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

method <- opt$method
chr_idx <- opt$chr
out_dir <- opt$outdir

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/enrich_sm")

geno_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/"
ldsc_dir <- paste0(geno_dir, "ldsc/")

pip_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/result_wald_label/", method, "/")

# read parameter file for fine-mapping
if (method == "jaxqtl"){
  prefix <- "nb"
  model <- "nb"
}else if (method == "tqtl"){
  prefix <- "tqtl"
  model <- "lm"
}else{
  print("none")
}

params <- read_tsv(paste0("../../finemap/code/", prefix, "_allegenes.tsv"), F)

colnames(params)[1:4] <- c("phenotype_id","chr","celltype","N")

pip_files <- list.files(pip_dir)
all_pip_summary <- data.frame()

params_onechr <- params %>% filter(chr == chr_idx)
chr_bim <- fread(paste0(geno_dir, "chr", chr_idx, ".bim"), header=F) %>% select(snp=V2)
chr_annot <- fread(paste0(ldsc_dir, "chr", chr_idx, ".annot.gz"), header=T) %>% 
  select(ends_with(".merge"))
chr_annot <- bind_cols(chr_bim, chr_annot)

for (idx in 1:nrow(params_onechr)){
  gene <- params_onechr$phenotype_id[[idx]]
  celltype <- params_onechr$celltype[[idx]]
  
  pattern <- paste0(gene, ".", celltype, ".", model, ".L10.*.tsv.gz")
  file <- pip_files[grepl(glob2rx(pattern), pip_files)]
  
  if (length(file) > 0){
    pip <- fread(paste0(pip_dir, file), header=T) %>% 
      left_join(chr_annot, by=c("snp")) %>% 
      gather(key = category, value = annot, ends_with(".merge")) %>% 
      separate(category, into=c("celltype", "category", "rm"), sep="\\.") %>% 
      select(-rm)
    
    pip_summary <- pip %>%
      group_by(category, celltype) %>% 
      summarize(n = sum(annot > 0),
                pip_mean = sum(pip*annot)/n,
                pip_or_mean = sum(pip/(1-pip) * annot)/n) %>% ungroup() %>% 
      rename(annot_celltype = celltype) %>% 
      mutate(phenotype_id = gene,
             celltype = celltype) 
    
    all_pip_summary <- bind_rows(all_pip_summary, pip_summary)
    print(chr_idx)
  }
}

all_pip_summary %>% fwrite(paste0(out_dir, "/", method, "_pip_annot.", chr_idx, ".tsv.gz"), sep="\t")

