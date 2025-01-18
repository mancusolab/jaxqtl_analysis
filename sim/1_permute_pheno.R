# permute single cell counts phenotype data and offset

library(dplyr)
library(data.table)

library(optparse)

option_list <- list(
  make_option(c("--file"), type="character", default=NULL, 
              help="phenotype file", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene name", metavar="character"),
  make_option(c("-s", "--seed"), type="integer", default=NULL, 
              help="Seed", metavar="number"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

df <- fread(opt$file, sep = "\t", header = TRUE)

gene_name <- opt$gene

set.seed(opt$seed)
perm_idx <- sample(1:nrow(df))

# permute gene expression and offset
perm_ct <- df[[gene_name]][perm_idx]
perm_offset <- df[["log_offset"]][perm_idx]

out <- df %>% select(individual:ct_pc2) %>% 
  mutate(log_offset = perm_offset,
         ct = perm_ct)

colnames(out)[colnames(out) == "ct"] <- gene_name

out %>% fwrite(opt$out, sep="\t")
