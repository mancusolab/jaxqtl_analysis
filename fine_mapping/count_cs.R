# count number of credible sets output from SuSiE
# Note: this contains all results from fine-mapping folder

library(tidyverse)
library(data.table)
library(optparse)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap")

option_list <- list(
  make_option(c("--indir"), type="character", default=NULL, 
              help="fine mapping result dir", metavar="character"),
  make_option(c("--outdir"), type="character", default=NULL, 
              help="output dir", metavar="character"),
  make_option(c("--model"), type="character", default=NULL, 
              help="model of eQTL, nb or lm", metavar="character"),
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

indir <- opt$indir  # indir=paste0("./result_wald_label/")
outdir <- opt$outdir # outdir="./result_summary/"
model <- opt$model

allfiles <- list.files(indir)

# summarize CS
suffix <- paste0(".", model, ".L10.estvarFALSE.wald.cs.tsv")
cs_files <- allfiles[grepl("cs.tsv$", allfiles)]

df_files <- tibble(file = cs_files) %>%
  mutate(file_tmp = gsub("score", "wald", file),
         to_split = gsub(suffix, "", file_tmp),
         method = ifelse(grepl("wald", file), "wald", "score")) %>%
  select(-file_tmp) %>%
  separate(to_split, into=c("phenotype_id", "celltype"), sep = "\\.")

outdf <- data.frame()
for (idx in 1:nrow(df_files)){
  onefile <- df_files$file[idx]
  df <- fread(paste0(indir, onefile), header = T,
              colClasses = c("integer", rep("numeric",3), "character")) %>%
    mutate(phenotype_id = df_files$phenotype_id[idx],
           celltype = df_files$celltype[idx],
           method = df_files$method[idx],
           filename = onefile)
  outdf <- bind_rows(outdf, df)
  print(idx)
}

outdf %>% write_tsv(paste0(outdir, model, "_CS.wald.label.tsv.gz"))

# summarize PIP max (marginal)
suffix <- paste0(".", model, ".L10.estvarFALSE.wald.tsv.gz")
pip_files <- allfiles[grepl(".tsv.gz$", allfiles)]

outdf <- tibble(file = pip_files) %>%
  mutate(file_tmp = gsub("score", "wald", file),
         to_split = gsub(suffix, "", file_tmp),
         method = ifelse(grepl("wald", file), "wald", "score")) %>%
  select(-file_tmp) %>%
  separate(to_split, into=c("phenotype_id", "celltype"), sep = "\\.") %>%
  mutate(pip_max = NA)

pip_cutoff <- seq(0.05, 0.95, 0.05)
pip_cutoff_df <- matrix(nrow=nrow(outdf), ncol=length(pip_cutoff))
for (idx in 1:nrow(outdf)){
  onefile <- outdf$file[idx]
  df <- fread(paste0(indir, onefile), header = T)

  pip_cutoff_df[idx,] <- sapply(pip_cutoff, function(x) mean(df$pip > x))
  outdf$pip_max[idx] <- max(df$pip)
  print(idx)
}

pip_cutoff_df <- as.data.frame(pip_cutoff_df)
colnames(pip_cutoff_df) <- paste0("pip_gt_", pip_cutoff)

# this is everything, need to subset to eGenes identified by jaxqtl linear score
bind_cols(outdf, pip_cutoff_df) %>% write_tsv(paste0(outdir, model, "_pip.wald.label.tsv.gz"))


# # summarize lambda
# suffix <- paste0(".", model, ".L10.wald.lambda.tsv")
# lamb_files <- allfiles[grepl(".lambda.tsv$", allfiles)]
# 
# outdf <- tibble(file = lamb_files) %>%
#   mutate(file_tmp = gsub("score", "wald", file),
#          to_split = gsub(suffix, "", file_tmp),
#          method = ifelse(grepl("wald", file), "wald", "score")) %>%
#   select(-file_tmp) %>%
#   separate(to_split, into=c("phenotype_id", "celltype"), sep = "\\.") %>%
#   mutate(lambda = NA)
# 
# for (idx in 1:nrow(outdf)){
#   onefile <- outdf$file[idx]
#   df <- fread(paste0(indir, onefile), header = T)
# 
#   outdf$lambda[idx] <- df$lambda
#   print(idx)
# }
# 
# outdf %>% write_tsv(paste0(outdir, model, "_lamda.wald.label.tsv.gz"))
