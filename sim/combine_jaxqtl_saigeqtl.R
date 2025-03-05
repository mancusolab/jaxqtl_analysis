# combine saigeqtl results
library(tidyverse)
library(glue)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--nsim"), type="integer", default=NULL, 
              help="number of repeats in each simulation", metavar="number"),
  make_option(c("--params"), type="character", default=NULL, 
              help="parameter file path", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

params_name <- opt$params # params_Feb2_25
nsim <- opt$nsim
fwer <- 0.05

params <- read_tsv(glue("./params/{params_name}"), col_names=F)

allres <- data.frame()
for (sim_idx in 1:nrow(params)){
  res_file <- glue("./output/{params_name}/sim{sim_idx}.tsv")
  tmp <- fread(res_file, header=T)
  allres <- bind_rows(allres, tmp)
}

allres <- allres %>% mutate(sim_idx = 1:nrow(.))

allres$rep <- NA
allres$rej_sqtl_score <- NA

for (sim_idx in allres$sim_idx){
  pval <- c()
  for (rep_idx in 1:nsim){
    res_file <- glue("./output/{params_name}/sim{sim_idx}.pheno{rep_idx}_sqtl_causal_snp")
    if (file.exists(res_file)){
      tmp <- fread(res_file, header=T)
      if (tmp$null_converged == TRUE){
        pval <- c(pval, tmp$p_value) 
      }
    }
  }
  pval <- pval[(!is.na(pval) & is.finite(pval))]
  allres$rep[sim_idx] <- length(pval)
  allres$rej_sqtl_score[sim_idx] <- mean(pval < fwer)
  print(sim_idx)
}

allres %>% write_tsv(glue("./output/{params_name}/{params_name}_allres.tsv"))

