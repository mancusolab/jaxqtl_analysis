library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/sim/NB_nocov")
null_alt <- "alt"
  
if (null_alt == "null"){
  params <- read_tsv(paste0("../../../code/sim/params_", null_alt, "_beta0"), F) %>%
    rename(alpha=X1, beta0=X2, maf=X3, threshold=X4, id=X5) %>% 
    mutate(id = paste0(id, "_", null_alt, ".tsv")) 
}else{
  params <- read_tsv(paste0("../../../code/sim/params_", null_alt, "_beta0"), F) %>%
    rename(alpha=X1, beta0=X2, beta=X3, maf=X4, threshold=X5, id=X6) %>% 
    mutate(id = paste0(id, "_", null_alt, ".tsv"))
}

headers <- read_tsv(paste0("sim1_", null_alt, ".tsv")) %>% colnames()

allres <- data.frame()
for (i in 1:nrow(params)){
  oneres <- fread(params$id[i])
  allres <- bind_rows(allres, oneres)
}

bind_cols(params, allres) %>% 
  fwrite(paste0("./sim_all_", null_alt, "_beta0_threshold.new.tsv"), sep = "\t")

# calculate adjusted power #

allres_null <- read_tsv("./sim_all_null_beta0.tsv")
allres_alt <- read_tsv("./sim_all_alt_beta0.tsv")
fwer <- 0.05

headers <- read_tsv(paste0("sim1_null.tsv")) %>% colnames()
params_null <- allres_null %>% select(alpha:id)
allres_null <- allres_null %>% select(rej_nb_wald:rej_pois_score)

correct_factor <- allres_null/fwer
correct_factor[correct_factor < 1] <- 1

correct_factor <- bind_cols(params_null, allres_null/fwer)

params_alt <- allres_alt %>% select(alpha:id)

params_alt <- params_alt %>% 
  left_join(correct_factor %>% select(-id), by = c("alpha", "maf")) %>% 
  mutate(id = gsub(".tsv", ".pval.tsv", id))

allres <- data.frame(rej_lm = rep(NA, nrow(params_alt)),
                     rej_nb_score = rep(NA, nrow(params_alt)), 
                     rej_pois_score = rep(NA, nrow(params_alt)))

for (i in 1:nrow(params_alt)){
  oneres <- fread(params_alt$id[i])
  pvals <- oneres$rej_lm[!is.na(oneres$rej_lm)]
  allres$rej_lm[i] <- mean(pvals < (fwer / params_alt$rej_lm[i]))
  
  pvals <- oneres$rej_nb_score[!is.na(oneres$rej_nb_score)]
  allres$rej_nb_score[i] <- mean(pvals < (fwer / params_alt$rej_nb_score[i]))
  
  pvals <- oneres$rej_pois_score[!is.na(oneres$rej_pois_score)]
  allres$rej_pois_score[i] <- mean(pvals < (fwer / params_alt$rej_pois_score[i]))
}


allres_alt <- allres_alt %>% 
  mutate(rej_lm = allres$rej_lm,
         rej_nb_score = allres$rej_nb_score,
         rej_pois_score = allres$rej_pois_score)

allres_alt %>% 
  write_tsv("./sim_all_alt_beta0_corrected.tsv")
