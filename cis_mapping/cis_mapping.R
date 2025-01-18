### run cis mapping in R
library(tidyverse)
library(data.table)
library(statmod)
library(broom)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/code/test_cis")

args <- commandArgs(trailingOnly=TRUE)
y_path <- args[1]
G_path <- args[2]
X_path <- args[3]
offset_path <- args[4]
nperm <- as.integer(args[5])

y <- fread(y_path,header=F) %>% pull(V1)
G <- fread(G_path,header=F) %>% as.data.frame()
X <- fread(X_path,header=F) %>% select(-V1) # remove intercept
offset_eta <- fread(offset_path, header=F) %>% pull(V1)

dat <- X %>% mutate(y = y, offset_eta = offset_eta)

run_map_score <- function(dat){
  pval_vec <- c()
  z_vec <- c()
  fit_null <- MASS::glm.nb(y~V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+offset(offset_eta), data=dat)
  
  for (j in 1:ncol(G)){
    z <-  glm.scoretest(fit_null, G[,j])
    pval_vec <- c(pval_vec, pnorm(abs(z), lower.tail = FALSE)*2)
    z_vec <- c(z_vec, z)
  } 
  min_p <- min(pval_vec)
  return(min_p)
}

run_map_wald <- function(dat){
  pval_vec <- c()
  for (j in 1:ncol(G)){
    dat <- dat %>% mutate(g = G[,j])
    fit <-  MASS::glm.nb(y~V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+g+offset(offset_eta), data=dat)
    pval_vec <- c(pval_vec, tidy(fit)$p.value['g'])
  } 
  min_p <- min(pval_vec)
  
  return(min_p)
}

# score test
start_time <- Sys.time()
print(paste0("Score test start: ", start_time))

score_time <- system.time({
  p_nom <- run_map_score(dat)
  n_obs <- nrow(dat)
  
  set.seed(1)
  pval_score <- c()
  
  for (idx in 1:nperm){
    perm_idx <- sample(1:n_obs, n_obs, replace=FALSE)
    dat_perm <- dat %>% 
      mutate(y = y[perm_idx], offset_eta = offset_eta[perm_idx])
    pval_score <- c(pval_score, run_map_score(dat_perm))
  }
})

end_time <- Sys.time()
print(paste0("Score test end: ", end_time))

print(paste0("Score test time: ", score_time['elapsed']))

###### wald test ######
start_time <- Sys.time()
print(paste0("wald test start: ", start_time))

wald_time <- system.time({
  p_nom <- run_map_wald(dat)
  n_obs <- nrow(dat)
  
  set.seed(1)
  pval_wald <- c()
  
  for (idx in 1:nperm){
    perm_idx <- sample(1:n_obs, n_obs, replace=FALSE)
    dat_perm <- dat %>% 
      mutate(y = y[perm_idx], offset_eta = offset_eta[perm_idx])
    pval_wald <- c(pval_wald, run_map_wald(dat_perm))
  }
})

end_time <- Sys.time()
print(paste0("wald test end: ", end_time))

print(paste0("Wald test time: ", wald_time['elapsed']))

tibble(score_time = score_time['elapsed'],
       wald_time = wald_time['elapsed']) %>% 
  write_tsv(paste0(gsub("_y", "", y_path), ".runtime"))
