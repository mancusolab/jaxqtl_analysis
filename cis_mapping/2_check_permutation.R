# check beta distribution and p value under permutation
library(tidyverse)
library(data.table)
library(qvalue)
library(statmod)

plots_dir <- "../result/cis/figures/"
celltype <- "NK"
res_dir <- paste0("../result/cis/celltype16/", celltype, "_perm")

# ENSG00000100181: converged
# ENSG00000254709: NA, not converge, pval all zeros
# ENSG00000100097: not NA, not converge
# ENSG00000235343: converged, but U shape

# ENSG00000254709: beta estimates NA
count <- fread("../data/test/y_NA.tsv",header=F) %>% pull(V1)
G <- fread("../data/test/G_NA.tsv",header=F) %>% as.data.frame()
X <- fread("../data/test/X.tsv",header=F) %>% select(-V1)
offset_eta <- fread("../data/test/offset.tsv", header=F) %>% pull(V1)

dat <- X %>% mutate(y = count)
dim(X); dim(G)

# find variants

# ENSG00000254709: lead SNP index 1484
fit_null <- glm(y~., family=poisson, data=dat, offset=offset_eta)
lead_idx <- 0
z_vec_score <- c()
z_vec_wald <- c()
z_wald_robust <- c()
p_vec <- c()

for (j in 1:ncol(G)){
  # score test
  z <-  glm.scoretest(fit_null, G[,j])
  p_vec <- c(p_vec, pnorm(abs(z), lower.tail = FALSE)*2)
  z_vec_score <- c(z_vec_score, z)
  
  # wald not robust
  fit_full <- glm(y~., family=poisson, data=dat %>% mutate(g=G[,j]), offset=offset_eta)
  res <- fit_full %>% broom::tidy() %>% filter(term == 'g')
  z_vec_wald <- c(z_vec_wald, res$statistic)
  
  # wald robust (Huber whites)
  se2 <- diag(sandwich::vcovHC(fit_full, type = "HC0"))['g']
  z <- fit_full$coefficients['g']/sqrt(se2)
  z_wald_robust <- c(z_wald_robust, z)
} 

plot(z_vec_wald, z_vec_score) # should be asymptotically similar
plot(z_vec_wald, z_wald_robust) # robust SE reduce Z

dat %>% 
  mutate(g = G[,1351]) %>% 
  ggplot(aes(x = as.factor(g), y = y/offset_eta)) + geom_boxplot()

# permute count
set.seed(123)
num_perm <- 20

max_z_score <- c()
max_z_wald_robust <- c()
max_z_wald <- c()

set.seed(1)
for (i in 1:num_perm){
  perm_idx <- sample(1:nrow(dat), nrow(dat))
  dat <- dat %>% mutate(y = count[perm_idx])
  z_tmp_score <- c()
  z_tmp_wald_robust <- c()
  z_tmp_wald <- c()
  
  fit_null <- glm(y~., family=poisson, data=dat, offset=offset_eta[perm_idx])
  
  for (j in 1:ncol(G)){
    # score test
    z <-  glm.scoretest(fit_null, G[,j])
    z_tmp_score <- c(z_tmp_score, z)
    
    # # wald robust
    # fit_full <- glm(y~., family=poisson, data=dat %>% mutate(g = G[,j]), offset=offset_eta[perm_idx])
    # se2 <- diag(sandwich::vcovHC(fit_full, type = "HC0"))['g']
    # z <- fit_full$coefficients['g']/sqrt(se2)
    # z_tmp_wald_robust <- c(z_tmp_wald_robust, z)
    # 
    # # wald not robust
    # res <- fit_full %>% broom::tidy() %>% filter(term == 'g')
    # z_tmp_wald <- c(z_tmp_wald, res$statistic)
  } 

  max_z_score <- c(max_z_score, max(abs(z_tmp_score)))
  # max_z_wald_robust <- c(max_z_wald_robust, max(abs(z_tmp_wald_robust)))
  # max_z_wald <- c(max_z_wald, max(abs(z_tmp_wald)))

  print(i)
}

hist(max_z_score)
p <- pnorm(abs(max_z_score), lower.tail = F)*2
hist(p)

dat %>% mutate(g = G[,2]) %>% 
  ggplot(aes(x = as.factor(g), y = y)) + geom_boxplot()


plot(z_tmp_wald, z_tmp_wald_robust)
plot(z_tmp_wald, z_tmp_score)

hist(count) # not zero inflated
hist(max_z_score)
hist(min_p_wald_robust)
hist(max_z_wald_robust)

plot(z_tmp_score, z_tmp_wald)

table(G[,77])
# perm_idx <- sample(1:nrow(dat), nrow(dat))
perm_idx <- 1:length(count)
dat <- dat %>% mutate(y = count[perm_idx],g = G[,77])

dat %>% 
  # filter(y < 100) %>%
  ggplot(aes(x = as.factor(g), y = y)) + geom_boxplot()
summary(glm(y~., family=poisson, data=dat, offset=offset_eta[perm_idx]))


summary(MASS::glm.nb(y~. + offset(offset_eta[perm_idx]), data=dat,
                     control=glm.control(maxit=1000)))


# 100 permutation for p_score_min

p = seq(0,1, length=100)
plot(p, dbeta(p, 0.6, 200), ylab="density", type ="l", col=4)
