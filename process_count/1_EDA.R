# distribution of count data
library(tidyverse)
library(data.table)

## distribution of NK cell
genelist <- read_tsv("./cis_mapping/NK_chr17_genelist.tsv", F) %>% pull(X1)
dat <- read_tsv("../data/OneK1K/pheno/NK_sum/natural_killer_cell.bed.gz")
intercept <- read_tsv("../result/mean_expr/natural_killer_cell.exp_intercept.allgenes.gz")
mu <- fread("../result/mean_expr/natural_killer_cell.mu.allgenes.gz")

# 23,419 expressed genes
all_sumstats <- tibble(mean_fitted = intercept$exp_intercept, 
                       var_fitted = apply(mu, 2, var),
                       dispersion = (var_fitted - mean_fitted)/mean_fitted^2,
                       expected_var = mean_fitted + dispersion * mean_fitted^2) %>% 
  filter(mean_fitted != 0)

summary(all_sumstats)

# adjust for library size
all_sumstats %>% 
  filter(mean_fitted < 0.09) %>%
  ggplot(aes(x = mean_fitted, y = var_fitted)) +
  geom_point()+
  geom_line(aes(y = expected_var), color = "blue", size = 0.5) + 
  # geom_abline(slope=1, intercept = 0)+
  theme_bw()+
  ggtitle("var(y) vs. mean(y), 23,417 expressed genes")


summary(all_sumstats$dispersion)

length(unique(dat$phenotype_id)) # total 27018 for NK cells 


# not adjust for library size
all_sumstats <- dat %>% 
  gather(key = iid, value = ct, `1000_1001`:`9_9`) %>% 
  group_by(phenotype_id) %>% 
  summarize(mean = mean(ct), var = var(ct), dispertion = (var - mean)/(mean^2)) %>% 
  filter(mean > 0)

med_disp <- median(all_sumstats$dispertion)

summary(all_sumstats$dispertion)
all_sumstats %>% 
  # filter(mean < 10000) %>%
  ggplot(aes(x = dispertion)) + geom_histogram()+ 
  geom_vline(xintercept = med_disp, color = "blue")

sumstats %>% 
  filter(mean < 10000) %>%
  ggplot(aes(x = mean, y = var)) +
  geom_point()+
  geom_abline(slope=1, intercept = 0)+
  theme_bw()+
  ggtitle("var(y) vs. mean(y), 27,017 genes")


dat143 <- dat %>% filter(phenotype_id %in% genelist)

