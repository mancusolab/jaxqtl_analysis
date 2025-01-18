# genes found by tqtl only (not at margin of qvalue threshold)
# 1657 tqtl only egenes (1609 converged, 1653 beta converged)
# remove unconverged: found 1605 eQTLs
NB_all_df %>% 
  inner_join(tqtl_only %>% 
               select(phenotype_id, celltype, chi2_tqtl=chi2, variant_id,
                      qval_tqtl=qval, beta1_tqtl=beta_shape1, beta2_tqtl=beta_shape2), 
             by=c("phenotype_id", "celltype")) %>% 
  filter(phenotype_id == "ENSG00000144820")
  ggplot(aes(x = chi2_tqtl, y = chi2, color = express_percent)) + geom_point() +
  # ggplot(aes(x = beta2_tqtl, y = beta_shape2, color = express_percent)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1)


# use example: ENSG00000256747 (chr 12)

library(tidyverse)
library(data.table)
library(ggbeeswarm)

options(nwarnings = 10000)

celltype <- "B_IN"
gene <- "ENSG00000144820"; # ENSG00000256747 (12); ENSG00000162520 (1); ENSG00000144820:896

raw_y <- fread(paste0("../data/OneK1K/pheno/celltype16/", celltype, ".bed.gz")) %>% 
  gather(key = iid, value = count, `1000_1001`:`9_9`) %>% 
  group_by(iid) %>% 
  mutate(offset = sum(count)) %>% 
  ungroup()

tmm_y <-fread(paste0("../data/OneK1K/pheno/celltype16/", celltype, ".tmm.bed.gz")) %>% 
  gather(key = iid, value = count, `1000_1001`:`9_9`) %>% 
  group_by(iid) %>% 
  mutate(offset = sum(count)) %>% 
  ungroup()

gene <- "ENSG00000182866"
tmm_y %>% filter(Geneid == gene) %>% ggplot(aes(x = count)) + geom_histogram()
raw_y %>% filter(Geneid == gene) %>% ggplot(aes(x = count/offset)) + geom_histogram()
  
raw_y %>% filter(Geneid == gene) %>% 
  mutate(rate = count/offset) %>% 
  left_join(tmm_y %>% filter(Geneid == gene) %>% 
              select(iid, Geneid, tmm_ct=count), by = c("iid", "Geneid")) %>% 
  # filter(tmm_ct < 1) %>%
  ggplot(aes(x = rate, y=tmm_ct)) + geom_point()

iid_order <- fread("../data/test/iid_order.tsv")

X <- fread(paste0("../data/test/", celltype, ".X.tsv.gz")) %>% select(-V1) %>% 
  mutate(iid = iid_order$iid) %>% 
  as_tibble()

# prs <- fread("../data/test/prs.tsv", header =TRUE)

# X <- X %>% left_join(prs, by="iid")

G <- fread(paste0("../data/test/", gene, ".G.tsv.gz")) %>% 
  mutate(iid = iid_order$iid) %>% 
  as_tibble()

snp_meta <- fread(paste0("../data/test/", gene, ".var.tsv.gz")) %>% 
  mutate(beta = NA, se = NA, converge = NA, loglik = NA)

y <- raw_y %>% filter(Geneid == gene)
dat <- X %>% left_join(y %>% select(count, iid, offset), by = "iid")


for (i in 1:(ncol(G)-1)){
  
  df <- dat %>% mutate(g = G %>% pull(paste0("V",i)),
                       rate = count/offset)
  df %>% ggplot(aes(x = as.factor(g), y = count)) + geom_quasirandom() + axis_theme+
    ggtitle("normalized counts")
  boxplot(rate ~ g, data=df %>% mutate(rate=count/offset))
  
  oneres <- MASS::glm.nb(count ~ V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + g+ offset(log_offset), 
                         data = df %>% mutate(log_offset=log(offset)), control=glm.control(maxit=100))
  
  tidyres <- broom::tidy(oneres)
  snp_meta$beta[i] <- tidyres %>% filter(term == "g") %>% pull(estimate)
  snp_meta$se[i] <- tidyres %>% filter(term == "g") %>% pull(std.error)
  snp_meta$loglik[i] <- -logLik(oneres)
  snp_meta$converge <- oneres$converged
  
  print(i)
  
  df <- df %>% select(-g)
}

View(snp_meta)
jqtl_pois %>% filter(phenotype_id == gene) %>% 
  mutate(beta_r = snp_meta$beta) %>% filter(abs(beta_r - beta_glm)>0.01) %>% 
  ggplot(aes(x = beta_glm, y = beta_r)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "red") + axis_theme

jqtl_pois %>% filter(phenotype_id == gene) %>% 
  mutate(beta_r = snp_meta$beta,
         beta_se = snp_meta$se,
         Z_r = snp_meta$beta/snp_meta$se,
         beta_diff = abs(beta_r - beta_glm)) %>% 
  ggplot(aes(x = beta_r, y = beta_glm)) +
  # ggplot(aes(x = Z_r, y = beta_glm/se_glm)) +
  geom_point(aes(color = se_glm), size=0.7) +
  scale_colour_gradient(low = "skyblue", high = "royalblue", na.value = NA) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")+
  axis_theme +
  ggtitle(paste0("Beta of cis-snps in ", gene, " using R and jaxqtl"))

ggsave2(paste0(plots_dir, "beta_glm_pois_R_jaxqtl.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)
