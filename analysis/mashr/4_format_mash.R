## prepare mash results

library(susieR);library(mashr)
thresh <- 0.05

# load("../result/mvsusie/result/mashr_Uk_fit_EZ_TRUE_complete_TRUE.RData")

# check Z score for selected eqtl
maxZ <- finemap_eqtl_Z %>%  # finemap_cslead_eqtl_Z
  mutate(Z = slope/slope_se) %>% 
  group_by(eqtl) %>%
  slice_max(abs(Z)) %>% 
  pull(Z)
mean(abs(maxZ)>4)

m_finemap_res <- m_finemap_highpip # m_finemap_cslead
SE <- data.highpip$rawdata$se; anyNA(SE) # data.cslead

# 335 eQTL-cell type has very large SE, occurs more often in CD4_SOX4, Plasma, and NK_R
apply(SE, 2, function(x) {sum(x>10)})
sum(as.vector(SE)>10)

# local false discovery rate
lfsr <- get_lfsr(m_finemap_res); dim(lfsr) # 8994 x 14 (no negative values)

pm.mash.Z <- get_pm(m_finemap_res); # Z score scale from EZ model
summary(as.vector(pm.mash.Z))
cor(as.vector(data.highpip$Z), as.vector(pm.mash.Z)) # largely consistent with input Z score

pm.mash.beta <- get_pm(m_finemap_res) * SE; # convert to beta scale
plot(as.vector(data.highpip$data$Bhat), as.vector(pm.mash.beta))

# SE scales with slope estimates
data.highpip$rawdata$beta[pm.mash.beta > 10]
data.highpip$rawdata$se[pm.mash.beta > 10]
plot(data.highpip$rawdata$beta[pm.mash.beta > 10], data.highpip$rawdata$se[pm.mash.beta > 10])

# call significant snp
sigmat  <- (lfsr<=thresh); anyNA(sigmat)
mash_sig_eqtls <- names(get_significant_results(m_finemap_res))
nsig    <-  rowSums(sigmat); anyNA(nsig)
mean(nsig>0); sum(nsig>0) # 97% eQTL are significant in at least 1 cell type

mean(nsig==1) # 4% in only one cell type
mean(nsig==14) # 40% in 14 cell type

nsig_mash <- sum(nsig>0); nsig_mash # 8764/8994 (97%) significant
neqtl <- nrow(lfsr);neqtl # 1038 with PIP >0.9


het.func <- function(normdat, threshold) {
  apply((normdat),1,function(x){sum(x > threshold)})
}

het.norm <- function(effectsize) {
  t(apply(effectsize,1,function(x){
    x/x[which.max(abs(x))]
  }))
}

shared_mag_ct <- het.func(het.norm(effectsize=pm.mash.beta[nsig>0, ]),threshold=0.5)
anyNA(shared_mag_ct)

# shared by sign across all cell types
shared_sign_df <- as_tibble(het.norm(effectsize=pm.mash.beta[nsig>0, ])>0) %>% 
  mutate(eqtl = rownames(pm.mash.beta[nsig>0, ])) %>% 
  select(eqtl, everything())

mean(rowSums(shared_sign_df[, 2:15]) == 14)
# shared by magnitude across all cell types
shared_mag_ct <- het.func(het.norm(effectsize=pm.mash.beta[nsig>0, ]),threshold=0.5)
anyNA(shared_mag_ct)

shared_mag <- het.norm(effectsize=pm.mash.beta[nsig>0, ]) > 0.5 # TRUE or FALSE
mean(rowSums(shared_mag)==14) # 15% 
mean(rowSums(shared_mag)==1) # 9% 

shared_mag_df <- as_tibble(shared_mag) %>% 
  mutate(eqtl = rownames(shared_mag),
         phenotype_id = gsub("_chr.*","",eqtl)) %>% 
  select(eqtl, everything())


# identify specific (==1) or shared (>1) by magnitude 
shared_mag_df_test <- shared_mag_df %>% select(-c(eqtl, phenotype_id))
shared_mag_df_test <- sigmat[nsig>0,] & shared_mag_df_test
shared_mag_df_test <- bind_cols(shared_mag_df %>% select(c(eqtl, phenotype_id)), 
                                shared_mag_df_test)

shared_mag_df_test %>% gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl) %>% 
  summarize(n_sig = sum(sig)) %>% 
  count(n_sig)

specific_eqtl_mash <- shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% ungroup() %>%
  filter(n_shared==1 & shared == TRUE) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  select(phenotype_id, eqtl, celltype, GeneSymbol)

shared_eqtl_mash <- shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% ungroup() %>%
  filter(n_shared>1 & shared == TRUE) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  distinct(phenotype_id, eqtl, GeneSymbol)

# calculate distance
dist_mash <- ss_df %>% filter(eqtl %in% specific_eqtl_mash$eqtl) %>% 
  distinct(eqtl, tss_distance) %>% pull(tss_distance) %>% abs()
dist_raw <- specific_eqtl_raw %>% distinct(eqtl, dist) %>% 
  pull(dist) %>% abs()
wilcox.test(dist_mash, dist_raw)
t.test(dist_mash, dist_raw)

summary(dist_raw)
summary(dist_mash)

# # create bed file
# snp_df %>% distinct(chr, pos_1, pos) %>% mutate(chr = paste0("chr", chr)) %>% 
#   write_tsv("./eqtl_pip0.5.bed", col_names = F)

# snp_df %>%
#   filter(eqtl %in% mash_sig_eqtls) %>%
#   distinct(chr, pos_1, pos) %>%
#   write_tsv("mash_sig_eqtls.bed", col_names = F)


