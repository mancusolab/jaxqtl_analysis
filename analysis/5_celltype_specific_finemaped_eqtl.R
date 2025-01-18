# examine cell type specificity 

library(corrplot)
library(ggplotify)
library(mvsusieR)
library(gplots)

# fine-map results
plots_dir <- "../result/cis/figures/celltype_specific/finemapped/"

finemap_indir <- "../result/finemap/result_wald/"
cols_names <- c("chr", "pos_1", "pos", "snp", "pip", "phenotype_id", "celltype")

load("../result/mvsusie/result/mashr_Uk_fit_EZ_TRUE.RData")

snp_annot <- read_tsv("../result/mvsusie/result/all_intersect.tsv", F) %>% 
  select(chr=X1, pos=X3, snplist=X7, celltype_annot=X8)

# take the strongest snp from each CS
nb_cs_leadsnp <- nb_cs %>%
  filter(celltype != "allcells") %>% 
  group_by(celltype, phenotype_id) %>%
  slice_max(pip, n=1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(eqtl = paste0(phenotype_id, "_", snp)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")
summary(nb_cs_leadsnp$pip)

# 9966 SNP-gene pairs
finemap_cslead_eqtl_Z <- finemap_cslead_eqtl_Z %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")

finemap_cslead_eqtl_Z %>% distinct(phenotype_id, snp)

nb_highpip <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap_pip0.5.tsv.gz"), 
                         col_names = cols_names) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  mutate(eqtl = paste0(phenotype_id,"_", snp)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, tss=end), by="phenotype_id") %>% 
  filter(celltype != "allcells") #%>% 
  #filter(pip >= 0.9)

length(unique(nb_highpip$eqtl)) # 1038 unique SNP-gene 
finemap_eqtl_Z %>% distinct(phenotype_id, snp) # 1038 unique SNP-gene

# finemap_cslead_eqtl_Z, finemap_eqtl_Z
# nb_cs_leadsnp, nb_highpip
snp_df <- nb_cs_leadsnp
ss_df <- finemap_cslead_eqtl_Z
suffix <- "cslead" # "pip0.9", "cslead"

specific_eqtl_raw <- snp_df %>% 
  add_count(eqtl) %>% 
  filter(n < 2)

specific_eqtl_raw %>% 
  distinct(chr,pos_1, pos) %>% 
  arrange(chr, pos) %>% 
  write_tsv("sp_snp_before_mash.bed", col_names = F)

# 48/530; 
specific_eqtl_raw %>% filter(celltype == "Mono_C") %>% distinct(eqtl)

# 84%
length(unique(specific_eqtl_raw$eqtl))/length(unique(snp_df$eqtl))

shared_eqtl_raw <- snp_df %>% 
  add_count(eqtl) %>% 
  filter(n > 1)
length(unique(shared_eqtl$eqtl))

shared_eqtl_raw %>% 
  distinct(chr,pos_1, pos) %>% 
  arrange(chr, pos) %>% 
  write_tsv("shared_snp_before_mash.bed", col_names = F)

sp_shared_intersect <- intersect(unique(specific_eqtl_raw$phenotype_id), 
                                 unique(shared_eqtl_raw$phenotype_id))

# 2% , 1% missing from finemap summary stats
n_eqtl <- length(unique(snp_df$eqtl));n_eqtl
((n_eqtl * 14) - nrow(ss_df))/nrow(ss_df)

# all eQTLs have >= 4 cell types observed
ss_df %>% count(eqtl) %>% pull(n) %>% table()

#### shared eGenes ####
n_specific <- length(unique(specific_eqtl_raw$phenotype_id)); n_specific

# all these studied genes are expressed in at least 2 cell types
p_sp_express <- allgenes_express %>% 
  filter(phenotype_id %in% specific_eqtl_raw$phenotype_id) %>% 
  group_by(phenotype_id) %>% 
  summarize(coverage_1per = sum(express_percent >= 0.01),
            coverage_10per = sum(express_percent >= 0.1)) %>% 
            #gtex_threshold = sum(tpm_0.1 >= 0.2 & express_percent_6 >= 0.2)) %>% 
  gather(key = thresh, value = ct, -phenotype_id) %>% 
  count(thresh, ct) %>% 
  mutate(prop = n/n_specific) %>% rename(threshold = thresh) %>% 
  add_row(threshold = "coverage_1per", ct = 0, n = 0, prop = 0) %>% # for share space with bar x = 0
  add_row(threshold = "coverage_10per", ct = 1, n = 0, prop = 0) %>% 
  ggplot(aes(x = ct, y = prop, fill=threshold)) + 
  geom_col(position = "dodge") + axis_theme +
  scale_x_continuous(breaks=seq(0,14, 1))+
  xlab("expressed cell types") + ylab("proportion of eGenes")+
  axis_theme + theme(legend.position = "bottom",
                     legend.key.size = unit("0.3", units = "cm"),
                     legend.text = element_text(size = 7)) 
p_sp_express
ggsave2(paste0(plots_dir, "specific_egenes_express.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)

#### shared eQTL between cell types ####

# plot shared eQTL between cell types
# 728 eQTLs found in only one cell type
# new 84%, 8326
snp_df %>% 
  count(phenotype_id, snp) %>% 
  count(n, name = "count_celltype") %>% 
  mutate(prop = count_celltype/sum(count_celltype))

# use different expression threshold
genes_1per <- allgenes_express %>% 
  filter(express_percent >= 0.01) %>% 
  count(phenotype_id) %>% 
  filter(n>13) %>% pull(phenotype_id)

genes_10per <- allgenes_express %>% 
  filter(express_percent >= 0.1) %>% 
  count(phenotype_id) %>% 
  filter(n>13) %>% pull(phenotype_id)

# genes_gtex <- allgenes_express %>% 
#   filter(tpm_0.1 >= 0.2 & express_percent_6 >= 0.2) %>% 
#   count(phenotype_id) %>% 
#   filter(n>13) %>% pull(phenotype_id)

# restrcit to genes expressed by some thresh across all 14 cell types
eqtl_cts <- snp_df %>% 
  count(phenotype_id, snp) %>% 
  mutate(thresh = "coverage_1per") %>% 
  bind_rows(snp_df %>% 
              filter(phenotype_id %in% genes_10per) %>% 
              count(phenotype_id, snp) %>% 
              mutate(thresh = "coverage_10per")) %>%
  count(thresh, n, name = "count_celltype") %>% 
  group_by(thresh) %>%
  mutate(n_eqtl = sum(count_celltype)) %>% 
  ungroup() %>% 
  mutate(prop = count_celltype / n_eqtl)
eqtl_cts
# 60-70% fine-mapped eQTLs are specific to one cell type
eqtl_cts %>% filter(n == 1) %>% mutate(ratio = prop/(1-prop))

p_eqtl_share <- eqtl_cts %>% rename(threshold=thresh) %>% 
  ggplot(aes(x = n, y = prop, fill = threshold)) +
  geom_bar(position = "dodge", stat="identity") + 
  scale_x_continuous(breaks=seq(1,14, 1))+
  xlab("number of cell types") + ylab("fraction of eQTLs")+
  axis_theme + theme(legend.position = "bottom",
                     legend.key.size = unit("0.3", units = "cm"),
                     legend.text = element_text(size = 7)) 

p_eqtl_share
legend <- get_legend(p_eqtl_share)
grid <- plot_grid(
  p_sp_express+theme(legend.position = "none"), 
  p_eqtl_share+theme(legend.position = "none"), 
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1,
  rel_widths = c(1, 1)
)

plot_grid(
  grid, legend, 
  align = 'h',
  nrow = 2,
  axis="b",
  rel_heights = c(1, 0.1)
)

ggsave2(paste0(plots_dir, "specific_express_bycutoff.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

n1_eqtl_cellcts <- specific_eqtl_raw %>% 
  count(celltype) %>% 
  left_join(cell_counts, by="celltype") %>% 
  mutate(prop = total_cell / sum(total_cell))

p_ct_spec <- n1_eqtl_cellcts %>% 
  ggplot(aes(x = prop, y = n)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  xlab("Cell type proportion") + ylab("count of cell-type-specific eQTLs")+
  axis_theme

plot_grid(
  p_eqtl_share, p_ct_spec,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1,
  rel_widths = c(1.5, 1)
)

ggsave2(paste0(plots_dir, "NB_shared_finemap_eqtl.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


#### specific eQTLs ~ number of cells ####
p_ct_spec
ggsave2(paste0(plots_dir, "specific_gene_total_cells.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

fit <- glm(n ~ prop, data=n1_eqtl_cellcts, family = "gaussian") # "poisson", "gaussian"
broom::tidy(fit)

# linear model R2: 91%
# poisson model R2: 87%
1-sum((fit$data$n - fit$fitted.values)^2)/sum((fit$data$n - mean(fit$data$n))^2)
fit$R

#### cor of specific gene across cell types ####
# show specific gene's expression in null cell type
specific_gene_list <- list()
for (cell in unique(specific_eqtl_raw$celltype)){
  specific_gene_list[[cell]] <- specific_eqtl_raw %>% filter(celltype == cell) %>% 
    distinct(phenotype_id) %>% pull(phenotype_id)
}

express_rate_mean_cor <- matrix(NA, nrow = 14, ncol = 14, dimnames=list(celltypes_grouped, celltypes_grouped))
express_rate_mean_test <- matrix(NA, nrow = 14, ncol = 14, dimnames=list(celltypes_grouped, celltypes_grouped))

for (i in 1:14){
  cell_row <- celltypes_grouped[i]
  gene_list <- specific_gene_list[[cell_row]]
  for (j in 1:14){
    cell_col <- celltypes_grouped[j]
    # test rate difference
    tmp <- allgenes_express_wide %>% filter(phenotype_id %in% gene_list)
    express_rate_mean_cor[i, j] <- cor(tmp[[cell_row]], tmp[[cell_col]], method = "spearman")
    
    express_rate_mean_test[i, j] <- wilcox.test(tmp[[cell_row]], tmp[[cell_col]], 
                                                alternative = "t") %>% 
      broom::tidy() %>% pull(p.value)
  }
}

# mean 87%, SD=6.8%
mean(c(express_rate_mean_cor[upper.tri(express_rate_mean_cor)],
       express_rate_mean_cor[lower.tri(express_rate_mean_cor)]), na.rm = T)
sd(c(express_rate_mean_cor[upper.tri(express_rate_mean_cor)],
       express_rate_mean_cor[lower.tri(express_rate_mean_cor)]), na.rm = T)

sum(express_rate_mean_test >= 0.5)
sum(express_rate_mean_cor<0, na.rm = T)

p <- as.grob(~heatmap.2(express_rate_mean_cor, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',
                        col = pal_red))

p <- as.grob(~corrplot(express_rate_mean_cor, method = c('color'), 
                               #type = "upper",
                               is.corr = F, 
                               col.lim = c(-1, 1),
                               col = pal_bry,
                               number.cex = 0.75,
                               tl.cex = 0.7,
                       diag = T,
                       tl.col = "black",
                       tl.srt = 45,
                       tl.offset = 0.6,
))

plot_grid(
  p,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="r",
  labels = c(""), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "specific_gene_express_cor.png"),
        width = onecol_width*2.1, height = onecol_height*2.1, units = 'cm', dpi = 300)

#### cochran ####

cochran <- function(effect, SE){
  wk <- 1/SE^2
  pooled_effect <- sum(effect * wk) / sum(wk)
  
  TS <- sum(wk * (effect - pooled_effect)^2)
  
  # follows chi2 (df = K-1)
  chi2_df <- length(wk) - 1
  pval <- pchisq(TS, df = chi2_df, lower.tail = FALSE)
  
  return(list(chi2 = TS, p = pval))
}

specific_eqtl$cochran_p <- NA
specific_eqtl$cochran_chi2 <- NA

for (i in 1:nrow(specific_eqtl)){
  gene <- specific_eqtl$phenotype_id[[i]]
  snp <- specific_eqtl$snp[[i]]
  
  tmp <- nb_cs_pip0.5_wald %>% 
    filter(eqtl == paste0(gene, "_", snp)) %>% 
    filter(!is.na(slope) & !is.na(slope_se))
  
  coch_res <- cochran(tmp$slope, tmp$slope_se)
  specific_eqtl$cochran_p[[i]] <- coch_res$p
  specific_eqtl$cochran_chi2[[i]] <- coch_res$chi2
}

# 433/728=59% show hetergeneous effect
nrow(specific_eqtl)
specific_eqtl %>% filter(cochran_p < 0.05/nrow(specific_eqtl))

# use all points to test distance ~ number of shared cell type
lm(n ~ abs(dist), data = snp_df %>%
     left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end)) %>% 
     mutate(dist = pos - tss) %>% 
     mutate(eqtl = paste0(phenotype_id, "_", snp)) %>% 
     count(eqtl, dist)) %>% broom::tidy()

p_share_dist <- snp_df %>%
  left_join(gene_lookup %>% select(phenotype_id, tss=end)) %>% 
  mutate(dist = pos - tss) %>% 
  count(eqtl, dist) %>% 
  group_by(n) %>% 
  summarize(mean = mean(abs(dist))) %>%
  ggplot(aes(x = n, y = mean)) + 
  geom_point() +
  scale_x_continuous(breaks=seq(1,14,1))+
  geom_smooth(method = "lm") + 
  ylab("mean distance to TSS") + xlab("number of shared cell types") + axis_theme
p_share_dist

#### correlation between effect size and gene expression across tissues ####
# remove SNPs when:
# 1) over half of tissue effect size in opposite direction of discovery cell type
# 2) any effect size in opposite direction and a magnitude > 0.5
# 3) extreme discovery effect size
# TODO: flip all effect size (*-1) if discovery effect size is negative


#### upset plot for shared eQTls ####
## 9966 unique SNP-egene, 4520 unique egenes
nb_cs %>%
  filter(celltype != "allcells") %>%
  left_join(nb_cs_ct %>%
              count(celltype, phenotype_id, name="cs_ct"),
            by=c("phenotype_id", "celltype")) %>%
  group_by(celltype, phenotype_id, cs_ct) %>%
  slice_max(pip, n=2, with_ties = TRUE) %>%
  ungroup() %>% # filter(cs_ct>1) %>% View
  distinct(phenotype_id, snp) # %>% 
 # write_tsv("cs_leadsnp.tsv")

nb_cs_ct %>% filter(celltype != "allcells") %>% count(celltype, phenotype_id)

tmp <- snp_df %>% 
  left_join(gene_lookup %>% select(phenotype_id, tss=end)) %>% 
  mutate(dist = pos - tss)

lt <- list()
for (cell in celltypes){
  lt[[cell]] <- unique(tmp$eqtl[which(tmp$celltype == cell)])
}

m <- make_comb_mat(lt)
# m <- make_comb_mat(lt, min_set_size = 50)
# m <- make_comb_mat(lt, min_set_size = 100, top_n_sets = 14)
m2 <- m[, comb_size(m) >= 20] # for nb_cs_pip0.9

p_upset <- grid.grabExpr(draw(UpSet(m2,
                                    comb_order = order(comb_size(m2), decreasing = TRUE),
                                    set_order = c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B",
                                                  "NK", "NK_R",
                                                  "B_Mem", "B_IN", "Plasma",
                                                  "Mono_C","Mono_NC", "DC"),
                                    top_annotation = upset_top_annotation(m2,
                                                                          axis_param = list(gp = gpar(fontsize = 5)),
                                                                          annotation_name_gp = gpar(fontsize = 5),
                                                                          height = unit(1, "cm")),
                                    right_annotation = upset_right_annotation(m2,
                                                                              axis_param = list(gp = gpar(fontsize = 5)),
                                                                              annotation_name_gp = gpar(fontsize = 5),
                                                                              width = unit(1, "cm")),
                                    row_names_gp = gpar(fontsize = 5))))

plot_grid(
  p_upset, p_ct_spec,
  align = 'h',
  hjust = 0.1, # -1
  nrow = 1,
  rel_widths=c(1.7, 1),
  labels = c("A", "B")
)

ggsave2(paste0(plots_dir, "upset_shared_finemap_eqtl_", suffix, ".png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

#### after mashr ####
library(susieR);library(mashr)
thresh <- 0.05
source("./mvsusie/normfuncs.R")

m_finemap_res <- m_finemap_cslead
SE <- data.cslead$rawdata$se

lfsr <- get_lfsr(m_finemap_res); dim(lfsr) # 9966 x 14 (no negative values)
pm.mash.Z <- get_pm(m_finemap_res); # Z score scale from EZ model
pm.mash.beta <- get_pm(m_finemap_res) * SE; # convert to beta scale
sigmat  <- (lfsr<=thresh); anyNA(sigmat)
nsig    <-  rowSums(sigmat); anyNA(nsig)
mean(nsig>1)
sum(nsig<2)/9658
sum(nsig==14)/9658

mean(abs(pm.mash.beta)>100, na.rm = T)
nsig_mash <- sum(nsig>0); nsig_mash # 939/1038 (90%) significant; 9658/9966 (97%)
neqtl <- nrow(lfsr);neqtl # 1038 with PIP >0.9

het.func <- function(normdat, threshold) {
  apply((normdat),1,function(x){sum(x > threshold, na.rm = TRUE)})
}

het.norm <- function(effectsize) {
  t(apply(effectsize,1,function(x){
    x/x[which.max(abs(x))]
  }))
}

shared_mag_ct <- het.func(het.norm(effectsize=pm.mash.beta[nsig>0, ]),threshold=0.5)
anyNA(shared_mag_ct)


shared_by_sign <- as_tibble(het.norm(effectsize=pm.mash.beta[nsig>0, ])>0) %>% 
  mutate(eqtl = rownames(pm.mash.beta[nsig>0, ])) %>% 
  select(eqtl, everything())

anyNA(shared_by_sign)
mean(rowSums(shared_by_sign[, -1], na.rm = T)==14)
mean(rowSums(shared_by_sign[, -1], na.rm = T)==1)
rowSums(shared_by_sign[,-1])

# shared by magnitude
# upset plot for shared magnitude
shared_mag <- het.norm(effectsize=pm.mash.beta[nsig>0, ]) > 0.5
shared_mag[is.na(shared_mag)] <- FALSE
mean(rowSums(shared_mag, na.rm = T)==14)

shared_mag_df <- as_tibble(shared_mag) %>% 
  mutate(eqtl = rownames(shared_mag),
         phenotype_id = gsub("_chr.*","",eqtl))

# shared_mag_df or shared_by_sign
tmp <- shared_mag_df %>% #mutate(eqtl = gsub("_b37$", "", eqtl)) %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(dist = abs(pos - tss)) %>% 
  gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl, dist) %>% 
  summarize(n_sig = sum(sig, na.rm = T)) %>% 
  ungroup() 

ggsave2(paste0(plots_dir, "shared_raw_tss.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

tmp %>% 
  filter(sig == TRUE) %>% 
  distinct(eqtl, n_sig, celltype) %>% 
  filter(n_sig == 1) %>% 
  pull(celltype) %>% table()

tmp %>% 
  ggplot(aes(x = as.factor(n_sig), y = abs(dist))) +
  geom_boxplot()

boxplot(as.factor(tmp$n_sig), abs(tmp$dist))
# shared by mag: p = 0.016 Poisson
# shared by direction: p = 8.54E-4
glm(dist ~ n_sig, data=tmp, family = "gaussian") %>% broom::tidy()
glm(n_sig ~ dist, data=tmp, family = "poisson") %>% broom::tidy()

glm(mean ~ n_sig, weights = wts, data=tmp %>% 
      group_by(n_sig) %>% summarize(mean = mean(dist),
                                    wts = 1/(var(dist))),
    family = "gaussian") %>% broom::tidy()

glm(n_sig ~ mean, weights = wts, data=tmp %>% 
      group_by(n_sig) %>% summarize(mean = mean(dist),
                                    wts = 1/(var(dist))),
    family = "poisson") %>% broom::tidy()

# closer to TSS is more likely to have shared magnitude; shared_by_sign; shared_mag_df
# shared_by_sign %>% 
p_share_sign_dist <- tmp %>% 
  group_by(n_sig) %>% summarize(mean = mean(dist),
                                wts = 1/(var(dist)),
                                n = n()) %>% 
  ungroup() %>% 
  mutate(mean = mean / 1000) %>%
  ggplot(aes(x = n_sig, y = mean)) + geom_point() + 
  geom_smooth(method="lm", mapping = aes(weight = wts)) +
  ylim(c(40, 80)) +
  axis_theme + xlab("number of cell types") + ylab("mean distance to TSS (kb)")+
  ggtitle("shared by sign")
p_share_sign_dist

p_share_mag_dist <- tmp %>% 
  group_by(n_sig) %>% summarize(mean = mean(dist),
                                wts = 1/(var(dist))) %>% 
  ungroup() %>% 
  mutate(mean = mean / 1000) %>%
  ggplot(aes(x = n_sig, y = mean)) + geom_point() + 
  geom_smooth(method="lm", mapping = aes(weight = wts)) +
  ylim(c(0, 80)) +
  axis_theme + xlab("number of cell types") + ylab("mean distance to TSS (kb)")+
  ggtitle("shared by magnitude")
p_share_mag_dist

tmp <- snp_df %>%
  count(eqtl, dist) %>% 
  group_by(n) %>% 
  summarize(mean_dist = mean(dist),
            wts = 1/var(dist),
            wts = ifelse(!is.finite(wts), 1, wts),
            n_ct = n()) %>% 
  ungroup()

glm(n ~ mean_dist, weights = wts, data=tmp,
    family = "poisson") %>% broom::tidy()

p_share_raw_dist <- tmp %>% 
  mutate(mean_dist = mean_dist / 1000) %>% 
  ggplot(aes(x = n, y = mean_dist)) + geom_point() +
  geom_smooth(aes(weight = wts), method = "lm") + 
  ggtitle("shared without mash") + ylim(c(0,80))+
  axis_theme

plot_grid(
  p_share_raw_dist, p_share_mag_dist,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  labels = c("A", "B")
)

ggsave2(paste0(plots_dir, "shared_sign_mag_tss.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# calculate distance bin

p_dist_percent_raw <- snp_df %>%
  count(eqtl, dist, name = "n_celltype") %>% 
  mutate(dist = ifelse(dist < 1, dist+1, dist),
         bins = cut(dist, breaks = c(0, 5000, 100000, 500000))) %>% 
  group_by(n_celltype) %>% 
  count(bins) %>% 
  mutate(n_sig_ct = sum(n)) %>% 
  ungroup() %>% 
  mutate(prop = n / n_sig_ct) %>% 
  ggplot(aes(fill=bins, y=prop, x=n_celltype)) +
  geom_bar(position="fill", stat="identity") + 
  scale_x_continuous(breaks = seq(from = 1, to = 14, by = 1))+
  ggtitle("raw estimates")+
  axis_theme


p_dist_percent_mash <- tmp %>% 
  mutate(dist = ifelse(dist < 1, dist+1, dist),
         bins = cut(dist, breaks = c(0, 5000, 100000, 500000))) %>% 
  group_by(n_sig) %>% 
  count(bins) %>% 
  mutate(n_sig_ct = sum(n)) %>% 
  ungroup() %>% 
  #filter(!is.na(bins)) %>% 
  mutate(prop = n / n_sig_ct) %>% 
  rename(n_celltype = n_sig) %>% 
  ggplot(aes(fill=bins, y=prop, x=n_celltype)) +
  geom_bar(position="fill", stat="identity") + 
  scale_x_continuous(breaks = seq(from = 1, to = 14, by = 1))+
  ggtitle("mash estimates")+
  axis_theme

legend <- get_legend(p_dist_percent_raw + theme(legend.position = "bottom"))
grid1 <- plot_grid(
  p_dist_percent_raw + theme(legend.position = "none"), 
  p_dist_percent_mash + theme(legend.position = "none"),
  align = 'h',
  hjust = 0, # -1
  nrow = 1,
  labels = c("A", "B")
)
plot_grid(
  grid1, legend,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.1)
)

ggsave2(paste0(plots_dir, "shared_mag_ct_dist_bins.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# significant eQTL
as_tibble(sigmat) %>% mutate(eqtl = rownames(sigmat)) %>% 
  gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl) %>% 
  summarize(shared_mash = sum(sig, na.rm = T)) %>%
  inner_join(snp_df %>% 
               mutate(eqtl = paste0(phenotype_id, "_", snp)) %>% 
               count(eqtl, name="shared")) %>% 
  gather(key = group, value = ct, shared_mash:shared) %>% 
  count(group, ct) %>% 
  mutate(prop = n/neqtl) %>% 
  add_row(group = "shared", ct = 13, n = 0, prop = 0) %>%
  filter(ct != 0) %>% 
  ggplot(aes(x = ct, y = prop, fill = group)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("shared"="grey", "shared_mash"="black"))+
  scale_x_continuous(breaks=seq(1,14, 1))+
  xlab("number of cell types") + ylab("fraction of eQTLs")+
  axis_theme+theme(legend.position = "bottom",
                   legend.key.size = unit("0.3", units = "cm"),
                   legend.text = element_text(size=6))

ggsave2(paste0(plots_dir, "shared_eqtl_sig_ct_before_after_mash.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

# 89/98 eQTL are specific eQTLs
shared_by_sign %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  group_by(eqtl) %>% 
  summarize(ct = sum(shared, na.rm = T)) %>%
  mutate(n_eqtl = n(), group = "shared_by_sign") %>% 
  bind_rows(snp_df %>% 
               count(eqtl, name = "ct") %>% 
               mutate(n_eqtl = n(), group="shared")) %>% 
  bind_rows(tibble(eqtl = names(shared_mag_ct), 
                   ct = shared_mag_ct,
                   n_eqtl = length(unique(names(shared_mag_ct))),
                   group = "shared_by_magnitude")) %>% 
  count(group, ct, n_eqtl) %>% 
  mutate(prop = n/n_eqtl) %>% View

p1 <- tibble(group=c("specific", "shared_by_sign", "shared_by_magnitude"),
       n_eqtl = c(length(unique(snp_df$eqtl)),
                  length(unique(shared_by_sign$eqtl)),
                  length(shared_mag_ct))) %>% 
  mutate(group = fct_relevel(group, "specific")) %>% 
  ggplot(aes(x = group, y = n_eqtl, fill = group)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("specific"="grey", "shared_by_sign"="red3", "shared_by_magnitude"="#4575B4"))+
  ylab("number of eQTLs") + xlab("")+
  axis_theme+theme(legend.position = "bottom",
                   legend.key.size = unit("0.3", units = "cm"),
                   legend.text = element_text(size=6),
                   axis.text.x = element_blank())

p2 <- shared_by_sign %>% gather(key = celltype, value = specific, B_IN:Plasma) %>% 
  group_by(eqtl) %>% 
  summarize(ct = sum(specific, na.rm = T)) %>%
  mutate(n_eqtl = n(), group = "shared_by_sign") %>% 
  bind_rows(snp_df %>% 
              count(eqtl, name = "ct") %>% 
              mutate(n_eqtl = n(), group="specific")) %>% 
  bind_rows(tibble(eqtl = names(shared_mag_ct), 
                   ct = shared_mag_ct,
                   n_eqtl = length(unique(names(shared_mag_ct))),
                   group = "shared_by_magnitude")) %>% 
  count(group, ct, n_eqtl) %>% 
  mutate(prop = n/n_eqtl) %>% 
  mutate(group = fct_relevel(group, "specific")) %>% 
  # add_row(group = "shared_by_sign", ct = 1, n = 0, prop = 0) %>% 
  # add_row(group = "specific", ct = 13, n = 0, prop = 0) %>% 
  ggplot(aes(x = ct, y = prop, fill = group)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("specific"="grey", "shared_by_sign"="red3", "shared_by_magnitude"="#4575B4"))+
  scale_x_continuous(breaks=seq(1,14, 1))+
  xlab("number of cell types") + ylab("fraction of eQTLs")+
  axis_theme+theme(legend.position = "bottom",
                   legend.key.size = unit("0.3", units = "cm"),
                   legend.text = element_text(size=6))
  
p2

legend <- get_legend(p1)
grid <- plot_grid(
  p1+theme(legend.position = "none"), 
  p2+theme(legend.position = "none"), 
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1,
  rel_widths = c(1, 2)
)

plot_grid(
  grid, legend, 
  align = 'h',
  nrow = 2,
  axis="b",
  rel_heights = c(1, 0.1)
)
ggsave2(paste0(plots_dir, "shared_eqtl_ct_before_after_mash.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

lt <- list()
for (cell in celltypes){
  lt[[cell]] <- shared_mag_df$eqtl[shared_mag_df[[cell]] == TRUE]
}

m <- make_comb_mat(lt)
m2 <- m[, comb_size(m) >= 40] # for nb_cs_pip0.9

p_upset <- grid.grabExpr(draw(UpSet(m2,
                                    comb_order = order(comb_size(m2), decreasing = TRUE),
                                    set_order = c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET", "CD8_S100B",
                                                  "NK", "NK_R",
                                                  "B_Mem", "B_IN", "Plasma",
                                                  "Mono_C","Mono_NC", "DC"),
                                    top_annotation = upset_top_annotation(m2,
                                                                          axis_param = list(gp = gpar(fontsize = 5)),
                                                                          annotation_name_gp = gpar(fontsize = 5),
                                                                          height = unit(1, "cm")),
                                    right_annotation = upset_right_annotation(m2,
                                                                              axis_param = list(gp = gpar(fontsize = 5)),
                                                                              annotation_name_gp = gpar(fontsize = 5),
                                                                              width = unit(1, "cm")),
                                    row_names_gp = gpar(fontsize = 5))))

plot_grid(
  p_upset,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1, 1), rel_widths=c(1, 1, 1)
)

ggsave2(paste0(plots_dir, "upset_shared_bymag_finemap_eqtl.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

## shared eGene (by magnitude)
# no intersection between specific_genes.tsv and shared_egenes.tsv
shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(eqtl, phenotype_id, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  summarize(n_shared = sum(shared)) %>% pull(n_shared) %>% table()
  filter(n_shared>1) %>% distinct(phenotype_id) %>% write_tsv("shared_genes.tsv", col_names = F)

sp_beta <- shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(eqtl, phenotype_id, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% 
  ungroup() %>% 
  filter(n_shared == 1 & shared == TRUE) # %>% distinct(eqtl)

sp_beta %>% filter(celltype == "Mono_C")
1844/9658 # 20%
967/9658  # 10%
## before and after eQTL effect
finemap_eqtl_Z

as_tibble(pm.mash.beta) %>% 
  mutate(eqtl = rownames(pm.mash.beta)) %>% 
  gather(key = celltype, value = "pm_beta", B_IN:Plasma) %>% 
  left_join(finemap_cslead_eqtl_Z %>% select(eqtl, celltype, slope, slope_se), by=c("celltype", "eqtl")) %>% 
  inner_join(sp_beta %>% select(eqtl, celltype), by=c("eqtl", "celltype")) %>% 
  filter(abs(pm_beta) > 10) %>% 
  ggplot(aes(x = slope, y = pm_beta)) + geom_point()

finemap_cslead_eqtl_Z %>% filter(slope_se > 100) %>% View
  pull(celltype) %>% table()

#### GSEA ####
library(fgsea)

gmt_symbols <- gmtPathways("../data/OneK1K/annotation/gsea/C5_GO/c5.go.mf.v2024.1.Hs.symbols.gmt")
gmt_symbols <- gmtPathways("../data/OneK1K/annotation/gsea/C7_immune/c7.all.v2024.1.Hs.symbols.gmt")
gmt_symbols <- gmtPathways("../data/OneK1K/annotation/gsea/C3_reg/c3.tft.v2024.1.Hs.symbols.gmt")

shared_genelist <- shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  summarize(n_shared = sum(shared)) %>% # pull(n_shared) %>% table()
  filter(n_shared>1) %>% distinct(phenotype_id) %>% 
  left_join(gene_lookup %>% dplyr::select(GeneSymbol, phenotype_id)) %>% 
  distinct(GeneSymbol) %>% pull(GeneSymbol)

sp_genelist <- shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  summarize(n_shared = sum(shared)) %>% # pull(n_shared) %>% table()
  filter(n_shared<2) %>% distinct(phenotype_id) %>% 
  left_join(gene_lookup %>% dplyr::select(GeneSymbol, phenotype_id)) %>% 
  distinct(GeneSymbol) %>% pull(GeneSymbol)

shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% ungroup() %>% 
  filter(n_shared>1 & shared==TRUE) %>% distinct(eqtl)

shared_genelist
sp_genelist

shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  summarize(n_shared = sum(shared)) %>% ungroup() %>% 
  filter(n_shared>1) %>% separate(eqtl, into=c("gene", "chr", "pos", "ref", "alt"),sep="_") %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos = as.integer(pos),
         pos_1 = pos - 1) %>% 
  distinct(chr,pos_1, pos) %>% 
  arrange(chr, pos) %>% 
  write_tsv("shared_snp_after_mash.bed", col_names = F)

background_genes <- read_tsv("background_genes.txt", F) %>% rename(phenotype_id=X1) %>% 
  left_join(gene_lookup %>% select(GeneSymbol, phenotype_id)) %>% 
  distinct(GeneSymbol) %>% pull(GeneSymbol)
gsea_res <- fora(pathways = gmt_symbols, genes = genelist,universe=background_genes,
                 minSize = 15, maxSize = 500)

gsea_res %>% arrange(padj) %>% arrange(pval) %>% View
  filter(padj < 0.05) %>% View

names(gmt_symbols)
gmt_symbols[["SOX3_TARGET_GENES"]]

#### after mash ATAC ####
# here we try to verify specific eQTLs have cell type matched annotations
# looks like many effect can be found in promoter or enhancer region
atac_sp_sh <- read_tsv("../result/mvsusie/result/all_intersect_atac", F) %>% 
  select(chr=X1, pos=X3, celltype_annot=X7, group=X8)

atac_sp_sh <- read_tsv("../result/mvsusie/result/all_intersect_epimap", F) %>% 
  select(chr=X1, pos=X3, celltype_annot=X7, reg=X8, group=X9)

table(atac_sp_sh$group)

shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% ungroup() %>% 
  separate(eqtl, into=c("gene", "chr", "pos", "ref", "alt"),sep="_") %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos = as.integer(pos)) %>% 
  filter(n_shared == 1 & shared == TRUE) %>% 
  left_join(atac_sp_sh %>% filter(group == "sp_snp_after"), by=c("chr", "pos")) %>% 
  filter(!is.na(celltype_annot)) %>% 
  arrange(chr, pos) %>% 
  View

specific_eqtl %>% 
  left_join(atac_sp_sh %>% filter(group == "sp_snp_before"), by=c("chr", "pos")) %>% 
  filter(!is.na(celltype_annot)) %>% View

#### TF ####
snp_df %>% 
  left_join(gene_lookup %>% select(phenotype_id, strand)) %>% 
  mutate(score = 100) %>% 
  distinct(chr, pos_1, pos, eqtl, score, strand) %>% 
  arrange(chr, pos_1) # %>% 
  # write_tsv("../result/mvsusie/result/cs_lead.bed", col_names = F)

# download from hTFtarget: http://bioinfo.life.hust.edu.cn/hTFtarget/#!/
tf_df <- read_tsv("../data/OneK1K/annotation/TF/TF-Target-information.txt")
table(tf_df$tissue)

# blood: "blood", "Blood", "peripheral blood", 
# cell type: "lymphocyte", "cell:primary human memory B cells"
tf_list <- c()
for (tf in unique(tf_df$tissue)){
  new <- unique(str_split_1(tf, ","))
  tf_list <- c(tf_list, new[!new %in% tf_list])
}

tf_df %>% 
  mutate(tf_tissue = ifelse(grepl("blood|Blood|peripheral blood", tissue), "blood",
                            ifelse(grepl("lymphocyte", tissue), "T_Bcell",
                                   ifelse(grepl("cell:primary human memory B cells", tissue), "Bcell", NA)))) %>% 
  filter(!is.na(tf_tissue)) %>% filter(tf_tissue == "Bcell")
  distinct(TF, target) #%>% 
  # filter(tf_tissue == "blood")

# less TF factor target at those genes
# count #TF for those eGenes
tf_df <- tf_df %>% 
  mutate(tf_tissue = ifelse(grepl("blood|Blood|peripheral blood", tissue), "blood",
                            ifelse(grepl("lymphocyte", tissue), "T_Bcell",
                                   ifelse(
                                     grepl("cell:primary human memory B cells", tissue), "Bcell",
                                     NA)))) %>% 
  filter(!is.na(tf_tissue)) # %>%
  distinct(TF, target)

# no mash
specific_eqtl <- specific_eqtl_mash
shared_eqtl <- shared_eqtl_mash
sp_genelist <- specific_eqtl %>% pull(GeneSymbol) %>% unique()
shared_genelist <- shared_eqtl%>% pull(GeneSymbol) %>% unique()

to_rm <- intersect(sp_genelist, shared_genelist);length(to_rm)
sp_genelist <- sp_genelist[!sp_genelist %in% to_rm];length(sp_genelist)
shared_genelist <- shared_genelist[!shared_genelist %in% to_rm];length(shared_genelist)
tf_df_tmp <- tf_df %>%
  count(target) %>%
  inner_join(bind_rows(tibble(target = sp_genelist, group = "specific"),
                       tibble(target = shared_genelist, group = "shared")))

tf_df_tmp %>% group_by(group) %>% summarise(mean = mean(n))

# p: 1.01E-31 vs. 0.58
wilcox.test(tf_df_tmp %>% filter(group == "specific") %>% pull(n),
            tf_df_tmp %>% filter(group == "shared") %>% pull(n)) %>% 
  broom::tidy()

p_tf_old <- tf_df_tmp %>% 
  ggplot(aes(x = group, y = n)) + geom_boxplot() +
  ylab("number of TF") + xlab("group of genes")+ggtitle("raw estimates")+
  axis_theme

plot_grid(p_tf_old, p_tf_mash, 
          nrow = 1,
          labels = c("A", "B"))

ggsave2(paste0(plots_dir, "egene_mash_num_tf.png"),
       width = full_width, height = full_height, units = 'cm', dpi = 300)

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

pm.mash.beta.df <- as_tibble(pm.mash.beta) %>% mutate(eqtl = rownames(pm.mash.beta)) %>% 
  gather(key = celltype, value = pm_beta, B_IN:Plasma) %>% 
  group_by(eqtl) %>% mutate(ct_complete=sum(!is.na(pm_beta))) %>% ungroup()

df <- specific_eqtl_mash %>% 
  inner_join(tfbs %>% select(eqtl, TF)) %>% 
  distinct() %>% mutate(meta_TE = NA, meta_seTE = NA)

for (i in 1:nrow(df)){
  symbol_i <- df$GeneSymbol[[i]]
  eqtl_i <- df$eqtl[[i]]
  tf_tmp <- tf_df %>% filter(target == symbol)
  
  if (nrow(tf_tmp) > 0){
    TE_vec <- c()
    seTE_vec <- c()
    for (tf in unique(tf_tmp$TF)){
      tmp <- pm.mash.beta.df %>% filter(eqtl == eqtl_i) %>% 
        left_join(allgenes_express %>% 
                    left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
                    filter(GeneSymbol == tf), by="celltype")
      fit <- lm(abs(pm_beta) ~ rate_mean, data=tmp) %>% broom::tidy()
      TE_vec <- c(TE_vec, fit$estimate[2])
      seTE_vec <- c(seTE_vec, fit$std.error[2])
    }
    meta_res <- meta_fix(TE_vec, seTE_vec)
    df$meta_TE[i] <- meta_res$meta_TE
    df$meta_seTE[i] <- meta_res$meta_seTE
  }
  print(i)
}
View(df)

## TF-marker
tf_marker <- read_tsv("../data/OneK1K/annotation/TF/TFMarker.txt") %>% 
  janitor::clean_names() %>% 
  filter(cell_type == "Normal cell")

tf_marker %>% count(cell_name) %>% View
table(tf_marker$cell_name)

## ENCODE-TF
tfbs <- read_tsv("../result/mvsusie/result/all_intersect_tf_motif", F) %>% 
  select(chr=X1, pos=X3, eqtl=X4, TF=X10, start=X8, end=X9) %>%
  separate(TF, into=c("TF", "suffix"), sep='_', remove = F) %>% 
  distinct(chr, pos, eqtl, TF)

specific_eqtl <- specific_eqtl_mash
shared_eqtl <- shared_eqtl_mash
sp_genelist <- specific_eqtl %>% pull(GeneSymbol) %>% unique()
shared_genelist <- shared_eqtl%>% pull(GeneSymbol) %>% unique()

to_rm <- intersect(sp_genelist, shared_genelist);length(to_rm)

tfbs_tmp <- bind_rows(specific_eqtl %>% mutate(group = "specific"),
                      shared_eqtl %>% mutate(group = "shared")) %>% 
  filter(!GeneSymbol %in% to_rm) %>% 
  inner_join(tfbs, by="eqtl") %>% 
  distinct(GeneSymbol, TF, group) %>% 
  group_by(group) %>% 
  count(GeneSymbol) %>% ungroup()

tfbs_tmp %>% group_by(group) %>% summarize(mean = mean(n),
                                           med = median(n))

# mash: P = 0.03, mean sp 1.98 vs. shared 2.28, med 1 vs. 2
# raw: P = 0.51, mean 2.16 vs. 1.95; med 1 vs 2
wilcox.test(tfbs_tmp %>% filter(group == "specific") %>% pull(n),
            tfbs_tmp %>% filter(group == "shared") %>% pull(n)) %>% 
  broom::tidy()

length(unique(tfbs_tmp$GeneSymbol))
p_tf_mash <- tfbs_tmp %>% #View
  ggplot(aes(x = group, y = n)) + geom_boxplot() + axis_theme

plot_grid(p_tf_old, p_tf_mash, 
          nrow = 1,
          labels = c("A", "B"))

ggsave2(paste0(plots_dir, "egene_mash_num_tfbs_encode.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# TF expression and eQTL effect size

tf_df %>% filter(target == "ABRACL")
tmp <- shared_eqtl_mash %>% 
  inner_join(tfbs, by=c("eqtl")) %>% 
  distinct()
tmp %>% View

# ENSG00000119408_chr9_126985858_T_C_b37, SPI1
for (i in 1:nrow(tmp)){
TF_gene <- tmp$TF[i] # "NR2C2"
which_eqtl <- tmp$eqtl[i] # "ENSG00000174791_chr11_66104057_T_A_b37"
as_tibble(pm.mash.beta) %>% mutate(eqtl = rownames(pm.mash.beta)) %>% 
  gather(key = celltype, value = pm_beta, B_IN:Plasma) %>% 
  filter(eqtl == which_eqtl) %>% 
  left_join(allgenes_express %>% 
              left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
              filter(GeneSymbol == TF_gene), by="celltype") %>% 
  ggplot(aes(x = rate_mean, y = pm_beta, color = celltype)) + geom_point() +
  geom_smooth(method = "lm", formula = y~poly(x,1)) +
  scale_color_manual(values = celltype_cols) + 
  geom_hline(yintercept = 0, linetype = "dashed")+
  axis_theme+ggtitle(paste0(which_eqtl, ":", TF_gene)) +
  theme(legend.key.size = unit(0.3, units = "cm"))

ggsave2(paste0(plots_dir, "shared/", which_eqtl, "_", TF_gene, ".png"),
        width = onecol_width*1.3, height = onecol_height, units = 'cm', dpi = 300)
}
## selection parameter

s_tmp <- bind_rows(specific_eqtl %>% mutate(group = "specific"),
                      shared_eqtl %>% mutate(group = "shared")) %>% 
  filter(!GeneSymbol %in% to_rm) %>% 
  left_join(pLI_score %>% select(phenotype_id, pLI, loeuf), by="phenotype_id") %>%
  left_join(eds_score %>% select(phenotype_id, EDS), by="phenotype_id") %>%
  left_join(gene_lookup %>% select(phenotype_id, gene_type), by = "phenotype_id") %>% 
  left_join(rvis_score, by=c("phenotype_id")) %>% 
  distinct(group, phenotype_id, pLI, loeuf, EDS, RVIS)

s_tmp %>% group_by(group) %>% summarize(pLI = mean(pLI, na.rm = T),
                                        eds = median(EDS, na.rm = T),
                                        loeuf = mean(loeuf, na.rm = T),
                                        RVIS = mean(RVIS, na.rm = T))

# raw: pLI 0.0956, EDS 0.505, loeuf 0.927, RVIS 0.0737
# mash: pLI 0.126, EDS 0.00697, loeuf 0.683, RVIS 0.0231
wilcox.test(s_tmp %>% filter(group == "specific") %>% pull(EDS),
            s_tmp %>% filter(group == "shared") %>% pull(EDS)) %>% 
  broom::tidy()

length(unique(tfbs_tmp$GeneSymbol))
p_tf_RVIS <- s_tmp %>% #View
  ggplot(aes(x = group, y = RVIS)) + geom_boxplot() + axis_theme

plot_grid(p_tf_pLI, p_tf_EDS, p_tf_loeuf, p_tf_RVIS,
          nrow = 1,
          labels = c("pLI", "EDS", "loeuf", "RVIS"),
          label_size = 8,
          hjust = 0)

ggsave2(paste0(plots_dir, "egene_raw_selection.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

## GxE
ge <- read_tsv("../result/mvsusie/result/cASE_GxE_rep.tab") %>% filter(significant == TRUE) %>% 
  filter(cell.type %in% c("CD4_Tcells", "dendritic cells", "dendriticCells",
                          "monocytes", "PBMC", "PBMCs", "whole blood"))
table(ge$cell.type)

ge_tmp <- bind_rows(specific_eqtl %>% mutate(group = "specific"),
                   shared_eqtl %>% mutate(group = "shared")) %>% 
  filter(!GeneSymbol %in% to_rm) %>% 
  distinct(group, phenotype_id) %>% 
  inner_join(ge, by=c("phenotype_id" = "ensg"))

ge_tmp %>% group_by(group) %>% summarize(n=n())

## ctcf motif
ctcf <- read_tsv("../result/mvsusie/result/all_intersect_ctcf_motif", F) %>% 
  select(chr=X1, pos=X3, eqtl=X4, start=X8, end=X9, celltype_ctcf=X13) %>% 
  distinct()

specific_eqtl_mash %>% 
  inner_join(ctcf %>% select(eqtl, celltype_ctcf), by="eqtl") %>% 
  distinct() %>% 
  View

#### pairwise sharing ####
## intersection of finemapped snps / union of finemapped snps
share_M <- matrix(0, nrow = 14, ncol = 14)

# shared_mag_df
shared_mag_df_tmp <- shared_mag_df %>% gather(key = celltype, value = share, B_IN:Plasma) %>% 
  filter(share == TRUE)
# nb_highpip, 

share_M <- create_jaccard_mat(shared_mag_df_tmp, "eqtl")
share_M <- share_M[match(celltypes_grouped, rownames(share_M)), match(celltypes_grouped, rownames(share_M))]

p_share_M <- as.grob(~corrplot(share_M, method = c('color'), 
                               type = "upper",
                               is.corr = F, 
                               col.lim = c(0, 1),
                               col = pal_blue,
                               #addCoef.col = "black",
                               number.cex = 0.75,
                               tl.cex = 0.7,
                               tl.col = "black",
                               tl.srt = 45
))

plot_grid(
  p_share_M,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="l",
  labels = c(""), label_size = 10, vjust=1
)
ggsave(paste0(plots_dir, "pairwise_eqtl_share_mag_mash.png"),
       width = onecol_width*2, height = onecol_height*2, units = 'cm', dpi = 300)

#### replication ####

# eQTL gen
joined <- snp_df %>%
  separate(snp, into=c("rm1", "rm2", "ref", "alt"), remove = F) %>% 
  mutate(variant_id = paste0(chr, ":", pos)) %>% 
  count(phenotype_id, snp, variant_id, ref, alt) %>% 
  left_join(eqtlgen %>% dplyr::select(phenotype_id, variant_id, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
            by=c("phenotype_id", "variant_id"))
mean(!is.na(joined$Zscore_eqtlgen))

pm_beta_df <- as.tibble(pm.mash.beta[nsig>0, ]) %>% 
  mutate(eqtl = rownames(pm.mash.beta[nsig>0, ])) %>% 
  gather(key = celltype, value = slope, B_IN:Plasma)

pm_beta_df %>% filter(abs(slope) > 10) %>% 
  left_join(finemap_cslead_eqtl_Z %>% select(celltype, eqtl, raw_slope=slope, raw_slope_se=slope_se)) %>% 
  add_count(eqtl) %>% 
  arrange(eqtl) %>% filter(n==1) %>% 
  View

NB_all_df %>% filter(phenotype_id == "ENSG00000099864")
nb_cs_pip0.5_wald %>% filter(phenotype_id == "ENSG00000099864") %>% View
shared_mag_df_tmp %>%
  group_by(eqtl) %>% 
  mutate(n = sum(share)) %>% ungroup() %>%
  left_join(pm_beta_df, by=c("celltype", "eqtl")) %>%
  filter(abs(slope) < 20) %>%
  ggplot(aes(x = as.factor(n), y = abs(slope))) + geom_boxplot()

eqtlgen %>% filter(variant_id == "11:6343486" 
                   & phenotype_id == "ENSG00000170955")

nb_highpip %>% 
  left_join(finemap_eqtl_Z %>% select(celltype, eqtl, slope)) %>% 
  group_by(eqtl) %>% 
  mutate(n = n()) %>% ungroup() %>% 
  # left_join(pm_beta_df, by=c("celltype", "eqtl")) %>% 
  # filter(abs(slope) < 20) %>% 
  ggplot(aes(x = as.factor(n), y = abs(slope))) + geom_boxplot()

joined <- shared_mag_df_tmp %>% 
  group_by(eqtl) %>% 
  summarize(n = sum(share)) %>% ungroup() %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_") %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         variant_id = paste0(chr, ":", pos)) %>% 
  left_join(eqtlgen %>% dplyr::select(phenotype_id, variant_id, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
            by=c("phenotype_id", "variant_id"))
  

aligns <- allele.qc(joined$ref,joined$alt,joined$OtherAllele, joined$AssessedAllele)
sum(aligns$keep)

# replicate at eQTL level: 
# NB: 14406/18907, 76%; 12967/18907, 69%
# linear score: 12765/16654, 75%; effect consistent: 11473/16654, 69%
# for ambiguous snps, sign_flip will be TRUE, so need to correct this
joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            rep = sum(!is.na(Zscore_eqtlgen))) %>% 
  ungroup() %>% 
  mutate(prop = rep / n_ct)
  # mutate(Zscore_eqtlgen = ifelse(sign_flip == TRUE, Zscore_eqtlgen * -1, Zscore_eqtlgen)) %>% 
  # filter(Zscore_eqtlgen * slope > 0) %>% #nrow() 
  View

prop.test(c(669, 6740), c(967, 8326)) %>% broom::tidy()

prop.test(c(399, 6178), c(559, 7455)) %>% broom::tidy()

# TODO: for GTEx
snp_hg38 <- read_tsv("./sp_snp_hg38", F) %>% 
  select(chr=X1, pos_38=X3, X4) %>% 
  mutate(X4 = gsub(".*:", "", X4),
         chr = as.integer(gsub("chr", "", chr))) %>% 
  separate(X4, into=c("rm", "pos_19"), sep="-") %>% select(-rm) %>% 
  mutate(pos_19 = as.integer(pos_19))
  
# shared_mag_df_tmp
joined <- specific_eqtl_mash %>% 
  # group_by(eqtl) %>%
  # summarize(n = sum(share)) %>% ungroup() %>%
  # filter(n==1) %>%
  separate(eqtl, into=c("phenotype_id", "chr", "pos_19", "ref", "alt"), sep="_", remove=F) %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos_19 = as.integer(pos_19)) %>% 
  left_join(snp_hg38, by=c("chr", "pos_19")) %>% 
  left_join(gtex %>% select(phenotype_id, chr, pos_38, gtex_ref=ref, gtex_alt=alt, slope_gtex), 
             by=c("phenotype_id", "chr", "pos_38"))

aligns <- allele.qc(joined$ref,joined$alt,joined$gtex_ref, joined$gtex_alt)
sum(aligns$strand_flip)
sum(!aligns$keep) # keep all snps

# replicate at eQTL level: 
# 0.38
joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip) %>% 
  filter(keep == TRUE) %>% # nrow()
  mutate(slope_gtex = ifelse((strand_flip == TRUE & alt != gtex_alt) | (sign_flip == TRUE & alt != gtex_alt), 
                             slope_gtex * -1, slope_gtex)) %>% 
  mutate(n = 1) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            rep = sum(!is.na(slope_gtex))) %>% 
  ungroup() %>% 
  mutate(prop = rep / n_ct)

prop.test(c(218, 4233), c(559, 8326)) %>% broom::tidy()


#### intersect with specific annotation ####
snp_annot

shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
  distinct(phenotype_id, eqtl, celltype, shared) %>% 
  group_by(phenotype_id, eqtl) %>% 
  mutate(n_shared = sum(shared)) %>% ungroup() %>% 
  filter(n_shared == 1 & shared == TRUE) %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_") %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos = as.integer(pos),
         pos_1 = pos - 1) %>% 
  left_join(snp_annot %>% filter(snplist == "sp_snp_after_mash"), by=c("chr", "pos")) %>% 
  filter(!is.na(celltype_annot))
  
specific_eqtl %>% 
  left_join(snp_annot %>% filter(snplist == "sp_snp_before_mash"), by=c("chr", "pos")) %>% 
  filter(!is.na(celltype_annot)) %>% distinct() %>% View
  