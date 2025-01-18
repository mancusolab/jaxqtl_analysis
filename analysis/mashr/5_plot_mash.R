# finemap

library(mvsusieR)
library(corrplot)
library(ggplotify)
load("../result/mvsusie/result/mashr_Uk_fit_EZ_TRUE_complete_TRUE.RData")

#### mixture weights ####
# m is fitted using weight on random set
# using cov estimated by top provide more weights to null, 
# which is as expected, since top eQTL not neccesarity catpure "true" eQTL effect
mixture_wts <- get_estimated_pi(m_finemap_highpip) # m_finemap_cslead
length(mixture_wts)
sort(mixture_wts)

U_k_wt_avg <- matrix(0, nrow=14, ncol=14)
cov_list <- c(U.ed.cslead, U.c)
for (i in names(cov_list)){
  U_k_wt_avg <- U_k_wt_avg + mixture_wts[i] * cov_list[[i]]
}

p_wts <- tibble(cov = names(mixture_wts),
                wts = mixture_wts) %>% 
  mutate(cov = fct_reorder(cov, desc(wts))) %>% 
  ggplot(aes(x = cov, y = wts)) +
  geom_col() +
  axis_theme +
  ylab("mixture weight") + xlab("covariance structure")+
  theme(axis.text.x = element_text(angle=35, hjust = 1))
p_wts

ggsave2(paste0(plots_dir, "Fig_mvsusie_weights.png"),
        width = onecol_width*1.7, height = onecol_height, units = 'cm', dpi = 300)


### largest component from data driven approach ###
U.ed$ED_PCA_1

# ED_empirical: 19%,  ED_tPCA: 8%, 
# cor plot looks similar
U_cor <- cov2cor(U_k_wt_avg)

U_cor <- U_cor[match(celltypes_grouped, rownames(U_cor)), match(celltypes_grouped, rownames(U_cor))]
mean(U_cor > 0) # all positive values

p_U_em <- as.grob(~corrplot(U_cor, method = c('color'), 
         #type = "upper",
         is.corr = T, 
         #col.lim = c(0,1),
         col = pal_bry,
         tl.cex = 0.7,
         mar = c(4,0,4,0),
         tl.offset = 0.5,
         tl.srt = 45,
         tl.col = "black"
))


plot_grid(
  p_U_em,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="r",
  labels = c(""), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "Fig_mvsusie_U_k_wt_avg.png"),
        width = onecol_width*2, height = onecol_height*2.2, units = 'cm', dpi = 300)

# plot Vest matrix
U_cor <- cov2cor(Vest)

U_cor <- U_cor[match(celltypes_grouped, rownames(U_cor)), match(celltypes_grouped, rownames(U_cor))]
mean(U_cor > 0) # all positive values

p_U_em <- as.grob(~corrplot(U_cor, method = c('color'), 
                            #type = "upper",
                            is.corr = T, 
                            #col.lim = c(0,1),
                            col = pal_br,
                            tl.cex = 0.7,
                            mar = c(4,0,4,0),
                            tl.offset = 0.5,
                            tl.srt = 45,
))


plot_grid(
  p_U_em,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="r",
  labels = c(""), label_size = 10, vjust=1
)
ggsave2(paste0(plots_dir, "Fig_mvsusie_resid_Vest.png"),
        width = onecol_width*2, height = onecol_height*2.2, units = 'cm', dpi = 300)


sum(U_norm)
U.ed$ED_empirical
p_pc4 <- as.grob(~corrplot(U_cor, method = c('color'), 
                           type = "upper",
                           is.corr = T, 
                           col = pal(200)
))
p_pc1 <- as.grob(~corrplot(cov2cor(U.ed$ED_tPCA), method = c('color'), 
                           type = "upper",
                           is.corr = T, 
                           col = pal(200)
))
p_pc2 <- as.grob(~corrplot(cov2cor(U.ed$ED_PCA_2), method = c('color'), 
                           type = "upper",
                           is.corr = T, 
                           col = pal(200)
))

grid_pc <- plot_grid(
  p_pc4, p_pc1, p_pc2,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="r",
  labels = c("A", "B", "C"), label_size = 10, vjust=1
)

plot_grid(
  p_wts, grid_pc,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  axis="r",
  labels = c("A"), label_size = 10, vjust=1
)

p_wts
ggsave2(paste0(plots_dir, "Fig_mvsusie_weights.png"),
        width = onecol_width*1.5, height = onecol_height, units = 'cm', dpi = 300)


grid_pc
ggsave2(paste0(plots_dir, "Fig_mvsusie_cov.png"),
        width = full_width*2.5, height = full_height*2, units = 'cm', dpi = 300)


#### effect sharing between cell types ####

# get lfsr
m_finemap_lfsr <- get_lfsr(m_finemap); dim(m_finemap_lfsr) # 9966 snps x 14 cell types
m_finemap_lfsr <- tibble::rownames_to_column(as.data.frame(m_finemap_lfsr), "snp")

#### significant snps ####
# take snps with lfsr < 0.05 in at least celltype

thresh <- 0.05 # set threshold for lfsr

sig_eqtls <- get_significant_results(m_finemap, thresh = thresh)
length(sig_eqtls) # 9609
length(sig_eqtls)/nrow(m_finemap_lfsr)

##### share by sign and/or magnitude #####

# see code: https://github.com/stephenslab/gtexresults/blob/master/analysis/SharingHist.Rmd
source("./mvsusie/normfuncs.R")

# prepare post estimates
lfsr <- get_lfsr(m_finemap); dim(lfsr) # 9966 x 14 (no negative values)
pm.mash.beta <- get_pm(m_finemap); dim(pm.mash.beta) # 9966 x 14

sigmat  <- (lfsr<=thresh)
nsig    <-  rowSums(sigmat)
mean(nsig > 0) # 96% are significant in at least one cell type

# overall sharing by sign: share the same sign with strongest effect: 92%
# het.norm computes the ratio of x/rowMax(x)
signall <- mean(het.norm(pm.mash.beta[nsig>0,])>0)
signall

# overall sharing by magnitude: 29.4%
magall <- mean(het.norm(pm.mash.beta[nsig>0,])>0.5)
magall

# share by magnitude (compared against the cell type with largest effect)
# get proportion of pairwise sharing between traits
share_0.5_sign <- get_pairwise_sharing(m_finemap, factor = 0.5, lfsr_thresh = thresh)
share_0.5_sign

which.min(share_0.5_sign) # NK and CD4_NC
rowMeans(share_0.5_sign)

Matrix::triu(share_0.5_sign)

sigmat <- (lfsr <= thresh)
nsig   <- rowSums(sigmat)
ct <- het.func(het.norm(effectsize=pm.mash.beta[nsig>0,]),threshold=0.5)

p_share_mag <- tibble::rownames_to_column(as.data.frame(ct), "eQTL") %>% 
  count(ct) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = ct, y = prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1, 14, 1))+
  ylim(0,1)+
  ylab("Proportion of lead fine-mapped eQTLs")+
  ggtitle("sharing by magnitude") +
  axis_theme

sign.func <- function(normeffectsize){
  apply(normeffectsize,1,function(x)(sum(x>0)))
}

sigmat <- (lfsr<=thresh)
nsig   <- rowSums(sigmat)
ct <- sign.func(het.norm(effectsize=pm.mash.beta[nsig>0,]))

p_share_sign <- tibble::rownames_to_column(as.data.frame(ct), "eQTL") %>% 
  count(ct) %>% 
  mutate(prop = n / sum(n)) %>% 
  ggplot(aes(x = ct, y = prop)) + geom_col() +
  scale_x_continuous(breaks = seq(1, 14, 1))+
  ylim(0, 1)+
  ylab("Proportion of lead fine-mapped eQTLs") +
  ggtitle("sharing by sign") +
  axis_theme

plot_grid(
  p_share_sign, p_share_mag,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="r",
  labels = c("A", "B"), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "mash_eqtl_share_overall.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

##### pairwise sharing #####

# share by magnitude
share_mag <- get_pairwise_sharing(m_finemap_cslead, factor = 0.5, lfsr_thresh = thresh)
share_mag
share_mag <- share_mag[match(celltypes_grouped, rownames(share_mag)), 
                       match(celltypes_grouped, rownames(share_mag))]

p_mag <- as.grob(~corrplot(share_mag, method = c('color'), 
         #type = "upper",
         is.corr = T, 
         col.lim = c(0,1),
         col = pal_bry,
         tl.cex = 0.8,
         #mar = c(4,0,4,0),
         #tl.offset = 1,
         tl.srt = 45,
))

share_sign <- get_pairwise_sharing(m_finemap_cslead, factor = 0, lfsr_thresh = thresh)
share_sign

share_sign <- share_sign[match(celltypes_grouped, rownames(share_sign)), 
                       match(celltypes_grouped, rownames(share_sign))]


p_sign <- as.grob(~corrplot(share_sign, method = c('color'), 
                           #type = "upper",
                           is.corr = T, 
                           col.lim = c(0,1),
                           col = pal,
                           tl.cex = 0.6,
                           mar = c(4,0,4,0),
                           tl.offset = 1,
                           tl.srt = 45,
))

plot_grid(
  p_mag,
  align = 'h',
  hjust = -1, # -1
  nrow = 1
)

ggsave2(paste0(plots_dir, "mash_eqtl_share_pairwise.png"),
        width = onecol_width*1.8, height = onecol_height*2.2, units = 'cm', dpi = 300)

#### plot shared by magnitude ####

lt <- list()
for (cell in celltypes){
  lt[[cell]] <- shared_mag_df$eqtl[shared_mag_df[[cell]] == TRUE]
}

m <- make_comb_mat(lt)
m2 <- m[, comb_size(m) >= 10] # for nb_cs_pip0.9

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


