## make plot for other analysis

#### eGene sharing ####

lt <- list()
for (cell in celltypes){
  lt[[cell]] <- unique(NB_egenes_df$phenotype_id[which(NB_egenes_df$celltype == cell)])
}

m <- make_comb_mat(lt)
m2 <- m[, comb_size(m) >= 30] # for nb_cs_pip0.9

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
  hjust = 0.1, # -1
  nrow = 1,
  rel_widths=c(1, 1)
)

ggsave2(paste0(plots_dir, "upset_shared_egenes.png"),
        width = onecol_width*1.5, height = onecol_height, units = 'cm', dpi = 300)

#### count shared/sp before and after mash ####
p1 <- tibble(group=c("specific", "shared_by_sign", "shared_by_magnitude"),
             n_eqtl = c(length(unique(snp_df$eqtl)),
                        length(unique(shared_sign_df$eqtl)),
                        length(unique(shared_mag_df$eqtl)))) %>% 
  mutate(group = fct_relevel(group, "specific")) %>% 
  ggplot(aes(x = group, y = n_eqtl, fill = group)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("specific"="grey", "shared_by_sign"="red3", "shared_by_magnitude"="#4575B4"))+
  ylab("number of eQTLs") + xlab("")+
  axis_theme+theme(legend.position = "bottom",
                   legend.key.size = unit("0.3", units = "cm"),
                   legend.text = element_text(size=6),
                   axis.text.x = element_blank())
p1
p2 <- shared_sign_df %>% gather(key = celltype, value = specific, B_IN:Plasma) %>% 
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

legend <- get_plot_component(p1, 'guide-box-bottom')
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

#### TSS dist analysis for shared and sp ####

mash_sp_dist <- specific_eqtl_mash %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(dist = abs(pos - tss)) %>% pull(dist)

raw_sp_dist <- specific_eqtl_raw$dist

median(mash_sp_dist); median(raw_sp_dist)
wilcox.test(mash_sp_dist, raw_sp_dist)

## plot mean TSS
snp_df_tmp <- snp_df %>%
  filter(eqtl %in% mash_sig_eqtls) %>% 
  mutate(dist = abs(pos - tss)) %>% 
  count(eqtl, dist) %>% 
  group_by(n) %>% 
  summarize(mean_dist = mean(dist),
            wts = 1/var(dist),
            wts = ifelse(!is.finite(wts), 1, wts),
            n_ct = n()) %>% 
  ungroup()

p_share_raw_dist <- snp_df_tmp %>% 
  mutate(mean_dist = mean_dist / 1000) %>% 
  ggplot(aes(x = n, y = mean_dist)) + geom_point() +
  geom_smooth(aes(weight = wts), method = "lm") + 
  ggtitle("shared without mash") + ylim(c(-20,110))+
  scale_x_continuous(breaks = seq(1,14))+
  axis_theme
p_share_raw_dist

glm(mean_dist ~ n, data = snp_df_tmp, weights = wts, family = "gaussian") %>% broom::tidy()

shared_mag_df_tmp <- shared_mag_df %>% #mutate(eqtl = gsub("_b37$", "", eqtl)) %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(dist = abs(pos - tss)) %>% 
  gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl, dist) %>% 
  summarize(n_sig = sum(sig, na.rm = T)) %>% 
  ungroup() 

p_share_mag_dist <- shared_mag_df_tmp %>% 
  group_by(n_sig) %>% summarize(mean = mean(dist),
                                wts = 1/(var(dist))) %>% 
  ungroup() %>% 
  mutate(mean = mean / 1000) %>%
  ggplot(aes(x = n_sig, y = mean)) + geom_point() + 
  geom_smooth(method="lm", mapping = aes(weight = wts)) +
  ylim(c(-20,110))+
  scale_x_continuous(breaks = seq(1,14))+
  axis_theme + xlab("number of cell types") + ylab("mean distance to TSS (kb)")+
  ggtitle("shared by magnitude")
p_share_mag_dist

glm(mean ~ n_sig, data = shared_mag_df_tmp %>% 
  group_by(n_sig) %>% summarize(mean = mean(dist),
                                wts = 1/(var(dist))) %>% 
  ungroup(), weights = wts, family = "gaussian") %>% broom::tidy()

shared_sign_df_tmp <- shared_sign_df %>% #mutate(eqtl = gsub("_b37$", "", eqtl)) %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(dist = abs(pos - tss)) %>% 
  gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl, dist) %>% 
  summarize(n_sig = sum(sig, na.rm = T)) %>% 
  ungroup() 

p_share_sign_dist <- shared_sign_df_tmp %>% 
  group_by(n_sig) %>% summarize(mean = mean(dist),
                                wts = 1/(var(dist)),
                                sd = sd(dist)) %>% 
  ungroup() %>% 
  mutate(mean = mean / 1000) %>%
  ggplot(aes(x = n_sig, y = mean)) + geom_point() +
  # geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd)) +
  geom_smooth(method="lm", mapping = aes(weight = wts)) +
  scale_x_continuous(breaks = seq(1,14))+
  ylim(c(-20,110))+
  axis_theme + xlab("number of cell types") + ylab("mean distance to TSS (kb)")+
  ggtitle("shared by sign")
p_share_sign_dist

glm(mean ~ n_sig, data = shared_sign_df_tmp %>% 
      group_by(n_sig) %>% summarize(mean = mean(dist),
                                    wts = 1/(var(dist))) %>% 
      ungroup(), weights = wts, family = "gaussian") %>% broom::tidy()

plot_grid(
  p_share_raw_dist, p_share_sign_dist, p_share_mag_dist,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  labels = c("A", "B", "C")
)

ggsave2(paste0(plots_dir, "shared_tss_mean.png"),
        width = full_width*1.2, height = full_height, units = 'cm', dpi = 300)

# compare distance
wilcox.test(specific_eqtl_mash %>% distinct(eqtl) %>% left_join(snp_df %>% distinct(eqtl, dist)) %>% pull(dist),
            shared_eqtl_mash %>% distinct(eqtl) %>% left_join(snp_df %>% distinct(eqtl, dist)) %>% pull(dist))

t.test(specific_eqtl_raw %>% distinct(eqtl, dist) %>% pull(dist),
            shared_eqtl_raw %>% distinct(eqtl, dist) %>% pull(dist))


willcox.test(shared_eqtl_mash %>% distinct(eqtl) %>% left_join(snp_df %>% distinct(eqtl, dist)) %>% pull(dist),
            shared_eqtl_raw %>% distinct(eqtl, dist) %>% pull(dist))

## boxplot of specific vs. shared

p1 <- snp_df_sig %>% count(eqtl, dist) %>% 
  mutate(group = ifelse(n == 1, "specific", "shared"),
         group = fct_relevel(group, "specific"),
         dist = dist / 1000) %>% 
  ggplot(aes(x = group, y = dist)) + geom_boxplot() + axis_theme + ylab("TSS dist (kb)")+
  ggtitle("simple counting")
p1

p2 <- shared_mag_df_tmp %>% 
  mutate(group = ifelse(n_sig == 1, "specific", "shared"),
         group = fct_relevel(group, "specific"),
         dist = dist / 1000) %>% 
  ggplot(aes(x = group, y = dist)) + geom_boxplot() + axis_theme + ylab("TSS dist (kb)")+
  ggtitle("mash share by magnitude")

p_tss_share_by_mag <- p2
plot_grid(
  p1, p2,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  labels = c("A", "B")
)

ggsave2(paste0(plots_dir, "shared_tss_mean_sp_shared.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

## plot TSS bin proportion 
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

p_dist_percent_mash <- shared_mag_df_tmp %>% 
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
  axis_theme + theme(legend.position = "bottom")

legend <- get_plot_component(p_dist_percent_mash, "guide-box-bottom")
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


#### TF analysis: gene-TF pairs ####

# download from hTFtarget: http://bioinfo.life.hust.edu.cn/hTFtarget/#!/
tf_df_raw <- read_tsv("../data/OneK1K/annotation/TF/TF-Target-information.txt")
table(tf_df_raw$tissue)

# extract TF-gene pairs observed in blood tissue/cells
# blood: "blood", "Blood", "peripheral blood", 
# cell type: "lymphocyte", "cell:primary human memory B cells"
tf_list <- c()
for (tf in unique(tf_df_raw$tissue)){
  new <- unique(str_split_1(tf, ","))
  tf_list <- c(tf_list, new[!new %in% tf_list])
}

tf_df_blood <- tf_df_raw %>% 
    filter(grepl("blood|Blood|peripheral blood|lymphocyte|cell:primary human memory B cells", tissue)) %>% 
    distinct(TF, target)
tf_df_t_b <- tf_df_raw %>% filter(grepl("lymphocyte|cell:primary human memory B cells", tissue)) %>% distinct(TF, target); dim(tf_df_t_b)
tf_df_none_t_b <- tf_df_blood %>% anti_join(tf_df_t_b); dim(tf_df_none_t_b)

tf_df_no_blood <- tf_df_raw %>% distinct(TF, target) %>% anti_join(tf_df_blood)
tf_df <- tf_df_blood
# extract genes with specific eQTL or shared eQTL and remove those genes with both
for (dat in c("raw", "mash")){
  if (dat == "raw"){
    specific_eqtl <- specific_eqtl_raw %>% filter(eqtl %in% mash_sig_eqtls)
    shared_eqtl <- shared_eqtl_raw %>% filter(eqtl %in% mash_sig_eqtls)
  }else{
    specific_eqtl <- specific_eqtl_mash
    shared_eqtl <- shared_eqtl_mash
  }
  
  sp_genelist <- unique(specific_eqtl$GeneSymbol)
  shared_genelist <- unique(shared_eqtl$GeneSymbol)
  
  to_rm <- intersect(sp_genelist, shared_genelist);length(to_rm)
  sp_genelist <- sp_genelist[!sp_genelist %in% to_rm];length(sp_genelist)
  shared_genelist <- shared_genelist[!shared_genelist %in% to_rm];length(shared_genelist)
  
  tf_df_tmp <- tf_df %>%
    count(target) %>%
    inner_join(bind_rows(tibble(target = sp_genelist, group = "specific"),
                         tibble(target = shared_genelist, group = "shared")))
  
  tf_df_tmp %>% group_by(group) %>% summarize(mean = mean(n),
                                              med = median(n)) %>% print()
  print(wilcox.test(tf_df_tmp %>% filter(group == "specific") %>% pull(n),
              tf_df_tmp %>% filter(group == "shared") %>% pull(n)) %>% 
    broom::tidy())
  
  title <- ifelse(dat == "raw", "simple counting", "mash estimates")
  p <- tf_df_tmp %>% 
    ggplot(aes(x = group, y = n)) + geom_boxplot() +
    ylab("number of TF") + xlab("group of genes")+
    ggtitle(title) +
    axis_theme + xlab("")
  
  assign(paste0("p_tf_", dat), p)
}

grid3 <- plot_grid(p_tf_raw, p_tf_mash, 
          nrow = 1)

plot_grid(
  grid1, grid2, grid3,
  align = 'h',
  hjust = 0, # -1
  vjust=1,
  nrow = 3,
  axis="b",
  labels = c("A: T/B blood cells", "B: no T/B blood cells", "C: background"), 
  label_size = 8
)

ggsave2(paste0(plots_dir, "egene_mash_num_tf.png"),
        width = full_width, height = full_height*2.2, units = 'cm', dpi = 300)

# TF is generally more expressed than other genes due to stoschasticity + searching time
# needed to bind TFBS and initiate transcription of important client genes (constantly expressed)
allgenes_express %>% left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  mutate(blood_tf = ifelse(GeneSymbol %in% (tf_df %>% distinct(TF) %>% pull(TF)), "yes", "no")) %>% 
  # filter(rate_mean < 0.001) %>% 
  ggplot(aes(x = blood_tf, y = express_percent)) + geom_boxplot() + facet_wrap(~celltype)

## ENCODE-TF
tfbs <- read_tsv("../result/mvsusie/result/all_intersect_tf_motif", F) %>% 
  select(chr=X1, pos=X3, eqtl=X4, TF=X10, start=X8, end=X9) %>%
  separate(TF, into=c("TF", "suffix"), sep='_', remove = F) %>% 
  distinct(chr, pos, eqtl, TF) %>% 
  filter(eqtl %in% snp_df$eqtl)

specific_eqtl <- specific_eqtl_raw
shared_eqtl <- shared_eqtl_raw
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
p_tf_raw <- tfbs_tmp %>% #View
  ggplot(aes(x = group, y = n)) + geom_boxplot() + axis_theme

plot_grid(p_tf_raw, p_tf_mash, 
          nrow = 1,
          labels = c("A", "B"))

ggsave2(paste0(plots_dir, "egene_mash_num_tfbs_encode.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

#### ATACseq overlap ####
atacseq_map

atac <- read_tsv("../result/mvsusie/result/all_intersect_atac",F)
colnames(atac) <- c("chr", "pos_1", "pos", "chr_atac", "start", "end", "atac_celltype")

specific_eqtl_mash <- specific_eqtl_mash %>% 
  left_join(snp_df %>% distinct(eqtl, chr, pos), by="eqtl")

# eQTL - cell type
snp_df_sig <- snp_df %>% filter(eqtl %in% mash_sig_eqtls)

atac_enrich <- function(yazar_ct, atac_celltype, snp_df_sig, sp_dat) {
  # specific eQTL
  n_spec <- sp_dat %>% distinct(eqtl) %>% nrow(); n_spec
  
  # eQTL within a cell type
  n_sig_snp_ct <- snp_df_sig %>%
    filter(celltype == yazar_ct) %>% 
    inner_join(atac %>% filter(atac_celltype == atac_ct), by=c("chr", "pos")) %>% nrow()
  
  # eQTL specific to a cell type
  n_spec_ct <- sp_dat %>% 
    filter(celltype == yazar_ct) %>% 
    inner_join(atac %>% filter(atac_celltype == atac_ct), by=c("chr", "pos")) %>% nrow()
  
  enrich <- n_spec_ct/n_sig_snp_ct/(n_spec/n_sig_snp)
  return(enrich)
}

n_sig_snp <- snp_df_sig %>% distinct(eqtl) %>% nrow(); n_sig_snp

table(snp_df_sig$celltype)
table(specific_eqtl_mash$celltype)
table(specific_eqtl_raw$celltype)

atac_intersect <- expand.grid(celltype = unique(snp_df$celltype), method = c("mash", "raw"))
atac_intersect <- atac_intersect %>% mutate(enrich = NA, enrich_se = NA)
for (i in 1:nrow(atac_intersect)){
  yazar_ct <- atac_intersect$celltype[i]
  atac_ct <-atacseq_map$annot_celltype[which(atacseq_map$yazar_celltype == yazar_ct)]; atac_ct
  method <- atac_intersect$method[i]
  
  if (method == "raw"){
    sp_dat <- specific_eqtl_raw
  }else{
    sp_dat <- specific_eqtl_mash
  }
  
  atac_intersect$enrich[i] <- atac_enrich(yazar_ct, atac_celltype, snp_df_sig, sp_dat)
  
  # bootstrap
  n_bt <- 1000
  set.seed(2024)
  enrich_vec <- c()
  for (j in 1:n_bt){
    snp_df_sig_perm <- snp_df_sig %>% mutate(celltype = sample(celltype, replace=T))
    sp_dat_perm <- snp_df_sig_perm %>% filter(eqtl %in% sp_dat$eqtl) %>% 
      distinct(eqtl, .keep_all = T)
    enrich_vec <- c(enrich_vec, atac_enrich(yazar_ct, atac_celltype, snp_df_sig_perm, sp_dat_perm))
  }
  
  atac_intersect$enrich_se[i] <- sd(enrich_vec, na.rm = T)
  
  print(i)
}

atac_intersect <- atac_intersect %>% mutate(enrich_lb = enrich - 1.96 * enrich_se,
                          enrich_ub = enrich + 1.96 * enrich_se)

forest_plot_theme(atac_intersect %>% 
                    mutate(enrich_lb = ifelse(enrich == 0 | enrich_lb < 0, 0, enrich_lb),
                           enrich_ub = ifelse(enrich == 0, 0, enrich_ub)) %>% 
                    mutate(celltype = fct_relevel(celltype, rev(celltype_order)),
                           method = ifelse(method == "raw", "simple counting", "mash")) %>% 
                    ggplot(aes(y = celltype, x = enrich, color = method)) +
                    geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = enrich_lb, xmax = enrich_ub),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    scale_color_manual(values = c("simple counting"="grey", "mash"="#4575B4"))+
                    xlab("Enrichment (95% CI)") + 
                    ylab(" ") + 
                    theme_bw()) +
  theme(legend.key.size = unit(0.3, units = "cm"))

ggsave2(paste0(plots_dir, "mash_raw_enrich_atac.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)


#### SNP2TF analysis ####

# show specific eqtl are enriched in TFBS
# match in database: 1786 vs. 136
mash_tfbs <- read_tsv("../result/mvsusie/result/SNP2TF/raw_sp_TF_enrich.txt") %>% 
  janitor::clean_names()
mash_match <- read_tsv("../result/mvsusie/result/SNP2TF/mash_sp_annot.txt") %>% 
  janitor::clean_names()
names(mash_tfbs)

mash_tfbs %>% filter(number_p_value < 0.05/1786)

#### RegulomeDB ####

db <- read_tsv(paste0(annot_dir, "RegulomeDB/ENCFF250UJY.tsv"))

joined_db_raw <- specific_eqtl_raw %>% 
  separate(eqtl, into=c("gene", "snp"), sep = "_", extra = "merge") %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp" = "variant_id")) %>% 
  inner_join(db, by=c("rsid"))

summary(joined_db_mash$probability_score)
summary(joined_db_raw$probability_score)

#### selection parameter ####
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
p_tf_EDS <- s_tmp %>% #View
  ggplot(aes(x = group, y = EDS)) + geom_boxplot() + axis_theme

plot_grid(p_tf_pLI, p_tf_EDS, p_tf_loeuf, p_tf_RVIS,
          nrow = 1,
          labels = c("pLI", "EDS", "loeuf", "RVIS"),
          label_size = 8,
          hjust = 0)

ggsave2(paste0(plots_dir, "egene_raw_selection.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

#### POPscore ####

popscore <- read_tsv("../data/OneK1K/annotation/popscore/PoPS_FullResults.txt.gz") %>% 
  filter(cohort == "PASS")
hist(popscore$pops_score)

table(popscore$trait)
bind_rows(specific_eqtl %>% mutate(group = "specific"),
          shared_eqtl %>% mutate(group = "shared")) %>% 
  filter(!GeneSymbol %in% to_rm) %>% 
  inner_join(popscore %>% 
               group_by(ensgid) %>% 
               summarize(mean_score = mean(pops_score)), by=c("phenotype_id" = "ensgid")) %>% 
  ggplot(aes(x = group, y = mean_score)) + geom_boxplot()
  


