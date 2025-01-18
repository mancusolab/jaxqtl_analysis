# look at gene expression across 14 cell types


# plot shared eQTL between cell types; 76%
snp_df %>% 
  count(eqtl, name = "n_eqtl") %>% 
  count(n_eqtl) %>% 
  mutate(prop = n / sum(n))

##### gene expression by threshold ####

# use different expression threshold
genes_10per <- allgenes_express %>% 
  filter(express_percent >= 0.1) %>% 
  count(phenotype_id) %>% 
  filter(n==14) %>% pull(phenotype_id)

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

# identified in one celltype: 83% (1% cutoff), 81% (10% cutoff)
as.data.frame(eqtl_cts)

p_eqtl_share <- eqtl_cts %>% rename(threshold=thresh) %>% 
  ggplot(aes(x = n, y = prop, fill = threshold)) +
  geom_bar(position = "dodge", stat="identity") + 
  scale_x_continuous(breaks=seq(1,14, 1))+
  xlab("number of cell types") + ylab("fraction of eQTLs")+
  axis_theme + theme(legend.position = "bottom",
                     legend.key.size = unit("0.3", units = "cm"),
                     legend.text = element_text(size = 7)) 

p_eqtl_share

n_specific <- length(unique(specific_eqtl_raw$phenotype_id)); n_specific
p_sp_express <- allgenes_express %>% 
  filter(phenotype_id %in% specific_eqtl_raw$phenotype_id) %>% 
  group_by(phenotype_id) %>% 
  summarize(coverage_1per = sum(express_percent >= 0.01),
            coverage_10per = sum(express_percent >= 0.1)) %>% 
  gather(key = thresh, value = ct, -phenotype_id) %>% 
  count(thresh, ct) %>% 
  mutate(prop = n/n_specific) %>% rename(threshold = thresh) %>% 
  ggplot(aes(x = ct, y = prop, fill=threshold)) + 
  geom_col(position = "dodge") + axis_theme +
  scale_x_continuous(breaks=seq(0,14, 1))+
  xlab("expressed cell types") + ylab("proportion of eGenes")+
  axis_theme + theme(legend.position = "bottom",
                     legend.key.size = unit("0.3", units = "cm"),
                     legend.text = element_text(size = 7)) 
p_sp_express

legend <- get_plot_component(p_sp_express, 'guide-box-bottom')
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


#### correlation of specific gene across cell types ####
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

#### gene expression after mash ####

summary(as.vector(SE[rownames(SE) %in% specific_eqtl_mash$eqtl,]))
summary(as.vector(SE[rownames(SE) %in% shared_eqtl_mash$eqtl,]))
wilcox.test(as.vector(SE[rownames(SE) %in% specific_eqtl_mash$eqtl,]),
as.vector(SE[rownames(SE) %in% shared_eqtl_mash$eqtl,]))

allgenes_express %>% 
  filter(phenotype_id %in% specific_eqtl_mash$phenotype_id) %>% 
  mutate(group = "specific") %>% 
  bind_rows(allgenes_express %>% 
              filter(phenotype_id %in% shared_eqtl_mash$phenotype_id) %>% 
              filter(!phenotype_id %in% specific_eqtl_mash$phenotype_id) %>% 
              mutate(group = "shared")) %>% 
  ggplot(aes(x = log(rate_mean), group = group, linetype = group)) +
  geom_density(aes(color=celltype), alpha=0.8) +
  facet_wrap(~celltype)+
  scale_color_manual(values = celltype_cols) + axis_theme
specific_eqtl_mash

ggsave2(paste0(plots_dir, "specific_shared_gene_express_density.png"),
        width = full_width*2, height = full_height*2, units = 'cm', dpi = 300)

