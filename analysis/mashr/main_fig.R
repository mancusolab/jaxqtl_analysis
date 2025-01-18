### main figure

p2 <- shared_sign_df %>% gather(key = celltype, value = specific, B_IN:Plasma) %>% 
  group_by(eqtl) %>% 
  summarize(ct = sum(specific, na.rm = T)) %>%
  mutate(n_eqtl = n(), group = "shared by sign") %>% 
  bind_rows(snp_df %>% 
              filter(eqtl %in% mash_sig_eqtls) %>% # keep mash sig SNPs
              count(eqtl, name = "ct") %>% 
              mutate(n_eqtl = n(), group="simple counting")) %>% 
  bind_rows(tibble(eqtl = names(shared_mag_ct), 
                   ct = shared_mag_ct,
                   n_eqtl = length(unique(names(shared_mag_ct))),
                   group = "shared by magnitude")) %>%
  count(group, ct, n_eqtl) %>% 
  mutate(prop = n/n_eqtl) %>% 
  mutate(group = fct_relevel(as.factor(group), c("simple counting"))) %>% 
  # filter(group!= "simple counting") %>% 
  # add_row(group = "shared_by_sign", ct = 1, n = 0, prop = 0) %>% 
  # add_row(group = "specific", ct = 13, n = 0, prop = 0) %>% 
  ggplot(aes(x = ct, y = prop, fill = group)) +
  geom_col(position = "dodge") + 
  scale_fill_manual(values = c("simple counting"="grey", "shared by sign"="#E69F00", "shared by magnitude"="#4575B4"))+
  scale_x_continuous(breaks=seq(1,14, 1))+
  xlab("number of cell types") + ylab("Fraction of fine-mapped sc-eQTLs")+
  axis_theme+theme(legend.position = "bottom",
                   legend.key.size = unit("0.2", units = "cm"), # use 0.3 for main fig
                   legend.text = element_text(size=7)) # use 8 for main fig

p2

plot_grid(
  p2, p_tss_share_by_mag,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  labels = c("A", "B")
)


ggsave2(paste0(plots_dir, "Fig4_A.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# TSS distance and specificity
shared_mag_df_tmp <- shared_mag_df %>% #mutate(eqtl = gsub("_b37$", "", eqtl)) %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(dist = abs(pos - tss)) %>% 
  gather(key = celltype, value = sig, B_IN:Plasma) %>% 
  group_by(eqtl, dist) %>% 
  summarize(n_sig = sum(sig, na.rm = T)) %>% 
  ungroup() 

snp_df_tmp <- snp_df %>%
  filter(eqtl %in% mash_sig_eqtls) %>% 
  mutate(dist = abs(pos - tss)) %>% 
  count(eqtl, dist)

p_share_raw_dist <- snp_df_tmp %>% 
  mutate(mean_dist = mean_dist / 1000) %>% 
  ggplot(aes(x = n, y = mean_dist)) + geom_point() +
  geom_smooth(aes(weight = wts), method = "lm") + 
  ggtitle("shared without mash") + ylim(c(-20,110))+
  scale_x_continuous(breaks = seq(1,14))+
  axis_theme

p_dist <- shared_mag_df_tmp %>% 
  rename("shared by magnitude" = "n_sig") %>% 
  left_join(snp_df_tmp %>% select(eqtl, `simple counting`=n)) %>% 
  gather(key = method, value = n_sig, c("shared by magnitude", "simple counting")) %>%
  mutate(group = as.factor(ifelse(n_sig == 1, "specific", "shared")),
         method = fct_relevel(method, "simple counting"),
         group = fct_relevel(group, "specific")) %>% 
  mutate(dist = dist / 1000) %>%
  ggplot(aes(x = group, y = dist, fill = method)) + 
  geom_violin(width=0.7,scale = "width") + 
  geom_boxplot(width=0.05, color="black", alpha=0.2, outlier.size = 0.5,
               position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("simple counting"="grey", "shared by sign"="#E69F00", "shared by magnitude"="#4575B4"))+
  axis_theme + xlab("") + ylab("Distance to TSS (kb)")+
  ggtitle("shared by magnitude")
p_dist

# share by magnitude
share_mag <- get_pairwise_sharing(m_finemap_res, factor = 0.5, lfsr_thresh = thresh)
share_mag
share_mag <- share_mag[match(celltypes_grouped, rownames(share_mag)), 
                       match(celltypes_grouped, rownames(share_mag))]

pal_blue <- colorRampPalette(c("white", "#4575B4"), space = "rgb")(200)
legend <- get_legend(p2)
legend <- get_plot_component(p2, 'guide-box-bottom', return_all = TRUE)
p_mag <- as.grob(~corrplot(share_mag, method = c('color'), 
                           #type = "upper",
                           is.corr = F, 
                           col.lim = c(0,1),
                           col = pal_blue,
                           tl.cex = 0.6,
                           # tl.cex = 0.3,
                           #tl.offset = 1,
                           tl.srt = 45,
                           cl.ratio = 0.2, 
                           cl.align = "c",
                           cl.length = 5,
                           cl.cex = 0.8,
                           cl.offset = 1,
                           tl.col = "black",
                           mar = c(4,0,4,0)
))

matrix.colors <- getGrob(p_mag, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
p_mag <- editGrob(p_mag,
               gPath("square"), grep = TRUE,
               gp = gpar(col = NA,
                         fill = NA))

# apply the saved colours to the underlying matrix grob
p_mag <- editGrob(p_mag,
               gPath("symbols-rect-1"), grep = TRUE,
               gp = gpar(fill = matrix.colors))

# convert the background fill from white to transparent, while we are at it
p_mag <- editGrob(p_mag,
               gPath("background"), grep = TRUE,
               gp = gpar(fill = NA))

grid1 <- plot_grid(
  p2 + theme(legend.position = "none"), 
  p_dist + theme(legend.position = "none") + ggtitle(""),
  #p_mag, 
  align = 'b',
  hjust = -1, # -1
  nrow = 1,
  axis="l",
  labels = c("A", "B"), 
  label_size = 10, vjust=1,
  rel_widths = c(1, 1)
)
legend <- get_plot_component(p2 + theme(legend.text=element_text(size=8)), 
                             'guide-box-bottom', return_all = TRUE)
plot_grid(
  grid1, 
  legend, 
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  axis="l",
  rel_heights = c(1, 0.1)
)

ggsave2(paste0(plots_dir, "Fig4.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

plot_grid(p_mag)
ggsave2(paste0(plots_dir, "Fig4_pairwise_share.png"),
        width = onecol_width*1.5, height = onecol_height*1.5, units = 'cm', dpi = 300)
