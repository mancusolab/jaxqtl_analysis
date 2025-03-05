# Figure: finemapping and enrichment analysis

# average within eGene
tmp <- bind_rows(nb_cs_pip0.5 %>% inner_join(both_hits %>% select(phenotype_id, celltype)) %>% mutate(model="NegBinom"),
                 lm_cs_pip0.5 %>% inner_join(both_hits %>% select(phenotype_id, celltype)) %>% mutate(model="Linear")) %>% 
  left_join(gtf %>% select(phenotype_id, chr, tss=end), by=c("chr", "phenotype_id")) %>% 
  mutate(tss_dist = pos - tss)

tmp %>% distinct(phenotype_id, snp, tss_dist, model) %>% group_by(model) %>% 
  summarize(mean = mean(tss_dist),
            median = median(tss_dist),
            min = min(tss_dist), max = max(tss_dist))

wilcox.test(abs(tss_dist) ~ model, data=tmp %>% 
              distinct(phenotype_id, snp, tss_dist, model))

# TSS
tss_breaks <- c(0, 5000, 100000, 5000000) # 

wilcox.test(tss_dist ~ model,
            data=tmp %>% 
              mutate(tss_dist = abs(tss_dist),
                     tss_dist = ifelse(tss_dist < 1, tss_dist+1, tss_dist), # put tss_dist=0 to first bin
                     bins = cut(tss_dist, 
                                breaks = tss_breaks)) %>% 
              distinct(model, phenotype_id, snp, tss_dist, bins) # %>% 
            # group_by(sc_bulk) %>% summarize(mean = mean(tss_dist), median = median(tss_dist))
)


# baseline enrichment
## NB sc vs. bulk
p_enrich_baseline <- forest_plot_theme(baseline_meta_celltype %>% 
                                         # filter(meta_OR < 500) %>% 
                                         rename(model=method) %>% 
                                         mutate(annot = gsub("_phastCons46way$", "", annot)) %>% 
                                         group_by(annot) %>% 
                                         mutate(jaxqtl_or = meta_OR[which(model == "NegBinom")]) %>% ungroup() %>% 
                                         mutate(annot = fct_reorder(annot, jaxqtl_or)) %>% 
                                         ggplot(aes(y = annot, x = meta_OR, color = model)) +
                                         geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                                         geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                        height = 0.25, 
                                                        position = position_dodge(width=0.5)) +
                                         geom_vline(xintercept = 1, color = "red", 
                                                    linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                         scale_color_manual(values = method_cols) +
                                         xlab("Odds Ratio (95% CI)") + 
                                         ylab(" ") + 
                                         theme_bw()+
                                         ggtitle("Baseline annotations"))
p_enrich_baseline


# check annot order are the same
atac_dat <- atac_meta_celltype %>% 
  mutate(celltype = fct_relevel(celltype, rev(celltype_order))) %>% 
  rename(model = method)

p_enrich_atac <- forest_plot_theme(atac_dat %>% 
                                     # filter(meta_OR<100) %>% 
                                     ggplot(aes(y = celltype, x = meta_OR, color = model)) +
                                     geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                                     geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                    height = 0.25, 
                                                    position = position_dodge(width=0.5)) +
                                     geom_vline(xintercept = 1, color = "red", 
                                                linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                     scale_color_manual(values = method_cols) +
                                     xlab("Odds Ratio (95% CI)") + 
                                     ylab(" ") + 
                                     theme_bw() +
                                     ggtitle("ATACseq"))
p_enrich_atac
## epimap
epimap_dat <- epimap_meta_celltype %>% 
  rename(model = method) %>% 
  mutate(category = ifelse(grepl("enhancer", annot), "EpiMap enhancer", "EpiMap promoter")) %>% 
  select(colnames(atac_dat), everything())

p_enrich_epimap <- forest_plot_theme(bind_rows(atac_dat %>% mutate(category = "ATACseq"),
                                               epimap_dat) %>% 
                                       filter(meta_OR < 100) %>% 
                                       mutate(celltype = fct_relevel(celltype, rev(celltype_order)),
                                              category = fct_relevel(category, c("ATACseq", "EpiMap promoter", "EpiMap enhancer"))) %>% 
                                       ggplot(aes(y = celltype, x = meta_OR, color = model)) +
                                       geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                                       geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                      height = 0.25, 
                                                      position = position_dodge(width=0.5)) +
                                       facet_grid(.~category)+
                                       geom_vline(xintercept = 1, color = "red", 
                                                  linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                       scale_color_manual(values = method_cols) +
                                       xlab("Odds Ratio (95% CI)") + 
                                       ylab(" ") + 
                                       theme_bw())

p_enrich_epimap

# SCENT result
scent_map
nb_allcell_df_tmp <- nb_enrich %>% select(celltype) %>% 
  left_join(scent_map, by=c("celltype"="yazar_celltype")) %>% 
  left_join(nb_allcell_scent %>% select(-celltype), by=c("annot_celltype"="ref_cell")) %>% 
  select(-annot_celltype)

p_enrich_scent <- forest_plot_theme(bind_rows(nb_enrich %>% mutate(model = "NegBinom"),
                                              lm_score_enrich %>% mutate(model = "Linear")) %>% 
                                      rename(enrich = mean_enrich) %>% 
                                      mutate(celltype = fct_relevel(celltype, rev(celltype_order))) %>% 
                                      ggplot(aes(y = celltype, x = enrich, color = model)) + 
                                      geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
                                      geom_errorbarh(aes(xmin = OR_L, xmax = OR_U),
                                                     height = 0.25,
                                                     position = position_dodge(0.5)) +
                                      geom_vline(xintercept = 1, color = "red", 
                                                 linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                      scale_color_manual(values = method_cols)+
                                      axis_theme + ylab("")+
                                      xlab("Odds Ratio (95% CI)") + 
                                      ggtitle("SCENT enhancer-gene"))
p_enrich_scent

# add meta result across 14 cell type vs. bulk result
atac_meta
epimap_meta
scent_meta

p_epimap_atac_meta <- forest_plot_theme(bind_rows(epimap_meta %>%
                                                    rename(model = method) %>% 
                                                    mutate(category = paste0("EpiMap ", category),
                                                           model = ifelse(model == "tensorqtl", "Linear",  "NegBinom")),
                                                  atac_meta %>% mutate(category = "ATACseq",
                                                                       model = ifelse(model == "tensorqtl", "Linear", model))) %>% 
                                                    select(model, category, everything()) %>% 
                                          mutate(category = fct_relevel(category, c("ATACseq", "EpiMap promoter", "EpiMap enhancer"))) %>% 
                                          ggplot(aes(y = model, x = meta_OR, color = model)) + 
                                          geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
                                          geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                         height = 0.25,
                                                         position = position_dodge(0.5)) +
                                          geom_vline(xintercept = 1, color = "red", 
                                                     linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                          scale_color_manual(values = method_cols)+
                                          facet_grid(.~category)+
                                          axis_theme + ylab("")+
                                          xlab("Odds Ratio (95% CI)")) + theme(strip.text = element_blank())

p_epimap_atac_meta

p_scent_meta <- scent_meta %>% 
  ggplot(aes(y = model, x = meta_OR, color = model)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = method_cols)+
  axis_theme + ylab("")+
  theme(legend.position = "bottom")+
  xlab("Odds Ratio (95% CI)")
p_scent_meta

grid1 <- plot_grid(
  p_enrich_epimap+theme(legend.position = "none",
                        plot.margin = unit(c(0,0,0,0), "cm")) + xlab("")+
    coord_fixed(ratio = 0.2, xlim = c(1, 3.5)),
  p_epimap_atac_meta+theme(legend.position = "none",
                           plot.margin = unit(c(0,0,0,0), "cm"))+
    xlab("Enrichment (95% CI)") +
    coord_fixed(ratio = 0.225, xlim = c(1, 3.5)),
  # xlim(1,4.7)
  align = 'v',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.3),
  # rel_heights = c(-20, -4),
  rel_widths=c(1, 1), axis="b",
  labels = c(""), label_size = 10, vjust=1
)
grid1

grid2 <- plot_grid(
  p_enrich_scent+theme(legend.position = "none") + xlab("") +
    scale_x_continuous(limits=c(0, 15), breaks=seq(0,15,1)),
  #coord_fixed(xlim = c(1,4.7), ratio=0.4),
  p_scent_meta+
    theme(legend.position = "none")+
    #coord_fixed(xlim = c(1,4.7), ratio=0.6)+
    scale_x_continuous(limits=c(0, 15), breaks=seq(0,15,1))+
    xlab("Enrichment (95% CI)"),
  align = 'v',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.3),
  # rel_heights = c(-20, -6),
  rel_widths=c(1, 1), axis="b",
  labels = c(""), label_size = 10, vjust=1
)

grid2

grid_row2 <- plot_grid(
  grid1,
  grid2,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  # rel_heights = c(1, 1.2, 1),
  rel_widths=c(1.5, 1), axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1
)
grid_row2

legend <- cowplot::get_plot_component(p_scent_meta, 
                                      'guide-box-bottom', return_all = TRUE)
legend <- get_legend(p_scent_meta)
plot_grid(
  grid_row2,
  legend,
  # legend,
  align = 'v',
  hjust = -1, # -1
  nrow = 2, ncol = 1,
  rel_heights = c(1, 0.05), rel_widths=c(1, 1)
)

ggsave2(paste0(plots_dir, "enrich_method.png"),
        width = full_width*1.5, height = full_height*1.7, units = 'cm', dpi = 300)

