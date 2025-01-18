# Figure: finemapping and enrichment analysis

# average within eGene
tmp <- nb_cs_pip0.5 %>% 
  # filter(pip >= 0.9) %>% 
  left_join(gtf %>% select(phenotype_id, chr, tss=end), by=c("chr", "phenotype_id")) %>% 
  mutate(tss_dist = pos - tss,
         # bins = (bins - 1000),
         sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

tmp %>% distinct(phenotype_id, snp, tss_dist, sc_bulk) %>% group_by(sc_bulk) %>% 
  summarize(mean = mean(tss_dist),
            median = median(tss_dist),
            min = min(tss_dist), max = max(tss_dist))

wilcox.test(abs(tss_dist) ~ sc_bulk, data=tmp %>% 
              distinct(phenotype_id, snp, tss_dist, sc_bulk))

# TSS
# p_tss <- tmp %>% distinct(sc_bulk, phenotype_id, snp, .keep_all = T) %>% 
#   ggplot(aes(x = tss_dist)) + geom_histogram(aes(fill = sc_bulk)) +
#   scale_fill_manual(values=sc_bulk_cols)+
#   facet_grid(.~sc_bulk)+
#   axis_theme
# p_tss

# TSS of lead SNPs, create bins of each size of 125000 bps
# count percent of variants fall into each bin per method, so that prop sum to 1 per method
# tss_breaks <- c(0, 5000, 50000, 100000, 5000000) # 
tss_breaks <- c(0, 5000, 100000, 500000) # 

wilcox.test(tss_dist ~ sc_bulk,
  data=tmp %>% 
  mutate(tss_dist = abs(tss_dist),
         tss_dist = ifelse(tss_dist < 1, tss_dist+1, tss_dist), # put tss_dist=0 to first bin
         bins = cut(tss_dist, 
                    breaks = tss_breaks)) %>% 
  distinct(sc_bulk, phenotype_id, snp, tss_dist, bins) # %>% 
  # group_by(sc_bulk) %>% summarize(mean = mean(tss_dist), median = median(tss_dist))
  )

tss_bins <- tmp %>% 
  mutate(tss_dist = abs(tss_dist),
         tss_dist = ifelse(tss_dist < 1, tss_dist+1, tss_dist), # put tss_dist=0 to first bin
         bins = cut(tss_dist, 
                    breaks = tss_breaks)) %>% 
  distinct(sc_bulk, phenotype_id, snp, tss_dist, bins) %>% # remove duplicate of cell types
  # count(sc_bulk)
  add_count(sc_bulk, name = "num_var_method") %>% 
  add_count(sc_bulk, bins, name = "num_var_bins") %>% 
  distinct(sc_bulk, bins, num_var_method, num_var_bins) %>% 
  mutate(prop = num_var_bins/num_var_method)
tss_bins %>% arrange(sc_bulk, bins)

test_prop <- tibble(bins = unique(tss_bins$bins), p = NA) %>% arrange(bins)
for (bin in test_prop$bins){
  tss_bins_tmp <- tss_bins %>% filter(bins == bin)
  test_res <- prop.test(x = c(tss_bins_tmp$num_var_bins[which(tss_bins_tmp$sc_bulk == "bulk-eQTL")], 
                              tss_bins_tmp$num_var_bins[which(tss_bins_tmp$sc_bulk == "sc-eQTL")]), 
                        n = c(tss_bins_tmp$num_var_method[which(tss_bins_tmp$sc_bulk == "bulk-eQTL")], 
                              tss_bins_tmp$num_var_method[which(tss_bins_tmp$sc_bulk == "sc-eQTL")])) 
  test_prop$p[which(test_prop$bins == bin)] <- broom::tidy(test_res) %>% pull(p.value)
  
}

test_prop

# TODO: add P values on top of each bin bars
library(superb)
p_tss <- tss_bins %>% 
  mutate(bins = case_when(bins == "(0,5e+03]" ~ "<5kb",
                          bins == "(1e+05,5e+05]" ~ ">100kb",
                          bins == "(5e+03,1e+05]" ~ "5kb-100kb")) %>% 
  ggplot(aes(x = factor(bins, levels = c("<5kb", "5kb-100kb", ">100kb")),
             y = prop, group=sc_bulk)) +
  geom_bar(aes(fill = sc_bulk), stat="identity", position = "dodge2") +
  scale_fill_manual(values=sc_bulk_cols)+
  showSignificance(c(0.8,1.2), 0.28, -0.01, "P = 6.71E-4") +
  showSignificance(c(2.8,3.2), 0.28, -0.01, "P = 1.07E-4") +
  xlab("Distance to TSS")+ylab("Fraction of eQTLs in each distance bin")+
  guides(fill=guide_legend(title="group"))+
  ggtitle("Distribution of fine-mapped eQTLs")+
  axis_theme

p_tss
ggsave2(paste0(plots_dir, "Fig5_C.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)


# note: combine the last two bins because B_Mem and Mono_C don't have PIP>0.8
p_pip <- tmp %>%  
  mutate(tss_dist = abs(tss_dist),
         tss_dist = ifelse(tss_dist < 1, tss_dist+1, tss_dist), # put tss_dist=0 to first bin
         bins = cut(tss_dist, 
                    breaks = tss_breaks, labels = 1:(length(tss_breaks)-1))) %>% 
  add_count(sc_bulk, celltype, name = "num_var_method") %>% 
  add_count(sc_bulk, celltype, bins, name = "num_var_bins") %>% 
  distinct(sc_bulk, celltype, bins, num_var_method, num_var_bins) %>% 
  mutate(prop = num_var_bins/num_var_method) %>% 
  group_by(sc_bulk, celltype) %>% 
  summarize(ratio_bin3_bin1 = prop[which(bins==3)]/prop[which(bins==1)],
            sum_prop = sum(prop)) %>% ungroup() %>% 
  mutate(celltype = fct_relevel(celltype, rev(c(celltype_order, "allcells")))) %>% 
  ggplot(aes(y = celltype, x = ratio_bin3_bin1, fill = celltype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=celltype_cols) +
  axis_theme
p_pip

tmp %>%  
  left_join(nb_cs_pip0.5_wald %>% select(celltype, phenotype_id, snp, slope, slope_se)) %>% 
  mutate(tss_dist = abs(tss_dist),
         tss_dist = ifelse(tss_dist < 1, tss_dist+1, tss_dist), # put tss_dist=0 to first bin
         bins = cut(tss_dist, 
                    breaks = tss_breaks, labels = 1:(length(tss_breaks)-1))) %>% 
  add_count(sc_bulk, celltype, name = "num_var_method") %>% 
  add_count(sc_bulk, celltype, bins, name = "num_var_bins") %>% 
  group_by(sc_bulk, celltype, bins) %>% 
  mutate(mean_slope = mean(abs(slope))) %>% ungroup() %>% 
  distinct(sc_bulk, celltype, bins, num_var_method, num_var_bins, mean_slope) %>% 
  mutate(prop = num_var_bins/num_var_method) %>% 
  group_by(sc_bulk, celltype) %>% 
  summarize(ratio_bin3_bin1 = prop[which(bins==3)]/prop[which(bins==1)],
            ratio_slope_bins3_bin1 = mean_slope[which(bins==3)]/mean_slope[which(bins==1)],
            sum_prop = sum(prop)) %>% ungroup() %>% 
  mutate(celltype = fct_relevel(celltype, rev(c(celltype_order, "allcells")))) %>% 
  ggplot(aes(y = celltype, x = ratio_bin3_bin1, fill = celltype)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values=celltype_cols) +
  axis_theme

  
# count number of finemapped eQTLs, at least in one cell type
# fine-mapped eQTLs is more distal than bulk (P=)
tmp %>% distinct(sc_bulk, phenotype_id, snp, .keep_all = T) %>% 
  group_by(sc_bulk) %>% 
  summarise(mean_dist = mean(abs(tss_dist)),
            med_dist = median(abs(tss_dist)),
            n = n())

wilcox.test(abs(tss_dist) ~ sc_bulk, data = tmp %>% distinct(sc_bulk, phenotype_id, snp, .keep_all = T)) %>% 
  broom::tidy()

# baseline enrichment
## NB sc vs. bulk
p_enrich_baseline <- forest_plot_theme(bind_rows(baseline_meta_celltype %>% 
                                                   filter(method == "NegBinom") %>% mutate(method = "sc-eQTL"),
                                                 nb_meta %>% filter(celltype == "allcells" & annot %in% baseline_annot_list) %>% 
                                                   select(annot, meta_OR, meta_lower, meta_upper, meta_p) %>% 
                                                   mutate(method = "bulk-eQTL")) %>% 
                                         # filter(meta_OR < 500) %>% 
                                         mutate(annot = gsub("_phastCons46way$", "", annot)) %>% 
                                         group_by(annot) %>% 
                                         mutate(jaxqtl_or = meta_OR[which(method == "sc-eQTL")]) %>% ungroup() %>% 
                                         mutate(annot = fct_reorder(annot, jaxqtl_or)) %>% 
                                         ggplot(aes(y = annot, x = meta_OR, color = method)) +
                                         geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                                         geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                        height = 0.25, 
                                                        position = position_dodge(width=0.5)) +
                                         geom_vline(xintercept = 1, color = "red", 
                                                    linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                         scale_color_manual(values = sc_bulk_cols) +
                                         xlab("Odds Ratio (95% CI)") + 
                                         ylab(" ")+ 
                                         theme_bw()+
                                         ggtitle("Baseline annotations")) + 
  theme(legend.position = c(.95, .5), # for stacked plot in main figure
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.5, units = "cm"))
p_enrich_baseline

## ATAC
# check annot order are the same
all.equal(atac_meta_celltype %>% filter(method == "NegBinom") %>% mutate(method = "sc-eQTL") %>% pull(annot),
          atac_meta_allcells %>% filter(method == "NegBinom") %>% pull(annot))

all.equal(atac_meta_celltype %>% filter(method == "NegBinom") %>% pull(annot), 
          atac_meta_allcells %>% filter(method == "NegBinom") %>% pull(annot))

atac_dat <- bind_rows(atac_meta_celltype %>% filter(method == "NegBinom") %>% mutate(method = celltype), 
          atac_meta_allcells %>% filter(method == "NegBinom") %>% 
            mutate(method = "bulk-eQTL",
                   # here checked the order is the same, need to match cell type
                   celltype = atac_meta_celltype %>% filter(method == "NegBinom") %>% pull(celltype))) %>% 
  mutate(celltype = fct_relevel(celltype, rev(celltype_order)))

p_enrich_atac <- forest_plot_theme(atac_dat %>% 
                                     # filter(meta_OR<100) %>% 
                    ggplot(aes(y = celltype, x = meta_OR, color = method)) +
                    geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    scale_color_manual(values = celltype_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    theme_bw() +
                    ggtitle("ATAC-seq"))
p_enrich_atac

## epimap
epimap_dat <- bind_rows(epimap_meta_celltype %>% 
                          filter(method == "NegBinom" & celltype != "allcells") %>% select(annot, starts_with("meta"), celltype) %>% 
                          mutate(sc_bulk = "sc-eQTL",
                                 method = celltype),
                        epimap_meta_celltype %>% 
                          filter(method == "NegBinom" & celltype == "allcells") %>% 
                          distinct(annot, .keep_all = T) %>% 
                          select(annot, starts_with("meta")) %>% 
                          right_join(epimap_meta_celltype %>% 
                                       filter(method == "NegBinom" & celltype != "allcells") %>% 
                                       select(annot, celltype), by="annot") %>% 
                          mutate(sc_bulk = "bulk-eQTL",
                                 method = "bulk-eQTL")) %>% 
  mutate(category = ifelse(grepl("enhancer", annot), "EpiMap enhancer", "EpiMap promoter")) %>% 
  select(colnames(atac_dat), everything())

allfunc_dat <- bind_rows(epimap_dat,
          bind_rows(nb_enrich %>% mutate(method = celltype),
          nb_allcell_df_tmp %>% mutate(method = "bulk-eQTL")) %>% 
  mutate(category = "SCENT enhancer-gene") %>% 
  rename(enrich = mean_enrich) %>% 
  filter(celltype != "allcells") %>% 
  mutate(annot = celltype, meta_OR = enrich, meta_lower = OR_L, meta_upper = OR_U,
         meta_p = NA, meta_se = se, category = "SCENT enhancer-gene",
         sc_bulk = method) %>% 
  select(names(epimap_dat)))

p_enrich_epimap <- forest_plot_theme(bind_rows(atac_dat %>% mutate(category = "ATAC-seq"),
                                               epimap_dat) %>%
                                       rename(`Cell type` = "method") %>% 
                                              # allfunc_dat) %>% 
                                       #filter(meta_OR < 100) %>% 
                                       mutate(celltype = fct_relevel(celltype, rev(celltype_order)),
                                              category = fct_relevel(category, c("ATAC-seq", "EpiMap promoter", "EpiMap enhancer"))) %>% 
                                       ggplot(aes(y = celltype, x = meta_OR, color = `Cell type`)) +
                                       geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                                       geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                                      height = 0.25, 
                                                      position = position_dodge(width=0.5)) +
                                       facet_grid(.~category)+
                                       geom_vline(xintercept = 1, color = "red", 
                                                  linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                       scale_color_manual(values = celltype_cols, labels = celltype_bquote) +
                                       scale_y_discrete(
                                         labels = celltype_bquote
                                       )+ 
                                       xlab("Odds Ratio (95% CI)") + 
                                       ylab(" ")) + guides(fill=guide_legend(title="Cell type"))
  
p_enrich_epimap

# SCENT result
scent_map
nb_allcell_df_tmp <- nb_enrich %>% select(celltype) %>% 
  left_join(scent_map, by=c("celltype"="yazar_celltype")) %>% 
  left_join(nb_allcell_scent %>% select(-celltype), by=c("annot_celltype"="ref_cell")) %>% 
  select(-annot_celltype)

p_enrich_scent <- forest_plot_theme(bind_rows(nb_enrich %>% mutate(method = celltype),
                                              nb_allcell_df_tmp %>% mutate(method = "bulk-eQTL")) %>% 
                                      mutate(category = "SCENT enhancer-gene") %>% 
                                      rename(enrich = mean_enrich) %>% 
                                      filter(celltype != "allcells") %>% 
                                      mutate(celltype = fct_relevel(celltype, rev(celltype_order))) %>% 
                                      ggplot(aes(y = celltype, x = enrich, color = method)) + 
                                      geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
                                      geom_errorbarh(aes(xmin = OR_L, xmax = OR_U),
                                                     height = 0.25,
                                                     position = position_dodge(0.5)) +
                                      geom_vline(xintercept = 1, color = "red", 
                                                 linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                      facet_grid(.~category)+
                                      scale_color_manual(values = celltype_cols)+
                                      axis_theme + ylab("")+
                                      xlab("Odds Ratio (95% CI)"))
p_enrich_scent

# add meta result across 14 cell type vs. bulk result
atac_meta
epimap_meta
scent_meta

p_epimap_atac_meta <- forest_plot_theme(bind_rows(epimap_meta %>% filter(method == "NegBinom") %>% 
            mutate(category = paste0("EpiMap ", category)),
          atac_meta %>% mutate(category = "ATAC-seq",
                               method = "NegBinom") %>% select(sc_bulk, category, everything()) %>%
          #scent_meta %>% mutate(category = "SCENT enhancer-gene", method = "NegBinom") %>%
            select(sc_bulk, category, method, everything())) %>% 
            mutate(category = fct_relevel(category, 
                                          c("ATAC-seq", "EpiMap promoter", "EpiMap enhancer", "SCENT enhancer-gene"))) %>% 
  ggplot(aes(y = sc_bulk, x = meta_OR, color = sc_bulk)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = sc_bulk_cols)+
  facet_grid(.~category)+
  axis_theme + ylab("")+
  xlab("Odds Ratio (95% CI)")) + theme(strip.text = element_blank())
  
p_scent_meta <- scent_meta %>% 
  ggplot(aes(y = sc_bulk, x = meta_OR, color = sc_bulk)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = sc_bulk_cols)+
  axis_theme + ylab("")+
  xlab("Odds Ratio (95% CI)")

grid_row1 <- plot_grid(
  p_enrich_baseline + theme(legend.key.size = unit(0.5, "cm"),
                            legend.text = element_text(size = 8),
                            legend.title = element_text(size = 8)),
  p_tss + theme(legend.position = "none"), 
  #p_pip + theme(legend.position = "none"),
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1, 1), rel_widths=c(1.2, 1, 1), # axis="r",
  labels = c("A", "B", "C"), label_size = 10, vjust=1
)
grid_row1

grid1 <- plot_grid(
  p_enrich_epimap+theme(legend.position = "none",
                        plot.margin = unit(c(0,0,0,0), "cm")) + xlab("") + xlim(0,5),
    #coord_fixed(ratio = 0.3, xlim = c(0, 5)),
  p_epimap_atac_meta+theme(legend.position = "none",
                           plot.margin = unit(c(0,0,0.2,0), "cm"))+
    xlab("Enrichment (95% CI)") + xlim(0,5),
    #coord_fixed(ratio = 0.3, xlim = c(0, 5)),
  align = 'v',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.25), # 0.3
  # rel_heights = c(-20, -10), # -20, -4
  rel_widths=c(1, 1), axis="b",
  labels = c("C"), label_size = 10, vjust=1
)
grid1

ggsave2(paste0(plots_dir, "Fig5_scent.png"),
        width = onecol_width, height = onecol_height*1.5, units = 'cm', dpi = 300)

ct_legend <- get_plot_component(p_enrich_epimap+theme(legend.position = "bottom",
                                                      legend.key.size = unit(0.1, units="cm"),
                                                      legend.text = element_text(size = 7),
                                                      legend.title = element_text(size = 7))+
                                  guides(color=guide_legend(title="Cell type", nrow=2)), "guide-box-bottom")
plot_grid(
  grid_row1, grid1,
  ct_legend,
  align = 'v',
  hjust = -1, # -1
  nrow = 3,
  rel_heights = c(1, 1.2, 0.1), rel_widths=c(1, 1, 1)
)

ggsave2(paste0(plots_dir, "Fig5_sc_bulk.png"),
        width = full_width*1.2, height = full_height*2.3, units = 'cm', dpi = 300)


# ATACseq only

plot_grid(
  p_enrich_atac+scale_x_continuous(limits=c(0, 4), breaks=seq(0,4,1))+
    theme(legend.position = "right") + xlab("")+
    theme(
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 5),  # Adjust legend text size
      legend.title = element_text(size = 5)  # Adjust legend title size
    )+guides(color=guide_legend(title=""))+
    ggtitle("Enrichment of open chromatin (ATAC-seq) in eQTLs"),
  atac_meta %>% mutate(category = "ATAC-seq",
                       method = "NegBinom") %>% 
  ggplot(aes(y = sc_bulk, x = meta_OR, color = sc_bulk)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = sc_bulk_cols)+
  axis_theme + ylab("")+
  xlab("Odds Ratio (95% CI)") + theme(strip.text = element_blank(),
                                      legend.position = "right")+
    scale_x_continuous(limits=c(0, 4), breaks=seq(0,4,1))+
    theme(
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 5),  # Adjust legend text size
      legend.title = element_text(size = 5)  # Adjust legend title size
    )+guides(color=guide_legend(title="")),
  align = 'v',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.3),
  # rel_heights = c(-20, -6),
  rel_widths=c(1, 1), axis="b",
  labels = c(""), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "Fig5_atac.png"),
        width = onecol_width*1.3, height = onecol_height*1.5, units = 'cm', dpi = 300)

