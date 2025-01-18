
#### pLI ####
p_pLI <- NB_all_df_s %>%
  group_by(phenotype_id) %>%
  mutate(group = ifelse(sum(is_egene) > 0, "eGene", "non-eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, group, pLI) %>%
  ggplot(aes(x = group, y = pLI)) +  
  geom_violin(width=1, aes(fill = group)) +
  geom_boxplot(width=0.05, color="black", alpha=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = c("eGene"="#D55E00", "non-eGene"="white"))+
  xlab("")+ylab("pLI score")+
  axis_theme; p_pLI

p_loeuf <- NB_all_df_s %>%
  group_by(phenotype_id) %>%
  mutate(group = ifelse(sum(is_egene) > 0, "eGene", "non-eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, group, loeuf) %>%
  ggplot(aes(x = group, y = loeuf)) + #geom_boxplot(outlier.size = 0.5) + 
  geom_violin(width=0.6, aes(fill = group)) +
  geom_boxplot(width=0.2, color="black", alpha=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = c("eGene"="#D55E00", "non-eGene"="white"))+
  xlab("")+ylab("LOEUF score")+
  axis_theme; p_loeuf

p_rvis <- NB_all_df_s %>%
  group_by(phenotype_id) %>%
  mutate(group = ifelse(sum(is_egene) > 0, "eGene", "non-eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, group, RVIS_percentile) %>%
  ggplot(aes(x = group, y = RVIS_percentile)) + #geom_boxplot(outlier.size = 0.5) + 
  geom_violin(width=0.4, aes(fill = group)) +
  geom_boxplot(width=0.2, color="black", alpha=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = c("eGene"="#D55E00", "non-eGene"="white"))+
  xlab("")+ylab("RVIS quantiles")+
  axis_theme; p_rvis

p_eds <- NB_all_df_s %>%
  group_by(phenotype_id) %>%
  mutate(group = ifelse(sum(is_egene) > 0, "eGene", "non-eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, group, EDS) %>%
  ggplot(aes(x = group, y = EDS)) + #geom_boxplot(outlier.size = 0.5) + 
  geom_violin(width=0.6, aes(fill = group)) +
  geom_boxplot(width=0.2, color="black", alpha=0.2, outlier.size = 0.5) +
  scale_fill_manual(values = c("eGene"="#D55E00", "non-eGene"="white"))+
  xlab("")+ylab("Enhancer domain score (EDS)")+
  axis_theme; p_eds

grid1_cons <- plot_grid(
  p_pLI + theme(legend.position = "none"),
  p_loeuf + theme(legend.position = "none"),
  p_rvis + theme(legend.position = "none"),
  p_eds + theme(legend.position = "none"),
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  #rel_heights = c(1, 1, 1), rel_widths=c(1, 1, 1, 1), # axis="r",
  labels = c("A"), 
  label_size = 10, vjust=1
)

legend_row <- get_plot_component(p_pLI + theme(legend.position = "bottom",
                                               legend.key.size = unit(0.3, units='cm'),
                                               legend.margin = margin(0, 0, 20, 0)), 
                                 "guide-box-bottom")

grid1_cons <- plot_grid(
  grid1_cons, legend_row,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 0.1), 
  rel_widths=c(1, 1), axis="b",
  labels = c("A"), 
  label_size = 10, vjust=1
)
grid1_cons

#### effect size and AF overall ####
p_effectsize <- NB_egenes_df %>%
  mutate(slope = slope_wald, group = "eGene") %>%
  filter(!is.na(pLI)) %>%
  mutate(af = ifelse(slope > 0, af, 1-af),
         pLI_high = ifelse(pLI >= 0.9, TRUE, FALSE)) %>% 
  ggplot(aes(x = af, y = abs(slope))) +
  geom_point(size=0.5, color = "#D55E00") +
  xlab("Allele Frequency")+ylab("|sc-eQTL effect size|")+
  geom_smooth(method = "loess", color = "blue") +
  axis_theme + theme(legend.position = "none"); p_effectsize


#### effect size density plot #### 

p_density <- NB_egenes_df %>%
  mutate(slope = slope_wald) %>%
  filter(!is.na(pLI)) %>%
  mutate(pLI = ifelse(pLI >= 0.9, "High", "Low")) %>% 
  ggplot(aes(x = abs(slope), fill = pLI, color = pLI)) +
  geom_density(alpha=0.8) + 
  scale_fill_manual(values = c("High"="orange", "Low"="grey")) +
  scale_color_manual(values = c("High"="orange", "Low"="grey")) +
  xlab("|sc-eQTL Effect size|") + axis_theme +
  theme(legend.position = c(.75, .85),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.3, units = "cm")); p_density


grid2_ef <- plot_grid(
  p_effectsize, p_density,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  #rel_heights = c(1, 1, 1), rel_widths=c(1, 1, 1, 1), # axis="r",
  labels = c("B", "C"), 
  label_size = 10, vjust=1
)
grid2_ef

plot_grid(
  grid1_cons, grid2_ef,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  #rel_heights = c(1, 1, 1), rel_widths=c(1, 1, 1, 1), # axis="r",
  #labels = c("A"), 
  label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "mainfig_selection.png"),
        width = full_width, height = full_height*1.5, units = 'cm', dpi = 300)

