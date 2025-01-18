# main figure
######## egene (A) ######## 
# order eGenes by total number of cells
celltype_bquote_0.01 <- c("CD4_NC" = bquote(CD4[NC]),
                     "CD4_ET*" = bquote(CD4[ET] ~ "*"),
                     "CD4_SOX4" = bquote(CD4[SOX4]),
                     "CD8_ET*" = bquote(CD8[ET] ~ "*"),
                     "CD8_NC" = bquote(CD8[NC]),
                     "CD8_S100B" = bquote(CD8[S100B]),
                     "NK*" = bquote(NK ~ "*"),
                     "NK_R" = bquote(NK[R]),
                     "Plasma*" = bquote(Plasma ~ "*"),
                     "B_Mem*" = bquote(B[Mem] ~ "*"),
                     "B_IN" = bquote(B[IN]),
                     "Mono_C*" = bquote(Mono[C] ~ "*"),
                     "Mono_NC*" = bquote(Mono[NC] ~ "*"),
                     "DC" = bquote(DC))

p_egene_0.01 <- egenes_df %>% 
  mutate(cell = ifelse(p_diff_nb_lm_score < 0.05/14, paste0(cell,"*"), cell),
         cell = fct_reorder(cell, desc(total_cell))) %>% 
  gather(key = model, value = eGenes, c(jaxqtl_nb, jaxqtl_lm_score, jaxqtl_pois)) %>%
  mutate(model = case_when(model == "jaxqtl_nb" ~ "negbinom", 
                            model == "jaxqtl_pois" ~ "Poisson",
                            model == "jaxqtl_lm_score" ~ "linear",
                            model == "tqtl" ~ "tensorqtl"),
         model = fct_relevel(model, "negbinom")) %>% 
  ggplot(aes(x = cell, y = eGenes, fill = model)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = method_cols) + xlab("")+
  scale_x_discrete(
    labels = celltype_bquote_0.01
  )+
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        # legend.position = c(.88, .95), # for single plot
        legend.position = c(.93, .95), # for stacked plot in main figure
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.3, units = "cm")
)
p_egene_0.01

ggsave2(paste0(plots_dir, "Fig5_A.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

celltype_bquote_0.1 <- c("CD4_NC*" = bquote(CD4[NC]~ "*"),
                          "CD4_ET" = bquote(CD4[ET]),
                          "CD4_SOX4" = bquote(CD4[SOX4]),
                          "CD8_ET*" = bquote(CD8[ET] ~ "*"),
                          "CD8_NC" = bquote(CD8[NC]),
                          "CD8_S100B" = bquote(CD8[S100B]),
                          "NK*" = bquote(NK ~ "*"),
                          "NK_R" = bquote(NK[R]),
                          "Plasma*" = bquote(Plasma ~ "*"),
                          "B_Mem" = bquote(B[Mem]),
                          "B_IN" = bquote(B[IN]),
                          "Mono_C" = bquote(Mono[C]),
                          "Mono_NC" = bquote(Mono[NC]),
                          "DC*" = bquote(DC~ "*"))

p_egene_0.1 <- egenes_df_0.1 %>% 
  left_join(egenes_df %>% select(cell, total_cell), by="cell") %>% 
  mutate(cell = ifelse(p_diff < 0.05/14, paste0(cell,"*"), cell),
         cell = fct_reorder(cell, desc(total_cell))) %>% 
  gather(key = software, value = eGenes, c(jaxqtl_nb, tqtl, saige_egenes)) %>%
  mutate(eGenes = as.numeric(eGenes),
         software = case_when(software == "jaxqtl_nb" ~ "jaxQTL", 
                              software == "saige_egenes" ~ "SAIGE-QTL",
                              software == "tqtl" ~ "tensorQTL"),
         software = fct_relevel(as.factor(software), c("jaxQTL", "SAIGE-QTL", "tensorQTL"))) %>% 
  ungroup() %>% 
  ggplot(aes(x = cell, y = eGenes, fill = software)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_x_discrete(
    labels = celltype_bquote_0.1
  )+
  scale_fill_manual(values = method_cols) +
  xlab("") +
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.3, units = "cm")
  )
p_egene_0.1

ggsave2(paste0(plots_dir, "Fig5_B.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


plot_grid(
  p_egene_0.01, 
  p_egene_0.1,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 1), rel_widths=c(1, 1), axis="r",
  labels = c("A", "B"), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "mainfig_egene.png"),
        width = full_width, height = full_height*1.5, units = 'cm', dpi = 300)

# jaxQTL-Linear and tensorQTL
p_egene_lm <- egenes_df_0.1 %>% 
  left_join(egenes_df %>% select(cell, total_cell), by="cell") %>% 
  mutate(cell = ifelse(p_diff < 0.05/14, paste0(cell,"*"), cell),
         cell = fct_reorder(cell, desc(total_cell))) %>% 
  gather(key = software, value = eGenes, c(jaxqtl_lm_score, tqtl)) %>%
  mutate(eGenes = as.numeric(eGenes),
         software = case_when(software == "jaxqtl_lm_score" ~ "jaxQTL-Linear", 
                              software == "tqtl" ~ "tensorQTL")) %>% 
  ungroup() %>% 
  ggplot(aes(x = cell, y = eGenes, fill = software)) +
  geom_bar(position = "dodge", stat="identity") +
  scale_fill_manual(values = method_cols) +
  xlab("cell type") +
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size = unit(0.3, units = "cm")
  )

lt <- list(tensorQTL = unique(paste0(tqtl_egenes_df_0.1$phenotype_id, "_", tqtl_egenes_df_0.1$celltype)),
           jaxQTL_Linear = unique(paste0(jqtl_lm_score_egenes_df_0.1$phenotype_id, "_", 
                                         jqtl_lm_score_egenes_df_0.1$celltype)))
m <- make_comb_mat(lt)
p_upset <- grid.grabExpr(draw(UpSet(m,
                                    top_annotation = upset_top_annotation(m,
                                                                          axis_param = list(gp = gpar(fontsize = 5)),
                                                                          annotation_name_gp = gpar(fontsize = 5),
                                                                          height = unit(1, "cm")),
                                    right_annotation = upset_right_annotation(m,
                                                                              axis_param = list(gp = gpar(fontsize = 5)),
                                                                              annotation_name_gp = gpar(fontsize = 5),
                                                                              width = unit(1, "cm")),
                                    row_names_gp = gpar(fontsize = 5))))


plot_grid(
  p_egene_lm, 
  p_upset,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1.5, 1), axis="r",
  labels = c("A", "B"), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "egene_jaxqtl_lm_tqtl.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


######## upset (B) ######## 
p_upset
lt <- list(NegBinom = unique(paste0(NB_egenes_df$phenotype_id, "_", NB_egenes_df$celltype)), 
           # tensorqtl = unique(paste0(tqtl_egenes_df$phenotype_id, "_", tqtl_egenes_df$celltype)),
           Linear = unique(paste0(jqtl_lm_score_egenes_df$phenotype_id, "_", jqtl_lm_score_egenes_df$celltype)),
           Poisson = unique(paste0(pois_egenes_df$phenotype_id, "_", pois_egenes_df$celltype)))
m <- make_comb_mat(lt)
p_upset <- grid.grabExpr(draw(UpSet(m,
                                    top_annotation = upset_top_annotation(m,
                                                                          axis_param = list(gp = gpar(fontsize = 5)),
                                                                          annotation_name_gp = gpar(fontsize = 5),
                                                                          height = unit(1, "cm")),
                                    right_annotation = upset_right_annotation(m,
                                                                              axis_param = list(gp = gpar(fontsize = 5)),
                                                                              annotation_name_gp = gpar(fontsize = 5),
                                                                              width = unit(1, "cm")),
                                    row_names_gp = gpar(fontsize = 5))))

######## compare metric (C) ######## 
metric <- "express_percent" # pLI, 
cols_select <- c("phenotype_id", "celltype", "count_mean", "express_percent", "pLI", "loeuf", "EDS")
egenes_bind <- bind_rows(both_hits %>% select(all_of(cols_select)) %>% mutate(category = "both"),
                         nb_only %>% select(all_of(cols_select)) %>% mutate(category = "NegBinom"),
                         tqtl_only %>% select(all_of(cols_select)) %>% mutate(category = "Linear")) %>% 
  mutate(category = fct_relevel(category, c("both", "Linear", "NegBinom")))

not_egenes <- NB_all_df %>% select(all_of(cols_select)) %>% 
  anti_join(egenes_bind, by=c("celltype", "phenotype_id"))

# egenes alpha: median 0.045, mean 0.258
# not egene alpha: median 0.021 mean 0.32
egenes_bind %>% distinct(celltype, phenotype_id, .keep_all = T) %>% summary()
not_egenes %>% distinct(celltype, phenotype_id) %>% 
  inner_join(NB_all_df, by=c("phenotype_id", "celltype")) %>%  summary()

wilcox.test(express_percent ~ group,
data=bind_rows(egenes_bind, not_egenes %>% mutate(category = "none")) %>% 
  mutate(group = ifelse(category == "none", "none", "has_egene"))) %>% broom::tidy()

bind_rows(egenes_bind, not_egenes %>% mutate(category = "none")) %>% 
  mutate(group = ifelse(category == "none", "none", "has_egene")) %>% 
  distinct(group, phenotype_id, .keep_all = T) %>% group_by(group) %>% 
  summarise(mean = mean(EDS, na.rm = T), median = median(EDS, na.rm=T))

# EDS: eGenes vs. not eGenes, p=8.54E-4
# louef: eGenes vs. not eGenes, 3.30E-11
wilcox.test(loeuf ~ category,
            data=bind_rows(egenes_bind %>% mutate(category = "eGenes"), 
                           not_egenes %>% mutate(category = "not eGenes")) %>% 
              distinct(category, phenotype_id, .keep_all = T)) %>% broom::tidy()


#######

p_express_percent <- bind_rows(egenes_bind) %>% 
  mutate(model = fct_relevel(category, c("both", "Linear", "NegBinom"))) %>% 
  ggplot(aes(x = model, y = express_percent, fill = model)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_manual(values = method_cols)+
  axis_theme + ylab("coverage")
p_express_percent

wilcox.test(express_percent ~ category,
            data=egenes_bind %>% filter(category %in% c("both", "NegBinom"))) %>% broom::tidy()

ggsave2(paste0(plots_dir, "egene_express_diff.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)

bind_rows(egenes_bind) %>% 
  mutate(model = fct_relevel(category, c("both", "Linear", "NegBinom"))) %>% 
  ggplot(aes(x = model, y = log(count_mean), fill = model)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_fill_manual(values = method_cols)+
  axis_theme + ylab("coverage")

wilcox.test(count_mean ~ category,
            data=egenes_bind %>% filter(category %in% c("Linear", "NegBinom"))) %>% broom::tidy()

# look at all eGenes by cell type
p_byct <- bind_rows(NB_egenes_df %>% select(phenotype_id, celltype, express_percent, count_mean) %>% mutate(model = "NegBinom"),
          jqtl_lm_score_egenes_df %>% select(phenotype_id, celltype, express_percent, count_mean) %>% mutate(model = "Linear")) %>% 
  mutate(celltype = fct_relevel(celltype, celltype_order)) %>% 
  ggplot(aes(x = celltype, y = express_percent, fill = model)) + geom_boxplot(outlier.size=0.5) +
  stat_summary(aes(group=model), fun.y=mean, position = position_dodge(0.75),
               geom="point", shape=18, size=2, color="yellow") +
  scale_fill_manual(values = method_cols) + axis_theme +
  theme(axis.text.x = element_text(angle=30, hjust = 1)); p_byct

p_cts <- bind_rows(NB_egenes_df %>% select(phenotype_id, celltype, express_percent, count_mean) %>% mutate(model = "NegBinom"),
                    jqtl_lm_score_egenes_df %>% select(phenotype_id, celltype, express_percent, count_mean) %>% mutate(model = "Linear")) %>% 
  mutate(celltype = fct_relevel(celltype, celltype_order)) %>% 
  ggplot(aes(x = model, y = express_percent, fill = model)) + geom_boxplot() +
  stat_summary(aes(group=model), fun.y=mean, position = position_dodge(0.75),
               geom="point", shape=18, size=2, color="yellow") +
  scale_fill_manual(values = method_cols) + axis_theme 

plot_grid(
  p_cts + theme(legend.position = "none"), 
  p_byct,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 3.5), axis="r",
  labels = c("A", "B"), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "egene_express_diff_byCT.png"),
        width = full_width*1.2, height = full_height, units = 'cm', dpi = 300)


plot_grid(
  p_upset, 
  p_express_percent + theme(legend.position = "none"),
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1.2, 1), axis="r",
  labels = c("A", "B"), label_size = 10, vjust=1
)

ggsave2(paste0(plots_dir, "egene_upset_express.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

for (CT in unique(NB_egenes_df$celltype)){
  fit <- wilcox.test(NB_egenes_df %>% filter(celltype == CT) %>% pull(express_percent),
              jqtl_lm_score_egenes_df %>% filter(celltype == CT) %>% pull(express_percent),
              "l")
  print(glue("{CT}: {fit$p.value}"))
}

######## replication (D) ######## 

rep_eqtl_df <- data.frame(model = c("NegBinom", "Linear"),
                     total = c(18907, 16654),
                     eQTLgen = c(14406, 12765),
                     GTEx = c(9225, 8198))
                     #Yazar2022 = c(6559, 7574))

p_replicate_eqtl <- rep_eqtl_df %>% 
  gather(key = "rep_data", value = "rep", eQTLgen, GTEx) %>% 
  mutate(rep=ifelse(rep > 1, rep/total, rep)) %>% 
  ggplot(aes(x = rep_data, y = rep, fill = model)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = method_cols) +
  geom_text(aes(# label = scales::percent(round(rep, 3)), 
    label = round(rep, 2),
                group=model), 
            position = position_dodge(width = 0.9),
            vjust=-0.5, hjust=0.5, size=2.5)+
  ylim(0,1) + axis_theme + xlab("replication data")
p_replicate_eqtl

p_replicate_egene <- rep_egene_df %>% 
  gather(key = "rep_data", value = "rep", eQTLgen, GTEx, Yazar2022) %>% 
  mutate(rep=ifelse(rep > 1, rep/total, rep)) %>% 
  ggplot(aes(x = rep_data, y = rep, fill = method)) + 
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values = method_cols) +
  geom_text(aes(# label = scales::percent(round(rep, 3)), 
    label = round(rep, 2),
    group=method), 
    position = position_dodge(width = 0.9),
    vjust=-0.5, hjust=0.5, size=2.5)+
  ylim(0,1) + axis_theme + ggtitle("Replication of eGenes")

p_replicate_egene

grid2 <- plot_grid(
  p_upset, 
  p_express_percent + theme(legend.position = "none"), 
  p_replicate_eqtl + theme(legend.position = "none"),
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1, 1), rel_widths=c(1, 1.1, 1), axis="r",
  labels = c("C", "D", "E"), label_size = 10, vjust=1
)

plot_grid(
  grid1,
  grid2, 
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  rel_heights = c(1, 1), rel_widths=c(1, 1),
  labels = c(""),  label_size = 10
)

ggsave2(paste0(plots_dir, "Fig4_new.png"),
        width = full_width, height = full_height*1.5, units = 'cm', dpi = 300)


