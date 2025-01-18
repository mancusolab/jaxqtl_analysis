# selection metrics

# pLI
# assume we can assign genes to three categories: 
# pLI closer to 1 indicates that the gene or transcript cannot tolerate protein truncating variation.
# pLI is the probability of belonging in the haploinsufficient class (heterozygous variants are not tolerated).
# cutoff: pLI >= 0.9 
# only for protein-coding gene

# LOUEF (Loss-of-function Observed / expected upper bound fraction)
# https://gnomad.broadinstitute.org/help/constraint#loeuf
# ob/expected ratio and their 90% CI that account for gene length.
# upper bound of this CI is a conservative estimate of the ratio
# Low LOEUF scores indicate strong selection against predicted loss-of-function (pLoF) 
# high LOEUF scores suggest a relatively higher tolerance to inactivation
# cutoff: LOEUF score < 0.6 
# only for protein-coding gene

# EDS (enhancer domain score)
# higher EDS means more redundant enhancer regulatory domain

# RVIS
# departure of common functional mutation from average level across genes
# RVIS < 0, purifying selection
# RVIS > 0, absence of purifying selection, likely balanced or positive selection
# lower score (lower percentile) more intolerant

#
NB_egenes_df_s <- NB_egenes_df %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type),
            by="phenotype_id") %>%
  left_join(rvis_score, by=c("phenotype_id"))

NB_all_df_s <- NB_all_df %>% left_join(NB_egenes_df %>% 
                          mutate(is_egene = 1) %>% 
                          select(phenotype_id, celltype, is_egene)) %>% 
  mutate(is_egene = replace_na(is_egene, 0)) %>% 
  left_join(rvis_score, by=c("phenotype_id")) %>% 
  left_join(gene_lookup %>% select(phenotype_id, gene_type), by=c("phenotype_id"))


nb_cs_ct_s <- nb_cs_ct %>% 
  left_join(pLI_score %>% select(phenotype_id, pLI, loeuf), by="phenotype_id") %>%
  left_join(eds_score %>% select(phenotype_id, EDS), by="phenotype_id") %>%
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type), by = "phenotype_id") %>% 
  left_join(rvis_score, by=c("phenotype_id")) %>% 
  filter(sc_bulk == "sc-eQTL")

#### pLI #### 
# majority is protein-coding genes

# gene expression
p_pli_express <- NB_all_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.01, 0.1, 0.3, 0.5, 0.9, 1))) %>% 
  filter(!is.na(pLI_bins)) %>% 
  group_by(phenotype_id, pLI_bins) %>% 
  mutate(mean_express = mean(express_percent, na.rm = T)) %>% ungroup() %>% 
  ggplot(aes(x = pLI_bins, y = mean_express)) + geom_boxplot() + 
  ylab("coverage") + axis_theme

p_pli_all <- NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>%
  distinct(phenotype_id, category, pLI) %>% 
  ggplot(aes(x = category, y = pLI)) + geom_boxplot() + axis_theme
p_pli_all

p_pli_cell <- NB_all_df_s %>% 
  filter(is_egene > 0) %>% 
  mutate(celltype = fct_relevel(as.factor(celltype), celltype_order)) %>% 
  ggplot(aes(x = celltype, y = pLI, fill = celltype)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = celltype_cols)+
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1), legend.position = "none")
p_pli_cell

p_pli_cs <- nb_cs_ct_s %>%
  group_by(phenotype_id, pLI) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = pLI, y = mean_cs)) +
  geom_point(size=0.6) +
  geom_smooth(method="lm") + 
  axis_theme + xlab("pLI")+
  ylab("average number of CS per gene")
p_pli_cs

# mean CS ~ pLI, p = 0.113, neg slope
lm(mean_cs ~ pLI, data = nb_cs_ct_s %>% 
     group_by(phenotype_id, pLI) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

# pLI ~ eGene (at least one cell type): P=5.37E-51, slope=-0.0981 
lm(pLI ~ category + mean_express, 
   data=NB_all_df_s %>%
  group_by(phenotype_id) %>%
  mutate(category = ifelse(sum(is_egene) > 0, 1, 0),
         mean_express = mean(express_percent, na.rm = T)) %>%
  distinct(phenotype_id, mean_express, category, pLI)) %>% 
  broom::tidy()

# shared eGenes ~ pLI
# GTEx: tissue shared eGenes are depleted of high-pLI genes
lm(pLI ~ n, data= NB_egenes_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.01, 0.1, 0.3, 0.5, 0.9, 1))) %>% 
  group_by(phenotype_id, pLI, pLI_bins) %>% 
  summarize(n = n()) %>% ungroup()) %>% broom::tidy()

NB_egenes_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.01, 0.1, 0.3, 0.5, 0.9, 1))) %>% 
  group_by(phenotype_id, pLI, pLI_bins) %>% 
  summarize(n = n()) %>% 
  # group_by(n) %>% summarize(mean_pLI = mean(pLI, na.rm = T)) %>% ungroup() %>% 
  ggplot(aes(x = as.factor(n), y = pLI)) +
  geom_boxplot() # + geom_smooth(method="lm")

pli_grid <- plot_grid(
  p_pli_cell, p_pli_all,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 0.3)
)
pli_grid

NB_egenes_df %>% distinct(celltype, phenotype_id, pLI) %>% 
  left_join(gtf, by = "phenotype_id") %>% 
  group_by(celltype) %>% 
  summarize(n_egenes = n(),
            pLI_0.9 = sum(pLI >= 0.9, na.rm = T)/n_egenes,
            mean_pLI = mean(pLI, na.rm = T)) %>% 
  arrange(pLI_0.9)

p_pli_prop <- NB_egenes_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.01, 0.1, 0.3, 0.5, 0.9, 1))) %>% 
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, pLI_bins, name = "bins_ct") %>% 
  distinct(celltype, pLI_bins, total_egene, bins_ct) %>% 
  mutate(prop = bins_ct/total_egene) %>% 
  filter(!is.na(pLI_bins)) %>% 
  ggplot(aes(x = pLI_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols) + axis_theme
p_pli_prop

nb_cs_ct_s %>%
  group_by(phenotype_id, pLI, loeuf, RVIS, RVIS_percentile, EDS) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = mean_cs, y = EDS)) + geom_point() + geom_smooth()

  
#### LOEUF #### 
# majority is protein-coding genes

NB_egenes_df %>% distinct(phenotype_id, loeuf) %>% 
  left_join(gtf, by = "phenotype_id") %>% 
  summarize(loeuf_0.6 = sum(loeuf > 0.6, na.rm = T)/n_protein_coding,
            mean_loeuf = mean(loeuf, na.rm = T))

# gene expression
p_loeuf_express <- NB_all_df_s %>% 
  mutate(loeuf_bins = cut(loeuf, breaks=c(0, 0.3, 0.6, 0.9, 2))) %>% 
  filter(!is.na(loeuf_bins)) %>% 
  group_by(phenotype_id, loeuf_bins) %>% 
  mutate(mean_express = mean(express_percent, na.rm = T)) %>% ungroup() %>% 
  ggplot(aes(x = loeuf_bins, y = express_percent)) + geom_boxplot() + 
  ylab("coverage") + axis_theme

p_loeuf_all <- NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>%
  distinct(phenotype_id, category, loeuf) %>% 
  ggplot(aes(x = category, y = loeuf)) + geom_boxplot() + axis_theme
p_loeuf_all


# loeuf ~ eGene (at least one cell type): P=1.45e-51, slope=0.115
lm(loeuf ~ category + mean_express, 
   data=NB_all_df_s %>%
     group_by(phenotype_id) %>%
     mutate(category = ifelse(sum(is_egene) > 0, 1, 0),
            mean_express = mean(express_percent, na.rm = T)) %>%
     distinct(phenotype_id, mean_express, category, loeuf)) %>% 
  broom::tidy()

# LOEUF by cell type
p_loeuf_cell <- NB_all_df_s %>% 
  filter(is_egene > 0) %>% 
  mutate(celltype = fct_relevel(as.factor(celltype), celltype_order)) %>% 
  ggplot(aes(x = celltype, y = loeuf, fill = celltype)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = celltype_cols)+
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1), legend.position = "none")
p_loeuf_cell

p_loeuf_cs <- nb_cs_ct_s %>%
  group_by(phenotype_id, loeuf) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = loeuf, y = mean_cs)) +
  geom_point(size=0.6) +
  geom_smooth(method="lm") + 
  axis_theme + xlab("LOEUF")+
  ylab("average number of CS per gene")
p_loeuf_cs

# mean CS ~ loeuf, p = 3.77E-5, pos slope
lm(mean_cs ~ loeuf, data = nb_cs_ct_s %>% 
     group_by(phenotype_id, loeuf) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

loeuf_grid <- plot_grid(
  p_loeuf_cell, p_loeuf_all,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 0.3)
)
loeuf_grid


p_loeuf_prop <- NB_egenes_df_s %>% 
  mutate(loeuf_bins = cut(loeuf, breaks=c(0, 0.3, 0.6, 0.9, 2))) %>% 
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, loeuf_bins, name = "bins_ct") %>% 
  distinct(celltype, loeuf_bins, total_egene, bins_ct) %>% 
  mutate(prop = bins_ct/total_egene) %>% 
  filter(!is.na(loeuf_bins)) %>% 
  ggplot(aes(x = loeuf_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols) + axis_theme
p_loeuf_prop

#### EDS #### 
# not only protein-coding genes

NB_all_df_s %>% filter(!is.na(loeuf)) %>% pull(gene_type) %>% table()

NB_egenes_df %>% distinct(phenotype_id, EDS) %>% 
  left_join(gtf, by = "phenotype_id") %>% 
  summarize(EDS = sum(EDS, na.rm = T)/n_protein_coding,
            mean_EDS = mean(EDS, na.rm = T),
            med_EDS = median(EDS, na.rm = T))

# gene expression

p_eds_express <- NB_all_df_s %>% 
  mutate(eds_bins = cut(EDS, breaks=c(0, 0.25, 0.50, 0.75, 1))) %>% 
  filter(!is.na(eds_bins)) %>% 
  group_by(phenotype_id, eds_bins) %>% 
  mutate(mean_express = mean(express_percent, na.rm = T)) %>% ungroup() %>% 
  ggplot(aes(x = eds_bins, y = express_percent)) + geom_boxplot() + 
  ylab("coverage") + axis_theme

p_eds_all <- NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>%
  distinct(phenotype_id, category, EDS) %>% 
  ggplot(aes(x = category, y = EDS)) + geom_boxplot() + axis_theme
p_eds_all

# EDS ~ eGene (at least one cell type): P=5.69E-6, slope= -0.0109
lm(EDS ~ category + mean_express, 
   data=NB_all_df_s %>%
     group_by(phenotype_id) %>%
     mutate(category = ifelse(sum(is_egene) > 0, 1, 0),
            mean_express = mean(express_percent, na.rm = T)) %>%
     distinct(phenotype_id, mean_express, category, EDS)) %>% 
  broom::tidy()

# EDS by celltype
p_eds_cell <- NB_all_df_s %>% 
  filter(is_egene > 0) %>% 
  mutate(celltype = fct_relevel(as.factor(celltype), celltype_order)) %>% 
  ggplot(aes(x = celltype, y = EDS, fill = celltype)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = celltype_cols)+
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1), legend.position = "none")
p_eds_cell

p_eds_cs <- nb_cs_ct_s %>%
  group_by(phenotype_id, EDS) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = EDS, y = mean_cs)) +
  geom_point(size=0.6) +
  geom_smooth(method="lm") + 
  axis_theme + xlab("EDS")+
  ylab("average number of CS per gene")
p_eds_cs

# mean CS ~ loeuf, p = 0.436, neg slope
lm(mean_cs ~ EDS, data = nb_cs_ct_s %>% 
     group_by(phenotype_id, EDS) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

eds_grid <- plot_grid(
  p_eds_cell, p_eds_all,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 0.3)
)
eds_grid


p_eds_prop <- NB_egenes_df_s %>% 
  mutate(eds_bins = cut(EDS, breaks=c(0, 0.25, 0.50, 0.75, 1))) %>% 
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, eds_bins, name = "bins_ct") %>% 
  distinct(celltype, eds_bins, total_egene, bins_ct) %>% 
  mutate(prop = bins_ct/total_egene) %>% 
  filter(!is.na(eds_bins)) %>% 
  ggplot(aes(x = eds_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols) + axis_theme
p_eds_prop

#### RVIS #### 
# majority is protein-coding genes
# gene expression
p_rvis_express <- NB_all_df_s %>% 
  mutate(rvis_percentile_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100)),
         rvis_bins = cut(RVIS, breaks=c(-5, -2, 0, 0.5, 30))) %>% 
  filter(!is.na(rvis_percentile_bins)) %>% 
  group_by(phenotype_id, rvis_percentile_bins) %>% 
  mutate(mean_express = mean(express_percent, na.rm = T)) %>% ungroup() %>% 
  ggplot(aes(x = rvis_percentile_bins, y = express_percent)) + geom_boxplot() + 
  ylab("coverage") + axis_theme

p_rvis_all <- NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% ungroup() %>% 
  distinct(phenotype_id, category, RVIS_percentile) %>% 
  ggplot(aes(x = category, y = RVIS_percentile)) + geom_boxplot() + axis_theme
p_rvis_all

NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% ungroup() %>% 
  distinct(phenotype_id, category, RVIS, RVIS_percentile) %>% 
  group_by(category) %>% summarise(mean = mean(RVIS, na.rm = T),
                                   mean_p = mean(RVIS_percentile, na.rm = T))

# EDS ~ eGene (at least one cell type): P=1.12E-18, slope= 0.51
lm(RVIS ~ category + mean_express, 
   data=NB_all_df_s %>%
     group_by(phenotype_id) %>%
     mutate(category = ifelse(sum(is_egene) > 0, 1, 0),
            mean_express = mean(express_percent, na.rm = T)) %>%
     distinct(phenotype_id, mean_express, category, RVIS)) %>% 
  broom::tidy()

p_rvis_cell <- NB_all_df_s %>% 
  filter(is_egene > 0) %>% 
  mutate(celltype = fct_relevel(as.factor(celltype), celltype_order)) %>% 
  ggplot(aes(x = celltype, y = RVIS_percentile, fill = celltype)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = celltype_cols)+
  axis_theme+
  theme(axis.text.x = element_text(angle=35, hjust = 1), legend.position = "none")
p_rvis_cell

# eGenes are more tolerant to functional variation
NB_all_df_s %>% 
  filter(is_egene > 0) %>% 
  distinct(celltype, phenotype_id, RVIS, RVIS_percentile) %>% 
  group_by(celltype) %>% summarise(mean = mean(RVIS, na.rm = T),
                                   mean_p = mean(RVIS_percentile > 75, na.rm = T),
                                   mean_p2 = mean(RVIS_percentile < 25, na.rm = T))

# RVIS and their percentile
NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% 
  ungroup() %>% 
  filter(category > 0) %>%
  distinct(category, phenotype_id, RVIS, RVIS_percentile) %>% 
  filter(!is.na(RVIS)) %>% 
  group_by(category) %>% 
  summarize(mean_rvis = mean(RVIS),
            mean_rvis_p = mean(RVIS_percentile))

# immune disorders have greater tolerance to functional variation
NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, category, RVIS, RVIS_percentile) %>% 
  filter(!is.na(RVIS_percentile)) %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100), labels = c(25, 50, 75, 100))) %>% 
  group_by(category, rvis_bins) %>% summarise(n = n()) %>% ungroup() %>% 
  group_by(category) %>% mutate(n_genes = sum(n)) %>% ungroup() %>% 
  ggplot(aes(x = rvis_bins, y = n/n_genes)) + 
  geom_bar(stat="identity", position = "dodge2", aes(fill = category))

# group by cell types
NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% 
  ungroup() %>% 
  filter(is_egene > 0) %>% 
  filter(!is.na(RVIS_percentile)) %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100), labels = c(25, 50, 75, 100))) %>% 
  group_by(celltype, rvis_bins) %>% summarise(n = n()) %>% ungroup() %>% 
  group_by(celltype) %>% mutate(n_genes = sum(n)) %>% ungroup() %>% View
  ggplot(aes(x = rvis_bins, y = n/n_genes)) + 
  geom_bar(stat="identity", position = "dodge2", aes(fill = celltype)) +
  facet_wrap(~celltype)

NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, is_egene, category, loeuf, RVIS, RVIS_percentile, pLI, EDS) %>% 
  filter(!is.na(RVIS_percentile)) %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100), labels = c(25, 50, 75, 100))) %>% 
  filter(is_egene > 0) %>% #group_by(rvis_bins) %>% summarize(mean = mean(loeuf, na.rm = T))
  ggplot(aes(x = rvis_bins, y = EDS))+geom_boxplot()

# count egenes based on bins of RVIS score
NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, "eGene", "not eGene")) %>% 
  ungroup() %>% 
  distinct(phenotype_id, is_egene, category, loeuf, RVIS, RVIS_percentile, pLI, EDS) %>% 
  filter(!is.na(RVIS_percentile)) %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100), labels = c(25, 50, 75, 100))) %>% 
  group_by(rvis_bins) %>% summarize(prop = mean(is_egene),
                                    n = n(),
                                    sd = sqrt(prop * (1-prop)/ n))
  ggplot(aes(x = rvis_bins, y = loeuf))+geom_boxplot()

plot(rvis_score$RVIS,rvis_score$RVIS_percentile)
p_rvis_cs <- nb_cs_ct_s %>%
  group_by(phenotype_id, RVIS_percentile) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = RVIS_percentile, y = mean_cs)) +
  geom_point(size=0.6) +
  geom_smooth(method="lm") + 
  axis_theme + xlab("RVIS percentile")+
  ylab("average number of CS per gene")
p_rvis_cs

# mean CS ~ RVIS, p = 2.05E-5, pos slope
lm(mean_cs ~ RVIS_percentile, data = nb_cs_ct_s %>% 
     group_by(phenotype_id, RVIS_percentile) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

p_rvis_prop <- NB_egenes_df_s %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100))) %>%
  # mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 50, 100))) %>% 
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, rvis_bins, name = "bins_ct") %>% 
  distinct(celltype, rvis_bins, total_egene, bins_ct) %>% 
  mutate(prop = bins_ct/total_egene) %>% 
  filter(!is.na(rvis_bins)) %>% 
  ggplot(aes(x = rvis_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols) + axis_theme +
  theme(legend.key.size =unit(0.2, "cm"))

p_rvis_prop

non_egenes <- NB_all_df_s %>% group_by(phenotype_id) %>% 
  mutate(egene = sum(is_egene > 0)) %>% ungroup() %>% 
  filter(egene < 1) %>% distinct(phenotype_id, RVIS_percentile) %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100)),
         total_egene = n(),
         celltype = "non-eGenes") %>% 
  add_count(rvis_bins, name = "bins_ct") %>% 
  distinct(celltype, rvis_bins, total_egene, bins_ct) 

NB_egenes_df_s %>% 
  mutate(rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100))) %>%
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, rvis_bins, name = "bins_ct") %>% 
  distinct(celltype, rvis_bins, total_egene, bins_ct) %>% 
  bind_rows(non_egenes) %>% 
  mutate(prop = bins_ct/total_egene,
         celltype = fct_relevel(celltype, celltype_order)) %>% 
  filter(!is.na(rvis_bins)) %>% 
  ggplot(aes(x = rvis_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols) +
  axis_theme +
  theme(legend.key.size =unit(0.2, "cm"))

ggsave2(paste0(plots_dir, "NB_egenes_RVIS_prop.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


rvis_grid <- plot_grid(
  p_rvis_cell, p_rvis_all,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 0.3)
)
rvis_grid

plot_grid(
  pli_grid, loeuf_grid, eds_grid, rvis_grid,
  align = 'h',
  hjust = -1, # -1
  nrow = 4,
  labels = c("A", "B", "C", "D"), label_size = 10
)

ggsave2(paste0(plots_dir, "NB_egenes_selection.png"),
        width = full_width, height = full_height*2.5, units = 'cm', dpi = 300)

plot_grid(
  p_pli_cs, p_loeuf_cs, p_eds_cs, p_rvis_cs,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  labels = c("A", "B", "C", "D"), label_size = 10
)

ggsave2(paste0(plots_dir, "NB_egenes_selection_cs.png"),
        width = full_width, height = full_height/1.3, units = 'cm', dpi = 300)

# grid of gene expression and each metric

plot_grid(
  p_pli_express, p_loeuf_express, p_rvis_express, p_eds_express,
  align = 'h',
  hjust = -1, # -1
  nrow = 2,
  labels = c("A", "B", "C", "D"), label_size = 10
)

ggsave2(paste0(plots_dir, "gene_coverage_selection.png"),
        width = full_width, height = full_height*1.2, units = 'cm', dpi = 300)

# metric and expression
# increased gene expression with increasing pLI, 
NB_all_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.2, 0.5, 0.7, 0.9, 1)),
         rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100)),
         loeuf_bins = cut(loeuf, breaks=c(0, 0.2, 0.5, 0.9, 2)),
         eds_bins = cut(EDS, breaks=c(0, 0.25, 0.50, 0.75, 1))) %>% 
  filter(!is.na(eds_bins)) %>% 
  # group_by(rvis_bins) %>% 
  # summarize(mean_expression = mean(express_percent),
  #                                   mean_rate = mean(rate_mean),
  #           n=n()) %>% ungroup() %>%
  ggplot(aes(x = eds_bins, y = express_percent)) + geom_boxplot() + geom_smooth()

NB_all_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.2, 0.5, 0.7, 0.9, 1)),
         rvis_bins = cut(RVIS_percentile, breaks=c(0, 25, 50, 75, 100)),
         loeuf_bins = cut(loeuf, breaks=c(0, 0.2, 0.5, 0.9, 2)),
         eds_bins = cut(EDS, breaks=c(0, 0.25, 0.50, 0.75, 1))) %>% 
  distinct(phenotype_id, RVIS_percentile, EDS, rvis_bins, eds_bins) %>% 
  filter(!is.na(eds_bins) & !is.na(rvis_bins)) %>% 
  ggplot(aes(x = eds_bins, y = RVIS_percentile)) + geom_boxplot() + geom_smooth()

NB_all_df_s %>% 
  distinct(phenotype_id, EDS, RVIS) %>% 
  ggplot(aes(x = EDS, y = RVIS)) + geom_point() + geom_smooth(method='lm')

#### test eGenes/non-eGene vs. metric ####
# since gene expression is correlated with these metrics
# test gene expression ~ metric
# pLI: 7.65e-125, slope = 0.166
# loeuf: 1.76e-86, slope = -0.119
# EDS: 9.07E-7, slope = -0.0883
# rvis: 4.13E-41, slope = -0.0381

lm(mean_express ~ EDS, 
    data=NB_all_df_s %>% 
      group_by(phenotype_id) %>% 
      mutate(mean_express = mean(express_percent, na.rm = T)) %>% 
      ungroup() %>% 
      distinct(phenotype_id, mean_express, pLI, loeuf, EDS, RVIS, RVIS_percentile)) %>% 
  broom::tidy()

# U shape EDS: p = 9.93E-70
# U shape rvis (%): p = 1.67E-145
fit <- lm(mean_express ~ RVIS_percentile + RVIS_percentile2, 
   data=NB_all_df_s %>% 
     group_by(phenotype_id) %>% 
     mutate(mean_express = mean(express_percent, na.rm = T)) %>% 
     ungroup() %>% 
     distinct(phenotype_id, mean_express, pLI, loeuf, EDS, RVIS, RVIS_percentile) %>% 
     mutate(EDS2 = EDS^2,
            RVIS_percentile2 = RVIS_percentile^2)) %>% summary()
fit
pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail=FALSE)

# pLI: 3.26E-49, slope = -0.762
# loeuf: 3.46e-50, slope = 0.638
# EDS: 4.64e-6, slope = -0.575
# rvis: 3.30e-18, slope = 0.184

glm(category ~ RVIS_percentile + mean_express, 
     family=binomial(link='logit'),
  data=NB_all_df_s %>% 
  group_by(phenotype_id) %>% 
  mutate(category = ifelse(sum(is_egene) > 0, 1, 0),
         mean_express = mean(express_percent, na.rm = T)) %>% 
  ungroup() %>% 
  distinct(phenotype_id, mean_express, category, pLI, loeuf, EDS, RVIS, RVIS_percentile)) %>% 
  broom::tidy()

NB_all_df_s %>% filter(celltype == "CD4_SOX4" & is_egene > 0) %>% 
  filter(!is.na(pLI)) %>% arrange(desc(pLI)) %>% View

NB_all_df_s %>% filter(is_egene > 0 & grepl("RPS", `HGNC gene`))


#### test eGenes/non-eGene vs. metric within each cell type ####

# logit: eGene/non-eGene ~ metric + GE
logit_res <- tibble(expand.grid(celltypes, c("pLI", "loeuf", "EDS", "RVIS", "RVIS_percentile"))) %>% 
  rename(celltype = Var1, metric = Var2) %>% 
  mutate(Zstat = NA, pval = NA, slope = NA, se = NA, se_bs = NA)
n_bs <- 1000

for (i in 1:nrow(logit_res)){
  cell <- logit_res$celltype[i]
  metric <- logit_res$metric[i]
  
  formula <- as.formula(paste0("is_egene ~ ", metric, " + express_percent"))
  cell_data <- NB_all_df_s %>% filter(celltype == cell)
  
  fit <- glm(formula, 
      family=binomial(link='logit'),
      data=cell_data) %>% 
    broom::tidy()
  logit_res$Zstat[i] <- fit$statistic[2]
  logit_res$pval[i] <- fit$p.value[2]
  logit_res$slope[i] <- fit$estimate[2]
  logit_res$se[i] <- fit$`std.error`[2]
  
  # bootstrap SE
  set.seed(i)
  slope_bt <- c()
  for (j in 1:n_bs){
    indices <- sample(1:nrow(cell_data), replace=TRUE)
    cell_data_bt <- cell_data[indices, ]
    fit <- glm(formula, 
               family=binomial(link='logit'),
               data=cell_data_bt) %>% 
      broom::tidy()
    slope_bt <- c(slope_bt, fit$estimate[2])
  }
  logit_res$se_bs[i] <- sd(slope_bt)
  print(i)
}

logit_res %>% ggplot(aes(x = se, y = se_bs)) + geom_point() +
  geom_abline() + axis_theme

plot(logit_res$se, logit_res$se_bs)

# linear result: metric ~ eGene + GE

lm_res <- tibble(expand.grid(celltypes, c("pLI", "loeuf", "EDS", "RVIS", "RVIS_percentile"))) %>% 
  rename(celltype = Var1, metric = Var2) %>% 
  mutate(Zstat = NA, pval = NA, slope = NA, se = NA, se_bs = NA)
n_bs <- 1000

for (i in 1:nrow(lm_res)){
  cell <- lm_res$celltype[i]
  metric <- lm_res$metric[i]
  
  formula <- as.formula(paste0( metric, " ~ is_egene + express_percent"))
  cell_data <- NB_all_df_s %>% filter(celltype == cell)
  
  fit <- lm(formula, data=cell_data) %>% 
    broom::tidy()
  lm_res$Zstat[i] <- fit$statistic[2]
  lm_res$pval[i] <- fit$p.value[2]
  lm_res$slope[i] <- fit$estimate[2]
  lm_res$se[i] <- fit$`std.error`[2]
  
  # bootstrap SE
  set.seed(i)
  slope_bt <- c()
  for (j in 1:n_bs){
    indices <- sample(1:nrow(cell_data), replace=TRUE)
    cell_data_bt <- cell_data[indices, ]
    fit <- lm(formula, data=cell_data_bt) %>%
      broom::tidy()
    slope_bt <- c(slope_bt, fit$estimate[2])
  }
  lm_res$se_bs[i] <- sd(slope_bt)
  print(i)
}

lm_res %>% ggplot(aes(x = se, y = se_bs)) + geom_point() +
  geom_abline() + axis_theme

# plot
celltype_order_pli <- logit_res %>% 
  mutate(celltype = as.character(celltype)) %>% 
  filter(metric == "pLI") %>% arrange(slope) %>% pull(celltype)

logit_res %>% 
  filter(metric != "RVIS_percentile") %>% 
  mutate(metric = case_when(metric == "loeuf" ~ "LOEUF",
                            .default = metric)) %>% 
  mutate(celltype = fct_relevel(celltype, rev(celltype_order)),
         metric = fct_relevel(as.factor(metric), c("pLI", "LOEUF", "RVIS", "EDS"))) %>% 
  ggplot(aes(y = celltype, x = slope)) + 
  geom_col(aes(fill = celltype)) +
  geom_errorbar(aes(xmin = slope-1.96*se, xmax = slope + 1.96*se),
                width = 0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = celltype_cols) +
  facet_grid(~metric, scales = "free_x") +
  # facet_grid(~metric) + 
  axis_theme + theme(legend.position = "none")
  
ggsave2(paste0(plots_dir, "logit_egene_celltype.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

ggsave2(paste0(plots_dir, "lm_metric_egene_celltype.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

#### effect size and AF by CT ####
NB_egenes_df %>%
  mutate(slope = slope_wald) %>%
# nb_cs_pip0.5_wald %>%
#   inner_join(nb_cs_pip0.5 %>% filter(celltype != "allcells") %>% distinct(celltype, phenotype_id, snp)) %>%
#   left_join(pLI_score, by = "phenotype_id") %>%
  filter(!is.na(pLI)) %>%
  mutate(af = ifelse(slope > 0, af, 1-af),
         pLI_high = ifelse(pLI >= 0.9, TRUE, FALSE)) %>% 
  ggplot(aes(x = af, y = abs(slope), group=celltype)) +
  geom_point(size=0.3) + facet_wrap(.~celltype, scales = "free_y", nrow = 4) + 
  xlab("(expression increasing) allele frequency")+
  geom_smooth(method = "loess") + scale_color_manual(values = c(`TRUE`="orange", `FALSE`="grey1", `NA`="grey"))

ggsave2(paste0(plots_dir, "effectsize_af_leadsnp.png"),
        width = full_width, height = full_height*2, units = 'cm', dpi = 300)

# eqtl effect density plot
NB_egenes_df %>%
  mutate(slope = slope_wald) %>%
  # nb_cs_pip0.5_wald %>%
  #   inner_join(nb_cs_pip0.5 %>% filter(celltype != "allcells") %>% distinct(celltype, phenotype_id, snp)) %>%
  #   left_join(pLI_score, by = "phenotype_id") %>%
  filter(!is.na(pLI)) %>%
  mutate(pLI_high = ifelse(pLI >= 0.9, TRUE, FALSE)) %>% 
  ggplot(aes(x = abs(slope), fill = pLI_high, color = pLI_high)) +
  geom_density(alpha=0.7) + facet_wrap(.~celltype, scales = "free_y", nrow = 4) + 
  scale_fill_manual(values = c(`TRUE`="orange", `FALSE`="grey1", `NA`="grey")) +
  scale_color_manual(values = c(`TRUE`="orange", `FALSE`="grey1", `NA`="grey"))

ggsave2(paste0(plots_dir, "effectsize_hist_pLI.png"),
        width = full_width, height = full_height*2, units = 'cm', dpi = 300)

# effect size of highly constrained genes are smaller than other genes
tmp <-  
  NB_egenes_df %>%
    mutate(slope = slope_wald) %>%
  # nb_cs_pip0.5_wald %>%
  # inner_join(nb_cs_pip0.5 %>% filter(celltype != "allcells") %>% distinct(celltype, phenotype_id, snp)) %>%
  # left_join(pLI_score, by = "phenotype_id") %>%
  filter(!is.na(pLI)) %>%
  mutate(pLI_high = ifelse(pLI >= 0.9, TRUE, FALSE))

lm(abs(slope) ~ af, data=tmp) %>% broom::tidy()

wilcox.test(abs(tmp$slope[tmp$pLI_high == TRUE]), abs(tmp$slope[tmp$pLI_high == FALSE]), "t") %>% broom::tidy()

nb_cs_pip0.5_wald %>%
  left_join(pLI_score, by = "phenotype_id") %>% 
  filter(!is.na(pLI)) %>% 
  group_by(pLI >= 0.9) %>% summarize(mean = mean(abs(slope)))

NB_egenes_df %>% 
  filter(!is.na(pLI)) %>% 
  group_by(pLI >= 0.9) %>% summarize(mean = mean(abs(slope)))

#### proportion of eGenes in each bins ####

NB_egenes_df_s %>% 
  mutate(pLI_bins = cut(pLI, breaks=c(0, 0.01, 0.1, 0.3, 0.5, 0.9, 1))) %>% 
  group_by(celltype) %>% 
  mutate(total_egene = n()) %>% ungroup() %>% 
  add_count(celltype, pLI_bins, name = "bins_ct") %>% 
  distinct(celltype, pLI_bins, total_egene, bins_ct) %>% 
  mutate(prop = bins_ct/total_egene) %>% 
  filter(!is.na(pLI_bins)) %>% 
  ggplot(aes(x = pLI_bins, y = prop, fill = celltype)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = celltype_cols)
  
#### selection between species ####
# entire CI > 0 -> rapidly evolving due to positive darwinian selection
# entrie CI < 0 ->  locus is subject to negative selection
non_neutral <- read_tsv("../data/OneK1K/annotation/Bustamante2005/Bustamante2005_result.txt") %>% 
  janitor::clean_names() %>% 
  mutate(above0 = ifelse(lower_ci_2_5_percent > 0, 1, 0),
         below0 = ifelse(upper_ci_97_5_percent < 0, 1, 0))
names(non_neutral)
View(non_neutral)
summary(non_neutral)

NB_all_df %>% left_join(NB_egenes_df %>% 
                          mutate(is_egene = 1) %>% 
                          select(phenotype_id, celltype, is_egene)) %>% 
  mutate(is_egene = replace_na(is_egene, 0)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type), by=c("phenotype_id")) %>% 
  filter(is_egene>0) %>% 
  left_join(non_neutral, by=c("GeneSymbol"="gene_name")) %>% 
  group_by(celltype) %>% 
  summarize(mean = mean(posterior_mean, na.rm = T),
            mean_pos = mean(above0, na.rm = T),
            mean_neg = mean(below0, na.rm = T))

tmp <- NB_all_df %>% left_join(NB_egenes_df %>% 
                          mutate(is_egene = 1) %>% 
                          select(phenotype_id, celltype, is_egene)) %>% 
  mutate(is_egene = replace_na(is_egene, 0)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type), by=c("phenotype_id")) %>% 
  filter(is_egene>0) %>%
  left_join(non_neutral, by=c("GeneSymbol" = "gene_name")) %>% 
  distinct(celltype, is_egene, GeneSymbol, .keep_all = T)

tmp %>% 
  group_by(celltype, is_egene) %>% 
  summarize(mean = mean(posterior_mean, na.rm = T),
            mean_ub = mean(upper_ci_97_5_percent, na.rm = T),
            mean_lb = mean(lower_ci_2_5_percent, na.rm = T)) %>% View

non_neutral_meta <- data.frame(celltype = celltypes, 
                               meta_OR=NA, meta_lower=NA, meta_upper=NA, meta_p=NA, meta_se=NA)

for (idx in seq_along(non_neutral_meta$celltype)){
  cell <- non_neutral_meta$celltype[idx]
  
  meta_dat <- tmp %>% 
    filter(celltype == cell) %>% 
    mutate(lower = lower_ci_2_5_percent,
           upper = upper_ci_97_5_percent, 
           TE = posterior_mean,
           seTE = NA)
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = phenotype_id,
                       data = meta_dat,
                       sm = "", 
                       method.tau = "PM",
                       fixed = TRUE,
                       random = FALSE,
                       title = "Enrichment (Pre-calculated)")
  
  metares <- summary(m.gen_bin)
  
  non_neutral_meta$meta_OR[[idx]] <- metares$TE.fixed
  non_neutral_meta$meta_lower[[idx]] <- metares$lower.fixed
  non_neutral_meta$meta_upper[[idx]] <- metares$upper.fixed
  non_neutral_meta$meta_p[[idx]] <- metares$pval.fixed
  non_neutral_meta$meta_se[[idx]] <- metares$seTE.fixed 
}

non_neutral_meta

#### positive selection genes ####
# total 360 genes, 186 are found in >1 study
pos_genes <- readxl::read_xls("../data/OneK1K/annotation/Barreiro2009/Supp_Table2.xls", skip = 3) %>% 
  janitor::clean_names() %>% 
  filter(!is.na(gene_symbol)) %>% 
  distinct(gene_symbol) %>% 
  mutate(is_pos = 1) # %>% filter(nb > 1)

pos_genes %>% add_count(gene_symbol) %>% filter(n>1)
pos_genes %>% filter(nb > 1)

tmp <- NB_all_df %>% left_join(NB_egenes_df %>% 
                                 mutate(is_egene = 1) %>% 
                                 select(phenotype_id, celltype, is_egene)) %>% 
  mutate(is_egene = replace_na(is_egene, 0)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type), by=c("phenotype_id")) %>% 
  left_join(pos_genes %>% select(gene_symbol, is_pos), by=c("GeneSymbol" = "gene_symbol")) %>% 
  mutate(is_pos = replace_na(is_pos, 0)) %>% 
  arrange(celltype, GeneSymbol, desc(is_egene)) %>%
  distinct(celltype, GeneSymbol, is_egene, .keep_all = T)


# fisher exact test
pos_summary <- tmp %>% 
  group_by(celltype, is_egene) %>% 
  summarize(mean = mean(is_pos),
            pos = sum(is_pos),
            n = n(),
            not_pos = n - pos) %>% 
  ungroup()

pos_summary
fisher_pos <- tibble(celltype = celltypes, 
                     fisher_p = NA, fisher_OR = NA, fisher_LB = NA, fisher_UB = NA)
for (i in seq_along(celltypes)){
  cell <- fisher_pos$celltype[i]
  
  mat <- matrix(c(pos_summary %>% filter(celltype == cell & is_egene > 0) %>% pull(pos),
           pos_summary %>% filter(celltype == cell & is_egene < 1) %>% pull(pos),
           pos_summary %>% filter(celltype == cell & is_egene > 0) %>% pull(not_pos),
           pos_summary %>% filter(celltype == cell & is_egene < 1) %>% pull(not_pos)),
           nrow=2)
  
  fit <- fisher.test(mat, alternative = "t") %>% broom::tidy()
  
  fisher_pos$fisher_p[i] <- fit$`p.value`
  fisher_pos$fisher_OR[i] <- fit$estimate
  fisher_pos$fisher_LB[i] <- fit$`conf.low`
  fisher_pos$fisher_UB[i] <- fit$`conf.high`
}

fisher_pos

fisher_pos %>% 
  mutate(celltype = fct_reorder(as.factor(celltype), fisher_OR)) %>% 
  ggplot(aes(y = celltype, x = fisher_OR, color = celltype)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = fisher_LB, xmax = fisher_UB),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = celltype_cols)+
  axis_theme + ylab("")+
  xlab("Odds Ratio (95% CI)") + theme(strip.text = element_blank())+
  theme(legend.key.size =unit(0.2, "cm"))


ggsave2(paste0(plots_dir, "egene_pos_select_genes.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)


# all eGenes
tmp %>% distinct(phenotype_id, is_egene, is_pos) %>% 
  group_by(is_egene) %>% 
  summarize(n = n(), pos = sum(is_pos))

mat <- matrix(c(120,298,6934-120, 21002-298),nrow=2)

chisq.test(tmp %>% distinct(phenotype_id, is_egene, is_pos) %>% pull(is_pos), 
           tmp %>% distinct(phenotype_id, is_egene, is_pos) %>% pull(is_egene))

fit <- fisher.test(mat, alternative = "t") %>% broom::tidy()

#### iHS score ####
# both large pos and neg values can be of interest
# how unusual it is to observe the haplotype around each SNP


#### BMM #### 

bmm <- read_tsv("../data/OneK1K/annotation/anc_dna/bmm_v2.7.tsv.gz")

# MYO16,rs9520848, pip
nb_cs %>% 
  filter(celltype != "allcells") %>% 
  left_join(geno_rsid %>% select(snp=variant_id, rsid), by="snp") %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>%
  inner_join(bmm, by=c("rsid"="SNP")) %>% 
  filter(pip >= 0.9) %>% 
  arrange(desc(abs(Z))) %>% 
  #filter(rsid == "rs9520848") %>% 
  View
 # ggplot(aes(x = pip, y = abs(Z))) + geom_point()

nb_cs_pip0.5_wald %>% filter(phenotype_id=="ENSG00000041515" & celltype == "CD4_NC") %>% View

#### primateAI score ####
primate_ai <- read_tsv("../data/OneK1K/annotation/primateAI/PrimateAI_Deep_Learning_VUS_Classification-51955905/135822688/ds.1a659d3c545f4dfb8d251f55a8c0abef/PrimateAI_scores_v0.2.tsv.gz",
                       skip=11)
primate_ai <- primate_ai %>% mutate(snp = paste0(chr, "_", pos, "_", ref, "_", alt, "_b37"))

nb_cs_pip0.5_primate_ai <- nb_cs_pip0.5 %>% 
  left_join(primate_ai %>% select(snp, primateDL_score))

nb_cs_pip0.5_primate_ai %>% filter(!is.na(primateDL_score))

# # Saving on object in RData format
# save(data1, file = "data.RData")
# # Save multiple objects
# save(data1, data2, file = "data.RData")
# # To load the data again
# load("data.RData")