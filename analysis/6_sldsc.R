# look at sldsc results
library(rcartocolor)
# baseline result:  tau and tau Z/P values
# baselineld: proportion of h2 and h2 enrichment

method <- "jaxqtl"
ldsc_indir <- "../result/finemap/sldsc_wald/"
res <- #read_tsv(paste0(ldsc_indir, method, "_sldsc_traits107_union.label.newselection.tsv.gz")) %>% 
  read_tsv(paste0(ldsc_indir, "sldsc_addunion/jaxqtl_sldsc_unioncells_allgtex_onek1k.tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  select(Category, annot, celltype, everything()) %>% 
  # filter(Category == "L2_1") %>% 
  filter(Category == "L2_2") %>% 
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = F),
         trait = ifelse(trait == "UKBiobank_Eczema", "UKBiobank.Eczema", trait)) # one-sided p, P(Z>z)

h2g <- #read_tsv(paste0(ldsc_indir, method, "_sldsc_traits107_h2_baseline.newselection.tsv.gz")) %>% # filter(grepl("UKBiobank_Eczema", filename))
  read_tsv(paste0(ldsc_indir, "sldsc_addunion/jaxqtl_sldsc_unioncells_allgtex_onek1k_h2.tsv.gz")) %>% 
  mutate(filename = gsub(".log$", "", filename),
         filename = gsub("UKBiobank_Eczema", "UKBiobank.Eczema", filename)) %>% 
  separate(filename, into=c("trait", "suffix"), sep=".baseline.") %>% 
  separate(suffix, into=c("annot", "celltype"), sep="\\.")

gtex_res <- read_tsv(glue("{ldsc_indir}/gtex_selection/gtex_brain_sldsc_selection_union_allgtex_onek1k.tsv.gz")) %>% 
  select(Category, annot, celltype, everything()) %>% 
  #filter(Category == "L2_1") %>% 
  filter(Category == "L2_2") %>% 
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = F)) # one-sided p, P(Z>z)

gtex_res_h2g <- read_tsv(glue("{ldsc_indir}/gtex_selection/gtex_brain_sldsc_selection_union_allgtex_onek1k_h2.tsv.gz")) %>% 
  mutate(filename = gsub(".log$", "", filename)) %>% 
  separate(filename, into=c("trait", "suffix"), sep=".baseline.") %>% 
  separate(suffix, into=c("annot", "celltype"), sep="\\.")

# 18 blood and immune related traits
res %>% distinct(trait) %>% View

immune_traits <- c("PASS.Type_1_Diabetes.Chiou2021",
                   # "PASS.Type_2_Diabetes.Xue2018",
                   # "PASS.Covid19_Infection.hg_v7",
                   "PASS.Rheumatoid_Arthritis.Ishigaki2022",
                   "PASS.Multiple_Sclerosis.IMSGC2019",
                   "PASS.Lupus.Bentham2015",
                   "PASS.Gout.Donertas2021",
                   "PASS.Inflammatory_Bowel_Disease.deLange2017",
                   "PASS.Celiac.Dubois2010",
                   "GBMI.Asthma",
                   "PASS.Ulcerative_Colitis.deLange2017",
                   "PASS.Crohns_Disease.deLange2017",
                   "PASS.Primary_Biliary_Cirrhosis.Cordell2015",
                   "UKB_460K.disease_HYPOTHYROIDISM_SELF_REP",
                   "PASS.Alzheimers_Disease.Wightman2021",
                   "UKBiobank.Eczema")
length(immune_traits)

blood_trait <- c(# "UKB_460K.blood_WHITE_COUNT",
                 # "UKB_460K.blood_RED_COUNT",
                 # "PASS.Reticulocyte_Count.Vuckovic2020",
                 # "UKB_460K.blood_RBC_DISTRIB_WIDTH", 
                 "PASS.Monocyte_Percentage.Vuckovic2020",
                 # "PASS.Eosinophil_Percentage.Vuckovic2020",
                 # "PASS.CHIP.Kessler2022",
                 # "PASS.Platelet_Distribution_Width.Vuckovic2020",
                 # "PASS.Platelet_Count.Vuckovic2020",
                 "PASS.Lymphocyte_Percentage.Vuckovic2020")

cancer_trait <- c("PANUKBB.Breast_Cancer_female",
                  # "PASS.Cancer_Esophageal.Gharahkhani2016",
                  "PASS.Cancer_Lung.McKay2017",
                  "PASS.Cancer_Prostate.Schumacher2018")

anc_dna <- c("anc_dna_2.7", "bmm_v2.7")
blood_immune <- c(blood_trait, immune_traits, anc_dna); length(blood_immune)
to_plot <- c("Lymphocyte_Percentage", 
             "Monocyte_Percentage",
             "Rheumatoid_Arthritis",
             "Inflammatory_Bowel_Disease",
             "Primary_Biliary_Cirrhosis",
             "Crohns_Disease",
             "Asthma")

trait_col <- c("immune" = "#009E73", "blood" = "#D55E00", "cancer" = "black")

# results adjusting for baselineLD annotation (for h2 and h2 enrichment)
# note: enrichment_p is not based on enrichment/enrichment_se
res_baselineld <- res %>% filter(grepl("baselineLD", filename))

# results adjusting for baseline annotation (for tau and tau Z/P)
# cell type-specific analysis use this tau_Z score to call significance
# calculate tau*
M <- 5961159
res_baseline <- res %>% filter(!grepl("baselineLD", filename)) %>% 
  left_join(h2g, by=c("trait", "annot", "celltype")) %>% 
  mutate(Coefficient_star = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient / h2g,
         Coefficient_star_std_error = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient_std_error / h2g,
         `Coefficient_star_z-score` = Coefficient_star / `Coefficient_star_std_error`)

res_clean <- res_baseline %>% 
  select(-c(starts_with("Prop"), starts_with("Enrichment"))) %>% 
  left_join(res_baselineld %>% select(c(starts_with("Prop"), starts_with("Enrichment"), trait, celltype)),
            by=c("trait", "celltype")) %>% 
  filter(trait %in% blood_immune) %>% 
  mutate(trait_group = ifelse(trait %in% blood_trait, "blood",
                              ifelse(trait %in% cancer_trait, "cancer", "immune")),
         trait_group = factor(trait_group, levels=c("blood", "immune", "cancer"))); dim(res_clean)

length(unique(res_clean$trait)) # 16 traits

##### Compare union cells and allcells  #####

# pertange of SNP
res_clean %>% filter(celltype %in% c("unioncells", "allcells")) %>% 
  distinct(celltype, `Prop._SNPs`)

# h2 and h2 enrichment
tmp <- res_clean %>% 
  filter(celltype == "allcells") %>% 
  filter(!trait %in% c("anc_dna", "bmm_v2.7")) %>% 
  rename(`Prop._h2_bulk` = `Prop._h2`,
         `Prop._SNPs_bulk`= `Prop._SNPs`,
         `Prop._h2_std_error_bulk` = `Prop._h2_std_error`,
         Enrichment_bulk = Enrichment,
         Enrichment_std_error_bulk = Enrichment_std_error,
         Enrichment_p_bulk = Enrichment_p) %>% 
  select(-c(celltype)) %>% 
  left_join(res_clean %>% 
              filter(celltype == "unioncells") %>% 
              select(trait, 
                     `Prop._h2_sc_union` = `Prop._h2`,
                     `Prop._SNPs_sc_union` = `Prop._SNPs`,
                     `Prop._h2_std_error_sc_union` = `Prop._h2_std_error`,
                     Enrichment_sc_union = Enrichment,
                     Enrichment_std_error_sc = Enrichment_std_error,
                     Enrichment_p_sc_union = Enrichment_p),
            by="trait") %>% 
  mutate(p_h2 = NA, h2_diff = NA)

wilcox.test(tmp$Prop._h2_sc_union,
            tmp$Prop._h2_bulk, "t")
broom::tidy(wilcox.test(tmp$Enrichment_sc_union,
            tmp$Enrichment_bulk, "t"))
wilcox.test(-log10(tmp$Enrichment_p_sc_union),
            -log10(tmp$Enrichment_p_bulk), "t")

meta_dat <- tmp %>% 
  mutate(TE = `Prop._h2_sc_union`, 
         lower = NA,
         upper = NA,
         seTE = `Prop._h2_std_error_sc_union`)

meta_dat <- tmp %>% 
  mutate(TE = `Prop._h2_bulk`, 
         lower = NA,
         upper = NA,
         seTE = `Prop._h2_std_error_bulk`)

# enrichment
meta_dat <- tmp %>% 
  mutate(TE = `Enrichment_bulk`, 
         lower = NA,
         upper = NA,
         seTE = `Enrichment_std_error_bulk`)
meta_dat <- tmp %>% 
  mutate(TE = `Enrichment_sc_union`, 
         lower = NA,
         upper = NA,
         seTE = `Enrichment_std_error_sc`)

m.gen_bin <- metagen(TE = TE,
                     seTE = seTE,
                     lower = lower,
                     upper = upper,
                     studlab = trait,
                     data = meta_dat,
                     sm = "",
                     # method.tau = "PM",
                     fixed = TRUE,
                     random = FALSE,
                     title = "Enrichment (Pre-calculated)")

summary(m.gen_bin)
m.gen_bin$seTE.fixed

# h2
# single cell eqtl: 0.0988 [0.0900; 0.1077];
# bulk eqtl:0.0609 [0.0533; 0.0685]

# enrichment
# sc_union: 2.8160 [2.5644; 3.0676]
# bulk: 2.7561 [2.4123; 3.0998]

# identify significantly different pairs
View(tmp)

for (idx in seq_along(tmp$trait)){
  which_trait <- tmp$trait[[idx]]; which_trait
  h2_diff <- tmp$Prop._h2_sc_union[[idx]] - tmp$Prop._h2_bulk[[idx]]
  h2_var <- (tmp$Prop._h2_std_error_sc_union[[idx]])^2 + (tmp$Prop._h2_std_error_bulk[[idx]])^2
  p <- pnorm(abs(h2_diff/sqrt(h2_var)), lower.tail = F) * 2
  tmp$p_h2[[idx]] <- p
  tmp$h2_diff[[idx]] <- h2_diff
}

tmp %>% mutate(percent_increase = (Prop._h2_sc_union - Prop._h2_bulk)/Prop._h2_bulk) %>% 
  arrange(desc(percent_increase)) %>% View
tmp %>% filter(p_h2 < 0.05/n())

tmp %>% filter(`Prop._h2_bulk` > `Prop._h2_sc_union`) %>% View

# only type I diabetes has h2_sc_union < h2_bulk
p_h2 <- tmp %>% filter(trait != "anc_dna_2.7") %>% 
  mutate(sig = ifelse(p_h2 < 0.05/n(), "yes", "no")) %>% 
  ggplot(aes(x = `Prop._h2_bulk`, y = `Prop._h2_sc_union`,color = trait_group)) +
  geom_point(size=0.8) +
  geom_errorbar(aes(ymin = `Prop._h2_sc_union` - 1 * `Prop._h2_std_error_sc_union`,
                    ymax = `Prop._h2_sc_union` + 1 * `Prop._h2_std_error_sc_union`),
                linewidth = 0.2) +
  geom_errorbarh(aes(xmin = `Prop._h2_bulk` - 1 * `Prop._h2_std_error_bulk`,
                     xmax = `Prop._h2_bulk` + 1 * `Prop._h2_std_error_bulk`),
                 linewidth = 0.2) +
  scale_color_manual(values = trait_col) +
  # scale_shape_manual(values = c("yes" = 17, "no" = 16)) +
  geom_abline(slope = 1, intercept = 0, color = "grey") + axis_theme +
  #ylim(0, 0.15)+xlim(0, 0.15)+
  ylim(0, 0.18)+xlim(0, 0.18)+
  theme(legend.position = "bottom")+
  ggtitle("Proportion of heritabily explained")+
  guides(color = guide_legend(title="trait group"))+
  ylab(bquote({h^2}~"(sc-eQTL_union)")) + xlab(bquote({h^2}~"(bulk-eQTL)"))+
  theme(axis.text.x = element_text(hjust = 1),
        legend.position = c(1, .45),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.key.size =unit(0.3, "cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 8)
  )

p_h2

ggsave2(paste0(plots_dir, "sldsc_107traits_h2_union.vs.bulk.png"),
        width = onecol_width/1.2, height = onecol_height*1.2, units = 'cm', dpi = 300)

# enrichment
p1 <- tmp %>% 
  ggplot(aes(x = Enrichment_bulk, y = Enrichment_sc_union, color = trait_group)) +
  geom_point(size=0.8) +
  scale_color_manual(values = trait_col) +
  geom_abline(slope = 1, intercept = 0, color = "grey") + axis_theme +
  xlim(0, 6) + ylim(0, 5)+
  ggtitle("Enrichment for sc-eQTL_union vs. bulk-eQTL")
p1

p2 <- tmp %>%
  ggplot(aes(x = -log10(Enrichment_p_bulk), y = -log10(Enrichment_p_sc_union), color = trait_group)) +
  geom_point(size=0.8) +
  scale_color_manual(values = trait_col) +
  geom_abline(slope = 1, intercept = 0, color = "grey") + axis_theme +
  xlim(0, 6) + ylim(0, 11)+
  # ggtitle("h2(C) for sc_union vs. bulk annotation")
  ggtitle("Enrichment p-values for sc-eQTL_union vs. bulk-eQTL")+
  guides(color = guide_legend(title="trait group"))
p2

plot_grid(
  p1 + theme(legend.position = "none"), p2,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(0.8, 1)
)

ggsave2(paste0(plots_dir, "sldsc_107traits_enrich_union.vs.bulk.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

########## cell type specific analysis ##########
# baseline annotation look at tau and tau_Z
tmp <- res_clean %>% 
  filter(celltype == "allcells") %>% 
  rename(Coefficient_star_bulk = Coefficient_star,
         `Coefficient_star_z-score_bulk` = `Coefficient_star_z-score`) %>% 
  select(-c(celltype)) %>% 
  left_join(res_clean %>% 
              filter(celltype == "unioncells") %>% 
              select(trait, 
                     Coefficient_star_sc_union = Coefficient_star,
                     `Coefficient_star_z-score_sc_union` = `Coefficient_star_z-score`),
            by="trait") 

######## find important cell type per trait ######## 
# pick cell type with highest enrichment for immune traits

trait_celltype <- res_clean %>% 
  filter(!trait %in% c("anc_dna", "bmm_v2.7")) %>% # filter(!is.na(Coefficient_p)) %>% 
  mutate(n_annot = n()) %>%
  mutate(qval = qvalue(Enrichment_p, fdr.level = 0.05, pi0=1)$qvalue) %>% 
  mutate(qval_sig = qval < 0.05,
         bonf_threshold = 0.05/n_annot,
         fdr_p_threshold = max(Enrichment_p[which(qval_sig == TRUE)])) %>% 
  separate(trait, into=c("rm1", "trait", "rm2"), sep="\\.") %>% 
  mutate(trait = gsub("_", " ", trait),
         trait = case_when(trait == "disease HYPOTHYROIDISM SELF REP" ~ "Hypothyroidsim",
                           .default = as.character(trait))) %>% 
  select(-c(rm1, rm2)) # %>% 
  # filter(!celltype %in% c("unioncells"))

trait_celltype %>% 
  rename(p=Coefficient_p) %>% 
  #filter(Enrichment > 0) %>%
  mutate(significant_enrich = ifelse(qval_sig == TRUE, "yes", "no"),
         celltype = case_when(celltype == "allcells" ~ "bulk-eQTL",
                              celltype == "unioncells" ~ "sc-eQTL_union",
                              .default = celltype),
         celltype = fct_relevel(celltype, c("bulk-eQTL", "sc-eQTL_union"))) %>%
  filter(qval_sig == TRUE) %>% 
  #filter(!celltype %in% c("sc-eQTL_union", "bulk-eQTL")) %>% 
  #group_by(trait) %>% filter(qval == min(qval)) %>% 
  select(trait, celltype, Coefficient_star, Coefficient_p=p) %>% 
  write_tsv("./sldsc_res_fig5.tsv")

p_celltype <- trait_celltype %>% 
  rename(P=Coefficient_p) %>% 
  #filter(Enrichment > 0) %>%
  mutate(significant_enrich = ifelse(qval_sig == TRUE, "yes", "no"),
         celltype = case_when(celltype == "allcells" ~ "bulk-eQTL",
                              celltype == "unioncells" ~ "sc-eQTL_union",
                              .default = celltype),
         celltype = fct_relevel(celltype, c("bulk-eQTL", "sc-eQTL_union"))) %>%
  filter(qval_sig == TRUE) %>% 
  mutate(trait = fct_relevel(trait, c("Multiple Sclerosis", "Asthma", "Lupus",
                                      "Eczema",
                                      "Rheumatoid Arthritis", "Hypothyroidsim", "Celiac", "Primary Biliary Cirrhosis",
                                      "Inflammatory Bowel Disease", "Crohns Disease", "Ulcerative Colitis"))) %>% 
  ggplot(aes(x = celltype, y = trait)) + 
  geom_point(aes(color = -log10(P), size = Coefficient_star, shape=significant_enrich)) + 
  scale_color_gradient(high = "red", low = "#fee8c8") +
  scale_shape_manual(values = c("yes" = 16, "no" = 1), guide=FALSE) +
  facet_grid(trait_group~., scales="free_y")+
  ggh4x::force_panelsizes(rows = c(0.3, 1))+
  scale_size_area(max_size = 3) +
  axis_theme + theme(axis.text.x = element_text(angle=35, hjust = 1, size=6),
                     axis.text.y = element_text(size=6),
                     plot.title = element_text(size=6),
                     legend.key.size =unit(0.2, "cm"),
                     legend.text = element_text(size=4),
                     #legend.title = element_text(size=6),
                     legend.position = "bottom") +
  xlab("") + ylab("")+
  scale_x_discrete(
    labels = celltype_bquote
  )+
  guides(size = guide_legend(title=bquote(tau^~"*"))) # +
  #guides(size = guide_legend(title="Enrichment"))

p_celltype

legend <- get_legend(p_celltype)
legend <- cowplot::get_plot_component(p_celltype + theme(legend.text = element_text(size=7)),
                                      'guide-box-bottom', return_all = TRUE)

grid1 <- plot_grid(
  p_h2 + ggtitle(""), 
  p_celltype + theme(legend.position = "none"), 
  align = 'h',
  hjust = -1,
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(1, 1.8), axis="b",
  labels = c("A", "B"), label_size = 10
)

plot_grid(
  grid1, legend,
  align = 'h',
  hjust = -1,
  nrow = 2,
  rel_heights = c(1, 0.1), rel_widths=c(1, 1), axis="b"
)

ggsave2(paste0(plots_dir, "Fig6.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

trait_celltype %>% 
  rename(p=Enrichment_p) %>% 
  # filter(Coefficient >= 0) %>%
  mutate(significant_enrich = ifelse(qval_sig == TRUE, "yes", "no"),
         celltype = case_when(celltype == "allcells" ~ "bulk-eQTL",
                              celltype == "unioncells" ~ "sc-eQTL_union",
                              .default = celltype),
         celltype = fct_relevel(celltype, c("bulk-eQTL", "sc-eQTL_union"))) %>%
  ggplot(aes(x = celltype, y = trait)) + 
  geom_point(aes(color = -log(p), size = Enrichment, shape=significant_enrich)) + 
  scale_color_gradient(high = "red", low = "#fee8c8") +
  scale_shape_manual(values = c("yes" = 16, "no" = 1), guide=FALSE) +
  facet_grid(trait_group~., scales="free_y")+
  ggh4x::force_panelsizes(rows = c(0.3, 1))+
  scale_size_area(max_size = 3) +
  axis_theme + theme(axis.text.x = element_text(angle=35, hjust = 1, size=6),
                     axis.text.y = element_text(size=6),
                     plot.title = element_text(size=6),
                     legend.key.size =unit(0.2, "cm"),
                     legend.text = element_text(size=4),
                     #legend.title = element_text(size=6),
                     legend.position = "bottom") +
  xlab("") + ylab("")+
  guides(color = guide_legend(title=bquote(tau^~"*")))

ggsave2(paste0(plots_dir, "sldsc_enrich.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

table(trait_celltype$celltype)

trait_celltype %>% filter(celltype == "unioncells") %>% pull(Coefficient_star) %>% summary()
trait_celltype %>% filter(celltype == "allcells") %>% pull(Coefficient_star) %>% summary()

# test overall difference of coef across traits: sc_union vs. bulk
wilcox.test(trait_celltype %>% filter(celltype == "unioncells") %>% pull(Coefficient_p),
            trait_celltype %>% filter(celltype == "allcells") %>% pull(Coefficient_p))
# look at tau to pick cell type
# qvalue is computed over all hypothesis tested and look at max p that pass FDR threhsold
res_baseline %>% 
  filter(!celltype %in% c("unioncells", "allcells")) %>% 
  filter(trait %in% immune_traits) %>%
  # group_by(trait) %>% 
  mutate(n_annot = n()) %>%
  # filter()
  mutate(qval = qvalue(Coefficient_p, fdr.level = 0.05)$qvalue) %>% 
  filter(Coefficient_p < 1/n_annot) %>% View
  filter(qval < 0.05) %>% View
  arrange(qval) %>% 
  select(trait, everything()) %>% View

allcat_res <- tibble(Category = unique(tmp$Category), 
                     meta_OR=NA, meta_lower=NA, meta_upper=NA, meta_p=NA)

for (idx in 1:nrow(allcat_res)){
  which_cat <- allcat_res$Category[[idx]]

  meta_dat <- tmp %>% 
    filter(celltype == "allcells" & Category == which_cat) %>% 
    mutate(TE = log(Enrichment), 
           lower = NA,
           upper = NA,
           seTE = Enrichment_std_error)
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = trait,
                       data = meta_dat,
                       sm = "OR",
                       method.tau = "PM",
                       fixed = FALSE,
                       random = TRUE,
                       title = "Enrichment (Pre-calculated)")
  
  metares <- summary(m.gen_bin)
  metares$lower.random
  metares$lower.random
  allcat_res$meta_OR[[idx]] <- exp(metares$TE.random)
  allcat_res$meta_lower[[idx]] <- exp(metares$lower.random)
  allcat_res$meta_upper[[idx]] <- exp(metares$upper.random)
  allcat_res$meta_p[[idx]] <- metares$pval.random
}

allcat_res_allcells <- allcat_res
allcat_res_unioncells <- allcat_res


metares <- summary(m.gen_bin); metares
metares$pval.random

# look at top cell type for fixed annotation
trait_celltype %>% filter(trait == "Eczema") %>% select(annot, celltype, qval) %>% 
  arrange(qval)

trait_celltype

##### ancient DNA analysis  #####
res_clean %>% filter(trait == "anc_dna")

res_anc <- res_baseline %>% 
  select(-c(starts_with("Prop"), starts_with("Enrichment"))) %>% 
  left_join(res_baselineld %>% select(c(starts_with("Prop"), starts_with("Enrichment"), trait, celltype)),
            by=c("trait", "celltype")) %>% 
  filter(trait %in% c("bmm_v2.7", "anc_dna_2.7")) # anc_dna_2.7

res_anc %>% arrange(desc(Enrichment)) %>% View
res_anc %>% arrange(desc(Coefficient_star)) %>% View
res_anc %>% arrange(desc(`Prop._h2`)) %>% View

# raw coeff p value
res_anc %>% 
  filter(!celltype %in% c("allcells", "unioncells")) %>% 
  mutate(celltype = fct_reorder(celltype, desc(Coefficient_p))) %>% 
  mutate(qval = qvalue(Coefficient_p, fdr.level=0.05, pi0 = 1)$qvalue) %>% 
  ggplot(aes(y = celltype, x = -log10(Coefficient_p))) + geom_col(aes(fill = celltype))+
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color="blue") +
  scale_fill_manual(values = celltype_cols)+
  ggtitle("-log10(Coefficient p-value) from S-LDSC")+
  axis_theme
  
# coeff Z score
p_coeff <- res_anc %>% 
  #filter(!celltype %in% c("unioncells")) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype)) %>% 
  mutate(qval = qvalue(Coefficient_p, fdr.level=0.05, pi0 = 1)$qvalue,
         sig = as.factor(qval < 0.05 & Coefficient > 0)) %>% 
  mutate(celltype = fct_reorder(celltype, `Coefficient_star`)) %>% 
  ggplot(aes(y = celltype, x = `Coefficient_star`, fill = sig)) + geom_col()+
  ggtitle("Standardized coefficient from S-LDSC")+xlab("standardized coefficient")+ylab("cell type")+
  scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "blue"))+
  axis_theme+theme(legend.key.size =unit(0.2, "cm"))+
  guides(fill=guide_legend(title="FDR<0.05"))

ggsave2(paste0(plots_dir, "anc_dna_percelltype_sldsc_coef.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)

# Enrichment
ct_order <- res_anc %>% 
  #filter(!celltype %in% c("unioncells")) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype)) %>% 
  mutate(qval = qvalue(Coefficient_p, fdr.level=0.05, pi0 = 1)$qvalue,
         sig = as.factor(qval < 0.05 & Coefficient > 0)) %>% 
  mutate(celltype = fct_reorder(celltype, `Coefficient_star`)) %>% pull(celltype)

p_enrich <- res_anc %>% 
  #filter(!celltype %in% c("unioncells")) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype)) %>% 
  mutate(qval = qvalue(Enrichment_p, fdr.level=0.05, pi0 = 1)$qvalue,
         sig = as.factor(qval < 0.05 & Coefficient > 0)) %>% # View
  mutate(celltype = fct_relevel(celltype, levels(ct_order))) %>% 
  ggplot(aes(y = celltype, x = Enrichment, fill = sig)) + geom_col()+
  ggtitle("Enrichment from S-LDSC")+
  scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "blue"))+
  axis_theme+theme(legend.key.size =unit(0.2, "cm"))+
  guides(fill=guide_legend(title="FDR<0.05"))
p_enrich

plot_grid(
  p_coeff + theme(legend.position = "none"), p_enrich,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(0.8, 1),
  labels = c("A", "B"), label_size = 10
)

ggsave2(paste0(plots_dir, "anc_dna_2.7_onek1k_sldsc_unioncells_allgtex_onek1k.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# joint model on top 5 cell types
res_anc_cond <- read_tsv(paste0(ldsc_indir, "anc_dna.baseline.cs.5celltypes.cond.results")) %>% 
  janitor::clean_names() %>% 
  filter(grepl("^L2_", category)) %>% 
  mutate(celltype = c("NK", "CD8_ET", "CD8_NC", "CD4_NC", "Mono_NC"))
res_anc_cond_ld <- read_tsv(paste0(ldsc_indir, "anc_dna.baselineLD.cs.5celltypes.cond.results")) %>% 
  janitor::clean_names() %>% 
  filter(grepl("^L2_", category)) %>% 
  mutate(celltype = c("NK", "CD8_ET", "CD8_NC", "CD4_NC", "Mono_NC"))

# coef
res_anc_cond %>% 
  mutate(celltype = fct_reorder(celltype, coefficient_z_score),
         coefficient_p = pnorm(coefficient_z_score, lower.tail = F)) %>% 
  ggplot(aes(y = celltype, x = coefficient_z_score)) + geom_col()+
  ggtitle("Coefficient Z score from S-LDSC")+
  axis_theme+theme(legend.key.size =unit(0.5, "cm"))


ggsave2(paste0(plots_dir, "anc_dna_top5_cond_sldsc.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

# enrich
res_anc_cond_ld %>% 
  mutate(celltype = fct_reorder(celltype, enrichment)) %>% 
  ggplot(aes(y = celltype, x = enrichment)) + 
  geom_col(aes(fill = -log10(enrichment_p)))+
  ggtitle("Enrichment from S-LDSC")+
  axis_theme+theme(legend.key.size =unit(0.5, "cm"))

ggsave2(paste0(plots_dir, "anc_dna_top5_cond_sldsc_enrich.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)


#### gtex selection ldsc ####
gtex_res

res %>% filter(!grepl("baselineLD", filename)) %>% 
  left_join(h2g, by=c("trait", "annot", "celltype")) %>% 
  mutate(Coefficient_star = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient / h2g,
         Coefficient_star_std_error = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient_std_error / h2g,
         `Coefficient_star_z-score` = Coefficient_star / `Coefficient_star_std_error`)

# coeff
tmp <- gtex_res %>% 
  filter(!grepl("baselineLD", filename)) %>% 
  left_join(gtex_res_h2g, by=c("trait", "annot", "celltype")) %>% 
  #filter(!celltype %in% c("unioncells")) %>% 
  mutate(Coefficient_star = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient / h2g,
         Coefficient_star_std_error = M * sqrt(Prop._SNPs*(1-Prop._SNPs)) * Coefficient_std_error / h2g,
         `Coefficient_star_z-score` = Coefficient_star / `Coefficient_star_std_error`) %>% 
  mutate(qval = qvalue(Coefficient_p, fdr.level=0.05, pi0 = 1)$qvalue,
         sig = as.factor(qval < 0.05 & Coefficient > 0)) %>% 
  mutate(celltype = fct_reorder(celltype, `Coefficient_star`)) 

p_coeff <- tmp %>% 
  ggplot(aes(y = celltype, x = `Coefficient_star`, fill = sig)) + geom_col()+
  ggtitle("Standardized coefficient")+xlab("standardized coefficient")+ylab("cell type")+
  scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "blue"))+
  axis_theme+theme(legend.key.size =unit(0.2, "cm"))+
  guides(fill=guide_legend(title="FDR<0.05"))

# Enrichment
p_enrich <- gtex_res %>% 
  filter(!grepl("baselineLD", filename)) %>% 
  #filter(!celltype %in% c("unioncells")) %>% 
  mutate(qval = qvalue(Enrichment_p, fdr.level=0.05, pi0 = 1)$qvalue,
         sig = as.factor(qval < 0.05 & Coefficient > 0)) %>% 
  right_join(tmp %>% select(celltype, Coefficient_star)) %>% 
  mutate(celltype = fct_reorder(celltype, `Coefficient_star`)) %>% 
  ggplot(aes(y = celltype, x = Enrichment, fill = sig)) + geom_col()+
  ggtitle("Enrichment")+
  scale_fill_manual(values = c(`FALSE` = "grey", `TRUE` = "blue"))+
  axis_theme+theme(legend.key.size =unit(0.2, "cm"))+
  guides(fill=guide_legend(title="FDR<0.05"))


plot_grid(
  p_coeff + theme(legend.position = "none"), p_enrich,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1), rel_widths=c(0.8, 1),
  labels = c("A", "B"), label_size = 10
)

ggsave2(paste0(plots_dir, "anc_dna_2.7_gtex_brain_sldsc_union_allgtex_onek1k.png"),
        width = full_width*1.3, height = full_height, units = 'cm', dpi = 300)


#### MAF ####

bmm <- read_tsv("../data/OneK1K/annotation/anc_dna/bmm_v2.3.sumstats.gz")
nb_cs_pip0.5 %>% 
  separate(snp, into=c("rm", "pos", "ref", "alt"), sep="_") %>%
  mutate(pos = as.integer(pos)) %>% 
  select(-rm) %>% 
  left_join(geno_rsid) %>% 
  left_join(bmm, by=c("rsid"="SNP", "ref"="A2", "alt"="A1")) %>% 
  filter(!is.na(Z)) %>% 
  mutate(Z = Z*-1) %>% 
  ggplot(aes(x = celltype, y = Z, fill=celltype))+geom_boxplot() + axis_theme+
  scale_fill_manual(values = celltype_cols)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.key.size =unit(0.3, "cm"))

ggsave2(paste0(plots_dir, "finemap_pip0.5_intersect_anc_dna.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)
