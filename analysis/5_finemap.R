# fine-map results

finemap_indir <- "../result/finemap/result_wald/"
cols_names <- c("chr", "pos_1", "pos", "snp", "pip", "phenotype_id", "celltype")

nb_cs <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap.cs.tsv.gz"), 
                  col_names = cols_names) %>% 
  mutate(celltype=gsub(".cs.bed", "", celltype)) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

# tqtl_cs <- read_tsv(paste0(finemap_indir, "tqtl/bed/allcelltype_finemap.cs.tsv.gz"))
lm_cs <- read_tsv(paste0(finemap_indir, "jaxqtl_lm_score/bed/allcelltype_finemap.cs.tsv.gz"),
                       col_names = cols_names) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  mutate(celltype=gsub(".cs.bed", "", celltype)) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

nb_cs_both <- nb_cs %>% inner_join(both_hits %>% select(phenotype_id, celltype))
lm_cs_both <- lm_cs %>% inner_join(both_hits %>% select(phenotype_id, celltype))

# nb_cs_both <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap_bothegenes_addtssbins.cs.tsv.gz"))
# tqtl_cs_both <- read_tsv(paste0(finemap_indir, "tqtl/bed/allcelltype_finemap_bothegenes_addtssbins.cs.tsv.gz"))

nb_cs_pip0.5 <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap_pip0.5.tsv.gz"), 
                         col_names = cols_names) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes))

lm_cs_pip0.5 <- read_tsv(paste0(finemap_indir, "jaxqtl_lm_score/bed/allcelltype_finemap_pip0.5.tsv.gz"), 
                         col_names = cols_names) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes))

sc_snps <- nb_cs_pip0.5 %>% filter(celltype != "allcells") %>% 
  #filter(pip >= 0.9) %>% 
  distinct(phenotype_id, snp)

nb_cs_pip0.5_wald <- read_tsv(paste0(jqtlcis_dir, "jaxqtl.finemap_eqtl_pip0.5_allest.nb.wald.tsv.gz")) %>% 
  distinct() %>% 
  filter(converged == TRUE) %>% 
  filter(celltype != "allcells") %>% 
  inner_join(sc_snps) %>% 
  mutate(eqtl = paste0(phenotype_id, "_", snp))

nb_cs_pip0.5_wald %>% write_tsv("eqtl_pip0.5.tsv.gz")

# ignore MHC genes

finemap_indir <- "../result/finemap/result_summary/"
method <- "wald" # "wald.all"
allgenes <- read_tsv("../result/finemap/result_wald/jaxqtl/nb_allegenes.tsv", F) %>% 
  filter(!X1 %in% c(MHC_genes, MAPT_genes))
dim(allgenes) # 24619 eGenes to map

# nb or lm, TODO: need to subset to eGenes
nb_cs_ct <- fread(paste0(finemap_indir, "NB_CS.", method, ".label.tsv.gz")) %>% 
  mutate(size = str_count(variable, ",")+1,
         sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  add_count(celltype, phenotype_id, name="num_cs")

lm_cs_ct <- fread(paste0(finemap_indir, "lm_CS.", method, ".label.tsv.gz")) %>% 
  mutate(size = str_count(variable, ",")+1,
         sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  add_count(celltype, phenotype_id, name="num_cs")

nb_cs_ct_both <- nb_cs_ct %>% inner_join(both_hits %>% select(phenotype_id, celltype))
lm_cs_ct_both <- lm_cs_ct %>% inner_join(both_hits %>% select(phenotype_id, celltype))


pip <- fread(paste0(finemap_indir, "NB_pip.", method,".label.tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "Bulk-eQTL", "sc-eQTL")) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)); 

lamb <- fread(paste0(finemap_indir, "NB_lambda.", method, ".label.tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "Bulk-eQTL", "sc-eQTL")) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes))

########## pre-check ########## 
# 18,247/18,281 have converged pip results
nrow(pip) # 18,247 / 18,281 succeed fine-mapping
NB_egenes_df %>% filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% nrow()
NB_egenes_df %>% filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% distinct(phenotype_id)
jaxqtl_allcell_egenes_df %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% distinct(phenotype_id) %>% nrow()

prop.test(x = c(6776, 6338), n=c(21183, 22878)) %>% broom::tidy()


allgenes %>% filter(X3 != "allcells") # 18,281 eGenes, 4520 unique genes
pip %>% filter(sc_bulk == "sc-eQTL") # 18,255 / 18,281
cs %>% filter(celltype != "allcells") # 14,979 has CS

tmp %>% select(cs, celltype, size, num_cs) %>% write_tsv("./cs.tsv.gz")

bind_rows(NB_egenes_df, jaxqtl_allcell_egenes_df) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  left_join(pip, by=c("phenotype_id", "celltype")) %>%
  mutate(susie_converge = as.factor(ifelse(is.na(pip_max), 0, 1))) %>% 
  filter(susie_converge == 0) %>% select(phenotype_id, celltype) %>% write_tsv("nb_notcvg.tsv",col_names = F)
  ggplot(aes(x = susie_converge, y = -log10(pval_nominal))) + 
  ggtitle("neglogp (max) ") +
  geom_boxplot(outlier.size = 0.6) + axis_theme

ggsave2(paste0(plots_dir, "NB_neglogp_cvg.png"),
          width = onecol_width/1.2, height = onecol_height, units = 'cm', dpi = 300)
  
# rerun 7 instances for no lamb value
dim(lamb)
allgenes %>% anti_join(lamb, by=c("phenotype_id", "celltype")) %>% write_tsv("./nb_nolamb.tsv",col_names = F)

# identify genes not converged,

lamb %>% 
  left_join(pip %>% select(phenotype_id, celltype, pip_max)) %>% 
  left_join(bind_rows(NB_egenes_df %>% select(phenotype_id, celltype, pval_nominal),
                      jaxqtl_allcell_egenes_df %>% select(phenotype_id, celltype, pval_nominal)), 
            by=c("phenotype_id", "celltype")) %>% 
  mutate(susie_converged = as.factor(ifelse(!is.na(pip_max), 1, 0))) %>% 
  filter(susie_converged == 0) %>% select(phenotype_id, celltype)
  ggplot(aes(x = susie_converged, y = -log10(pval_nominal))) + geom_boxplot(outlier.size = 0.5) +
  ggtitle("Discrepancy measure in obs_Z vs. exp_Z for eGenes") +
  axis_theme

ggsave2(paste0(plots_dir, "NB_lambda_cvg.png"),
        width = onecol_width/1.2, height = onecol_height, units = 'cm', dpi = 300)

# 4840 eGenes found by both bulk and sc
genes_bulk_sc <- unique(NB_geneset)[unique(NB_geneset) %in% unique(jaxqtl_allcell_egenes_df$phenotype_id)]
length(genes_bulk_sc)

########## analysis ########## 
#### PIP ####

# compare TSS on those genes fine-mapped by both
nb_cs %>% inner_join(both_hits %>% select(celltype, phenotype_id)) %>% 
  left_join(gtf %>% select(tss=end, phenotype_id), by=c("phenotype_id")) %>% 
  mutate(tss_dist = abs(pos - tss)) %>% pull(pip) %>% summary()

lm_cs %>% inner_join(both_hits %>% select(celltype, phenotype_id)) %>% 
  left_join(gtf %>% select(tss=end, phenotype_id), by=c("phenotype_id")) %>% 
  mutate(tss_dist = abs(pos - tss)) %>% pull(pip) %>% summary()

wilcox.test(pip ~ model, data=nb_cs %>% select(phenotype_id, celltype, nb_pip = pip, snp) %>% filter(celltype!="allcells") %>% 
  inner_join(lm_cs %>% select(phenotype_id, celltype, lm_pip = pip, snp) %>% filter(celltype!="allcells"),
             by=c("phenotype_id", "celltype", "snp")) %>% 
    gather(key = "model", value = "pip", c(nb_pip, lm_pip)))

nb_cs %>% select(phenotype_id, celltype, nb_pip = pip, snp) %>% filter(celltype!="allcells") %>% 
     inner_join(lm_cs %>% select(phenotype_id, celltype, lm_pip = pip, snp) %>% filter(celltype!="allcells"),
                by=c("phenotype_id", "celltype", "snp")) %>% 
  gather(key = "model", value = "pip", c(nb_pip, lm_pip)) %>% 
  ggplot(aes(x = model, y = pip)) + geom_boxplot()

#### number of CS ####

#### CS size and number of CS ####

# 12978 eGenes (cell type)
# mean number of CS: 1.15
# CS size: mean 32.24, median 16
nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>% summary()

nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>% 
  count(num_cs) %>% mutate(frac = n/sum(n))

nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% 
  distinct(phenotype_id, celltype, num_cs, size) %>% 
  summary()

## selection constraints
nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>%
  left_join(pLI_score, by=c("phenotype_id")) %>% 
  filter(!is.na(loeuf)) %>% 
  #group_by(pLI > 0.9) %>%
  group_by(loeuf < 0.6) %>% 
  count(num_cs) %>%
  mutate(total = sum(n),
         frac = n / total) %>%
  ggplot(aes(x = as.factor(num_cs), y = frac, fill = `loeuf < 0.6`)) + geom_col(position = "dodge") +
  scale_fill_manual(values = c(`TRUE`="blue", `FALSE`="grey")) + axis_theme

# eGenes with pLI > 0.9 have smaller number of CS than eGenes with higher pLI
# similar conclusion if using loeuf score

cs_low_constrain <- nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>%
  left_join(pLI_score, by=c("phenotype_id")) %>% 
  filter(loeuf > 0.6) %>% pull(num_cs)
cs_high_constrain <- nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>%
  left_join(pLI_score, by=c("phenotype_id")) %>% 
  filter(loeuf <= 0.6) %>% pull(num_cs)

mean(cs_low_constrain); mean(cs_high_constrain)
wilcox.test(cs_low_constrain, cs_high_constrain, alternative = "g") %>% broom::tidy()

cs_pli <- tibble(CT = unique(nb_cs_ct$celltype), slope = NA, slope_se = NA, pval = NA) %>% filter(CT != "allcells")
for (i in 1:nrow(cs_pli)){
  fit <- glm(num_cs ~ high_pLI,  
      data=nb_cs_ct %>% filter(sc_bulk == "sc-eQTL") %>% distinct(phenotype_id, celltype, num_cs) %>%
        left_join(gene_lookup %>% distinct(phenotype_id, gene_name), by="phenotype_id") %>%
        left_join(pLI_score, by=c("phenotype_id", "gene_name"="gene")) %>% 
        filter(celltype == cs_pli$CT[i]) %>% 
        filter(!is.na(pLI)) %>% 
        mutate(high_pLI = pLI >= 0.9), family="quasipoisson") %>% broom::tidy()
  cs_pli$slope[i] <- fit$estimate[2]
  cs_pli$slope_se[i] <- fit$`std.error`[2]
  cs_pli$pval[i] <- fit$p.value[2]
}
cs_pli

summary(fit)
# cs size comparing 
nb_cs_ct %>% 
  distinct(celltype, phenotype_id, .keep_all = T) %>% 
  group_by(sc_bulk) %>% 
  summarize(mean_cs = mean(num_cs),
            mean_cs_size = mean(size),
            med_cs_size = median(size))

# NB model has more CS: mean 1.167 vs. 1.124, 4.40E-6
v1 <- nb_cs_ct_both %>% distinct(phenotype_id, celltype, num_cs) %>% 
  filter(celltype != "allcells") %>% pull(num_cs); summary(v1)
v2 <- lm_cs_ct_both %>% distinct(phenotype_id, celltype, num_cs) %>% 
  filter(celltype != "allcells") %>% pull(num_cs); summary(v2)

wilcox.test(v1, v2, "t")

# NB model has less size in CS: mean 30.19, vs. 31.46, 3.11E-3
v1 <- nb_cs_ct_both %>% distinct(phenotype_id, celltype, num_cs, size) %>% 
  filter(celltype != "allcells") %>% pull(size); summary(v1)
v2 <- lm_cs_ct_both %>% distinct(phenotype_id, celltype, num_cs, size) %>% 
  filter(celltype != "allcells") %>% pull(size); summary(v2)

wilcox.test(v1, v2, "t")

nb_cs_ct %>% distinct(celltype, phenotype_id, .keep_all = T) %>% 
  filter(!celltype %in% c("allcells")) %>% 
  group_by(num_cs) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>% 
  mutate(percent = n_genes/sum(n_genes))

lm_cs_ct %>% distinct(celltype, phenotype_id, .keep_all = T) %>% 
  filter(!celltype %in% c("allcells")) %>% 
  group_by(num_cs) %>% 
  summarize(n_genes = n()) %>% 
  ungroup() %>% 
  mutate(percent = n_genes/sum(n_genes))

p1 <- nb_cs_ct %>% 
  filter(!celltype %in% c("allcells")) %>%
  ggplot(aes(x = num_cs)) + 
  geom_histogram() +
  ggtitle("Number of CS per eGene in a cell type") +
  axis_theme

p2 <- nb_cs_ct %>% 
  filter(!celltype %in% c("allcells")) %>%
  ggplot(aes(x = size)) + 
  geom_histogram() +
  ggtitle("Size of CS per eGene in a cell type") +
  axis_theme

# CS size and purity (average correlation), threshold min purity is 0.5^2
p3 <- cs %>% filter(sc_bulk == "sc-eQTL") %>% 
  ggplot(aes(x = cs_avg_r2)) + geom_histogram()+
  ggtitle("average r2 in CS") + 
  axis_theme

plot_grid(
  p1, p2, p3,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  rel_heights = c(1, 1, 1), rel_widths=c(1, 1, 1)
)

ggsave2(paste0(plots_dir, "NB_cisgenes_cs_num_size_cor.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

# number of CS across cell type
nb_cs_ct %>% 
  distinct(celltype, phenotype_id, num_cs) %>% 
  left_join(cell_counts, by=c("celltype"="cell_label")) %>% 
  mutate(mean_ct = ifelse(is.na(mean_ct), 10000, mean_ct),
         celltype = fct_reorder(celltype, desc(mean_ct))) %>% 
  ggplot(aes(x = celltype, y = num_cs)) +
  geom_violin() +
  # scale_fill_manual(values = celltype_cols)+
  ggtitle("Size of CS per eGene in a cell type") +
  axis_theme +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggsave2(paste0(plots_dir, "NB_cisgenes_num_cs.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

tmp <- nb_cs_ct %>% 
  distinct(celltype, phenotype_id, num_cs, sc_bulk) %>% 
  filter(sc_bulk == "sc-eQTL") %>% 
  left_join(cell_counts %>% mutate(prop = total_cell/sum(total_cell)), by="celltype") %>% 
  group_by(celltype, prop) %>% 
  summarize(mean_num_cs = mean(num_cs)) %>% 
  ungroup()

# num_cs ~ mean_ct: use poisson model, P=1.49E-5
glm(num_cs ~ prop,
    family = "poisson",
    data=nb_cs_ct %>% 
     distinct(celltype, phenotype_id, num_cs, sc_bulk) %>% 
     filter(sc_bulk == "sc-eQTL") %>% 
     left_join(cell_counts %>% mutate(prop = total_cell/sum(total_cell)), by="celltype")) %>% 
  broom::tidy()

p1 <- tmp %>% 
  ggplot(aes(x = prop, y = mean_num_cs)) + 
  geom_point(aes(color = celltype)) + 
  geom_smooth(method="lm") +
  xlab("cell type proportion") + ylab("average number of CS") +
  scale_color_manual(values = celltype_cols)+
  axis_theme+
  theme(
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 5),  # Adjust legend text size
    legend.title = element_text(size = 5)  # Adjust legend title size
  )+guides(color=guide_legend(title="cell type"))
p1
ggsave2(paste0(plots_dir, "NB_cisgenes_num_cs_ct.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

# size of CS
nb_cs_ct %>% 
  left_join(cell_counts, by=c("celltype"="cell_label")) %>% 
  mutate(mean_ct = ifelse(is.na(mean_ct), 10000, mean_ct),
         celltype = fct_reorder(celltype, desc(mean_ct))) %>% 
  ggplot(aes(x = celltype, y = log10(size))) +
  geom_violin() +
  ggtitle("Size of CS per eGene in a cell type") +
  axis_theme +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggsave2(paste0(plots_dir, "NB_cisgenes_cs_size.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

tmp <- nb_cs_ct %>% 
  left_join(cell_counts, by=c("celltype"="cell_label")) %>% 
  filter(sc_bulk == "sc-eQTL") %>% 
  group_by(celltype, total_cell) %>% 
  summarize(mean_size = mean(size),
            med_size = median(size)) %>% 
  ungroup()
tmp

# CS size ~ cell counts no relationship
lm(mean_size ~ total_cell, data = tmp) %>% broom::tidy()

tmp %>% 
  ggplot(aes(x = mean_ct, y = mean_size)) + geom_point() + geom_smooth(method="lm") +
  axis_theme

ggsave2(paste0(plots_dir, "NB_cisgenes_cs_size_mean_ct.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

#### number of CS vs. pLI or EDS ####

# number of CS proportional to pLI or EDS?
tmp <- nb_cs_ct %>% 
  left_join(pLI_score %>% select(phenotype_id, pLI, loeuf), by="phenotype_id") %>%
  left_join(eds_score %>% select(phenotype_id, EDS), by="phenotype_id") %>%
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by = "phenotype_id") %>% 
  filter(sc_bulk == "sc-eQTL")

tmp %>%
  # group_by(phenotype_id, EDS) %>%
  group_by(phenotype_id, loeuf, pLI) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup()

tmp %>% group_by(phenotype_id, loeuf) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = loeuf, y = mean_cs)) + geom_point() +
  geom_smooth(method = 'lm')

# no diff
lm(mean_cs ~ pLI, data = tmp %>% group_by(phenotype_id, pLI) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

# P=3.77E-5
lm(mean_cs ~ loeuf, data = tmp %>% group_by(phenotype_id, loeuf) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

# no diff
lm(mean_cs ~ EDS, data = tmp %>% group_by(phenotype_id, EDS) %>%
     summarize(mean_cs = mean(num_cs)) %>%
     ungroup()) %>% broom::tidy()

tmp %>%
  group_by(phenotype_id, loeuf) %>%
  summarize(mean_cs = mean(num_cs)) %>%
  ungroup() %>% 
  ggplot(aes(x = loeuf, y = mean_cs)) +
  geom_point(size=0.6) +
  geom_smooth(method="lm") + 
  axis_theme + xlab("LOEUF score")+
  ggtitle("#CS of eGenes and LOEUF score")+
  ylab("average number of CS per gene")

ggsave2(paste0(plots_dir, "NB_cisgenes_num_cs_selection.png"),
        width = onecol_width/1.2, height = onecol_height, units = 'cm', dpi = 300)

# # for Steven's grant: show how many eQTLs are shared with other cell types
# tmp <- nb_cs_pip0.5 %>% 
#   filter(celltype != "allcells") %>% 
#   add_count(phenotype_id, snp) %>% 
#   mutate(eqtl = paste0(phenotype_id, "_",snp)) %>% 
#   left_join(geno_rsid %>% select(snp=variant_id, rsid)) %>% 
#   left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, end))
# 
# tmp %>% 
#   distinct(celltype, phenotype_id,  snp, n) %>% 
#   mutate(shared = ifelse(n > 1, 1, 0)) %>% 
#   group_by(celltype) %>% 
#   summarize(n = n(),
#             n_shared = sum(shared)) %>% ungroup() %>% 
#   left_join(ind_celltype_ct %>% group_by(celltype) %>% summarize(ct = sum(counts))) %>%
#   mutate(fraction = n_shared / n) %>% 
#   arrange(ct)

#### heatmap of finemapped eQTL across cell types ####

nb_cs_pip0.5_wald %>% distinct() %>% 
  write_tsv("./nb_cs_pip0.5_wald.tsv.gz")

nb_cs_pip0.5_wald %>% 
  mutate(wts = (1/slope_se)^2) %>% 
  ggplot(aes(x = abs(tss_distance), y = abs(slope/slope_se))) + 
  geom_point(size=0.5) + 
  geom_smooth(method = "lm", aes(weight = wts)) +
  xlab("distance to TSS") + ylab("|effect size of finemapped sc-eQTL|")+
  axis_theme
ggsave2(paste0(plots_dir, "slope_tss.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

lm(abs_slope ~ abs(tss_distance), 
   data=nb_cs_pip0.5_wald %>% 
     mutate(abs_slope = abs(slope)), weights = (1/slope_se)^2) %>% 
  broom::tidy()

M <- nb_cs_pip0.5_wald %>% 
  distinct(phenotype_id, snp, celltype, slope, slope_se) %>% 
  mutate(Z = slope/slope_se) %>% 
  select(phenotype_id, snp, Z, celltype) %>% 
  spread(key = celltype, value = Z) %>% 
  select(B_IN:Plasma) %>% 
  as.matrix(); dim(M)

# 2698 x 14 cell types
sum(is.na(M))
colnames(M)
# M[is.na(M)] <- 0

k <- which(is.na(M), arr.ind=TRUE)
M[k] <- rowMeans(M, na.rm=TRUE)[k[,1]]

M <- t(scale(t(M)))


library(viridis)
colha <- ComplexHeatmap::HeatmapAnnotation(df = tibble(group = c("Bcell", "Bcell", "Tcell", "Tcell", "Tcell", "Tcell", "Tcell",
                                                                      "Tcell", "Myeloid", "Monocyte", "Monocyte", "NK", "NK", "Bcell")),
                                           col = list(group = c("Bcell" = "yellow", 
                                                              "Tcell" = "purple",
                                                              "Myeloid" = "red",
                                                              "Monocyte" = "orange",
                                                              "NK"="lightgreen")
                                           ), 
                                           which = "column")

ht <- ComplexHeatmap::Heatmap(M, column_title = "cell type",
                              row_title = "sc-eQTL", 
                              name = "scaled sc-eQTL Z-score",
                              # clustering_distance_rows = "pearson",
                              # clustering_method_rows = "ward.D2",
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "complete",
                              clustering_distance_columns = "euclidean",
                              clustering_method_columns = "complete",
                              col = magma(10),
                              top_annotation = colha,
                              show_column_names = TRUE,
                              show_row_names = FALSE)

draw(ht)

nb_cs_pip0.5 %>% 
  select("phenotype_id", "celltype", "snp") %>% 
  write_tsv("./finemap_eqtl_pip0.5.tsv")


#### factorgo analysis on shared eQTL between cell types ####

# write out input data for factorgo
nb_cs_pip0.5_wald %>% 
  filter(converged == TRUE) %>% 
  distinct(phenotype_id, snp, phenotype_id, celltype, slope, slope_se) %>% 
  mutate(Z = slope/slope_se) %>% 
  select(phenotype_id, snp, Z, celltype) %>% 
  spread(key = celltype, value = Z) %>% 
  unite("id", c(phenotype_id, snp), sep = "_") %>% 
  mutate_all(~replace_na(., 0)) %>% 
  select(-"allcells") %>% 
  write_tsv("./analysis/factorgo/Z.tsv.gz")

# prepare factorgo data
nb_cs_pip0.5_wald %>% 
  filter(converged == TRUE) %>% 
  distinct(phenotype_id, snp, phenotype_id, celltype, slope, slope_se) %>% 
  mutate(Z = slope/slope_se) %>% 
  select(phenotype_id, snp, Z, celltype) %>% 
  spread(key = celltype, value = Z) %>% 
  unite("id", c(phenotype_id, snp), sep = "_") %>% 
  mutate_all(~replace_na(., 0)) %>% 
  select(-"allcells") %>% 
  write_tsv("./analysis/factorgo/Z.tsv.gz")


df <- read_tsv("./analysis/factorgo/Z.tsv.gz")
snps <- df$id
traits <- colnames(df)[-1]

Zm <- read_tsv("./analysis/factorgo/finemap_eqtl_k5.Zm.tsv.gz", F, col_types = "n")
Wm <- read_tsv("./analysis/factorgo/finemap_eqtl_k5.Wm.tsv.gz", F, col_types = "n")

Zm %>% 
  mutate(celltype = traits,
         celltype = fct_reorder(celltype, X1)) %>% 
  ggplot(aes(x = X1, y = celltype)) + geom_col(aes(fill = celltype))+
  scale_fill_manual(values = celltype_cols) + axis_theme


Zm %>% 
  mutate(celltype = traits,
         celltype = fct_reorder(celltype, X1)) %>% 
  ggplot(aes(x = X1, y = X2)) + geom_point(aes(color = celltype))+
  scale_color_manual(values = celltype_cols)

Wm %>% 
  mutate(snp = snps) %>% 
  separate(snp, into=c("phenotype_id", "snp"), sep="_") %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  arrange(desc(abs(X1))) %>% View

########## share with others ########## 
# send finemapping result to Camelia
# to fine eQTL in CS across multiple cell types

nb_cs %>% group_by(celltype, phenotype_id) %>% summarize(sum = sum(pip)) %>% ungroup() %>% 
  summary()

nb_cs %>% filter(celltype != "allcells") %>% 
  select(chr, pos, snp, pip, phenotype_id, celltype) %>% 
  write_tsv("./nb_cs.tsv.gz")

nb_cs %>% add_count(pheno)

########## selection scan ########## 
bmm <- read_tsv("../../../Steven/shared/Ancient DNA Selection/bmm_v2.7.tsv.gz")

nb_cs_pip0.5 %>% 
  filter(pip >= 0.9) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp" = "variant_id")) %>% 
  left_join(bmm %>% select(rsid=SNP, Z_selection = Z), by="rsid") %>% 
  filter(!is.na(Z_selection)) %>% arrange(desc(abs(Z_selection))) %>% View
  # filter(celltype != "allcells") %>% 
  ggplot(aes(x = celltype, y = abs(Z_selection))) + geom_boxplot()
