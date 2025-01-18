########## Case study: RA ########## 
# Ishigaki2022: MA analysis on RA

# read coloc analysis result
RA_coloc <- read_tsv("../result/coloc/RA/allres") %>% 
  arrange(PP.H4.abf) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, tss=end), by=c("gene"="phenotype_id")) %>% 
  left_join(geno_rsid, by=c("hit1" = "variant_id")) %>% rename(rsid_eqtl = rsid) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid_gwas=rsid), by=c("hit2" = "variant_id")) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype),
         tss_dist = pos - tss) %>% 
  rename(phenotype_id = gene)

coloc_snps <- read_tsv("../result/coloc/RA/jaxqtl.cis_qtl_pairs.nb.wald.tsv.gz") %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id") %>% 
  filter(converged == TRUE) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype))

# coloc results from yazar
yazar_coloc <- readxl::read_xlsx("../ref/yazar2022_tables_s6_to_s19/science.abf3041_tables_s6_to_s19.xlsx", 
                               sheet ="Table S19", skip = 2)
table(yazar_coloc$Condition)

yazar_coloc %>% 
  filter(Condition == "RA" & SNP.PP.H4 >= 0.9) %>% 
  filter(`Gene ID` == "IL6ST") %>% View

RA_coloc %>% left_join(nb_wald %>% select(celltype, snp, phenotype_id, slope, slope_se, pval_nominal),
                       by=c("hit1"="snp", "phenotype_id", "celltype")) %>% 
  filter(celltype != "bulk") %>% 
  filter(is.na(slope)) %>% View

View(RA_coloc)
hist(RA_coloc$PP.H4.abf)

### find gene ###
# fine a gene that is not found by linear model and not by bulk

RA_coloc_bulk <- RA_coloc %>% left_join(nb_wald %>% select(celltype, snp, phenotype_id, slope, slope_se, pval_nominal),
                       by=c("hit1"="snp", "phenotype_id", "celltype")) %>% 
  filter(PP.H4.abf >= 0.9) %>% 
  filter(celltype == "bulk")

RA_coloc %>% filter(PP.H4.abf >= 0.9) %>% View
  anti_join(RA_coloc_bulk) %>% 
  anti_join(jqtl_lm_score_egenes_df %>% 
              distinct(phenotype_id, celltype), by=c("phenotype_id", "celltype")) %>% View

gene <- "ENSG00000134352"
jqtl_lm_score_egenes_df %>% 
  #distinct(phenotype_id, celltype) %>% 
  filter(phenotype_id == gene) %>% View

NB_egenes_df %>% 
  #distinct(phenotype_id, celltype) %>% 
  filter(phenotype_id == gene) %>% View

RA_coloc %>% 
  #filter(celltype == "bulk") %>% 
  filter(phenotype_id == gene) %>% View

RA_coloc_bulk %>% filter(phenotype_id == "ENSG00000101017")

jaxqtl_allcell_egenes_df %>% filter(phenotype_id == "ENSG00000101017")
jaxqtl_linear_allcell_egenes_df %>% filter(phenotype_id == "ENSG00000101017")

NB_egenes_df %>% distinct(phenotype_id, celltype, pval_nominal) %>% 
  anti_join(jaxqtl_allcell_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
  anti_join(jqtl_lm_score_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
  anti_join(jaxqtl_linear_allcell_egenes_df %>% distinct(phenotype_id, celltype)) %>% 
  arrange(pval_nominal) %>% View

# 43 gene-cell type, 15 genes were only found in a single cell type
RA_coloc %>% filter(PP.H4.abf >= 0.9) %>% filter(celltype!="bulk") %>% 
  distinct(GeneSymbol, celltype, phenotype_id) %>% View

RA_coloc %>% filter(PP.H4.abf >= 0.9) %>% filter(celltype!="bulk") %>% 
  add_count(GeneSymbol) %>% filter(n>1) %>% 
  arrange(GeneSymbol) %>% 
  View

RA_coloc %>% filter(PP.H4.abf >= 0.9) %>% 
  left_join(coloc_snps %>% dplyr::select(celltype, phenotype_id, snp, slope, slope_se),
            by=c("celltype", "phenotype_id", "hit1"="snp")) %>%
  arrange(GeneSymbol, desc(abs(tss_dist))) %>%
  add_count(hit1, phenotype_id) %>% View

# distant enhancer
# CMC1: ENSG00000187118

# PIP around 1, rs7731626
RA_coloc %>% filter(GeneSymbol == "IL6ST" & PP.H4.abf >= 0.9) %>% View
RA_coloc %>% filter(GeneSymbol == "ANKRD55" & PP.H4.abf >= 0.9) %>% View


coloc_snps %>% filter(snp == "chr5_55444683_G_A_b37") %>% View
gwas %>% filter(snp == "chr5_55444683_G_A_b37")

# PIP is small
RA_coloc %>% filter(GeneSymbol == "RBM39" & PP.H4.abf >= 0.9) %>% View

# PIP is 0.964, 0.655 (CD8)
RA_coloc %>% filter(GeneSymbol == "SESN3" & PP.H4.abf >= 0.9) %>% View


nb_cs_pip0.5 %>% filter(snp == "chr4_80894682_A_G_b37")

coloc_snps %>% filter(snp == "chr3_27799034_C_A_b37" & GeneSymbol == "CMC1") %>% View

# rs7731626
RA_coloc %>% filter(rsid_eqtl == "rs7731626" & PP.H4.abf >= 0.9) %>% 
  left_join(coloc_snps %>% select(celltype, phenotype_id, snp, slope, slope_se),
            by=c("celltype", "phenotype_id", "hit1"="snp")) %>%
  View

# effect size
NB_egenes_df %>% filter(phenotype_id %in% c("ENSG00000134352", "ENSG00000164512")) %>% 
  mutate(rate = exp(slope)-1, 
         rate_lb = exp(slope-1.96*slope_se)-1,
         rate_ub = exp(slope+1.96*slope_se)-1) %>% 
  arrange(phenotype_id) %>% View

# total 22 genes overlap
coloc_egenes <- RA_coloc %>% 
  filter(PP.H4.abf >= 0.9) %>% 
  distinct(GeneSymbol) %>% pull(GeneSymbol)
      
sc_egenes <- RA_coloc %>% 
  filter(PP.H4.abf >= 0.9) %>% 
  group_by(GeneSymbol) %>% mutate(found_bulk = sum(celltype == "bulk")) %>% ungroup() %>% 
  distinct(GeneSymbol, found_bulk) %>% 
  filter(found_bulk < 1) %>% 
  pull(GeneSymbol)


# plot 
expand.grid(GeneSymbol = RA_coloc %>% 
              filter(PP.H4.abf >= 0.9) %>% 
              distinct(GeneSymbol) %>% pull(GeneSymbol),
            celltype = c(celltypes, "bulk")) %>% 
  left_join(RA_coloc %>% 
              filter(PP.H4.abf >= 0.9) %>% 
              distinct(GeneSymbol, celltype) %>% 
              mutate(coloc = 1)) %>% 
  mutate(coloc = replace_na(coloc, 0),
         coloc = as.factor(coloc),
         celltype = fct_relevel(celltype, c(celltype_order, "bulk")),
         GeneSymbol = fct_relevel(GeneSymbol, sc_egenes)) %>% 
  ggplot(aes(y = GeneSymbol, x = celltype, fill = coloc)) +
  geom_tile(color = "black")+
  scale_fill_manual(values = c(`0`="white", `1`="cyan3"))+
  axis_theme_font6+
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        legend.margin = margin(6, 6, 6, 6)
  )
  
ggsave2(paste0(plots_dir, "Fig6_gene_coloc.png"),
        width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)


# snp: chr-pos-ref-alt

# eQTL or sQTL
RA_coloc_blueprint <- readxl::read_excel("../data/gwas/Ishigaki2022/Tables.xlsx", 
                              sheet = "ST-9", skip = 4) %>% 
  janitor::clean_names() %>% 
  select(snp_gwas = rs_id_4, snp_gtex = rs_id_7, gene_name, pp, cell)

RA_coloc_gtex <- readxl::read_excel("../data/gwas/Ishigaki2022/Tables.xlsx", 
                                         sheet = "ST-10", skip = 4) %>% 
  janitor::clean_names() %>% 
  select(snp_gwas = rs_id_4, snp_gtex = rs_id_7, gene_name, pp, tissue)

View(RA_snps)
nrow(RA_snps)

# 14 out of 35 GWAS (pip>0.5) are found in eQTL
# variant_id_7 is hg19
RA_snps %>% filter(pip_gwas>=0.5) %>% View

nb_cs_pip0.5 %>%
  left_join(geno_rsid %>% select(snp=variant_id, rsid), by="snp") %>% 
  filter(!is.na(rsid) & celltype != "allcells") %>% 
  inner_join(RA_snps %>% filter(pip_gwas>=0.5), by="rsid") %>% 
  inner_join(gene_lookup %>% select(phenotype_id, GeneSymbol, end)) %>% 
  mutate(tss_dist = pos - end) %>% 
  left_join(nb_cs_pip0.5_wald %>% 
              distinct(phenotype_id, snp, slope, slope_se, celltype, pval_nominal), 
            by=c("celltype", "phenotype_id", "snp")) %>% 
  filter(!is.na(gwas)) %>% arrange(snp) %>% 
  select(snp, rsid, phenotype_id, celltype, gene_name, GeneSymbol, slope, slope_se, tss_dist) %>% 
  # distinct(snp) %>% 
  View

#### popscore ####

popscore <- read_tsv('../data/OneK1K/annotation/popscore/PoPS_FullResults.txt.gz') %>% 
  filter(trait == "RA") %>% 
  distinct(ensgid, gene, pops_score)
gene_lookup %>% filter(GeneSymbol %in% c("IL6ST", "ANKRD55")) %>% View

popscore %>% arrange(desc(pops_score))
mean(popscore$pops_score < 0.743)
mean(popscore$pops_score < 0.13)
popscore %>% 
  filter(ensgid %in% c("ENSG00000134352", "ENSG00000164512"))


#### fine gene ####

# linear model sc-eqtl
jqtl_lm_score_egenes_df %>% 
  distinct(phenotype_id, celltype) %>% 
  filter(phenotype_id == "ENSG00000101017")

# linear model bulk
jaxqtl_linear_allcell_egenes_df

jaxqtl_allcell_egenes_df

# tried: 
NB_egenes_df %>% 
  distinct(phenotype_id, snp, celltype, tss_distance) %>% 
  add_count(phenotype_id, name = "num_CT") %>% 
  inner_join(nb_cs_pip0.5 %>% select(phenotype_id, snp, celltype, pip)) %>% 
  anti_join(jaxqtl_allcell_egenes_df %>% select(phenotype_id, celltype)) %>% 
  anti_join(jaxqtl_linear_allcell_egenes_df %>% select(phenotype_id, celltype)) %>% 
  anti_join(jqtl_lm_score_egenes_df %>% select(phenotype_id, celltype)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, end), by="phenotype_id") %>%
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  arrange(desc(abs(tss_distance))) %>% View
  #filter(snp == "chr12_2423857_A_G_b37")
  #distinct(phenotype_id, snp) %>% write_tsv("eqtl_nb_only.tsv")

NB_egenes_df %>% filter(phenotype_id == "ENSG00000163687") %>% pull(celltype)
