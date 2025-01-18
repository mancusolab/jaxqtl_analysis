## read in data needed for mashr analysis
library(corrplot)
library(ggplotify)
library(mvsusieR)
library(gplots)

load("../result/mvsusie/result/mashr_Uk_fit_EZ_TRUE_complete_TRUE.RData")

# fine-map results
plots_dir <- "../result/cis/figures/celltype_specific/finemapped/"

finemap_indir <- "../result/finemap/result_wald/"
cols_names <- c("chr", "pos_1", "pos", "snp", "pip", "phenotype_id", "celltype")

nb_cs <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap.cs.tsv.gz"), 
                  col_names = cols_names) %>% 
  mutate(celltype=gsub(".cs.bed", "", celltype)) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

# take the strongest snp from each CS and keep the complete cases
nb_cs_leadsnp <- nb_cs %>%
  filter(celltype != "allcells") %>% 
  group_by(celltype, phenotype_id) %>%
  slice_max(pip, n=1, with_ties = FALSE) %>%
  ungroup() %>% 
  mutate(eqtl = paste0(phenotype_id, "_", snp)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")
summary(nb_cs_leadsnp$pip)
length(unique(nb_cs_leadsnp$eqtl))

# 2698 SNP with PIP >= 0.5
nb_cs_pip0.5 <- read_tsv(paste0(finemap_indir, "jaxqtl/bed/allcelltype_finemap_pip0.5.tsv.gz"), 
                         col_names = cols_names) %>% 
  filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
  filter(celltype != "allcells") %>% 
  mutate(eqtl = paste0(phenotype_id, "_", snp)) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")
nb_cs_pip0.5 %>% distinct(eqtl)

# 9966 SNP-gene pairs; 8994 complete cases
length(unique(finemap_cslead_eqtl_Z$eqtl))
finemap_cslead_eqtl_Z <- finemap_cslead_eqtl_Z %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")

finemap_eqtl_Z <- finemap_eqtl_Z %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by="phenotype_id")
finemap_eqtl_Z %>% distinct(eqtl)
# finemap_cslead_eqtl_Z, finemap_eqtl_Z

# nb_cs_leadsnp, nb_highpip
nb_highpip <- nb_cs_pip0.5
snp_df <- nb_highpip %>% filter(eqtl %in% finemap_eqtl_Z$eqtl) %>% 
  left_join(gene_lookup %>% select(phenotype_id, tss=end)) %>% 
  mutate(dist = abs(pos - tss))
ss_df <- finemap_eqtl_Z # finemap_cslead_eqtl_Z
suffix <- "pip0.5" # "pip0.5" "pip0.9", "cslead"

specific_eqtl_raw <- snp_df %>% 
  filter(eqtl %in% mash_sig_eqtls) %>% 
  add_count(eqtl, name = "n_celltype") %>% filter(n_celltype == 1)
shared_eqtl_raw <- snp_df %>% 
  filter(eqtl %in% mash_sig_eqtls) %>% 
  add_count(eqtl, name = "n_celltype") %>% filter(n_celltype > 1)

# 83% unique(n==1), 17% shared (>1)
length(unique(specific_eqtl_raw$eqtl))/length(unique(snp_df$eqtl))
length(unique(shared_eqtl_raw$eqtl))/length(unique(snp_df$eqtl))

sp_shared_intersect <- intersect(unique(specific_eqtl_raw$phenotype_id), 
                                 unique(shared_eqtl_raw$phenotype_id))

# 2% , 1% missing from finemap summary stats
n_eqtl <- length(unique(snp_df$eqtl));n_eqtl
