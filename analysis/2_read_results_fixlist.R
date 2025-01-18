# read results

# restrict to genelist used in Yazar et.al.
# both yazar and SAIGEQTL analyzed genes expressed >10% individuals
new_threshold <- 0.1
yazar_genelist <- readxl::read_excel("../data/OneK1K/yazar_result/science.abf3041_tables_s6_to_s19.xlsx", 
                   sheet = 3, skip = 2) %>% 
  gather(key = cell, value = gene) %>% 
  mutate(cell = gsub(" ", "_", cell)) %>% 
  rename(celltype = cell) %>% 
  drop_na()

yazar_genelist %>% distinct(celltype, gene) # 156,687 gene-celltype tested

yazar_genelist %>% count(celltype)
  
yazar_genelist <- yazar_genelist %>% 
  distinct(celltype, gene) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by=c("gene" = "GeneSymbol")) %>% 
  distinct(celltype, gene, phenotype_id)

jqtl_colclass <- c("character", "numeric", "numeric", "character", 
                   rep("numeric", 6), "logical", rep("numeric", 7))
tqtl_colclass <- c("character", rep("numeric", 5), "character", rep("numeric", 10))


egenes_df_0.1 <- tibble(cell = celltypes, gene_pass = NA, 
                    jaxqtl_nb = NA, jaxqtl_lm_score = NA, jaxqtl_lm_wald = NA, jaxqtl_pois = NA, tqtl = NA,
                    rate_mean = NA)

NB_geneset_0.1 <- c()
tqtl_geneset_0.1 <- c()
lm_score_geneset_0.1 <- c()

NB_egenes_df_0.1 <- c()
tqtl_egenes_df_0.1 <- c()
jqtl_lm_score_egenes_df_0.1 <- c()

# all gene result that pass threshold
NB_all_df_0.1 <- c()
tqtl_all_df_0.1 <- c()

# threshold 0.01
for (i in 1:length(celltypes)){
  cell_type <- celltypes[i]
  
  # 6340 genes
  yazar_list <- yazar_genelist %>% 
    filter(celltype == cell_type) %>% 
    distinct(phenotype_id) %>% 
    pull(phenotype_id)
  
  gene_summary <- read_tsv(paste0(genemeta_dir, cell_type, "_gene_summary.tsv.gz")) %>% 
    rename(phenotype_id = Geneid)
  
  genes_pass <- gene_summary %>% 
    filter(express_percent>=new_threshold) %>% 
    filter(phenotype_id %in% yazar_list) %>% 
    pull(phenotype_id)
  
  egenes_df_0.1$gene_pass[i] <- length(genes_pass)
  
  jaxqtl_nb_res <- fread(paste0(jqtlcis_dir, "jaxqtl_allres_cis_score.newbeta.nb.", cell_type, ".tsv.gz")) %>% 
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>%
    filter(phenotype_id %in% genes_pass) %>%
    mutate(celltype = cell_type,
           beta_shape1 = as.numeric(beta_shape1), beta_shape2 = as.numeric(beta_shape2), 
           qval = qvalue(pval_beta, fdr.level=fdr_val)$qvalue) %>% 
    left_join(gene_summary, by = "phenotype_id")
  
  jaxqtl_lm_score_res <- fread(paste0(jqtlcis_dir, "jaxqtl_allres_cis_score.newbeta.lm.", cell_type, ".tsv.gz")) %>%
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>%
    filter(phenotype_id %in% genes_pass) %>%
    mutate(celltype = cell_type,
           beta_shape1 = as.numeric(beta_shape1), beta_shape2 = as.numeric(beta_shape2), 
           qval = qvalue(pval_beta, fdr.level=fdr_val)$qvalue) %>% 
    left_join(gene_summary, by = "phenotype_id")
  
  
  NB_all_df_0.1 <- bind_rows(NB_all_df_0.1,
                         jaxqtl_nb_res)
  
  jaxqtl_cisgenes_nb <- jaxqtl_nb_res %>% filter(qval < fdr_val) %>% pull(phenotype_id) %>% unique()
  jaxqtl_cisgenes_lm_score <- jaxqtl_lm_score_res %>% filter(qval < fdr_val) %>% pull(phenotype_id) %>% unique()
  
  # combine egenes results
  NB_egenes_df_0.1 <- bind_rows(NB_egenes_df_0.1, 
                            jaxqtl_nb_res %>% 
                              filter(phenotype_id %in% jaxqtl_cisgenes_nb))
  jqtl_lm_score_egenes_df_0.1 <- bind_rows(jqtl_lm_score_egenes_df_0.1,
                                       jaxqtl_lm_score_res %>%
                                         filter(phenotype_id %in% jaxqtl_cisgenes_lm_score))
  
  
  NB_geneset_0.1 <- c(NB_geneset_0.1, jaxqtl_cisgenes_nb)
  lm_score_geneset_0.1 <- c(lm_score_geneset_0.1, jaxqtl_cisgenes_lm_score)

  tqtl_res <- fread(paste0(tqtlcis_dir, "tqtl_allres.cis.", cell_type, ".tsv.gz")) %>% 
    mutate(celltype = cell_type) %>% 
    filter(phenotype_id %in% genes_pass) %>%
    mutate(qval = qvalue(pval_beta, fdr.level=fdr_val)$qvalue) %>% 
    left_join(gene_summary, by="phenotype_id")
  
  tqtl_cisgenes <- tqtl_res %>% filter(qval < fdr_val) %>% pull(phenotype_id) %>% unique()
  tqtl_geneset_0.1 <- c(tqtl_geneset_0.1, tqtl_cisgenes)
  
  # combine egenes results
  tqtl_egenes_df_0.1 <- bind_rows(tqtl_egenes_df_0.1, tqtl_res %>% 
                                filter(phenotype_id %in% tqtl_cisgenes))
  
  tqtl_all_df_0.1 <- bind_rows(tqtl_all_df_0.1, tqtl_res)
  
  egenes_df_0.1$jaxqtl_nb[i] <- length(jaxqtl_cisgenes_nb)
  egenes_df_0.1$tqtl[i] <- length(tqtl_cisgenes)
  egenes_df_0.1$jaxqtl_lm_score[i] <- length(jaxqtl_cisgenes_lm_score)
}

egenes_df_0.1

# read saigeqtl results
saigeqtl_genes_filter <- readxl::read_excel("../ref/saigeqtl/media-1.xlsx", sheet = 1, skip = 1) %>% 
  janitor::clean_names()
saigeqtl_egenes <- readxl::read_excel("../ref/saigeqtl/media-1.xlsx", sheet = 6, skip = 2) %>% 
  janitor::clean_names() %>% 
  slice(1:28) %>% 
  left_join(saigeqtl_genes_filter, by="cell_type") %>% 
  mutate(method = ifelse(method == "TensorQTL", "tensorQTL (Zhou et.al)", method),
         number_of_e_genes_fdr_0_05 = as.integer(number_of_e_genes_fdr_0_05))


egenes_df_0.1 <- egenes_df_0.1 %>% mutate(p_diff = NA) %>% 
  left_join(saigeqtl_egenes %>% 
              filter(method == "SAIGEQTL") %>% 
              select(cell=cell_type, 
                     saige_pass=number_of_genes_expressed_in_10_percent_of_donors,
                     saige_egenes=number_of_e_genes_fdr_0_05),
            by="cell")
for (i in 1:nrow(egenes_df_0.1)){
  egenes_df_0.1$p_diff[i] <- prop.test(x = c(egenes_df_0.1$jaxqtl_nb[i], egenes_df_0.1$saige_egenes[i]), 
                                           n = c(egenes_df_0.1$gene_pass[i], egenes_df_0.1$saige_pass[i])) %>% 
    broom::tidy() %>% pull(p.value)
}

# pick NK, Mono_C, Plasma
egenes_df_0.1 %>% filter(p_diff < 0.05/14)
