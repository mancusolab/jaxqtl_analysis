# read results

jqtl_colclass <- c("character", "numeric", "numeric", "character", 
                   rep("numeric", 6), "logical", rep("numeric", 7))
tqtl_colclass <- c("character", rep("numeric", 5), "character", rep("numeric", 10))


egenes_df <- tibble(cell = celltypes, gene_pass = NA, 
             jaxqtl_nb = NA, jaxqtl_lm_score = NA, jaxqtl_lm_wald = NA, jaxqtl_pois = NA, tqtl = NA,
             jaxqtl_nb_pi0 = NA, tqtl_pi0 = NA, rate_mean = NA)

NB_geneset <- c()
lm_score_geneset <- c()
lm_wald_geneset <- c()
pois_geneset <- c()
tqtl_geneset <- c()

NB_egenes_df <- c()
pois_egenes_df <- c()
jqtl_lm_score_egenes_df <- c()
jqtl_lm_wald_egenes_df <- c()
tqtl_egenes_df <- c()

# all gene result that pass threshold
NB_all_df <- c()
lm_score_all_df <- c()
lm_wald_all_df <- c()
tqtl_all_df <- c()
pois_all_df <- c()

# threshold 0.01
for (i in 1:length(celltypes)){
  cell_type <- celltypes[i]
  
  gene_summary <- read_tsv(paste0(genemeta_dir, cell_type, "_gene_summary.tsv.gz")) %>% 
    rename(phenotype_id = Geneid)
  
  genes_pass <- gene_summary %>% 
    filter(express_percent>=threshold) %>% pull(phenotype_id)
  
  egenes_df$gene_pass[i] <- length(genes_pass)
  
  jaxqtl_cisgenes_nb <- fread(paste0(jqtlcis_dir, 
                                     "jaxqtl_allres_cis_score.newbeta.nb.threshold", threshold,".",
                                     cell_type,".cisgenes.tsv"), header = F) %>% pull(V1)
  
  jaxqtl_cisgenes_lm_score <- fread(paste0(jqtlcis_dir,
                                           "jaxqtl_allres_cis_score.newbeta.lm.threshold", threshold,".",
                                           cell_type,".cisgenes.tsv"), header = F) %>% pull(V1)
  
  jaxqtl_cisgenes_lm_wald <- fread(paste0(jqtlcis_dir,
                                          "jaxqtl_allres_cis_wald.newbeta.lm.threshold", threshold,".",
                                          cell_type,".cisgenes.tsv"), header = F) %>% pull(V1)
  
  jaxqtl_cisgenes_pois <- fread(paste0(jqtlcis_dir,
                                       "jaxqtl_allres_cis_score.newbeta.pois.threshold", threshold,".",
                                       cell_type,".cisgenes.tsv"), header = F) %>% pull(V1)
  
  jaxqtl_nb_res <- fread(paste0(jqtlcis_dir, "jaxqtl_allres_cis_score.newbeta.nb.", cell_type, ".tsv.gz")) %>% 
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>%
    filter(phenotype_id %in% genes_pass) %>%
    mutate(celltype = cell_type,
           beta_shape1 = as.numeric(beta_shape1), beta_shape2 = as.numeric(beta_shape2), 
           qval = qvalue(pval_beta, fdr.level=fdr_val)$qvalue) %>% 
    left_join(gene_summary, by = "phenotype_id")
  
  jaxqtl_pois_res <- fread(paste0(jqtlcis_dir, "jaxqtl_allres_cis_score.newbeta.pois.", cell_type, ".tsv.gz")) %>%
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
  
  jaxqtl_lm_wald_res <- fread(paste0(jqtlcis_dir, "jaxqtl_allres_cis_wald.newbeta.lm.", cell_type, ".tsv.gz")) %>%
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>%
    filter(phenotype_id %in% genes_pass) %>%
    mutate(celltype = cell_type) %>%
    left_join(gene_summary, by = "phenotype_id")
  
  NB_all_df <- bind_rows(NB_all_df,
                         jaxqtl_nb_res)
  
  lm_score_all_df <- bind_rows(lm_score_all_df,
                               jaxqtl_lm_score_res)
  
  lm_wald_all_df <- bind_rows(lm_wald_all_df,
                              jaxqtl_lm_wald_res)
  # combine egenes results
  NB_egenes_df <- bind_rows(NB_egenes_df, 
                            jaxqtl_nb_res %>% 
                              filter(phenotype_id %in% jaxqtl_cisgenes_nb))

  pois_egenes_df <- bind_rows(pois_egenes_df,
                              jaxqtl_pois_res %>%
                                filter(phenotype_id %in% jaxqtl_cisgenes_pois))
  
  jqtl_lm_score_egenes_df <- bind_rows(jqtl_lm_score_egenes_df,
                                       jaxqtl_lm_score_res %>%
                                         filter(phenotype_id %in% jaxqtl_cisgenes_lm_score))
  
  jqtl_lm_wald_egenes_df <- bind_rows(jqtl_lm_wald_egenes_df,
                                      jaxqtl_lm_wald_res %>%
                                        filter(phenotype_id %in% jaxqtl_cisgenes_lm_wald))
  
  NB_geneset <- c(NB_geneset, jaxqtl_cisgenes_nb)
  lm_score_geneset <- c(lm_score_geneset, jaxqtl_cisgenes_lm_score)
  lm_wald_geneset <- c(lm_wald_geneset, jaxqtl_cisgenes_lm_wald)
  pois_geneset <- c(pois_geneset, jaxqtl_cisgenes_pois)
  
  tqtl_cisgenes <- fread(paste0(tqtlcis_dir, 
                                "tqtl_allres.", cell_type, ".threshold", threshold, ".cisgenes.tsv"),
                         header = F) %>% pull(V1)
  tqtl_geneset <- c(tqtl_geneset, tqtl_cisgenes)
  
  tqtl_res <- fread(paste0(tqtlcis_dir, "tqtl_allres.cis.", cell_type, ".tsv.gz")) %>% 
    mutate(celltype = cell_type) %>% 
    filter(phenotype_id %in% genes_pass) %>%
    mutate(qval = qvalue(pval_beta, fdr.level=fdr_val)$qvalue) %>% 
    left_join(gene_summary, by="phenotype_id")
  
  # combine egenes results
  tqtl_egenes_df <- bind_rows(tqtl_egenes_df, tqtl_res %>% 
                                filter(phenotype_id %in% tqtl_cisgenes))
  
  tqtl_all_df <- bind_rows(tqtl_all_df, tqtl_res)
  
  egenes_df$jaxqtl_nb[i] <- length(jaxqtl_cisgenes_nb)
  egenes_df$tqtl[i] <- length(tqtl_cisgenes)
  egenes_df$jaxqtl_lm_score[i] <- length(jaxqtl_cisgenes_lm_score)
  egenes_df$jaxqtl_lm_wald[i] <- length(jaxqtl_cisgenes_lm_wald)
  egenes_df$jaxqtl_pois[i] <- length(jaxqtl_cisgenes_pois)
  
  egenes_df$jaxqtl_nb_pi0[i] <- qvalue(jaxqtl_nb_res$pval_beta, fdr.level = fdr_val)$pi0
  egenes_df$tqtl_pi0[i] <- qvalue(tqtl_res$pval_beta, fdr.level = fdr_val)$pi0
  egenes_df$rate_mean[i] <- mean(gene_summary %>% filter(phenotype_id %in% genes_pass) %>% pull(rate_mean))
  
}

dim(tqtl_all_df); dim(NB_all_df)
dim(tqtl_egenes_df); dim(NB_egenes_df)

# allcell bulk eQTL results
gene_summary <- read_tsv(paste0(genemeta_dir, "allcells_gene_summary.tsv.gz")) %>% 
  rename(phenotype_id = Geneid)
genes_pass <- gene_summary %>% filter(express_percent>=threshold) %>% pull(phenotype_id)

jaxqtl_cisgenes_nb <- fread(paste0(jqtlcis_dir, 
                                   "jaxqtl_allres_cis_score.newbeta.nb.threshold",
                                   threshold,".allcells.cisgenes.tsv"), header = F) %>% pull(V1)

jaxqtl_allcell_egenes_df <- fread(paste0(jqtlcis_dir, 
                                         "jaxqtl_allres_cis_score.newbeta.nb.allcells.tsv.gz")) %>% 
  filter(phenotype_id %in% jaxqtl_cisgenes_nb) %>%
  mutate(celltype = "allcells", snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into = c("chrom", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>% 
  mutate(pos = as.integer(pos), chrom = as.integer(chrom)) %>% 
  left_join(gene_summary, by = "phenotype_id") %>% 
  left_join(geno_rsid, by=c("chr", "pos", "alt", "ref"))

jaxqtl_allcell_egenes_df <- jaxqtl_allcell_egenes_df %>% left_join(
  fread(paste0(jqtlcis_dir, 
               "jaxqtl_allres_cis_score.newbeta.nb.threshold0.01.allcells.leadsnp.tsv.gz")) %>% 
    select(phenotype_id, qval))

jaxqtl_cisgenes_linear <- fread(paste0(jqtlcis_dir, 
                                   "jaxqtl_allres_cis_score.newbeta.lm.threshold",
                                   threshold,".allcells.cisgenes.tsv"), header = F) %>% pull(V1)
jaxqtl_linear_allcell_egenes_df <- fread(paste0(jqtlcis_dir, 
                                         "jaxqtl_allres_cis_score.newbeta.lm.allcells.tsv.gz")) %>% 
  filter(phenotype_id %in% jaxqtl_cisgenes_linear) %>%
  mutate(celltype = "allcells", snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into = c("chrom", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>% 
  mutate(pos = as.integer(pos), chrom = as.integer(chrom)) %>% 
  left_join(gene_summary, by = "phenotype_id") %>% 
  left_join(geno_rsid, by=c("chr", "pos", "alt", "ref"))

tqtl_cisgenes <- fread(paste0(tqtlcis_dir, 
                              "tqtl_allres.allcells.threshold", threshold, ".cisgenes.tsv"), 
                       header = F) %>% pull(V1)

tqtl_allcell_egenes_df <- fread(paste0(tqtlcis_dir, "tqtl_allres.cis.allcells.tsv.gz")) %>% 
  filter(phenotype_id %in% tqtl_cisgenes) %>%
  mutate(celltype = "allcells", snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into = c("chrom", "pos", "ref", "alt"), sep = "_", remove = FALSE) %>% 
  mutate(pos = as.integer(pos), chrom = as.integer(chrom)) %>% 
  left_join(gene_summary, by = "phenotype_id") %>% 
  left_join(geno_rsid, by=c("chr", "pos", "alt", "ref"))

# find shared eGenes, eSNPs
# replace p value 0 with machine min
machine_min <- .Machine$double.xmin

NB_egenes_df <- as_tibble(NB_egenes_df) %>% 
  mutate(snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into = c("chrom", "pos", "ref", "alt"), sep = "_") %>% 
  mutate(variant_id = paste0(chr, ":", pos), pos = as.integer(pos)) %>% 
  mutate(Z = slope/slope_se, # score test Z
         chi2 = Z^2,
         pval_beta = ifelse(pval_beta==0, machine_min, pval_beta),
         pval_nominal = ifelse(pval_nominal==0, machine_min, pval_nominal)) %>% 
  left_join(cellmeta, by="celltype")

jqtl_lm_score_egenes_df <- as_tibble(jqtl_lm_score_egenes_df) %>% 
  mutate(snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into = c("chrom", "pos", "ref", "alt"), sep = "_") %>% 
  mutate(variant_id = paste0(chr, ":", pos), pos = as.integer(pos)) %>% 
  mutate(Z = slope/slope_se, # score test Z
         chi2 = Z^2,
         pval_beta = ifelse(pval_beta==0, machine_min, pval_beta),
         pval_nominal = ifelse(pval_nominal==0, machine_min, pval_nominal)) %>% 
  left_join(cellmeta, by="celltype")

# wald results for eGene-eSNP pair in each cell type
nb_wald <- fread(paste0(jqtlcis_dir, "jaxqtl.cis_qtl_pairs.nb.wald.egenes.tsv.gz")) %>%
  # filter(converged == TRUE) %>%
  rename(chr = chrom, ref = a1, alt = a0) %>%
  left_join(geno_rsid, by=c("chr", "pos", "alt", "ref"))

jqtl_lm_wald <- fread(paste0(jqtlcis_dir, "jaxqtl.cis_qtl_pairs.lm.wald.egenes.tsv.gz")) %>%
  filter(converged == TRUE) %>%
  rename(chr = chrom, ref = a1, alt = a0) %>%
  left_join(geno_rsid, by=c("chr", "pos", "alt", "ref"))

# for input of liftover script pylift_Nick.py: CHR    POS   REF   ALT   SNP
# nb_wald %>% distinct(chr, pos, ref, alt, rsid) %>% mutate(rsid=gsub("rs", "", rsid)) %>% 
#   write_tsv(paste0(jqtlcis_dir, "NB_cis_leadsnps.rsid.tsv.gz"))

# liftover using UCSC browder
nb_hg38 <- read_tsv(paste0(jqtlcis_dir, "/allcisgenes.hg38.bed"), col_names = F) %>% 
  rename(chr=X1, pos_38=X3) %>% 
  separate(X4, into=c("rm", "pos_19"), sep="-") %>% 
  mutate(pos_19 = as.integer(pos_19),
         chr = as.integer(gsub("chr", "", chr))) %>% 
  select(-c(rm, X2, X5)) %>% 
  distinct()

jqtl_lm_score_hg38 <- read_tsv(paste0(jqtlcis_dir, "/jqtl_lm_score_allcisgenes.hg38.bed"), col_names = F) %>% 
  rename(chr=X1, pos_38=X3) %>% 
  separate(X4, into=c("rm", "pos_19"), sep="-") %>% 
  mutate(pos_19 = as.integer(pos_19),
         chr = as.integer(gsub("chr", "", chr))) %>% 
  select(-c(rm, X2, X5)) %>% 
  distinct()

# write bed file of snps for assembly conversion
# tqtl_egenes_df %>% distinct(variant_id) %>% mutate(variant_id = gsub("chr|_b37","", variant_id)) %>% 
#   separate(variant_id, into=c("chr", "pos", "ref", "alt"),sep="_") %>% 
#   mutate(chr=paste0("chr", chr), pos_1 = as.integer(pos) - 1) %>% 
#   select(chr, pos_1, pos) %>% 
#   write_tsv(paste0(tqtlcis_dir, "allcisgenes.hg19.bed"), col_names = F)

tqtl_hg38 <- read_tsv(paste0(tqtlcis_dir, "allcisgenes.hg38.bed"), col_names = F) %>% 
  rename(chr=X1, pos_38=X3) %>% 
  separate(X4, into=c("rm", "pos_19"), sep="-") %>% 
  mutate(pos_19 = as.integer(pos_19),
         chr = as.integer(gsub("chr", "", chr))) %>% 
  select(-c(rm, X2, X5)) %>% 
  distinct()

# 41 not mapped to hg38
NB_egenes_df <- NB_egenes_df %>% 
  left_join(nb_wald %>% select(celltype, phenotype_id, snp, rsid,
                               slope_wald=slope, slope_se_wald=slope_se),
            by = c("celltype", "phenotype_id", "snp")) %>% 
  left_join(nb_hg38, by = c("chr", "pos"="pos_19")) %>% 
  left_join(eds_score, by="phenotype_id") %>% 
  left_join(pLI_score, by="phenotype_id") %>% 
  mutate(chrom = as.integer(chrom))

# here
jqtl_lm_score_egenes_df <- jqtl_lm_score_egenes_df %>% 
  left_join(jqtl_lm_wald %>% select(celltype, phenotype_id, snp, rsid,
                               slope_wald=slope, slope_se_wald=slope_se),
            by = c("celltype", "phenotype_id", "snp")) %>% 
  left_join(jqtl_lm_score_hg38, by = c("chr", "pos"="pos_19")) %>% 
  left_join(eds_score, by="phenotype_id") %>% 
  left_join(pLI_score, by="phenotype_id") %>% 
  mutate(chrom = as.integer(chrom))

# use the raw Z statistics (t-stat really here)
# 43 not mapped to hg38
tqtl_egenes_df <- as_tibble(tqtl_egenes_df) %>% 
  mutate(snp = variant_id) %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into=c("chr", "pos", "ref", "alt"), sep="_") %>% 
  left_join(geno_rsid %>% mutate(chr=as.character(chr),
                                 pos=as.character(pos)), 
            by=c("chr", "pos", "alt", "ref")) %>% 
  mutate(variant_id = paste0(chr, ":", pos),
         Z = slope/slope_se, # t-stat actually
         pval_beta = ifelse(pval_beta==0, machine_min, pval_beta),
         pval_nominal = ifelse(pval_nominal==0, machine_min, pval_nominal),
         chi2 = Z^2) %>% 
  left_join(NB_all_df %>% 
              select(phenotype_id, celltype, alpha_cov), 
            by = c("phenotype_id", "celltype")) %>% 
  left_join(eds_score, by="phenotype_id") %>% 
  left_join(pLI_score, by="phenotype_id") %>% 
  left_join(cellmeta, by="celltype") %>% 
  mutate(chr = as.integer(chr), pos = as.integer(pos)) %>% 
  left_join(tqtl_hg38, by = c("chr", "pos"="pos_19"))


# score test chi2
NB_all_df <- NB_all_df %>% 
  mutate(Z = slope/slope_se,
         pval_beta = ifelse(pval_beta==0, machine_min, pval_beta),
         pval_nominal = ifelse(pval_nominal==0, machine_min, pval_nominal),
         chi2 = Z^2) %>% 
  left_join(eds_score, by="phenotype_id") %>% 
  left_join(pLI_score, by="phenotype_id")


tqtl_all_df <- tqtl_all_df %>% 
  mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
  separate(variant_id, into=c("chr", "pos", "ref", "alt"), sep="_") %>% 
  mutate(variant_id = paste0(chr, ":", pos),
         Z = slope/slope_se,
         pval_beta = ifelse(pval_beta==0, machine_min, pval_beta),
         pval_nominal = ifelse(pval_nominal==0, machine_min, pval_nominal),
         chi2 = Z^2) %>% 
  left_join(eds_score, by="phenotype_id") %>% 
  left_join(pLI_score, by="phenotype_id") %>% 
  mutate(chr = as.integer(chr))


# summarize compare tqtl and NB, egenes-celltype pair
length(NB_geneset); length(tqtl_geneset) # 18907 vs. 16387
sum(!pois_geneset %in% NB_geneset) # poisson egenes are mostly found by NB model

# share eGene-cell type
# both_hits <- NB_egenes_df %>% 
#   rename(Z_jqtl = Z, pval_beta_jqtl = pval_beta, chi2_jqtl = chi2,
#          variant_id_jqtl = variant_id, tss_distance_jqtl = tss_distance) %>% 
#   inner_join(tqtl_egenes_df %>% 
#                select(phenotype_id, celltype, Z_lm=Z, pval_beta_tqtl=pval_beta, chi2_lm=chi2,
#                       variant_id_tqtl = variant_id, tss_distance_tqtl = start_distance), 
#              by = c("phenotype_id", "celltype"))

both_hits <- NB_egenes_df %>% 
  rename(Z_jqtl = Z, pval_beta_jqtl = pval_beta, chi2_jqtl = chi2,
         variant_id_jqtl = variant_id, tss_distance_jqtl = tss_distance) %>% 
  inner_join(jqtl_lm_score_egenes_df %>% 
               select(phenotype_id, celltype, Z_lm=Z, pval_beta_tqtl=pval_beta, chi2_lm=chi2,
                      variant_id_tqtl = variant_id, tss_distance_tqtl = tss_distance), 
             by = c("phenotype_id", "celltype"))

nrow(both_hits)

# share eGene and eSNP
# both_hits_samevar <- NB_egenes_df %>% 
#   rename(Z_jqtl = Z, pval_beta_jqtl = pval_beta, chi2_jqtl = chi2, slope_jqtl = slope, slope_jaxqtl_wald = slope_wald) %>% 
#   inner_join(tqtl_egenes_df %>% 
#                select(phenotype_id, variant_id, celltype, Z_lm=Z, pval_beta_tqtl=pval_beta, chi2_lm=chi2,
#                       slope_lm = slope, ref, alt, rsid, snp), 
#              by = c("phenotype_id", "celltype", "snp"))

both_hits_samevar <- NB_egenes_df %>% 
  rename(Z_jqtl = Z, pval_beta_jqtl = pval_beta, chi2_jqtl = chi2, slope_jqtl = slope, slope_jaxqtl_wald = slope_wald) %>% 
  inner_join(jqtl_lm_score_egenes_df %>% 
               select(phenotype_id, variant_id, celltype, Z_lm=Z, pval_beta_tqtl=pval_beta, chi2_lm=chi2,
                      slope_lm = slope, ref, alt, rsid, snp), 
             by = c("phenotype_id", "celltype", "snp"))


nrow(both_hits_samevar)
nrow(both_hits_samevar)/nrow(both_hits) # 65%

nrow(both_hits)/nrow(NB_egenes_df) # 78%

# genes only found by NB
# nb_only <- NB_egenes_df %>% 
#   anti_join(tqtl_egenes_df %>% select(phenotype_id, celltype), 
#             by = c("phenotype_id", "celltype"))

nb_only <- NB_egenes_df %>% 
  anti_join(jqtl_lm_score_egenes_df %>% select(phenotype_id, celltype), 
            by = c("phenotype_id", "celltype"))
dim(nb_only)

# genes only found by lm
# tqtl_only <- tqtl_egenes_df %>% 
#   anti_join(NB_egenes_df %>% select(phenotype_id, celltype), 
#             by = c("phenotype_id", "celltype"))

tqtl_only <- jqtl_lm_score_egenes_df %>% 
  anti_join(NB_egenes_df %>% select(phenotype_id, celltype), 
            by = c("phenotype_id", "celltype"))
dim(tqtl_only)


egenes_df <- egenes_df %>% 
  left_join(cell_counts, by = c("cell" = "celltype")) 

# note: multiple ensemble id can point to the same gene symbol
# tqtl_egenes_df %>% left_join(gene_lookup, by=c("phenotype_id")) %>% 
#   mutate(GeneSymbol = gsub("-","\\.", GeneSymbol)) %>% 
#   # add_count(celltype, GeneSymbol) %>% 
#   # filter(n > 1) %>% arrange(GeneSymbol) %>%  View
#   distinct(celltype, chr.x, GeneSymbol)


# confirmed that for linear model, score test stats always smaller than wald test
tqtl_egenes_df %>% mutate(Z_lm=slope/slope_se) %>% 
  select(celltype, phenotype_id, snp, Z_lm, true_df, pval_beta_tqtl=pval_beta) %>% 
  inner_join(jqtl_lm_score_egenes_df %>% 
               mutate(z_jqtl_lm=slope/slope_se) %>% 
               select(celltype, phenotype_id, snp=variant_id, z_jqtl_lm, true_nc, pval_beta_jqtl=pval_beta))


# threshold 0.1

egenes_df <- egenes_df %>% mutate(p_diff_nb_tqtl = NA, p_diff_nb_pois = NA, p_diff_nb_lm_score = NA)
for (i in 1:nrow(egenes_df)){
  egenes_df$p_diff_nb_tqtl[i] <- prop.test(x = c(egenes_df$jaxqtl_nb[i], egenes_df$tqtl[i]), 
                                           n = c(egenes_df$gene_pass[i], egenes_df$gene_pass[i])) %>% 
    broom::tidy() %>% pull(p.value)
  
  egenes_df$p_diff_nb_pois[i] <- prop.test(x = c(egenes_df$jaxqtl_nb[i], egenes_df$jaxqtl_pois[i]), 
                                           n = c(egenes_df$gene_pass[i], egenes_df$gene_pass[i])) %>% 
    broom::tidy() %>% pull(p.value)
  
  egenes_df$p_diff_nb_lm_score[i] <- prop.test(x = c(egenes_df$jaxqtl_nb[i], egenes_df$jaxqtl_lm_score[i]), 
                                               n = c(egenes_df$gene_pass[i], egenes_df$gene_pass[i])) %>% 
    broom::tidy() %>% pull(p.value)
}


# gather genes expression level across cell types
allgenes_express <- data.frame()
for (i in 1:length(celltypes)){
  cell_type <- celltypes[i]
  
  gene_summary <- read_tsv(paste0(genemeta_dir, cell_type, "_gene_summary.tsv.gz")) %>% 
    rename(phenotype_id = Geneid)
  
  allgenes_express <- bind_rows(allgenes_express, 
                                gene_summary %>% 
                                  select(phenotype_id, express_percent, 
                                         count_sum, tpm_0.1, express_percent_6,
                                         rate_mean, count_mean) %>% 
                                  mutate(celltype = cell_type))
}

summary(allgenes_express)

# put cell types on columns
allgenes_express_wide <- allgenes_express %>% select(phenotype_id, celltype, rate_mean) %>% 
  pivot_wider(
    names_from = celltype,
    values_from = rate_mean
  )

allgenes_express_per_wide <- allgenes_express %>% select(phenotype_id, celltype, express_percent) %>% 
  pivot_wider(
    names_from = celltype,
    values_from = express_percent
  )

allgenes_express_per_wide %>% 
  gather(key = celltype, value = express_percent, CD4_ET:CD4_SOX4) %>% 
  group_by(phenotype_id) %>% 
  summarize(ct = sum(express_percent > 0.01)) %>% 
  filter(ct > 0) %>% 
  distinct(phenotype_id) %>% write_tsv("./background_genes.txt",col_names = F)
