library(readxl)
library(glue)
outdir <- "../result/cis/tables"

atacseq_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "ATACseq_chiou")

# eGene results for :
# Table 1: jaxqtl-negbinom, linear, Poisson, bulk (sample-coverage > 1%)
# Table 2: jaxqtl-negbinom, tensorQTL (sample-coverage > 10%)
# Put all summary statistics results on zenodo
nb_tab <- NB_egenes_df %>% 
  select(cell_type=celltype, 
         phenotype_id, chrom, variant_id=snp, rsid, tss_distance, 
         ma_count, af, 
         slope=slope_wald, slope_se=slope_se_wald, pval_nominal, qvalue=qval,
         beta_shape1, beta_shape2,alpha_disp=alpha_cov) %>% 
  mutate(model = "negbinom")

lm_tab <- jqtl_lm_score_egenes_df %>% 
  select(cell_type=celltype, 
         phenotype_id, chrom, variant_id=snp, rsid, tss_distance, 
         ma_count, af, 
         slope=slope_wald, slope_se=slope_se_wald, pval_nominal, qvalue=qval,
         beta_shape1, beta_shape2, alpha_disp=alpha_cov) %>% 
  mutate(model = "linear")

pois_wald <- read_tsv("../result/cis/celltype16_new/all_celltype/jaxqtl.cis_qtl_pairs.pois.wald.egenes.tsv.gz")

pois_tab <- pois_egenes_df %>% 
  select(-c(slope, slope_se)) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by="variant_id") %>% 
  left_join(pois_wald %>% select(celltype, phenotype_id, snp, slope, slope_se), 
            by=c("celltype", "phenotype_id", "variant_id"="snp")) %>% 
  select(cell_type=celltype, 
         phenotype_id, chrom, variant_id, rsid, tss_distance, 
         ma_count, af, 
         slope, slope_se, pval_nominal, qvalue=qval,
         beta_shape1, beta_shape2, alpha_disp=alpha_cov) %>% 
  mutate(model = "Poisson")

# write out eqtl list to pull slope, slope_se
# pois_egenes_df %>% select(phenotype_id, celltype, variant_id) %>% write_tsv("../result/cis/celltype16_new/all_celltype/pois_score_allcisgenes.tsv.gz")
# jaxqtl_allcell_egenes_df %>% 
#   mutate(celltype="allcells") %>% 
#   select(phenotype_id, celltype, variant_id=variant_id.y) %>% write_tsv("../result/cis/celltype16_new_fixalpha/all_celltype/nb_allcells_score_allcisgenes.tsv.gz")

eqtl_tab <- bind_rows(nb_tab, lm_tab, pois_tab) %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  select(phenotype_id, GeneSymbol, everything())
table(eqtl_tab$model)

# for write big tables, use openxlsx to skip JAVA
openxlsx::write.xlsx(x = eqtl_tab, 
           file = glue("{outdir}/eGene_results.xlsx"))

# for genes expression > 10%
# report number of eGenes

egene_0.1_tab <- egenes_df_0.1 %>% 
  left_join(egenes_df %>% select(cell, total_cell), by="cell") %>% 
  mutate(cell = fct_reorder(cell, desc(total_cell))) %>% 
  gather(key = software, value = eGenes, c(jaxqtl_nb, tqtl, saige_egenes)) %>%
  mutate(eGenes = as.numeric(eGenes),
         software = case_when(software == "jaxqtl_nb" ~ "jaxQTL", 
                              software == "saige_egenes" ~ "SAIGE-QTL",
                              software == "tqtl" ~ "tensorQTL"),
         software = fct_relevel(as.factor(software), c("jaxQTL", "SAIGE-QTL", "tensorQTL"))) %>% 
  select(cell_type=cell, `number_of_gene_pass_10%`="gene_pass", software, eGenes)
View(egene_0.1_tab)
  
xlsx::write.xlsx(x=egene_0.1_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
           sheetName = "egene_10", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

# fine-mapping results:
# PIP, CS results values for variant-gene-CT triplets
# meta-enrichment results for each annotation
# put all fine-mapping results on zenodo
tss_breaks <- c(0, 5000, 100000, 500000) # 
pip0.5_tab <- nb_cs_pip0.5 %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol)) %>% 
  select(-pos_1) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp"="variant_id")) %>% 
  left_join(gtf %>% select(phenotype_id, chr, tss=end), by=c("chr", "phenotype_id")) %>% 
  mutate(celltype = ifelse(celltype == "allcells", "bulk", celltype),
         tss_dist = pos - tss) %>% 
  mutate(tss_dist_abs = abs(tss_dist),
         tss_dist_abs = ifelse(tss_dist_abs < 1, tss_dist_abs+1, tss_dist_abs), # put tss_dist=0 to first bin
         bins = cut(tss_dist_abs, 
                    breaks = tss_breaks)) %>% 
  mutate(bins = case_when(bins == "(0,5e+03]" ~ "<5kb",
                          bins == "(1e+05,5e+05]" ~ ">100kb",
                          bins == "(5e+03,1e+05]" ~ "5kb-100kb")) %>% 
  select(cell_type=celltype, phenotype_id, GeneSymbol, chr, pos, variant_id=snp, rsid, pip, tss_distance=tss_dist, bins)

xlsx::write.xlsx(x=pip0.5_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
                 sheetName = "finemap_eqtl_pip0.5", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)


## enrichment meta: baseline, EpiMap, ATAC-seq
baseline_meta_tab <- bind_rows(baseline_meta_celltype %>% 
            filter(method == "NegBinom") %>% mutate(method = "sc-eQTL"),
          nb_meta %>% filter(celltype == "allcells" & annot %in% baseline_annot_list) %>% 
            select(annot, meta_OR, meta_lower, meta_upper, meta_p, meta_se) %>% 
            mutate(method = "bulk-eQTL")) %>% 
  mutate(annot = gsub("_phastCons46way$", "", annot)) %>% 
  select(method, annotation=annot, meta_OR, meta_OR_lb = meta_lower, meta_OR_ub = meta_upper)

xlsx::write.xlsx(x=baseline_meta_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
                 sheetName = "baseline_meta", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)

atac_epimap_meta_tab <- bind_rows(atac_dat %>% mutate(category = "ATAC-seq"),
          epimap_dat) %>%
  rename(`Cell type` = "method") %>% 
  mutate(celltype = fct_relevel(celltype, rev(celltype_order)),
         category = fct_relevel(category, 
                                c("ATAC-seq", "EpiMap promoter", "EpiMap enhancer"))) %>% 
  mutate(annot = gsub("ATACseq", "ATAC-seq", annot),
         annot = gsub(".enhancer.merge|.promoter.merge", "", annot)) %>% 
  select(cell_type=celltype, sc_bulk, annotation=annot, category,
         meta_OR, meta_OR_lb = meta_lower, meta_OR_ub = meta_upper)

xlsx::write.xlsx(x=atac_epimap_meta_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
                 sheetName = "atac_epimap_meta", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)

## enrichment meta: SCENT
scent_tab <- nb_enrich %>% mutate(method = celltype) %>% 
   filter(celltype != "allcells") %>% 
  mutate(category = case_when(celltype %in% c("B_IN", "B_Mem", "Plasma") ~ "Bcell",
                              celltype %in% c("CD4_ET", "CD4_SOX4", "CD4_NC", "CD8_S100B", "CD8_ET", "CD8_NC",
                                              "NK", "NK_R") ~ "Tnk",
                              celltype %in% c("Mono_C", "Mono_NC", "DC") ~ "Myeloid")) %>% 
  rename(enrich = mean_enrich) %>% 
  bind_rows(nb_allcell_scent %>% 
              mutate(celltype = "bulk") %>% 
              select(celltype, enrich=mean_enrich, se, OR_L, OR_U, category=ref_cell)) %>% 
  select(cell_type=celltype, annot_cell_type=category, enrich, enrich_lb = OR_L, enrich_ub = OR_U)

xlsx::write.xlsx(x=scent_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
                 sheetName = "scent_meta", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)

# mashr results for 2012 eQTLs

mashr_pip0.5_tab <- snp_df %>% 
  mutate(tss=pos-tss) %>% 
  filter(eqtl %in% mash_sig_eqtls) %>% 
  #select(eqtl, celltype) %>% 
  add_count(eqtl, name="number_of_celltype") %>% 
  select(eqtl, cell_type=celltype, phenotype_id, GeneSymbol, chr, pos, tss_distance=tss, pip)

xlsx::write.xlsx(x=mashr_pip0.5_tab, 
                 file = glue("{outdir}/mashr_pip0.5.xlsx"), 
                 sheetName = "mashr_pip0.5", 
                 col.names = TRUE, row.names = TRUE, append = FALSE)

mashr_tab <- as.data.frame(pm.mash.beta[nsig>0, ]) %>% rownames_to_column(var = "eqtl") %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=FALSE) %>% 
  mutate(pos = as.integer(pos)) %>% 
  left_join(gene_lookup %>% dplyr::select(phenotype_id, tss=end), by="phenotype_id") %>% 
  mutate(tss_distance = pos - tss) %>% 
  gather(key = celltype, value = `posterior_effect_size`, B_IN:Plasma) %>% 
  left_join(nb_cs_pip0.5_wald %>% select(eqtl, celltype, raw_effect_size=slope))

mashr_tab <- mashr_tab %>% 
  select(eqtl, tss_distance, celltype, raw_effect_size, posterior_effect_size) %>% 
  left_join(shared_sign_df %>% gather(key = celltype, value = `shared by sign`, B_IN:Plasma), by=c("celltype", "eqtl")) %>% 
  left_join(shared_mag_df %>% 
              select(-phenotype_id) %>% 
              gather(key = celltype, value = `shared by magnitude`, B_IN:Plasma), by=c("celltype", "eqtl"))
View(mashr_tab)

xlsx::write.xlsx(x=mashr_tab, 
                 file = glue("{outdir}/jaxqtl_supp_tables_allres.xlsx"), 
                 sheetName = "mashr_pip0.5_effect_sharing", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)


# SLDSC results

h2_tab <- tmp %>% filter(trait != "bmm_v2.7") %>% 
  separate(trait, into=c("rm1", "trait", "rm2"), sep="\\.") %>% 
  mutate(trait = gsub("_", " ", trait),
         trait = case_when(trait == "disease HYPOTHYROIDISM SELF REP" ~ "Hypothyroidsim",
                           .default = as.character(trait))) %>% 
  select(trait, starts_with("Prop"), starts_with("Enrichment"),
         trait_group)

xlsx::write.xlsx(x=h2_tab, 
                 file = glue("{outdir}/sldsc.xlsx"), 
                 sheetName = "sldsc_h2", 
                 col.names = TRUE, row.names = TRUE, append = FALSE)

ct_tab <- trait_celltype %>% 
  rename(p=Coefficient_p) %>% 
  mutate(significant_enrich = ifelse(qval_sig == TRUE, "yes", "no"),
         celltype = case_when(celltype == "allcells" ~ "bulk-eQTL",
                              celltype == "unioncells" ~ "sc-eQTL_union",
                              .default = celltype),
         celltype = fct_relevel(celltype, c("bulk-eQTL", "sc-eQTL_union"))) %>%
  select(trait, celltype, Prop._SNPs, Coefficient_star, Coefficient_p=p, qvalue=qval)

xlsx::write.xlsx(x=ct_tab, 
                 file = glue("{outdir}/sldsc.xlsx"), 
                 sheetName = "sldsc_CT", 
                 col.names = TRUE, row.names = TRUE, append = TRUE)
