# integrate with GWAS

gwas_catalog <- fread("../data/gwas/gwas_catalog_v1.0.2-associations_e111_r2024-04-16.tsv") %>% 
  janitor::clean_names()

IBD <- c("Inflammatory bowel disease")
select_traits <- c("Inflammatory bowel disease",
                   "Multiple sclerosis", "Multiple sclerosis and HDL levels (pleiotropy)",
                   "Type 1 diabetes",
                   "Crohn's disease", "Crohn's disease-related phenotypes",
                    "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis) (pleiotropy)",
                   "Systemic lupus erythematosus", "Systemic lupus erythematosus (MTAG)",
                   "Rheumatoid arthritis", "Rheumatoid arthritis (ACPA-positive)", "Rheumatoid arthritis (ACPA-negative)",
                   "Ankylosing spondylitis")

gwas_catalog %>% # filter(grepl("erythematosus", disease_trait)) %>% distinct(pubmedid, .keep_all = T) %>% View
  filter(disease_trait %in% select_traits) %>% # distinct(pubmedid, .keep_all = T) %>% 
  filter(pvalue_mlog > -log(5e-8)) %>% 
  select(snps, disease_trait) %>% 
  distinct(snps, .keep_all = T) %>% 
  right_join(tqtl_egenes_df, by=c("snps"="rsid")) %>% 
  filter(!is.na(disease_trait)) %>% View
