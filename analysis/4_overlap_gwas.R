# overlap with gwas

library(tidyverse)
library(data.table)

# gwas catalog
gwas_cat <- fread("../data/GWAS/gwas_catalog_v1.0-associations_e111_r2024-03-28.tsv") %>% 
  janitor::clean_names()

names(gwas_cat)
head(gwas_cat) %>% View

summary(as.numeric(gwas_cat$p_value))

gwas_cat %>% distinct(disease_trait, .keep_all = T) %>% 
  select(disease_trait, everything()) %>% 
  filter(grepl("sclerosis", disease_trait, ignore.case=TRUE)) %>% View

# choose phenotypes
pheno <- c("Multiple sclerosis", "Type 1 diabetes", "Inflammatory bowel disease",
           "Crohn's disease", "Systemic lupus erythematosus", "Rheumatoid arthritis",
           "Ankylosing spondylitis")

gwas_cat %>% #filter(disease_trait %in% pheno) %>% 
  filter(grepl("diabetes|sclerosis|bowel|Crohn|lupus|Rheumatoid|Ankylosing", disease_trait, ignore.case=TRUE)) %>% 
  # distinct(disease_trait)
  mutate(pval = exp(-as.numeric(pvalue_mlog))) %>% 
  filter(pval < 5e-8) %>% 
  distinct(snps, .keep_all = T) %>% 
  inner_join(NB_egenes_df %>% distinct(rsid, celltype), by=c("snps"="rsid")) %>%
  select(snps, celltype, everything()) %>% 
  View

NB_egenes_df %>% View
