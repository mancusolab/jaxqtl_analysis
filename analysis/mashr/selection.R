nb_cs_leadsnp

bmm <- read_tsv("../data/OneK1K/annotation/anc_dna/bmm_v2.7.tsv.gz")

finemap_cslead_eqtl_Z %>% 
  # filter(pip>=0.9) %>% 
  # filter(abs(slope) < 10) %>% 
  ggplot(aes(x = af, y = abs(slope))) + geom_point(size=0.5)+
  # geom_smooth(method="loess")+
  facet_wrap(~celltype) + axis_theme

ggsave2(paste0(plots_dir, "cslead_complete_eqtl_effect_AF.png"),
        width = full_width, height = full_height*2, units = 'cm', dpi = 300)

nb_cs_leadsnp %>% 
  left_join(nb_wald %>% select(phenotype_id,snp,celltype,slope)) %>% 
  #filter(pip >= 0.9) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp"="variant_id")) %>% 
  #group_by(celltype) %>% 
  #slice_max(pip, n = 1, with_ties = F) %>% 
  inner_join(bmm, by=c("rsid"="SNP")) %>% #View
  slice_max(abs(Z), n=5) %>% View
  ggplot(aes(x = celltype, y = Z)) + geom_boxplot()

shared_mag_df %>% gather(key = celltype, value = shared, B_IN:Plasma) %>% 
    distinct(phenotype_id, eqtl, celltype, shared) %>% 
    group_by(phenotype_id, eqtl) %>% 
    summarize(n_shared = sum(shared)) %>% ungroup() %>% 
  separate(eqtl, into=c("gene", "snp"), sep = "_", extra = "merge") %>% 
  distinct(snp, n_shared) %>% add_count(snp, name="n_snp") %>% filter(n_snp < 2) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp"="variant_id")) %>% 
  inner_join(bmm, by=c("rsid"="SNP")) %>% 
  ggplot(aes(x = n_shared, y = abs(Z))) + geom_point()

snp_df %>% count(phenotype_id, snp, name = "n_celltype") %>% 
  add_count(snp, name="n_snp") %>% filter(n_snp < 2) %>% 
  distinct(snp, n_celltype) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp"="variant_id")) %>% 
  inner_join(bmm, by=c("rsid"="SNP")) %>% 
  ggplot(aes(x = n_celltype, y = abs(Z))) + geom_point()
