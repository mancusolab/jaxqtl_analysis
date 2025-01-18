## output snplist
specific_eqtl_mash %>% 
  separate(eqtl, into=c("gene", "chr", "pos", "ref", "alt"),sep="_") %>% 
  mutate(pos_1 = as.integer(pos) - 1,
         id = paste0(chr, ":", pos_1, "-", pos)) %>% 
  distinct(id) %>% 
  write_tsv("../result/mvsusie/result/mash_sp_snplist", col_names = F)


#### write rsid ####
specific_eqtl_mash %>% 
  separate(eqtl, into=c("rm", "variant_id"), sep="_", extra = "merge") %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by="variant_id") %>% 
  filter(rsid %in% c("rs1800759", "rs13257737", "rs17186084", "rs13279795", "rs62116613"))
  distinct(rsid) %>% 
  write_tsv("../result/mvsusie/result/mash_sp_rsid", col_names = F)

specific_eqtl_raw %>% 
  distinct(snp) %>% 
  left_join(geno_rsid %>% select(variant_id, rsid), by=c("snp" = "variant_id")) %>% 
  distinct(rsid) %>% 
  write_tsv("../result/mvsusie/result/raw_sp_rsid", col_names = F)

