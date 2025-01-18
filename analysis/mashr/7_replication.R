## replication 

#### eQTLGen ####
joined <- snp_df %>%
  filter(eqtl %in% mash_sig_eqtls) %>%
  separate(snp, into=c("rm1", "rm2", "ref", "alt"), remove = F) %>% 
  mutate(variant_id = paste0(chr, ":", pos)) %>% 
  count(phenotype_id, snp, variant_id, ref, alt) %>% 
  left_join(eqtlgen %>% dplyr::select(phenotype_id, variant_id, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
            by=c("phenotype_id", "variant_id"))

shared_mag_df_tmp <- shared_mag_df %>% gather(key = celltype, value = share, B_IN:Plasma) %>% 
  filter(share == TRUE)
joined <- shared_mag_df_tmp %>% 
  group_by(eqtl) %>% 
  summarize(n = sum(share)) %>% ungroup() %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_") %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         variant_id = paste0(chr, ":", pos)) %>% 
  left_join(eqtlgen %>% dplyr::select(phenotype_id, variant_id, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
            by=c("phenotype_id", "variant_id"))

aligns <- allele.qc(joined$ref,joined$alt,joined$OtherAllele, joined$AssessedAllele)

joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            found = sum(!is.na(Zscore_eqtlgen))) %>% 
  ungroup() 

p_mash <- joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            found = sum(!is.na(Zscore_eqtlgen))) %>% 
  ungroup() %>% 
  mutate(group=ifelse(n == 1, "specific", "shared")) %>% 
  group_by(group) %>% 
  summarize(found = sum(found),
            n_ct = sum(n_ct)) %>% 
  mutate(rep = found / n_ct,
         not_rep = 1 - rep) %>% 
  gather(key = eQTLGen, value = frac, rep:not_rep) %>% 
  ggplot(aes(fill=eQTLGen, y=frac, x=group)) +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("rep"="#17BECF", "not_rep"="grey"))+
  ggtitle("mash estimates")+xlab("")+ylab("fraction of eQTLs")+
  axis_theme+theme(legend.position = "bottom")

prop.test(c(113, 1148), c(189, 1488)) %>% broom::tidy()

legend <- get_plot_component(p_mash, 'guide-box-bottom')
legend <- get_legend(p_mash)

grid <- plot_grid(
  p_raw+theme(legend.position = "none"), 
  p_mash+theme(legend.position = "none"), 
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1,
  rel_widths = c(1, 1)
)

plot_grid(
  grid, legend, 
  align = 'h',
  nrow = 2,
  axis="b",
  rel_heights = c(1, 0.1)
)
ggsave2(paste0(plots_dir, "mash_raw_replication_eqtlgen.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


#### GTEx ####
# 
snp_hg38 <- read_tsv("./eqtl_pip0.5.hg38.bed", F) %>%  # eqtl_pip0.5.hg38.bed, sp_snp_hg38
  select(chr=X1, pos_38=X3, X4) %>% 
  mutate(X4 = gsub(".*:", "", X4),
         chr = as.integer(gsub("chr", "", chr))) %>% 
  separate(X4, into=c("rm", "pos_19"), sep="-") %>% select(-rm) %>% 
  mutate(pos_19 = as.integer(pos_19))

joined <- snp_df %>% filter(eqtl %in% mash_sig_eqtls) %>% 
  count(eqtl, phenotype_id, chr, pos) %>% 
  separate(eqtl, into=c("rm1", "rm2", "rm3", "ref", "alt"), sep="_", remove=F) %>% 
  left_join(snp_hg38, by=c("chr", "pos"="pos_19")) %>% 
  left_join(gtex %>% select(phenotype_id, chr, pos_38, gtex_ref=ref, gtex_alt=alt, slope_gtex), 
            by=c("phenotype_id", "chr", "pos_38"))

# shared_mag_df_tmp
joined <- shared_mag_df_tmp %>% 
  group_by(eqtl) %>% 
  summarize(n = sum(share)) %>% ungroup() %>% 
  separate(eqtl, into=c("phenotype_id", "chr", "pos", "ref", "alt"), sep="_", remove=F) %>% 
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos_19 = as.integer(pos)) %>% 
  left_join(snp_hg38, by=c("chr", "pos_19")) %>% 
  left_join(gtex %>% select(phenotype_id, chr, pos_38, gtex_ref=ref, gtex_alt=alt, slope_gtex), 
            by=c("phenotype_id", "chr", "pos_38"))

aligns <- allele.qc(joined$ref,joined$alt,joined$gtex_ref, joined$gtex_alt)
sum(aligns$strand_flip)
sum(!aligns$keep) # keep all snps

joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip) %>% 
  filter(keep == TRUE) %>% # nrow()
  mutate(slope_gtex = ifelse((strand_flip == TRUE & alt != gtex_alt) | (sign_flip == TRUE & alt != gtex_alt), 
                             slope_gtex * -1, slope_gtex)) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            found = sum(!is.na(slope_gtex))) %>% 
  ungroup() 

prop.test(c(63, 699), c(189, 1488)) %>% broom::tidy()

# replicate at eQTL level: 
# 0.38
p_raw <- joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip) %>% 
  filter(keep == TRUE) %>% # nrow()
  mutate(slope_gtex = ifelse((strand_flip == TRUE & alt != gtex_alt) | (sign_flip == TRUE & alt != gtex_alt), 
                             slope_gtex * -1, slope_gtex)) %>% 
  group_by(n) %>%
  summarize(n_ct = n(),
            found = sum(!is.na(slope_gtex))) %>% 
  ungroup() %>% 
  mutate(group=ifelse(n == 1, "specific", "shared")) %>% 
  group_by(group) %>% 
  summarize(found = sum(found),
            n_ct = sum(n_ct)) %>% 
  mutate(rep = found / n_ct,
         not_rep = 1 - rep) %>% 
  gather(key = GTEx, value = frac, rep:not_rep) %>% 
  ggplot(aes(fill=GTEx, y=frac, x=group)) +
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("rep"="#17BECF", "not_rep"="grey"))+
  ggtitle("raw estimates")+xlab("")+ylab("fraction of eQTLs")+
  axis_theme+theme(legend.position = "bottom")

legend <- get_legend(p_mash)

grid <- plot_grid(
  p_raw+theme(legend.position = "none"), 
  p_mash+theme(legend.position = "none"), 
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="b",
  labels = c("A", "B"), label_size = 10, vjust=1,
  rel_widths = c(1, 1)
)

plot_grid(
  grid, legend, 
  align = 'h',
  nrow = 2,
  axis="b",
  rel_heights = c(1, 0.1)
)
ggsave2(paste0(plots_dir, "mash_raw_replication_gtex.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)
