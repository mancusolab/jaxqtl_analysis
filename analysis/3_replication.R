# !!! Warning:
# careful with rsid of eqtls from catalogues (outdated or merged with another SNP)

##### bulk PBMC OneK1K findings #####
celltype_prop <- ind_celltype_ct %>% 
  group_by(individual) %>% 
  mutate(n = sum(counts)) %>% ungroup() %>% 
  mutate(prop = counts / n) %>% 
  group_by(celltype) %>% summarise(mean_prop = mean(prop)) %>% ungroup()

# replicate eQTL: 4788 / 18907 = 0.253
NB_egenes_df %>% select(phenotype_id, snp, celltype) %>% inner_join(jaxqtl_allcell_egenes_df %>% select(phenotype_id, snp), by=c("phenotype_id", "snp")) %>% nrow()

bind_rows(NB_egenes_df %>% select(phenotype_id, snp, celltype) %>% inner_join(jaxqtl_allcell_egenes_df %>% select(phenotype_id, snp), by=c("phenotype_id", "snp")) %>% mutate(success = 1),
          NB_egenes_df %>% select(phenotype_id, snp, celltype) %>% anti_join(jaxqtl_allcell_egenes_df %>% select(phenotype_id, snp), by=c("phenotype_id", "snp")) %>% mutate(success = 0)) %>% 
  left_join(celltype_prop, by="celltype") %>%
  group_by(celltype, mean_prop) %>% 
  summarize(rep_rate = mean(success)) %>% 
  ggplot(aes(x = mean_prop, y = rep_rate)) + geom_point()


##### eqtlgen findings #####

# external findings
# whole blood eqtls from multiple studies (FDR < 0.05); biallelic SNPs, distinct snp
# multi-ancestry meta-analysis
eqtlgen <- fread("../data/OneK1K/eQTLs/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>% 
  mutate(variant_id = paste0(SNPChr, ":", SNPPos)) %>% 
  dplyr::rename(phenotype_id = Gene)

eqtlgen %>% filter(GeneSymbol == "CD40") %>% arrange(desc(abs(Zscore))) %>% View
eqtlgen %>% filter(GeneSymbol == "IL6ST" & SNPChr == 5 & SNPPos == 55444683)
head(eqtlgen)
dim(eqtlgen)

# tqtl_egenes_df, jqtl_lm_score_egenes_df, NB_egenes_df, nb_
res <- jqtl_lm_score_egenes_df; nrow(res)

# 14408 join by position 
joined <- res %>%
  # distinct(phenotype_id, snp, variant_id, .keep_all = T) %>% 
  left_join(eqtlgen %>% select(phenotype_id, variant_id, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
             by=c("phenotype_id", "variant_id"))

aligns <- allele.qc(joined$ref,joined$alt,joined$OtherAllele, joined$AssessedAllele)
sum(aligns$keep)

# replicate at eQTL level: 
# NB: 14406/18907, 76%; 12967/18907, 69%
# linear score: 12765/16654, 75%; effect consistent: 11473/16654, 69%
# for ambiguous snps, sign_flip will be TRUE, so need to correct this
joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  filter(keep == TRUE & !is.na(Zscore_eqtlgen)) %>% #nrow()
  mutate(Zscore_eqtlgen = ifelse(sign_flip == TRUE, Zscore_eqtlgen * -1, Zscore_eqtlgen)) %>% 
  filter(Zscore_eqtlgen * slope > 0) %>% #nrow() 
  View
prop.test(c(14406, 12765), c(18907, 16654)) %>% broom::tidy()
prop.test(c(12967, 11473), c(18907, 16654)) %>% broom::tidy()

joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  filter(keep == TRUE) %>% 
  mutate(Zscore_eqtlgen = ifelse(sign_flip == TRUE, -1 * Zscore_eqtlgen, Zscore_eqtlgen)) %>% 
  ggplot(aes(x = Z, y = Zscore_eqtlgen)) + geom_point() + axis_theme

joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip,
                  sign_flip = ifelse(ref == OtherAllele & alt == AssessedAllele, FALSE, sign_flip)) %>% 
  filter(keep == TRUE) %>% 
  mutate(Zscore_eqtlgen = ifelse(sign_flip == TRUE, -1 * Zscore_eqtlgen, Zscore_eqtlgen)) %>% 
  filter(!is.na(Zscore_eqtlgen)) %>% 
  group_by(phenotype_id, snp) %>%
  summarize(n_celltype = n(),
            n_concordant = sum(Zscore_eqtlgen * Z > 0),
            n_discordant = sum(Zscore_eqtlgen * Z < 0)) %>%
  ggplot(aes(x = n_celltype, y = n_concordant)) + geom_point()


# replicate at eGene level: 
# NB: ;6010/6934, 17191/18907
# tqtl: 5491/6282; 14908/16390
NB_egenes_df %>% distinct(celltype, phenotype_id) %>% 
  inner_join(eqtlgen %>% distinct(phenotype_id), by="phenotype_id") %>% nrow()
nrow(NB_egenes_df)

# rs2246199, ENSG00000066379
tmp <- res %>% 
  left_join(bind_rows(success1, success2) %>% 
              select(celltype, phenotype_id, variant_id) %>% 
              mutate(replicate = 1)) %>% 
  mutate(replicate = ifelse(is.na(replicate), 0, replicate)) %>% 
  add_count(phenotype_id, variant_id, ref, alt) %>% 
  group_by(n) %>% 
  summarize(rep_rate = mean(replicate)) %>% 
  ungroup()
  
broom::tidy(lm(rep_rate ~ n, data = tmp))

tmp %>% 
  ggplot(aes(x = n, y = rep_rate)) + 
  geom_point() + geom_smooth(method = "lm") + axis_theme +
  ggtitle("mean replication rate of eQTLs vs. its number of \n cell types")

ggsave2(paste0(plots_dir, "NB_rep_rate_ncelltypes.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

##### GTEx findings #####

# GTEx v8 EUR, gene specific threshold
gtex <- fread("../data/OneK1K/eQTLs/GTEx/GTEx_Analysis_v8_eQTL_EUR/eqtls/Whole_Blood.v8.EUR.signif_pairs.txt.gz") %>%
#gtex <- fread("../data/OneK1K/eQTLs/GTEx/GTEx_Analysis_v8_eQTL_EUR/eqtls/Whole_Blood.v8.EUR.egenes.txt.gz") %>%
  mutate(phenotype_id = gsub("\\..*", "", phenotype_id),
         variant_id = gsub("_b38", "", variant_id)) %>% 
  separate(variant_id, into=c("chr", "pos_38", "ref", "alt"), sep="_") %>%
  mutate(chr = as.integer(gsub("chr", "", chr)),
         pos_38 = as.integer(pos_38)) %>% 
  select(phenotype_id:alt, slope_gtex = slope) %>%
  filter(!is.na(chr))

gtex %>% filter(phenotype_id == "ENSG00000101017")
gtex %>% filter(phenotype_id == "ENSG00000134352")

res <- NB_egenes_df # NB_egenes_df, jqtl_lm_score_egenes_df

joined <- res %>%
  inner_join(gtex %>% select(phenotype_id, chr, pos_38, gtex_ref=ref, gtex_alt=alt, slope_gtex), 
             by=c("phenotype_id", "chr", "pos_38"))

aligns <- allele.qc(joined$ref,joined$alt,joined$gtex_ref, joined$gtex_alt)
sum(aligns$strand_flip)
sum(!aligns$keep) # keep all snps

# replicate at eQTL level: 
# NB: 9225/18907, 49%, effect consistent: 8607/18907, 46%
# tqtl: 8198/16654, 50%; effect consistent: 7640/16654, 46%
joined %>% mutate(keep = aligns$keep,
                  sign_flip = aligns$sign_flip,
                  strand_flip = aligns$strand_flip) %>% 
  filter(keep == TRUE) %>% # nrow()
  mutate(slope_gtex = ifelse((strand_flip == TRUE & alt != gtex_alt) | (sign_flip == TRUE & alt != gtex_alt), 
                             slope_gtex * -1, slope_gtex)) %>% 
  filter((slope * slope_gtex)>0) %>% nrow()

prop.test(c(9225, 8198), c(18907, 16654)) %>% broom::tidy()
prop.test(c(8607, 7640), c(18907, 16654)) %>% broom::tidy()

# replicate eGene:
# 17118/18907; 6030/6934
# 14802/16390 ; 5474/6282
NB_egenes_df %>% distinct(celltype, phenotype_id) %>% 
  inner_join(gtex %>% distinct(phenotype_id), by="phenotype_id")

##### eQTL catalogue replication #####

lm_rep <- read_tsv(glue("{jqtlcis_dir}/jaxqtl_lm_leadsnp.CT_rep.tsv"))
nb_rep <- read_tsv(glue("{jqtlcis_dir}/jaxqtl_nb_leadsnp.CT_rep.tsv"))

ct_rep <- lm_rep %>% rename(`jaxQTL-Linear` = rate, CT = group) %>% 
  mutate(`jaxQTL-NegBinom` = nb_rep$rate)
ct_rep

ct_rep %>% 
  gather(key = model, value = replication, `jaxQTL-Linear`:`jaxQTL-NegBinom`) %>% 
  mutate(CT = fct_reorder(CT, desc(replication))) %>% 
  ggplot(aes(x = CT, y = replication, fill = model)) + geom_col(position = "dodge") +
  scale_fill_manual(values = method_cols) +
  xlab("Cell type group")+
  axis_theme+
  theme(legend.position = "bottom", axis.text.x = element_text(angle=30, hjust = 1),
        legend.key.size = unit(0.2, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
  
ggsave2(paste0(plots_dir, "sc_replication.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

##### yazar findings #####
# count genes test in each cell type
yazar_genes_pass <- readxl::read_excel("../data/OneK1K/yazar_result/science.abf3041_tables_s6_to_s19.xlsx", 
                            sheet = 3, skip = 2) %>% 
  gather(key = cell, value = gene) %>% 
  mutate(cell = gsub(" ", "_", cell)) %>% 
  drop_na() %>% 
  count(cell, name = "gene_pass") %>% 
  rename(celltype=cell)
yazar_genes_pass %>% arrange(gene_pass)

# alternative allele is the effect allele
yazar <- readxl::read_excel("../data/OneK1K/yazar_result/science.abf3041_tables_s6_to_s19.xlsx", 
                            sheet = 5, skip = 2) %>% 
  janitor::clean_names() %>% 
  filter(fdr < 0.05) %>% 
  rename(phenotype_id = gene_ensembl_id) %>% 
  mutate(cell_type = gsub(" ", "_", cell_type))
head(yazar)

yazar %>% filter(gene_id == "IL6ST")
res <- NB_egenes_df # jqtl_lm_score_egenes_df # NB_egenes_df, # NB_egenes_df_0.1

joined <- res %>%
  inner_join(yazar %>% 
               filter(e_snp_rank == "eSNP1") %>%
               select(celltype=cell_type, phenotype_id, chr=chromosome, pos=position, snp_assessed_allele, rho_correlation_coefficient), 
             by=c("celltype", "phenotype_id", "chr", "pos", "alt"="snp_assessed_allele"))
dim(joined)

joined %>% rename(slope_NegBinom=slope_wald) %>%
  ggplot(aes(x = slope_NegBinom, y = rho_correlation_coefficient)) + 
  geom_point(size=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + axis_theme

ggsave2(paste0(plots_dir, "NB_rep_onek1k_eSNP1.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)

# replicate at eQTL (eSNP) level: 
# NB: 6336/18907, 34%
# linear: 7366/16654, 44%

# replicate eGene:
# 12566/18907; 66%
# 12545/16654: 75%

#NB_egenes_df %>%
jqtl_lm_score_egenes_df %>%
  inner_join(yazar %>% distinct(phenotype_id, cell_type), by=c("phenotype_id", "celltype"="cell_type")) %>% 
  nrow()


############ cell type specific replication ############

# DICE: FACS-sorted immune cells, Adj. P value < 0.05, raw P value < 0.0001 and TPM > 1.0
dice_map <- read_tsv("../data/OneK1K/eQTLs/DICE/celltype_map.tsv")

allsuccess <- tibble(lineage = unique(cellmeta$lineage), eqtls = NA, found = NA,
                     external_eqtls = NA)
files <- dice_map %>% filter(yazar_celltype == cell) %>% pull(files)

success <- data.frame()
egene_res <- NB_egenes_df

for (idx in seq_along(allsuccess$lineage)){
  cell <- allsuccess$lineage[idx]
  print(cell)
  files <- dice_map %>% filter(yazar_celltype == cell) %>% pull(files)
  
  if (length(files)>0){
    tmp <- data.frame()
    for (file in files){
      dice <- fread(paste0("../data/OneK1K/eQTLs/DICE/", file), skip="#CHROM") %>% 
        separate(INFO, into=c("phenotype_id", "GeneSymbol", "Pval", "Beta"), sep=";") %>% 
        mutate(phenotype_id = gsub("Gene=", "", phenotype_id),
               Pval = as.numeric(gsub("Pvalue=", "", Pval)))
      tmp <- bind_rows(tmp, dice)
    }
    res <- egene_res %>% filter(lineage == cell)
    
    allsuccess$eqtls[idx] <- nrow(res)
    
    success <- res %>%
      inner_join(tmp %>% distinct(ID, phenotype_id), by=c("phenotype_id", "rsid"="ID"))
    
    allsuccess$found[idx] <- nrow(success) 
    allsuccess$external_eqtls[idx] <- nrow(tmp %>% distinct(ID, phenotype_id))
  }
}

allsuccess_nb <- allsuccess
allsuccess_tqtl <- allsuccess

allsuccess_nb %>% mutate(rep_rate_external = found/external_eqtls,
                         rep_rate = found/eqtls)
allsuccess_tqtl %>% mutate(rep_rate_external = found/external_eqtls,
                           rep_rate = found/eqtls)


###### Nathan 2022 (hg38) ####

rep_celltype <- c("CD4_memory", "CD8_T", "T_memory")
indir <- "../data/OneK1K/eQTLs/sceqtl/Nathan_2022/"
rep_res <- data.frame(celltype=rep_celltype, yazar_celltype=c("CD4", "CD8", NA), num_eqtl=NA, found = NA)

res <- tqtl_egenes_df
for (i in 1:2){
  nathan_celltype <- rep_res$celltype[i]
  match_celltype <- rep_res$yazar_celltype[i]
  nathan <- fread(paste0(indir, "Nathan2022Nature_", nathan_celltype, ".sig_qtl.tsv"), header=T)
  
  tmp <- res %>% filter(lineage %in% c(match_celltype))
  success <- tmp %>% 
    left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by = "phenotype_id") %>% 
    inner_join(nathan, by=c("GeneSymbol"="geneName", "rsid"="variantId"))
  
  rep_res$num_eqtl[i] <- nrow(tmp)
  rep_res$found[i] <- nrow(success)
}
rep_res

nathan <- fread(paste0(indir, "Nathan2022Nature_T_memory.sig_qtl.tsv"), header=T)
tmp <- res %>% filter(lineage %in% c("CD4", "CD8"))
success <- tmp %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by = "phenotype_id") %>% 
  inner_join(nathan, by=c("GeneSymbol"="geneName", "rsid"="variantId"))

rep_res$num_eqtl[3] <- nrow(tmp)
rep_res$found[3] <- nrow(success)
rep_res %>% mutate(rep_rate = found/num_eqtl)


###### sceqtlgen (hg19) #######

rep_celltype <- c("b-cells", "c-mono", "dend", "mono","nc-mono","nk","t_cd4","t_cd8")

rep_res <- data.frame(celltype=c("t_cd4", "nk", "t_cd4", "t_cd8", "t_cd8", "b-cells","t_cd8", "b-cells",
                                 "nk", "nc-mono", "mono", "dend", "b-cells", "t_cd4"), 
                      yazar_celltype=celltypes, num_eqtl=NA, found = NA)

res <- tqtl_egenes_df
for (i in 1:nrow(rep_res)){
  sceqtlgen_celltype <- rep_res$celltype[i]
  match_celltype <- rep_res$yazar_celltype[i]
  sceqtlgen <- fread(paste0("../data/OneK1K/eQTLs/sceqtlgen/genome-wide-eQTL_", sceqtlgen_celltype, ".txt")) %>% 
    filter(FDR < 0.05) %>% 
    rename(phenotype_id = GeneName) %>% 
    mutate(variant_id = paste0(SNPChr, ":", SNPChrPos))
  nrow(sceqtlgen)
  
  tmp <- res %>% filter(celltype == match_celltype); nrow(tmp)
  # success <- tmp %>% 
  #   inner_join(sceqtlgen, by=c("phenotype_id", "variant_id"))
  success <- tmp %>% 
    inner_join(sceqtlgen %>% distinct(phenotype_id), by=c("phenotype_id"))
  rep_res$num_eqtl[i] <- nrow(tmp)
  rep_res$found[i] <- nrow(success)
  nrow(success)
}

rep_res_tqtl <- rep_res
View(rep_res)

# sceqtl base

indir <- "../data/OneK1K/eQTLs/sceqtl/"
map <- readxl::read_excel(paste0(indir, "rep_map.xlsx"))

rep_res <- data.frame(yazar_celltype=celltypes, num_eqtl=NA, found = NA)

res <- yazar %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by=c("phenotype_id")) %>% 
  rename(celltype = cell_type, rsid=snp)

for (i in seq_along(celltypes)){
  cell <- rep_res$yazar_celltype[i]
  map_tmp <- map %>% filter(grepl(cell, yazar_celltype))
  print(cell)
  ext_eqtl <- data.frame()
  
  if (nrow(map_tmp) > 0){
    for (file in map_tmp$filename){
      onedf <- fread(paste0(indir, file), header=T) %>% mutate(chrom = as.character(chrom)) %>% 
        filter(FDR < 0.05)
      print(nrow(onedf))
      if (nrow(onedf) > 0){
        ext_eqtl <- bind_rows(ext_eqtl, onedf)
        print(paste0("found eQTLs ", nrow(res %>% filter(celltype == cell) %>% inner_join(ext_eqtl %>% distinct(geneName, variantId), by=c("GeneSymbol" = "geneName","rsid"="variantId")))))
      }
    }
    ext_eqtl <- ext_eqtl %>% distinct(geneName, variantId)
    
    res_tmp <- res %>% filter(celltype == cell)
    found <- nrow(res_tmp %>% inner_join(ext_eqtl, by=c("GeneSymbol" = "geneName","rsid"="variantId")))
    
    rep_res$num_eqtl[i] <- nrow(res_tmp)
    rep_res$found[i] <- found
  }
  
}

sum(rep_res$found)
sum(rep_res$found)/sum(rep_res$num_eqtl) # 22.9% vs. 24.7% vs. 20.45%


# build specific genes (warning: remove these genes at interpretation)
# ref: https://www-cell-com.libproxy1.usc.edu/ajhg/fulltext/S0002-9297(24)00168-X
build_specific_genes <- readxl::read_excel("../data/others/build_specific_genes_Ungar24.xlsx", sheet = 1, skip = 1)
table(build_specific_genes$annotation_specific)

h19_genes <- build_specific_genes %>% filter(annotation_specific=="hg19_spec_gene") %>% 
  distinct(gene_id) %>% pull(gene_id)

nb_only %>% 
  mutate(is_hg19_specific=ifelse(phenotype_id %in% h19_genes, "yes", "no")) %>% 
  group_by(is_hg19_specific) %>% 
  summarize(mean_express = mean(express_percent))

both_hits %>% filter(phenotype_id %in% h19_genes)
