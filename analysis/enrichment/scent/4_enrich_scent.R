# scent on fine-mapped eQTL

# indir <- "../result/cis/celltype16_tensorqtl_new/all_celltype/"
indir <- "../result/finemap/result_wald/jaxqtl/bed/"

egenes_res <- NB_egenes_df

# order plot by the order NB egenes findings
celltype_order <- NB_egenes_df %>% count(celltype) %>% rename(n_egenes=n)

# gene overlap with annot
gene_annot <- read_tsv("../result/enrich/egenes_cellspecific.overlap.subtract.PLS.new.tsv.gz")

allfiles <- list.files(indir)

scent_celltypes <- c("Bcell", "Tnk", "Myeloid")

sc_df <- data.frame()
for (which_cell in scent_celltypes){
  cell <- paste0(which_cell, ".subtract.PLS.ENCODE")

  sc_res <- allfiles[endsWith(allfiles, cell) & !startsWith(allfiles, "allcells")]
  print(sc_res)
  
  for (file in sc_res){
    empty <- file.size(paste0(indir, file)) == 0L
    if (!empty){
      df <- fread(paste0(indir, file), header = F) %>% 
        left_join(gene_lookup, by = c("V10" = "phenotype_id")) %>% 
        filter(V4 == GeneSymbol) %>% 
        mutate(celltype = sub("\\..*", "", file),
               ref_cell = which_cell)
      sc_df <- bind_rows(sc_df, df) 
      print(nrow(df))
    }
  }
}

dim(sc_df)
sc_df %>% distinct(celltype, V4, V7, V8)

tmp <- sc_df %>% mutate(chr = as.integer(chr)) %>% 
  left_join(gene_annot, by = c("GeneSymbol", "ref_cell", "chr")) %>% 
  mutate(enrich = num_var/hits) %>%
  add_count(celltype) %>% rename(n_hits = n) %>% 
  select(phenotype_id, GeneSymbol, enrich, celltype, n_hits) %>% 
  right_join(egenes_res, by=c("phenotype_id", "celltype")) %>% 
  mutate(enrich = ifelse(is.na(enrich), 0, enrich))

# 132 eQTL has match 
sum(tmp$enrich > 0)
tmp %>% filter(enrich > 0) %>% pull(celltype) %>% table()

tmp %>% group_by(celltype) %>% summarize(mean_enrich=mean(enrich))

n_bootstrap <- 1000
set.seed(2024)
se <- tibble(celltype = unique(tmp$celltype), se = NA, OR_L = NA, OR_U = NA);se

for (i in seq_along(se$celltype)){
  v_enrich <- tmp %>% filter(celltype == se$celltype[i]) %>% pull(enrich)

  bootstrap_mean <- replicate(n_bootstrap, mean(sample(v_enrich, replace = T)))
  se$se[i] <- sqrt(sum((bootstrap_mean - mean(bootstrap_mean))^2)/n_bootstrap)
  
  se$OR_L[i] <- quantile(bootstrap_mean, 0.05)
  se$OR_U[i] <- quantile(bootstrap_mean, 0.95)
}
se

tqtl_enrich <- tmp %>% 
  group_by(celltype) %>% 
  summarize(mean_enrich = mean(enrich), n_hits=mean(n_hits, na.rm = T), n_egenes = n()) %>% 
  ungroup() %>% 
  left_join(se, by="celltype")

OR_L_df <- tqtl_enrich %>% select(celltype, tqtl_enrich=OR_L) %>% 
  inner_join(nb_enrich %>% select(celltype, nb_enrich=OR_L), by="celltype") %>% 
  gather(key = method, value = OR_L, c(tqtl_enrich, nb_enrich))

OR_U_df <- tqtl_enrich %>% select(celltype, tqtl_enrich=OR_U) %>% 
  inner_join(nb_enrich %>% select(celltype, nb_enrich=OR_U), by="celltype") %>% 
  gather(key = method, value = OR_U, c(tqtl_enrich, nb_enrich))

tqtl_enrich %>% select(tqtl_enrich=mean_enrich, tqtl_n_hits = n_hits, tqtl_se=se, celltype) %>%
  inner_join(nb_enrich %>% select(nb_enrich=mean_enrich, nb_n_hist = n_hits, nb_se=se,celltype),
             by="celltype") %>%
  left_join(celltype_order, by="celltype") %>% 
  gather(key = method, value = enrich, c(tqtl_enrich, nb_enrich)) %>%
  left_join(OR_L_df, by = c("method", "celltype")) %>% 
  left_join(OR_U_df, by = c("method", "celltype")) %>% 
  mutate(method = ifelse(method == "nb_enrich", "NegBinom", "tensorqtl")) %>% 
  filter(method == "NegBinom") %>%
  mutate(celltype = fct_reorder(celltype, desc(enrich))) %>% 
  ggplot(aes(x = celltype, y = enrich, fill = method)) + 
  geom_bar(position = "dodge", stat="identity") +
  geom_errorbar(aes(x=celltype, ymin=OR_L, ymax=OR_U, group=method), 
                width=0.2, colour="orange", alpha=0.9, size=0.5, position=position_dodge(0.9))+
  geom_abline(slope=0, intercept = 1, linetype="dashed") + axis_theme +
  theme(axis.text.x = element_text(size=8, angle=45, hjust = 1)) +
  ggtitle("Enrichment of SCENT enhancer-gene")

ggsave2(paste0(plots_dir, "enrich/NB_cisgenes_scent_enrich_allcelltypes.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)

summary(nb_enrich$mean_enrich)
summary(tqtl_enrich$mean_enrich)

# number of eQTL has hits
bind_rows(nb_enrich %>% mutate(method = "jaxqtl_NB"),
          tqtl_enrich %>% mutate(method = "tqtl")) %>% 
  mutate(celltype = fct_reorder(celltype, desc(n_egenes))) %>%
  ggplot(aes(x = celltype, y = n_hits/n_egenes, fill = method)) + 
  geom_bar(position = "dodge", stat="identity") +
  geom_abline(slope=0, intercept = 1, linetype="dashed") + axis_theme +
  theme(axis.text.x = element_text(size=8, angle=45, hjust = 1)) +
  ggtitle("Number of eQTLs hits in SCENT enhancer-gene")

ggsave2(paste0(plots_dir, "jaxqtl_tqtl_cisgenes_scent_eqtlhits_allcelltypes.png"),
        width = full_width, height = full_height, units = 'cm', dpi = 300)


# meta analysis
nb_cs_both
t.test(nb_cs_both$pip, tqtl_cs_both$pip)

lm(pip ~ abs(tss_dist) + method, 
   data=bind_rows(nb_cs_both %>% mutate(method = "NegBinom"),
          tqtl_cs_both %>% mutate(method = "tensorqtl")) %>% 
     mutate(tss_dist = abs(pos - tss))) %>% 
  broom::tidy()

meta_dat <- nb_enrich %>% 
  mutate(TE = log(mean_enrich), 
         lower = log(OR_L),
         upper = log(OR_U),
         seTE = NA)
summary(nb_enrich$mean_enrich)

m.gen_bin <- metagen(TE = TE,
                     seTE = seTE,
                     lower = lower,
                     upper = upper,
                     studlab = celltype,
                     data = meta_dat,
                     sm = "OR",
                     method.tau = "PM",
                     fixed = FALSE,
                     random = TRUE,
                     title = "Enrichment (Pre-calculated)")
metares <- summary(m.gen_bin); metares
metares$pval.random
# new results
# tqtl: 2.2438 [1.6088; 3.1293] 4.76 < 0.0001
# jqtl: 2.1920 [1.6678; 2.8811]
