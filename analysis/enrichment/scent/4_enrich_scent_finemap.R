# scent on fine-mapped eQTL

# indir <- "../result/cis/celltype16_tensorqtl_new/all_celltype/"
compare_method <- FALSE
suffix <- "cs" # cs, PIP0.5, PIP0.8
pip_threshold <- 0

# for tqtl: set method = "jqtl"
for (method in c("jaxqtl", "jaxqtl_lm_score")){
  
  indir <- paste0("../result/finemap/result_wald/", method, "/bed/")
  if (method == "jaxqtl"){
    egenes_res <- NB_egenes_df
    prefix <- "nb"
  }else if (method == "tqtl"){
    egenes_res <- tqtl_egenes_df
    prefix <- "tqtl"
  }else if (method == "jaxqtl_lm_score"){
    egenes_res <- jqtl_lm_score_egenes_df
    prefix <- "jaxqtl_lm_score"
  }

  # gene overlap with annot for all eGenes
  if (method == "tqtl"){
    gene_annot <- read_tsv("../result/enrich/egenes_cellspecific.overlap.subtract.PLS.new.tsv.gz")
    gene_annot_bulk <- read_tsv("../result/enrich/egenes_bulk.overlap.subtract.PLS.new.tsv.gz")
    gene_annot <- bind_rows(gene_annot, gene_annot_bulk) %>% distinct()
  }else{
    gene_annot <- read_tsv("../result/enrich/egenes_allmodels.overlap.subtract.PLS.tsv.gz")
  }
  
  # finemapped eGenes: count number of variants in CS, or count number of variant > PIP_threshold
  if (pip_threshold > 0){
    gene_finemap <- read_tsv(paste0(indir, "/allcelltype_finemap_pip", pip_threshold, ".tsv.gz"), F) %>% 
      mutate(celltype = gsub(".cs.bed", "", X7)) %>% 
      rename(phenotype_id = X6, pip = X5, snp = X4) %>% 
      distinct(celltype, phenotype_id, snp, pip) %>% 
      group_by(celltype, phenotype_id) %>% 
      summarize(n_causal = sum(pip >= pip_threshold)) %>% 
      ungroup() 
  }else{
    gene_finemap <- read_tsv(paste0(indir, "/allcelltype_finemap.cs.tsv.gz"), F) %>% 
      mutate(celltype = gsub(".cs.bed", "", X7)) %>% 
      rename(phenotype_id = X6, pip = X5, snp = X4) %>% 
      distinct(celltype, phenotype_id, snp, pip) %>% 
      group_by(celltype, phenotype_id) %>% 
      summarize(n_causal = sum(pip >= pip_threshold)) %>% 
      ungroup() 
  }
  table(gene_finemap$celltype)
  
  # keep genes that are finemapped (i.e., remove MHC and MAPT genes)
  allegenes <- read_tsv(paste0(indir, prefix, "_allegenes.tsv"), F) %>% 
    select(phenotype_id = X1, celltype = X3) %>% 
    filter(!phenotype_id %in% c(MHC_genes, MAPT_genes))
  
  if (compare_method){
    allegenes <- allegenes %>%  filter(celltype != "allcells") %>% 
      inner_join(both_hits %>% select(celltype, phenotype_id), by=c("celltype", "phenotype_id"))
  }
  
  allfiles <- list.files(paste0(indir, "scent/"))
  
  scent_celltypes <- c("Bcell", "Tnk", "Myeloid")
  
  sc_df <- data.frame()
  for (which_cell in scent_celltypes){
    cell <- paste0(which_cell, ".subtract.PLS.ENCODE")
    
    sc_res <- allfiles[endsWith(allfiles, cell) & !startsWith(allfiles, "allcells")]
    print(sc_res)
    sc_res <- sc_res[grepl(suffix, sc_res)]
    
    for (file in sc_res){
      empty <- file.size(paste0(indir, "scent/", file)) == 0L
      if (!empty){
        df <- fread(paste0(indir,  "scent/", file), header = F) %>% rename(phenotype_id=V11) %>% 
          left_join(gene_lookup, by = c("phenotype_id")) %>% 
          filter(V4 == GeneSymbol) %>% 
          mutate(celltype = sub("\\..*", "", file),
                 ref_cell = which_cell) %>% 
          left_join(gene_finemap, by=c("celltype", "phenotype_id"))
        
        if (compare_method){
          df <- df %>% inner_join(both_hits %>% select(celltype, phenotype_id), 
                                  by=c("celltype", "phenotype_id"))
        }
        sc_df <- bind_rows(sc_df, df) 
        print(nrow(df))
      }
    }
  }
  
  dim(sc_df)
  
  sc_df <- sc_df %>% 
    mutate(chr = as.integer(chr)) %>% 
    add_count(celltype, phenotype_id, name="causal_hit") %>% 
    left_join(gene_annot, by = c("GeneSymbol", "ref_cell", "chr", "phenotype_id")) %>% 
    mutate(enrich = (causal_hit / n_causal) / (hits / num_var)) %>% 
    select(phenotype_id, GeneSymbol, enrich, celltype) %>% 
    right_join(allegenes, by=c("phenotype_id", "celltype")) %>% 
    mutate(enrich = ifelse(is.na(enrich) | !is.finite(enrich), 0, enrich))
  
  sum(sc_df$enrich > 0)

  # by cell type
  n_bootstrap <- 1000
  set.seed(2024)
  se <- tibble(celltype = unique(sc_df$celltype), se = NA, OR_L = NA, OR_U = NA)
  
  for (i in seq_along(se$celltype)){
    v_enrich <- sc_df %>% filter(celltype == se$celltype[i]) %>% pull(enrich)
    
    bootstrap_mean <- replicate(n_bootstrap, mean(sample(v_enrich, replace = T)))
    se$se[i] <- sqrt(sum((bootstrap_mean - mean(bootstrap_mean))^2)/n_bootstrap)
    
    se$OR_L[i] <- quantile(bootstrap_mean, 0.025)
    se$OR_U[i] <- quantile(bootstrap_mean, 0.975)
  }
  se
  
  sc_df <- sc_df %>% 
    group_by(celltype) %>% 
    summarize(mean_enrich = mean(enrich)) %>% 
    ungroup() %>% 
    left_join(se, by="celltype")
  
  
  if (method == "jaxqtl"){
    assign("nb_enrich", sc_df)
  }else if (method == "tqtl"){
    assign("tqtl_enrich", sc_df)
  }else if (method == "jaxqtl_lm_score"){
    assign("lm_score_enrich", sc_df)
  }
}

# allcells results
nb_enrich
lm_score_enrich

for (method in c("jaxqtl", "jaxqtl_lm_score")){
  
  indir <- paste0("../result/finemap/result_wald/", method, "/bed/")
  if (method == "jaxqtl"){
    egenes_res <- NB_egenes_df
    prefix <- "nb"
  }else if (method == "tqtl"){
    egenes_res <- tqtl_egenes_df
    prefix <- "tqtl"
  }else if (method == "jaxqtl_lm_score"){
    egenes_res <- jqtl_lm_score_egenes_df
    prefix <- "jaxqtl_lm_score"
  }
  
  # gene overlap with annot for all eGenes
  if (method == "tqtl"){
    gene_annot <- read_tsv("../result/enrich/egenes_cellspecific.overlap.subtract.PLS.new.tsv.gz")
    gene_annot_bulk <- read_tsv("../result/enrich/egenes_bulk.overlap.subtract.PLS.new.tsv.gz")
    gene_annot <- bind_rows(gene_annot, gene_annot_bulk) %>% distinct()
  }else{
    gene_annot <- read_tsv("../result/enrich/egenes_allmodels.overlap.subtract.PLS.tsv.gz")
  }
  
  # finemapped eGenes
  if (pip_threshold > 0){
    gene_finemap <- read_tsv(paste0(indir, "/allcelltype_finemap_pip", pip_threshold, ".tsv.gz"), F) %>% 
      mutate(celltype = gsub(".cs.bed", "", X7)) %>% 
      rename(phenotype_id = X6, pip = X5, snp = X4) %>% 
      distinct(celltype, phenotype_id, snp, pip) %>% 
      group_by(celltype, phenotype_id) %>% 
      summarize(n_causal = sum(pip >= pip_threshold)) %>% 
      ungroup() %>% 
      filter(celltype == "allcells")
  }else{
  gene_finemap <- read_tsv(paste0(indir, "/allcelltype_finemap.cs.tsv.gz"), F) %>% 
    mutate(celltype = gsub(".cs.bed", "", X7)) %>% 
    rename(phenotype_id = X6, pip = X5, snp = X4) %>% 
    distinct(celltype, phenotype_id, snp, pip) %>% 
    group_by(celltype, phenotype_id) %>% 
    summarize(n_causal = sum(pip >= pip_threshold)) %>% 
    ungroup() %>% 
    filter(celltype == "allcells")
  }
  table(gene_finemap$celltype)
  
  # all egenes (14228)
  
  # keep genes that are finemapped (i.e., remove MHC and MAPT genes)
  allegenes <- read_tsv(paste0(indir, prefix, "_allegenes.tsv"), F) %>% 
    select(phenotype_id = X1, celltype = X3) %>% 
    filter(!phenotype_id %in% c(MHC_genes, MAPT_genes)) %>% 
    filter(celltype == "allcells")
  
  allfiles <- list.files(paste0(indir, "scent/"))
  allfiles <- allfiles[grepl(suffix, allfiles)]

  scent_celltypes <- c("Bcell", "Tnk", "Myeloid")
  
  # all cell results (don't need to subset to both egenes)
  allcell_df <- data.frame()
  for (which_cell in scent_celltypes){
    cell <- paste0(which_cell, ".subtract.PLS.ENCODE")
    
    sc_res <- allfiles[endsWith(allfiles, cell) & startsWith(allfiles, "allcells")]
    print(sc_res)
    
    for (file in sc_res){
      empty <- file.size(paste0(indir, "scent/", file)) == 0L
      if (!empty){
        df <- fread(paste0(indir, "scent/", file), header = F) %>% rename(phenotype_id=V11) %>% 
          left_join(gene_lookup, by = c("phenotype_id")) %>% 
          filter(V4 == GeneSymbol) %>% 
          mutate(celltype = sub("\\..*", "", file),
                 ref_cell = which_cell) %>% 
          left_join(gene_finemap, by=c("celltype", "phenotype_id"))
        
        allcell_df <- bind_rows(allcell_df, df) 
        print(nrow(df))
      }
    }
  }
  
  allcell_df
  
  allcell_df <- allcell_df %>% 
    add_count(phenotype_id, name="causal_hit") %>% 
    left_join(gene_annot, by = c("GeneSymbol", "ref_cell", "phenotype_id")) %>% 
    mutate(enrich = (causal_hit / n_causal) / (hits / num_var),
           enrich = ifelse(!is.finite(enrich) | is.na(enrich), 0, enrich)) %>% 
    select(phenotype_id, GeneSymbol, enrich, celltype, ref_cell) 
  
  # by cell type
  n_bootstrap <- 1000
  set.seed(2024)
  se <- tibble(ref_cell = unique(allcell_df$ref_cell), se = NA, OR_L = NA, OR_U = NA, mean_enrich = NA)
  
  for (i in seq_along(se$ref_cell)){
    
    v_enrich <- allcell_df %>% filter(ref_cell == se$ref_cell[i]) %>% 
      right_join(allegenes, by=c("phenotype_id", "celltype")) %>% 
      mutate(enrich = ifelse(is.na(enrich) | !is.finite(enrich), 0, enrich)) %>% 
      pull(enrich)
    
    bootstrap_mean <- replicate(n_bootstrap, mean(sample(v_enrich, replace = T)))
    se$se[i] <- sqrt(sum((bootstrap_mean - mean(bootstrap_mean))^2)/n_bootstrap)
    
    se$OR_L[i] <- quantile(bootstrap_mean, 0.025)
    se$OR_U[i] <- quantile(bootstrap_mean, 0.975)
    se$mean_enrich[i] <- mean(v_enrich)
  }
  se
  
  allcell_df <- se %>% mutate(celltype = "allcells") %>% 
    select(celltype, mean_enrich, se, OR_L, OR_U, ref_cell)
  
  
  if (method == "jaxqtl"){
    assign("nb_allcell_scent", allcell_df)
  }else if (method == "tqtl"){
    assign("tqtl_allcell_scent", allcell_df)
  }else if (method == "jaxqtl_lm_score"){
    assign("lm_score_allcell_scent", allcell_df)
  }
}

sc_df %>% filter(GeneSymbol == "ANKRD55")
sc_df %>% filter(GeneSymbol == "IL6ST")

nb_enrich
nb_allcell_scent
lm_score_allcell_scent


OR_L_df <- lm_score_enrich %>% select(celltype, lm_score_enrich=OR_L) %>% 
  inner_join(nb_enrich %>% select(celltype, nb_enrich=OR_L), by="celltype") %>% 
  gather(key = method, value = OR_L, c(lm_score_enrich, nb_enrich))

OR_U_df <- lm_score_enrich %>% select(celltype, lm_score_enrich=OR_U) %>% 
  inner_join(nb_enrich %>% select(celltype, nb_enrich=OR_U), by="celltype") %>% 
  gather(key = method, value = OR_U, c(lm_score_enrich, nb_enrich))

summary(nb_enrich$mean_enrich)
summary(lm_score_enrich$mean_enrich)

nb_scent %>% distinct(celltype, gene_name) %>% nrow()
nb_scent %>% distinct(celltype, phenotype_id) %>% nrow()

# meta analysis
nb_enrich
scent_meta <- data.frame(sc_bulk = c("sc-eQTL", "bulk-eQTL"),
                        meta_OR=NA, meta_lower=NA, meta_upper=NA, seTE=NA)

# nb_enrich and nb_allcell_scent_tmp
for (idx in 1:2){
  method <- scent_meta$sc_bulk[idx]; method

  if (method == "sc-eQTL"){
    meta_dat <- nb_enrich %>% filter(celltype != "allcells")
  }else{
    meta_dat <- nb_allcell_scent
  }
  
  # meta analyze over OR and SE of OR via bootstraps
  meta_dat <- meta_dat %>% 
    mutate(TE = mean_enrich, 
           lower = NA,
           upper = NA,
           seTE = se); meta_dat
  
  summary(nb_enrich$mean_enrich)
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = celltype,
                       data = meta_dat,
                       sm = "",
                       # method.tau = "PM",
                       fixed = TRUE,
                       random = FALSE,
                       title = "Enrichment (Pre-calculated)")
  
  scent_meta$meta_OR[idx] <- m.gen_bin$TE.fixed
  scent_meta$meta_lower[idx] <- m.gen_bin$lower.fixed
  scent_meta$meta_upper[idx] <- m.gen_bin$upper.fixed
  scent_meta$seTE[idx] <- m.gen_bin$seTE.fixed
}

scent_meta
# new results
# sc: 1.942365 [1.762477; 2.122254]
# bulk: 1.389736 [1.146958; 1.632515]

stats <- (scent_meta$meta_OR[1] - scent_meta$meta_OR[2])/sqrt(sum(scent_meta$seTE^2))
pnorm(abs(stats), lower.tail = F)*2


nb_enrich
View(lm_score_enrich)


### compare NegBinom with linear
# method wise comparison
scent_meta <- data.frame(model = c("NegBinom", "Linear"),
                         meta_OR=NA, meta_lower=NA, meta_upper=NA, seTE=NA)


# nb_enrich and nb_allcell_scent_tmp
for (idx in 1:2){
  method <- scent_meta$model[idx]; method
  
  if (method == "NegBinom"){
    meta_dat <- nb_enrich %>% 
      filter(celltype != "allcells")
  }else{
    meta_dat <- lm_score_enrich %>% 
      filter(celltype != "allcells")
  }
  
  # meta analyze over OR and SE of OR via bootstraps
  meta_dat <- meta_dat %>% 
    mutate(TE = mean_enrich, 
           lower = NA,
           upper = NA,
           seTE = se); meta_dat
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = celltype,
                       data = meta_dat,
                       sm = "",
                       # method.tau = "PM",
                       fixed = TRUE,
                       random = FALSE,
                       title = "Enrichment (Pre-calculated)")
  
  scent_meta$meta_OR[idx] <- m.gen_bin$TE.fixed
  scent_meta$meta_lower[idx] <- m.gen_bin$lower.fixed
  scent_meta$meta_upper[idx] <- m.gen_bin$upper.fixed
  scent_meta$seTE[idx] <- m.gen_bin$seTE.fixed
}

scent_meta
# new results
# sc: 1.942365 [1.762477; 2.122254]
# bulk: 1.389736 [1.146958; 1.632515]

stats <- (scent_meta$meta_OR[1] - scent_meta$meta_OR[2])/sqrt(sum(scent_meta$seTE^2))
pnorm(abs(stats), lower.tail = F)*2


# plot

p_enrich_scent <- forest_plot_theme(bind_rows(nb_enrich %>% mutate(method = celltype),
                                              nb_allcell_df_tmp %>% mutate(method = "bulk-eQTL")) %>% 
                                      mutate(category = "SCENT enhancer-gene") %>% 
                                      rename(enrich = mean_enrich) %>% 
                                      filter(celltype != "allcells") %>% 
                                      mutate(celltype = fct_relevel(celltype, rev(celltype_order))) %>% 
                                      ggplot(aes(y = celltype, x = enrich, color = method)) + 
                                      geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
                                      geom_errorbarh(aes(xmin = OR_L, xmax = OR_U),
                                                     height = 0.25,
                                                     position = position_dodge(0.5)) +
                                      geom_vline(xintercept = 1, color = "red", 
                                                 linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                                      facet_grid(.~category)+
                                      scale_color_manual(values = celltype_cols)+
                                      axis_theme + ylab("")+
                                      xlab("Odds Ratio (95% CI)"))
p_enrich_scent
p_scent_meta <- scent_meta %>% 
  ggplot(aes(y = sc_bulk, x = meta_OR, color = sc_bulk)) + 
  geom_point(shape = 18, size = 2, position = position_dodge(0.5)) +  
  geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                 height = 0.25,
                 position = position_dodge(0.5)) +
  geom_vline(xintercept = 1, color = "red", 
             linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
  scale_color_manual(values = sc_bulk_cols)+
  axis_theme + ylab("")+
  xlab("Odds Ratio (95% CI)")

ct_legend <- get_plot_component(p_enrich_scent+theme(legend.position = "right",
                                                      legend.key.size = unit(0.3, units="cm"),
                                                      legend.text = element_text(size = 5),
                                                      legend.title = element_text(size = 5))+
                                  guides(color=guide_legend(title="Cell type")), "guide-box-right")

row <- plot_grid(
  p_enrich_scent+theme(legend.position = "none",
                        plot.margin = unit(c(0,0,0,0), "cm")) + xlab("") + xlim(0,7),
    #coord_fixed(ratio = 1, xlim = c(0,7)),
  p_scent_meta+theme(legend.position = "none",
                           plot.margin = unit(c(0,0,0.2,0), "cm"))+
    xlab("Enrichment (95% CI)") + xlim(0,7)+
    coord_fixed(ratio = 0.45, xlim = c(0,6.9)),
  align = 'v',
  hjust = -1, # -1
  nrow = 3,
  rel_heights = c(1, 0.25, 0.1), # 0.3
  # rel_heights = c(-20, -10), # -20, -4
  rel_widths=c(1, 1), axis="b",
  label_size = 10, vjust=1
)

plot_grid(row, ct_legend,
          align = 'v',
          hjust = -1, # -1
          nrow = 1,
          rel_heights = c(1, 1), # 0.3
          rel_widths=c(1, 0.2), axis="b",
          label_size = 10, vjust=1)

ggsave2(paste0(plots_dir, "Fig5_scent.png"),
        width = onecol_width, height = onecol_height*1.5, units = 'cm', dpi = 300)
