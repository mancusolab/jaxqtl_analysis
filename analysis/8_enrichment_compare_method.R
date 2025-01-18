
library(tidyverse)
library(data.table)
library(meta)

indir <- "../result/enrich/enrich_sm/"
annot_dir <- "../../jaxqtl_project/data/OneK1K/annotation/"

flag <- "logit.pip.fix.jqtl_linear"
method_suffix <- "" # "_cov"
enrich_method <- "logit" # enrichment, logit

nb <- fread(paste0(indir, "jaxqtl_annot_res.logit.pip.tsv.gz")) %>% 
  filter(model_converged == TRUE) %>% 
  mutate(slope_se = as.numeric(slope_se),
         sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

nb_meta <- fread(paste0(indir, "jaxqtl", method_suffix, "_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(nb_meta)
nb_both_meta <- fread(paste0(indir, "jaxqtl", method_suffix, "_both_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(nb_both_meta)
nb_only_meta <- fread(paste0(indir, "jaxqtl", method_suffix, "_only_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(nb_only_meta)

lm_score <- fread(paste0(indir, "linear_score_annot_res.logit.pip.tsv.gz")) %>% 
  filter(model_converged == TRUE) %>% 
  mutate(slope_se = as.numeric(slope_se),
         sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

lm_score_meta <- fread(paste0(indir, "linear_score", method_suffix, "_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(lm_score_meta)
lm_score_both_meta <- fread(paste0(indir, "linear_score", method_suffix, "_both_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(lm_score_both_meta)
lm_score_only_meta <- fread(paste0(indir, "linear_score", method_suffix, "_only_enrich_meta.", flag, ".tsv.gz")) %>% 
  mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL")) %>% 
  filter(effect == enrich_method) %>% select(-effect); anyNA(lm_score_only_meta)

# tqtl <- fread(paste0(indir, "tqtl_annot_res.logit.pip.tsv.gz")) %>% 
#   filter(model_converged == TRUE) %>% 
#   mutate(slope_se = as.numeric(slope_se),
#          sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))
# tqtl_meta <- fread(paste0(indir, "tqtl", method_suffix,"_enrich_meta.",flag,".tsv.gz")) %>% 
#   mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))
# tqtl_both_meta <- fread(paste0(indir, "tqtl", method_suffix, "_both_enrich_meta.", flag, ".tsv.gz")) %>% 
#   mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))
# tqtl_only_meta <- fread(paste0(indir, "tqtl", method_suffix, "_only_enrich_meta.", flag, ".tsv.gz")) %>% 
#   mutate(sc_bulk = ifelse(celltype == "allcells", "bulk-eQTL", "sc-eQTL"))

# check enrichment from logistic method and enrichment method
nb %>%
  filter(model_converged == TRUE) %>% 
  ggplot(aes(x = slope, y = enrichment)) + geom_point()

####### baseline annotation #######
# H3K9ac: active promoter
# H3K4me3: activate promoter
# H3K4me1: active enhancer
# H3K27ac: active enhancer

baseline_annot_list <- c("Coding_UCSC", 
                         "Promoter_UCSC","TFBS_ENCODE","UTR_3_UCSC", "UTR_5_UCSC","TSS_Hoffman",
                         "DGF_ENCODE", "CTCF_Hoffman",
                         "DHS_Trynka", 
                         "H3K9ac_Trynka", "H3K4me1_Trynka", "H3K27ac_PGC2", "H3K4me3_Trynka",
                         "SuperEnhancer_Hnisz", "Enhancer_Andersson",
                         #"Ancient_Sequence_Age_Human_Enhancer",
                         #"Ancient_Sequence_Age_Human_Promoter",
                         "Conserved_Primate_phastCons46way",
                         "Conserved_LindbladToh",
                         "Conserved_Mammal_phastCons46way",
                         "Conserved_Vertebrate_phastCons46way",
                         "Repressed_Hoffman")
baseline_extend_annot_list <- paste0(baseline_annot_list, ".extend.500")

nb_meta %>% distinct(annot) %>% View
nb_meta %>% filter(annot %in% baseline_annot_list)

# meta across cell types
baseline_meta_celltype <- data.frame()

compare_method <- TRUE
for (method in c("NegBinom", "tensorqtl")){
  if (method == "NegBinom"){
    if (compare_method){
      res <- nb_both_meta
    }else{
      res <- nb_meta 
    }
  }else{
    if (compare_method){
      res <- lm_score_both_meta # tqtl_both_meta
    }else{
      res <- lm_score_meta # tqtl_meta 
    }
  }
  
  tmp <- res %>%  
    filter(annot %in% baseline_annot_list) %>%
    distinct(annot) %>% 
    mutate(meta_OR=NA, meta_lower=NA, meta_upper=NA, meta_p=NA)
  
  for (idx in seq_along(tmp$annot)){
    which_annot <- tmp$annot[[idx]]
    print(idx)
    
    # meta-analyzed those 14 cell types
    if (enrich_method == "logit"){
      meta_dat <- res %>% 
        filter(celltype != "allcells" & annot == which_annot) %>% 
        mutate(TE = log(meta_OR), 
               lower = NA,
               upper = NA,
               seTE = meta_se) 
      sm_opt = "OR"
    }else{
      meta_dat <- res %>% 
        filter(celltype != "allcells" & annot == which_annot) %>% 
        mutate(TE = meta_OR, 
               lower = NA,
               upper = NA,
               seTE = meta_se)
      sm_opt = ""
    }
    
    if (nrow(meta_dat) > 0){
      m.gen_bin <- metagen(TE = TE,
                           seTE = seTE,
                           lower = lower,
                           upper = upper,
                           studlab = celltype,
                           data = meta_dat,
                           sm = sm_opt,
                           # method.tau = "PM",
                           fixed = TRUE,
                           random = FALSE,
                           title = "Enrichment (Pre-calculated)")
      
      metares <- summary(m.gen_bin)
      tmp$meta_OR[[idx]] <- exp(metares$TE.fixed)
      tmp$meta_lower[[idx]] <- exp(metares$lower.fixed)
      tmp$meta_upper[[idx]] <- exp(metares$upper.fixed)
      tmp$meta_p[[idx]] <- metares$pval.fixed
    }
  }
  baseline_meta_celltype <- bind_rows(baseline_meta_celltype,
                                      tmp %>% mutate(method = method))
}

View(baseline_meta_celltype %>% arrange(desc(meta_OR)))

forest_plot_theme(baseline_meta_celltype %>% 
                    filter(meta_OR < 1000) %>% 
                    mutate(annot = gsub("_phastCons46way$", "", annot)) %>% 
                    group_by(annot) %>% 
                    mutate(jaxqtl_or = meta_OR[which(method == "NegBinom")]) %>% ungroup() %>% 
                    mutate(annot = fct_reorder(annot, jaxqtl_or)) %>% 
                    ggplot(aes(y = annot, x = meta_OR, color = method)) +
                    geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    scale_color_manual(values = method_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    theme_bw()+
                    ggtitle("Baseline annotations")) + 
  theme(axis.text.x = element_text(angle=35, hjust = 1),
        legend.position = c(.95, .35),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)
  )

ggsave2(paste0(plots_dir, "/enrich_sm/",
               "jaxqtl_tqtl_bothegenes_enrich_baseline.png"),
        width = onecol_width*1.2, height = onecol_height*1.8, units = 'cm', dpi = 300)


####### ATACseq annotation #######

atacseq_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "ATACseq_chiou")

unique(atacseq_map$annot_celltype) # 10
# atacseq_map <- bind_rows(atacseq_map, data.frame(yazar_celltype="allcells", annot_celltype="PBMC", aggregate="bulk"))

atac_meta_celltype <- data.frame() # cell type
atac_meta_allcells <- data.frame() # allcells result

for (method in c("NegBinom", "tensorqtl")){
  if (method == "NegBinom"){
    if (compare_method){
      res <- nb_both_meta
    }else{
      res <- nb_meta 
    }
  }else{
    if (compare_method){
      res <- lm_score_both_meta # tqtl_both_meta
    }else{
      res <- lm_score_meta # tqtl_meta 
    }
  }
  
  tmp <- data.frame()
  tmp_allcells <- data.frame()
  for (idx in seq_along(atacseq_map$yazar_celltype)){
    yazar_cell <- atacseq_map$yazar_celltype[[idx]]
    atac_cell <- atacseq_map$annot_celltype[which(atacseq_map$yazar_celltype == yazar_cell)]
    
    tmp <- bind_rows(tmp, 
                     res %>% filter(celltype == yazar_cell & annot == paste0("ATACseq_", atac_cell)))
    
    tmp_allcells <- bind_rows(tmp_allcells,
                              res %>% filter(celltype == "allcells" & annot == paste0("ATACseq_", atac_cell)))
  }
  atac_meta_celltype <- bind_rows(atac_meta_celltype,
                                  tmp %>% mutate(method = method))
  atac_meta_allcells <- bind_rows(atac_meta_allcells,
                                  tmp_allcells %>% mutate(method = method))
}

View(atac_meta_celltype)
atac_meta_allcells

atac_meta <- data.frame(model = c("NegBinom", "tensorqtl"),
                        meta_OR=NA, meta_lower=NA, meta_upper=NA, seTE=NA)

for (idx in c("NegBinom", "tensorqtl")){

  if (enrich_method == "logit"){
    meta_dat <- atac_meta_celltype %>% 
      filter(method == idx) %>% 
      distinct(annot, .keep_all = T) %>% 
      mutate(TE = log(meta_OR), 
             lower = NA,
             upper = NA,
             seTE = meta_se)
    sm_opt = "OR"
  }else{
    meta_dat <- atac_meta_celltype %>% 
      filter(method == idx) %>% 
      distinct(annot, .keep_all = T) %>% 
      mutate(TE = meta_OR, 
             lower = NA,
             upper = NA,
             seTE = meta_se)
    sm_opt = ""
  }
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = celltype,
                       data = meta_dat,
                       sm = sm_opt,
                       # method.tau = "PM",
                       fixed = TRUE,
                       random = FALSE,
                       title = "Enrichment (Pre-calculated)")
  atac_meta$meta_OR[atac_meta$model == idx] <- exp(m.gen_bin$TE.fixed)
  atac_meta$meta_lower[atac_meta$model == idx] <- exp(m.gen_bin$lower.fixed)
  atac_meta$meta_upper[atac_meta$model == idx] <- exp(m.gen_bin$upper.fixed)
  atac_meta$seTE[atac_meta$model == idx] <- m.gen_bin$seTE.fixed
}

atac_meta
# meta:
# not different between two models: 0.74
stats <- (log(atac_meta$meta_OR[1]) - log(atac_meta$meta_OR[2]))/sqrt(sum(atac_meta$seTE^2))
pnorm(abs(stats), lower.tail = F) * 2

forest_plot_theme(atac_meta %>% 
                    mutate(model = fct_reorder(model, meta_OR)) %>% 
                    ggplot(aes(y = model, x = meta_OR, color = model)) +
                    geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    # scale_color_manual(values = method_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    ggtitle("Enrichment of enhancer or promoter"))

forest_plot_theme(atac_meta_celltype %>% 
                    group_by(celltype) %>% 
                    mutate(jaxqtl_or = meta_OR[which(method == "NegBinom")]) %>% ungroup() %>% 
                    mutate(celltype = fct_reorder(celltype, jaxqtl_or)) %>% 
                    ggplot(aes(y = celltype, x = meta_OR, color = method)) +
                    geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    scale_color_manual(values = method_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    theme_bw() +
                    ggtitle("cell type-specific ATACseq"))

ggsave2(paste0(plots_dir, "/enrich_sm/",
               "jaxqtl_tqtl_allegenes_enrich_atac.png"),
        width = onecol_width*1.4, height = onecol_height*1.8, units = 'cm', dpi = 300)

####### Epimap annotation ########

epimap_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "EpiMap")

epimap_map <- bind_rows(epimap_map, 
                        data.frame(yazar_celltype="allcells", annot_celltype="PBMC", aggregate="bulk"))

epimap_meta_celltype <- data.frame()

for (method in c("NegBinom", "tensorqtl")){
  if (method == "NegBinom"){
    res <- nb_meta
  }else{
    res <- lm_score_meta # tqtl_meta
  }
  
  for (category in c("promoter", "enhancer")){
    tmp <- data.frame()
    for (idx in seq_along(epimap_map$yazar_celltype)){
      yazar_cell <- epimap_map$yazar_celltype[[idx]]; print(yazar_cell)
      epimap_cell <- epimap_map$annot_celltype[which(epimap_map$yazar_celltype == yazar_cell)]; print(epimap_cell)
      
      tmp <- bind_rows(tmp,
                       res %>% 
                         filter(celltype == yazar_cell & annot == paste0(epimap_cell, ".", category, ".merge")))
      
      # find allcells matched with annotation
      tmp <- bind_rows(tmp,
                       res %>% 
                         filter(celltype == "allcells" & annot == paste0(epimap_cell, ".", category, ".merge")))
      
    }
    epimap_meta_celltype <- bind_rows(epimap_meta_celltype,
                                      tmp %>% mutate(method=method, category=category))
  }
}

epimap_meta_celltype <- epimap_meta_celltype %>% filter(!grepl("^PBMC", annot) & celltype!="allcells")
table(epimap_meta_celltype$celltype)

View(epimap_meta_celltype)

# # meta across celltypes: all cell type specific versus bulk

epimap_meta <- expand.grid(category = c("enhancer", "promoter"),
                           method = c("NegBinom", "tensorqtl")) %>% 
  mutate(meta_OR=NA, meta_lower=NA, meta_upper=NA, seTE=NA)

for (idx in seq_along(epimap_meta$method)){
  which_annot <- epimap_meta$category[[idx]]; 
  which_method <- epimap_meta$method[[idx]]; 
  
  meta_dat <- epimap_meta_celltype %>% 
    filter(method == which_method & celltype != "allcells" & grepl(paste0(which_annot, ".merge"), annot)) %>% 
    mutate(TE = log(meta_OR), 
           lower = NA,
           upper = NA,
           seTE = meta_se)
  
  if (nrow(meta_dat) > 0){
    m.gen_bin <- metagen(TE = TE,
                         seTE = seTE,
                         lower = lower,
                         upper = upper,
                         studlab = celltype,
                         data = meta_dat,
                         sm = "OR",
                         method.tau = "PM",
                         fixed = TRUE,
                         random = FALSE,
                         title = "Enrichment (Pre-calculated)")
    
    metares <- summary(m.gen_bin)
    epimap_meta$meta_OR[[idx]] <- exp(metares$TE.fixed)
    epimap_meta$meta_lower[[idx]] <- exp(metares$lower.fixed)
    epimap_meta$meta_upper[[idx]] <- exp(metares$upper.fixed)
    epimap_meta$seTE[[idx]] <- metares$seTE.fixed
  }
}

epimap_meta

epimap_meta_cat <- epimap_meta %>% filter(category == "promoter" & method == "NegBinom"); epimap_meta_cat
stats <- (log(epimap_meta_cat$meta_OR[1]) - log(epimap_meta_cat$meta_OR[2]))/sqrt(sum(epimap_meta_cat$seTE^2))
pnorm(abs(stats), lower.tail = F) * 2

# forest plot
forest_plot_theme(epimap_meta %>% 
                    # filter(method == "NegBinom") %>% 
                    mutate(method = fct_reorder(method, meta_OR)) %>% 
                    ggplot(aes(y = method, x = meta_OR, color = category)) +
                    geom_point(shape = 18, size = 1.5, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    # scale_color_manual(values = method_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    ggtitle("Enrichment of enhancer or promoter"))


ggsave2(paste0(plots_dir, "/enrich_sm/",
               "jaxqtl_epimap_sc_bulk.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)


####### Check pip and annotation (EpiMap) ########

pip_annot <- fread(paste0(indir, "jaxqtl_pip_annot.tsv.gz"))


