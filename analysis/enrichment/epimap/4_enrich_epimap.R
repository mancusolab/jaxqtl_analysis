# enrich analysis for qtlresults

# use one sample to check robustness of results
# enrich_file <- "allcelltype_newbeta_epimap_onesample_new.tsv.gz"
enrich_file <- "allcelltype_newbeta_epimap_new.tsv.gz"

enrich_res <- fread(paste0(enrich_dir, enrich_file)) %>% filter(!is.na(pval)) %>% 
  mutate(method = gsub("_new_fixalpha", "", method),
         method = ifelse(method == "celltype16", "NegBinom", "tensorqtl"),
         annot = gsub(".merge", "", annot),
         snplist = gsub("allcells", "Bulk-eQTL", snplist)) %>% 
  separate(annot, into=c("annot_celltype", "annot"), sep="\\.") 
table(enrich_res$method)
colnames(enrich_res)[12:14] <- c("OR_U", "OR", "OR_L") # flip OR_L and OR_U

match_annot <- epimap_map %>% 
  rename(celltype=yazar_celltype, annot=annot_celltype) %>% 
  select(celltype, annot)
match_annot <- bind_rows(match_annot, data.frame(celltype="Bulk-eQTL", annot="PBMC"))

# compare lead SNP for egenes found in methods (same eGenes)
snplist_opt = "all"
dat <- data.frame()
for (i in 1:nrow(match_annot)){
  celltype <- match_annot$celltype[i]
  if (snplist_opt == "private"){
    extract_list <- c(paste0(celltype, ".newbeta.nb.threshold0.01.leadSNP.jqtl.only"),
                      paste0(celltype, ".newbeta.nb.threshold0.01.leadSNP.tqtl.only"))
  }else if (snplist_opt == "both"){
    extract_list <- c(paste0(celltype, ".newbeta.nb.threshold0.01.leadSNP.both"),
                      paste0(celltype, ".threshold0.01.leadSNP.both"))
  }else{
    extract_list <- c(paste0(celltype, ".newbeta.nb.threshold0.01.leadSNP"),
                      paste0(celltype, ".threshold0.01.leadSNP"))
  }
  df <- enrich_res %>% 
    filter((snplist %in% extract_list) & annot_celltype == match_annot$annot[i]) %>% 
    mutate(celltype = celltype)
  
  dat <- bind_rows(dat, df)
}
View(dat)

dat %>% write_tsv("./analysis/enrichment/epimap/enrich_epimap.tsv")

forest_plot_theme(
  dat %>% 
    filter(method == "NegBinom") %>% 
    group_by(celltype) %>% 
    mutate(enhancer_or = OR[which(annot == "enhancer")]) %>% 
    ungroup() %>% 
    mutate(celltype = fct_reorder(celltype, enhancer_or)) %>% 
    ggplot(aes(y = celltype, x = OR, color = annot)) +
    geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
    geom_errorbarh(aes(xmin = OR_L, xmax = OR_U), 
                   height = 0.25, 
                   position = position_dodge(width=0.5)) +
    geom_vline(xintercept = 1, color = "red", 
               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
    # scale_color_manual(values = method_cols) +
    xlab("Odds Ratio (95% CI)") + 
    ylab(" ") + 
    theme_bw() +
    ggtitle("Enrichment of enhancer or promoter")
)

ggsave2(paste0(plots_dir, "/enrich/",
               "NB_", snplist_opt ,"_leadsnp_enrich_epimap.png"),
        width = onecol_width*1.1, height = onecol_height*1.5, units = 'cm', dpi = 300)

# fisher's test
p <- dat %>% filter(method == "jaxqtl" & OR >= 1) %>% pull(pval); length(p)
chisq <- sum(-2 * log(p)); chisq
pchisq(chisq, df = 2 * length(p), lower.tail = F)


# meta analysis

# meta analysis
allmeta_jqtl <- data.frame(annot = unique(dat$annot), 
                           meta_OR = NA, meta_lower = NA, meta_upper = NA, meta_p = NA)
allmeta_tqtl <- allmeta_jqtl

# meta analysis
for (annot in unique(dat$annot)){
  res_jqtl <- meta_annot(dat %>% filter(!grepl("Bulk-eQTL", snplist)), annot, "NegBinom")
  res_tqtl <- meta_annot(dat %>% filter(!grepl("Bulk-eQTL", snplist)), annot, "tensorqtl")
  
  allmeta_jqtl$meta_OR[which(allmeta_jqtl$annot == annot)] <- exp(res_jqtl$TE.random)
  allmeta_jqtl$meta_lower[which(allmeta_jqtl$annot == annot)] <- exp(res_jqtl$lower.random)
  allmeta_jqtl$meta_upper[which(allmeta_jqtl$annot == annot)] <- exp(res_jqtl$upper.random)
  allmeta_jqtl$meta_p[which(allmeta_jqtl$annot == annot)] <- res_jqtl$pval.random
  
  allmeta_tqtl$meta_OR[which(allmeta_tqtl$annot == annot)] <- exp(res_tqtl$TE.random)
  allmeta_tqtl$meta_lower[which(allmeta_tqtl$annot == annot)] <- exp(res_tqtl$lower.random)
  allmeta_tqtl$meta_upper[which(allmeta_tqtl$annot == annot)] <- exp(res_tqtl$upper.random)
  allmeta_tqtl$meta_p[which(allmeta_tqtl$annot == annot)] <- res_tqtl$pval.random
}

allmeta_tqtl <- allmeta_tqtl %>% mutate(method = "tensorqtl")
allmeta_jqtl <- allmeta_jqtl %>% mutate(method = "NegBinom")

which_method <- "NegBinom"
allmeta <- bind_rows(allmeta_tqtl %>% select(annot, method, meta_OR), 
                     allmeta_jqtl %>% select(annot, method, meta_OR),
                     dat %>% filter(grepl("Bulk-eQTL", snplist) & method == which_method) %>% 
                       mutate(method = "Bulk-eQTL") %>% 
                       select(annot, method, meta_OR=OR)) %>% 
  mutate(meta_lower = c(allmeta_tqtl$meta_lower, allmeta_jqtl$meta_lower,
                        dat %>% filter(grepl("Bulk-eQTL", snplist) & method == which_method) %>% pull(OR_L)),
         meta_upper = c(allmeta_tqtl$meta_upper, allmeta_jqtl$meta_upper,
                        dat %>% filter(grepl("Bulk-eQTL", snplist) & method == which_method) %>% pull(OR_U)),
         meta_p = c(allmeta_tqtl$meta_p, allmeta_jqtl$meta_p,
                    dat %>% filter(grepl("Bulk-eQTL", snplist) & method == which_method) %>% pull(pval)))

allmeta
forest_plot_theme(allmeta %>% 
                    filter(method %in% c(which_method, "Bulk-eQTL")) %>% 
                    mutate(method = ifelse(method == which_method, paste0("sc-eQTL"), "Bulk-eQTL")) %>% 
                    ggplot(aes(y = method, x = meta_OR, color = annot)) +
                    geom_point(shape = 18, size = 2, position = position_dodge(width=0.5)) +  
                    geom_errorbarh(aes(xmin = meta_lower, xmax = meta_upper),
                                   height = 0.25, 
                                   position = position_dodge(width=0.5)) +
                    geom_vline(xintercept = 1, color = "red", 
                               linetype = "dashed", cex = 1, alpha = 0.5, linewidth=0.6) +
                    # scale_color_manual(values = method_cols) +
                    xlab("Odds Ratio (95% CI)") + 
                    ylab(" ") + 
                    theme_bw() +
                    ggtitle("Enrichment of enhancer or promoter"))

ggsave2(paste0(plots_dir, "/enrich/",
               which_method, "_leadsnp_", snplist_opt, "_enrich_epimap.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)


