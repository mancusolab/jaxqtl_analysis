# calculate enrichment of per annot per cell type
library(tidyverse)
library(data.table)
library(meta)
library(optparse)

option_list <- list(
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="jaxqtl or tqtl [default= %default]", metavar="character"),
  make_option(c("-f", "--flag"), type="character", default=NULL, 
              help="name flag [default= %default]", metavar="character"),
  make_option(c("--metamodel"), type="character", default=NULL, 
              help="fix or random", metavar="character"),
  make_option(c("-c", "--combineannot"), action = "store_true", default = FALSE,
              help = "Print extra output [default]"),
  make_option(c("-o", "--outdir"), type="character", default="./enrich_meta_res", 
              help="out directory [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

method <- opt$method  # jaxqtl, tqtl, jaxqtl_lm_score
flag <- opt$flag # "logit.pip" or "logit.cs"
outdir <- opt$outdir
combine_annot <- opt$combineannot

meta_model <- opt$metamodel
if (meta_model == "fix"){
  fixed_flag = TRUE
  random_flag = FALSE
  meta_suffix = ".fix.jqtl_linear"
}else if (meta_model == "random"){
  fixed_flag = TRUE
  random_flag = TRUE
  meta_suffix = ".random.jqtl_linear"
}

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/enrich/enrich_sm")

# nb_allegenes.tsv, tqtl_allegenes.tsv, jaxqtl_lm_score_allegenes.tsv
genelist_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/code/"

# removed from finemapping
blacklist_egenes <- read_tsv(paste0(genelist_dir, "blacklist.tsv"), F) %>% pull(X1)

# read egenes list and remove non-finemapped genes
nb_egenes <- read_tsv(paste0(genelist_dir, "nb_allegenes.tsv"), F) %>% 
  select(phenotype_id=X1, celltype=X3) %>% 
  filter(!phenotype_id %in% blacklist_egenes)
lm_egenes <- read_tsv(paste0(genelist_dir, "jaxqtl_lm_score_allegenes.tsv"), F) %>% 
  select(phenotype_id=X1, celltype=X3) %>% 
  filter(!phenotype_id %in% blacklist_egenes)

# create only and both egenes (exclude all cells)
both_egenes_list <- nb_egenes %>% inner_join(lm_egenes, by = c("phenotype_id", "celltype")) %>% 
  filter(celltype != "allcells")
if (method == "jaxqtl"){
  only_egenes_list <- nb_egenes %>% anti_join(lm_egenes, by = c("phenotype_id", "celltype")) %>% 
    filter(celltype != "allcells")
}else{
  only_egenes_list <- lm_egenes %>% anti_join(nb_egenes, by = c("phenotype_id", "celltype")) %>% 
    filter(celltype != "allcells")
}

# combine annotation results
colclass <- c(rep("character", 2), rep("integer", 2), rep("numeric", 3), "logical", "integer", rep("numeric", 2))

if (combine_annot){
  allfiles <- list.files(method)
  allfiles <- allfiles[grepl(flag, allfiles)]
  
  res <- data.frame()
  for (file in allfiles){
    df <- read_tsv(paste0(method, "/", file), col_types="cciinnnlinn")
    res <- bind_rows(res, df %>% mutate(filename=file))
    print(file)
  }
  
  res <- res %>% 
    mutate(filename = gsub(paste0(".", flag, ".annot_res.gz$"), "", filename)) %>% 
    separate(filename, into=c("rm", "celltype"), sep="\\.") %>% 
    select(-rm)
  
  # inner join with egene list
  if (method == "jaxqtl"){
    res <- res %>% inner_join(nb_egenes, by=c("phenotype_id", "celltype"))
  }else{
    res <- res %>% inner_join(lm_egenes, by=c("phenotype_id", "celltype"))
  }
  
  res %>% fwrite(paste0(outdir, "/", method, "_annot_res.", flag, ".tsv.gz"), sep="\t")
}

res <- read_tsv(paste0(outdir, "/", method, "_annot_res.", flag, ".tsv.gz"), 
                col_types="cciinnnlinnc") %>% 
  filter(model_converged == TRUE)


# for each cell type, run meta-analysis over logit.pip results from each gene
run_meta <- function(enrich_res, enrich_col, enrich_se_col, sm_opt){
  
  meta_res <- enrich_res %>% 
    distinct(celltype, annot) %>% 
    mutate(meta_OR=NA, meta_lower=NA, meta_upper=NA, meta_p=NA, meta_se=NA)
  
  for (idx in 1:nrow(meta_res)){
    which_cell <- meta_res$celltype[[idx]]
    which_annot <- meta_res$annot[[idx]]
    print(idx)
    
    meta_dat <- enrich_res %>% 
      filter(celltype == which_cell & annot == which_annot) %>% 
      mutate(lower = NA,
             upper = NA)
    
    meta_dat$TE = meta_dat[[enrich_col]] # effect estimate
    meta_dat$seTE = meta_dat[[enrich_se_col]] # SE
    
    # remove NA
    meta_dat <- meta_dat %>% filter(!is.na(TE) & !is.na(seTE))
    
    if (nrow(meta_dat) > 0){
      tryCatch({
        if (fixed_flag == TRUE){
          # manually calculate fixed effect results
          wts <- 1/(meta_dat$seTE^2)
          sum_wts <- sum(wts)
          TE_fixed <- sum(meta_dat$TE * wts)/sum_wts
          TE_fixed_se <- sqrt(1/sum_wts)
          
          meta_res$meta_OR[[idx]] <- exp(TE_fixed)  # exp(metares$TE.fixed)
          meta_res$meta_lower[[idx]] <- exp(TE_fixed - 1.96*TE_fixed_se) # exp(metares$lower.fixed)
          meta_res$meta_upper[[idx]] <- exp(TE_fixed + 1.96*TE_fixed_se) # exp(metares$upper.fixed)
          meta_res$meta_p[[idx]] <- pnorm(abs(TE_fixed/TE_fixed_se), lower.tail = F)*2 # metares$pval.fixed
          meta_res$meta_se[[idx]] <- TE_fixed_se  # metares$seTE.fixed
        }else{
          # run metagen for random effect model
          m.gen_bin <- metagen(TE = TE,
                               seTE = seTE,
                               lower = lower,
                               upper = upper,
                               studlab = phenotype_id,
                               data = meta_dat,
                               sm = sm_opt, # "OR" or ""
                               method.tau = "PM",
                               fixed = TRUE,
                               random = random_flag,
                               title = "Enrichment (Pre-calculated)")
          
          metares <- summary(m.gen_bin)
          
          meta_res$meta_OR[[idx]] <- exp(metares$TE.random)
          meta_res$meta_lower[[idx]] <- exp(metares$lower.random)
          meta_res$meta_upper[[idx]] <- exp(metares$upper.random)
          meta_res$meta_p[[idx]] <- metares$pval.random
          meta_res$meta_se[[idx]] <- metares$seTE.random 
        }
      },
      error = function(err){
        print(paste0("Error occurs for ", idx))
      })
    }
  }
  return(meta_res)
}

# run meta on egenes found by method specific egenes
tmp <- res %>%
  inner_join(only_egenes_list, by=c("celltype", "phenotype_id"))

out_slope <- run_meta(tmp, "slope", "slope_se", "OR")
out_enrich <- run_meta(tmp, "enrichment", "enrichment_se", "")

bind_rows(out_slope %>% mutate(effect = "logit"),
          out_enrich %>% mutate(effect = "enrichment")) %>% 
  write_tsv(paste0(outdir, "/", method, "_only_enrich_meta.", flag, meta_suffix, ".tsv.gz"))

# run meta on egenes found by both methods
tmp <- res %>%
  inner_join(both_egenes_list, by=c("celltype", "phenotype_id"))

out_slope <- run_meta(tmp, "slope", "slope_se", "OR")
out_enrich <- run_meta(tmp, "enrichment", "enrichment_se", "")

bind_rows(out_slope %>% mutate(effect = "logit"),
          out_enrich %>% mutate(effect = "enrichment")) %>% 
  write_tsv(paste0(outdir, "/", method, "_both_enrich_meta.", flag, meta_suffix, ".tsv.gz"))

# all egenes per method
out_slope <- run_meta(res, "slope", "slope_se", "OR")
out_enrich <- run_meta(res, "enrichment", "enrichment_se", "")

bind_rows(out_slope %>% mutate(effect = "logit"),
          out_enrich %>% mutate(effect = "enrichment")) %>% 
  write_tsv(paste0(outdir, "/", method, "_enrich_meta.", flag, meta_suffix, ".tsv.gz"))
