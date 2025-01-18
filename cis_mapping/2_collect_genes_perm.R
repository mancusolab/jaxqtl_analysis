### gather results across chromosomes

library(tidyverse)
library(data.table)
library(qvalue)
library(optparse)

model <- "lm" # nb, lm, pois
test <- "score"
threshold <- 0.01
fdr_val <- 0.05

jqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha_perm/all_celltype/"
tqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl_new/all_celltype/"

jqtl_colclass <- c("character", "numeric", "numeric", "character", 
                   rep("numeric", 6), "logical", rep("numeric", 7))
tqtl_colclass <- c("character", rep("numeric", 5), "character", rep("numeric", 10))

gene_sum_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/"

# specify cell types
allcelltypes <- c("CD4_NC", "B_IN", "Plasma")

for (cell_type in allcelltypes){
  print(cell_type)
  
  allres <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_", test, ".newbeta.", model, ".", cell_type,".tsv.gz"), 
                  header = TRUE, colClasses = jqtl_colclass) %>% 
    mutate(pval_beta = as.numeric(pval_beta),
           beta_shape1 = as.numeric(beta_shape1),
           beta_shape2 = as.numeric(beta_shape2))
  
  gene_summary <- fread(paste0(gene_sum_dir, cell_type, "_gene_summary.tsv.gz"), header=TRUE)
  
  genes_pass <- gene_summary %>% filter(express_percent >= threshold) %>% pull(Geneid)
  
  # write out jaxqtl cis_genes
  print(summary(allres$pval_beta))
  jaxqtl_cisgenes <- allres %>% 
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>% 
    filter(phenotype_id %in% genes_pass)
  
  print(paste0("Number of egenes:", nrow(jaxqtl_cisgenes)))
  
  jaxqtl_cisgenes %>%
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_", test, ".newbeta.", model, ".threshold", threshold, ".", cell_type, ".leadsnp.tsv.gz"))
}


## gather cis results tensorqtl
if (opt$method == "tqtl"){
  for (cell_type in allcelltypes){
    print(cell_type)
    
    allres <- fread(paste0(tqtl_dir, "tqtl_allres.cis.", cell_type,".tsv.gz"), 
                    header = TRUE,
                    colClasses = tqtl_colclass)
    
    gene_summary <- fread(paste0(gene_sum_dir, cell_type, "_gene_summary.tsv.gz"), header=TRUE)
    
    genes_pass <- gene_summary %>% filter(express_percent >= threshold) %>% pull(Geneid)
    
    # write out tqtl cis_genes
    tqtl_cisgenes <- allres %>% 
      # mutate(pval_beta = ifelse(pval_beta > 1 | pval_beta < 0, pbeta(pval_nominal, beta_shape1, beta_shape2), pval_beta)) %>% 
      filter(!is.na(pval_beta)) %>% 
      filter(phenotype_id %in% genes_pass) %>% 
      mutate(qval = qvalue(pval_beta, fdr.level = fdr_val)$qvalue) %>% 
      filter(qval < fdr_val) 
    print(paste0("Number of egenes:", nrow(tqtl_cisgenes)))
    
    tqtl_cisgenes %>% 
      write_tsv(paste0(tqtl_dir, "/tqtl_allres.", cell_type, ".threshold", threshold, ".leadsnp.tsv.gz"))
    
    tqtl_cisgenes %>% 
      select(phenotype_id) %>% 
      write_tsv(paste0(tqtl_dir, "/tqtl_allres.", cell_type, ".threshold", threshold, ".cisgenes.tsv"), col_names = F)
    
  }
}
