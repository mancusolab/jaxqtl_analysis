### gather results across chromosomes

library(tidyverse)
library(data.table)
library(qvalue)

args <- commandArgs(trailingOnly=TRUE)
model <- args[1] # nb, lm, pois
no_tqtl <- args[2] # notqtl
threshold <- as.numeric(args[3]) # 0.01, 020

jqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16/all_celltype/"
tqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl/all_celltype/"

jqtl_colclass <- c("character", "numeric", "numeric", "character",
                   rep("numeric", 10), "logical")
tqtl_colclass <- c("character", rep("numeric", 5), "character", rep("numeric", 10))

gene_sum_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16/metadata/"
tss_strand_path <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/Homo_sapiens.GRCh37.87.strand.bed.tsv.gz"

# threshold <- 0.01 # threshold for jaxqtl (and tqtl)
fdr_val <- 0.05

# get 14 cell types
allcelltypes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv", F) %>%
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_")) %>% 
  pull(celltype)

tss_strand <- fread(tss_strand_path) %>% 
  rename(phenotype_id=gene_id)

for (cell_type in allcelltypes){
  print(cell_type)
  
  allres <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_score.", model, ".", cell_type,".tsv.gz"), 
                  header = TRUE, colClasses = jqtl_colclass) %>% 
    mutate(pval_beta = as.numeric(pval_beta),
           beta_shape1 = as.numeric(beta_shape1),
           beta_shape2 = as.numeric(beta_shape2))
  
  gene_summary <- fread(paste0(gene_sum_dir, cell_type, "_gene_summary.tsv.gz"), header=TRUE)
  
  genes_pass <- gene_summary %>% filter(express_percent >= threshold) %>% pull(Geneid)
  
  # write out genes not converged
  jaxqtl_genes_notcvg <- allres %>% 
    filter(cov_model_converged == FALSE | beta_converged < 1) %>% 
    filter(phenotype_id %in% genes_pass)
  print(paste0("Number of genes not converged:", nrow(jaxqtl_genes_notcvg)))
  
  gene_summary %>% 
    filter(Geneid %in% jaxqtl_genes_notcvg$phenotype_id) %>% 
    write_tsv(paste0(jqtl_dir, "jaxqtl_allres_cis_score.", model, ".threshold", threshold, ".", cell_type, ".genes.fail.tsv.gz"))
  
  # write out jaxqtl cis_genes
  jaxqtl_cisgenes <- allres %>% 
    mutate(pval_beta = ifelse(pval_beta > 1 | pval_beta < 0, pbeta(pval_nominal, beta_shape1, beta_shape2), pval_beta)) %>% 
    filter(!is.na(pval_beta) & cov_model_converged == TRUE & beta_converged > 0) %>% 
    filter(phenotype_id %in% genes_pass) %>% 
    mutate(qval = qvalue(pval_beta, fdr.level = fdr_val)$qvalue) %>% 
    filter(qval < fdr_val) 
  print(paste0("Number of egenes:", nrow(jaxqtl_cisgenes)))
  
  jaxqtl_cisgenes %>% 
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_score.", model, ".threshold", threshold, ".", cell_type, ".leadsnp.tsv.gz"))
  
  jaxqtl_cisgenes %>% 
    select(phenotype_id) %>% 
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_score.", model, ".threshold", threshold, ".", cell_type, ".cisgenes.tsv"), col_names = F)
  
  # call cisegenes on pval_nominal
  jaxqtl_cisgenes_nom <- allres %>% 
    filter(cov_model_converged == TRUE & pval_nominal >= 0 & pval_nominal <= 1) %>% 
    filter(phenotype_id %in% genes_pass) %>% 
    mutate(qval = qvalue(pval_nominal, fdr.level = fdr_val, pi0=1)$qvalue) %>% 
    filter(qval < fdr_val) 
  print(paste0("Number of egenes on pval_nom:", nrow(jaxqtl_cisgenes_nom)))
  
  jaxqtl_cisgenes_nom %>% 
    select(phenotype_id) %>% 
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_score.", model, ".threshold", threshold, ".", cell_type, ".pvalnom.cisgenes.tsv"), col_names = F)
  
  # write jaxqtl lead SNPs for qtl.bed
  qtl_bed <- jaxqtl_cisgenes %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(chrom, variant_id, phenotype_id, strand) %>%
    separate(variant_id, into=c("chrom", "pos")) %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.numeric(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  qtl_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".", model, ".threshold", threshold, ".leadSNP.bed"), col_names = F)
  
  # write tss file for qtltools (all background enter for QTL mapping)
  tss_bed <- gene_summary %>%
    select(phenotype_id = Geneid, chr) %>% 
    mutate(chr = as.character(chr)) %>% 
    filter(phenotype_id %in% genes_pass) %>%
    left_join(tss_strand, by=c("phenotype_id", "chr")) %>%
    mutate(variant_id = "SNP") %>%
    select(chr, start, end, phenotype_id, variant_id, strand)
  
  tss_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".", model, ".threshold", threshold, ".genespass.tss.bed"), col_names = F)
  
}


## gather cis results tensorqtl
if (no_tqtl != "notqtl"){
  for (cell_type in allcelltypes){
    print(cell_type)
    
    allres <- fread(paste0(tqtl_dir, "tqtl_allres.cis.", cell_type,".tsv.gz"), 
                    header = TRUE,
                    colClasses = tqtl_colclass)
    
    gene_summary <- fread(paste0(gene_sum_dir, cell_type, "_gene_summary.tsv.gz"), header=TRUE)
    
    genes_pass <- gene_summary %>% filter(express_percent >= threshold) %>% pull(Geneid)
    
    
    # write out tqtl cis_genes
    tqtl_cisgenes <- allres %>% 
      mutate(pval_beta = ifelse(pval_beta > 1 | pval_beta < 0, pbeta(pval_nominal, beta_shape1, beta_shape2), pval_beta)) %>% 
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
    
    # nominal
    tqtl_cisgenes_nom <- allres %>% 
      filter(!is.na(pval_nominal) & pval_nominal >= 0 & pval_nominal <= 1) %>% 
      filter(phenotype_id %in% genes_pass) %>% 
      mutate(qval = qvalue(pval_nominal, fdr.level = fdr_val, pi0=1)$qvalue) %>% 
      filter(qval < fdr_val) 
    print(paste0("Number of egenes:", nrow(tqtl_cisgenes_nom)))
    
    tqtl_cisgenes_nom %>% 
      select(phenotype_id) %>% 
      write_tsv(paste0(tqtl_dir, "/tqtl_allres.", cell_type, ".threshold", threshold, ".pvalnom.cisgenes.tsv"), col_names = F)
    
    # write tqtl lead SNPs for qtl.bed
    qtl_bed <- tqtl_cisgenes %>% 
      left_join(tss_strand, by="phenotype_id") %>% 
      select(variant_id, phenotype_id, strand) %>% 
      separate(variant_id, into=c("chrom", "pos")) %>% 
      mutate(snp = paste0(chrom, "_", pos), 
             pos_1 = as.numeric(pos) - 1) %>% 
      select(chrom, pos_1, pos, snp, phenotype_id, strand)
    
    qtl_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".threshold", threshold, ".leadSNP.bed"), col_names = F)
    
    # write tss bed file for qtltools (shared with tqtl)
    tss_bed <- gene_summary %>%
      select(phenotype_id = Geneid, chr) %>% 
      mutate(chr = as.character(chr)) %>% 
      filter(phenotype_id %in% genes_pass) %>%
      left_join(tss_strand, by=c("phenotype_id", "chr")) %>%
      mutate(variant_id = "SNP") %>%
      select(chr, start, end, phenotype_id, variant_id, strand)
    
    tss_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".threshold", threshold, ".genespass.tss.bed"), col_names = F)
  }
}

# call genes and lead SNP distinct to either jaxqtl or tensorqtl
# Note: cis genes directory changed
for (cell_type in allcelltypes){
  print(cell_type)
  jaxqtl_res <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_score.", model, ".", cell_type,".tsv.gz"), header = TRUE)
  jaxqtl_genes <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_score.", model, ".threshold", threshold, ".", 
                               cell_type, ".cisgenes.tsv"), header = FALSE) %>% pull(V1)
  
  tqtl_res <- fread(paste0(tqtl_dir, "tqtl_allres.cis.", cell_type,".tsv.gz"), header = TRUE)
  tqtl_genes <- fread(paste0(tqtl_dir, "tqtl_allres.", cell_type, ".threshold", threshold, ".cisgenes.tsv"), header = FALSE) %>% pull(V1)
  
  # find gene set distinct to each
  genes_jaxqtl_only <- jaxqtl_genes[!jaxqtl_genes %in% tqtl_genes]; print(length(genes_jaxqtl_only))
  genes_tqtl_only <- tqtl_genes[!tqtl_genes %in% jaxqtl_genes]; print(length(genes_tqtl_only))
  genes_both <- jaxqtl_genes[jaxqtl_genes %in% tqtl_genes]; print(length(genes_both))
  
  # write qtl for genes only found by jaxqtl
  jqtl_bed <- jaxqtl_res %>%
    filter(phenotype_id %in% genes_jaxqtl_only) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(chrom, variant_id, phenotype_id, strand) %>%
    separate(variant_id, into=c("chrom", "pos")) %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.numeric(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  jqtl_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".", model, ".threshold", threshold, ".leadSNP.jqtl.only.bed"), col_names = F)
  
  # write qtl for genes only found by tqtl
  tqtl_bed <- tqtl_res %>%
    filter(phenotype_id %in% genes_tqtl_only) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(variant_id, phenotype_id, strand) %>%
    separate(variant_id, into=c("chrom", "pos")) %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.numeric(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  tqtl_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".", model, ".threshold", threshold, ".leadSNP.tqtl.only.bed"), col_names = F)
  
  # test both in tqtl result (could do this in jaxqtl)
  tqtl_bed <- tqtl_res %>%
    filter(phenotype_id %in% genes_both) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(variant_id, phenotype_id, strand) %>%
    separate(variant_id, into=c("chrom", "pos")) %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.numeric(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  tqtl_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".", model, ".threshold", threshold, ".leadSNP.both.bed"), col_names = F)
}


