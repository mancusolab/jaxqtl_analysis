### gather results across chromosomes

library(tidyverse)
library(data.table)
library(qvalue)
library(optparse)

option_list <- list(
  make_option(c("-m", "--model"), type="character", default=NULL, 
              help="which jaxqtl model, [nb, lm, pois]", metavar="character"),
  make_option(c("--method"),  type="character",
              help = "Process for tensorqtl results", metavar="character"),
  make_option(c("--test"),  type="character",
              help = "wald or score", metavar="character"),
  make_option(c("-t", "--threshold"), type="numeric", default=NULL, 
              help="threshold for filtering genes", metavar="number"),
  make_option(c("-c", "--celltype"), type="character", default=NULL, 
              help="which cell type (process for one) ", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

model <- opt$model # nb, lm, pois
test <- opt$test
threshold <- as.numeric(opt$threshold) # 0.01, 020
fdr_val <- 0.05

jqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/all_celltype/"
tqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl_new/all_celltype/"

jqtl_colclass <- c("character", "numeric", "numeric", "character", 
                   rep("numeric", 6), "logical", rep("numeric", 7))
tqtl_colclass <- c("character", rep("numeric", 5), "character", rep("numeric", 10))

gene_sum_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/"
tss_strand_path <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/Homo_sapiens.GRCh37.82.strand.bed.tsv.gz"

# get 14 cell types
if (is.null(opt$celltype)){
  allcelltypes <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv", F) %>%
    rename(celltype = X1) %>% 
    mutate(celltype = str_replace(celltype, " ", "_")) %>% 
    pull(celltype) 
}else{
  allcelltypes <- c(opt$celltype)
}

tss_strand <- fread(tss_strand_path) %>% 
  rename(phenotype_id=gene_id)

for (cell_type in allcelltypes){
  print(cell_type)
  
  allres <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_", test, ".newbeta.", model, ".", cell_type,".tsv.gz"), 
                  header = TRUE, colClasses = jqtl_colclass) %>% 
    mutate(pval_beta = as.numeric(pval_beta),
           beta_shape1 = as.numeric(beta_shape1),
           beta_shape2 = as.numeric(beta_shape2))
  
  # remove ill-behaved estimates of alpha with bounded values
  # if (model == "nb"){
  #   allres <- allres %>% filter(alpha_cov > 1e-8 & alpha_cov < 1e10)
  # }
  
  gene_summary <- fread(paste0(gene_sum_dir, cell_type, "_gene_summary.tsv.gz"), header=TRUE)
  
  genes_pass <- gene_summary %>% filter(express_percent >= threshold) %>% pull(Geneid)
  
  # write out jaxqtl cis_genes
  print(summary(allres$pval_beta))
  jaxqtl_cisgenes <- allres %>% 
    # mutate(pval_beta = ifelse(pval_beta > 1 | pval_beta < 0, pbeta(pval_nominal, beta_shape1, beta_shape2), pval_beta)) %>% 
    filter(!is.na(pval_beta) & model_converged == TRUE & beta_converged > 0 & opt_status == TRUE) %>% 
    filter(phenotype_id %in% genes_pass) %>% 
    mutate(qval = qvalue(pval_beta, fdr.level = fdr_val)$qvalue) %>% 
    filter(qval < fdr_val) 
  print(paste0("Number of egenes:", nrow(jaxqtl_cisgenes)))
  
  allres %>% 
    filter(is.na(pval_beta) | model_converged == FALSE | beta_converged < 1 | opt_status == FALSE) %>% 
    filter(phenotype_id %in% genes_pass) %>% 
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_", test, ".newbeta.", model, ".threshold", threshold, ".", cell_type, ".failed.tsv.gz"))
  
  jaxqtl_cisgenes %>%
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_", test, ".newbeta.", model, ".threshold", threshold, ".", cell_type, ".leadsnp.tsv.gz"))

  jaxqtl_cisgenes %>%
    select(phenotype_id) %>%
    write_tsv(paste0(jqtl_dir, "/jaxqtl_allres_cis_", test, ".newbeta.", model, ".threshold", threshold, ".", cell_type, ".cisgenes.tsv"), col_names = F)

  # write jaxqtl lead SNPs for qtl.bed
  qtl_bed <- jaxqtl_cisgenes %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(chrom, variant_id, phenotype_id, strand) %>%
    mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>%
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep = "_") %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.integer(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)

  qtl_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".leadSNP.bed"), col_names = F)

  # write tss file for qtltools (all background enter for QTL mapping)
  tss_bed <- gene_summary %>%
    select(phenotype_id = Geneid, chr) %>%
    mutate(chr = as.character(chr)) %>%
    filter(phenotype_id %in% genes_pass) %>%
    left_join(tss_strand, by=c("phenotype_id", "chr")) %>%
    mutate(variant_id = "SNP") %>%
    select(chr, start, end, phenotype_id, variant_id, strand)

  tss_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".genespass.tss.bed"), col_names = F)

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
    
    # write tqtl lead SNPs for qtl.bed
    qtl_bed <- tqtl_cisgenes %>% 
      left_join(tss_strand, by="phenotype_id") %>% 
      select(variant_id, phenotype_id, strand) %>% 
      mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
      separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep="_") %>% 
      mutate(snp = paste0(chrom, "_", pos), 
             pos_1 = as.integer(pos) - 1) %>% 
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
  jaxqtl_res <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_score.newbeta.", model, ".", cell_type,".tsv.gz"), header = TRUE)
  jaxqtl_genes <- fread(paste0(jqtl_dir, "jaxqtl_allres_cis_score.newbeta.", model, ".threshold", threshold, ".", 
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
    mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep = "_") %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.integer(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  jqtl_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".leadSNP.jqtl.only.bed"), col_names = F)
  
  # write qtl for genes only found by tqtl
  tqtl_bed <- tqtl_res %>%
    filter(phenotype_id %in% genes_tqtl_only) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(variant_id, phenotype_id, strand) %>%
    mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep = "_") %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.integer(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  tqtl_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".leadSNP.tqtl.only.bed"), col_names = F)
  
  # test both in tqtl result
  tqtl_bed <- tqtl_res %>%
    filter(phenotype_id %in% genes_both) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(variant_id, phenotype_id, strand) %>%
    mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep="_") %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.integer(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  tqtl_bed %>% write_tsv(paste0(tqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".leadSNP.both.bed"), col_names = F)
  
  jqtl_bed <- jaxqtl_res %>%
    filter(phenotype_id %in% genes_both) %>%
    left_join(tss_strand, by="phenotype_id") %>%
    select(variant_id, phenotype_id, strand) %>%
    mutate(variant_id = gsub("chr|_b37", "", variant_id)) %>% 
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep="_") %>%
    mutate(snp = paste0(chrom, "_", pos),
           pos_1 = as.integer(pos) - 1) %>%
    select(chrom, pos_1, pos, snp, phenotype_id, strand)
  
  jqtl_bed %>% write_tsv(paste0(jqtl_dir, cell_type, ".newbeta.", model, ".threshold", threshold, ".leadSNP.both.bed"), col_names = F)
}
