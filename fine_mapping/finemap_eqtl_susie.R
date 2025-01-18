# fine-map for each genotype
# use r4sc env

library(optparse)
library(tidyverse)
library(data.table)
library(susieR)
library(arrow)

option_list <- list(
  make_option(c("--celltype"), type="character", default=NULL, 
              help="cell type", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene phenotype id", metavar="character"),
  make_option(c("--chr"), type="character", default=NULL, 
              help="chromosome of gene", metavar="character"),
  make_option(c("--model"), type="character", default=NULL, 
              help="model of eQTL", metavar="character"),
  make_option(c("--software"), type="character", default=NULL, 
              help="model of eQTL", metavar="character"),
  make_option(c("--method"), type="character", default=NULL, 
              help="method for hypothesis test, wald or score", metavar="character"),
  make_option(c("--sampleN"), type="integer", default=NULL, 
              help="sample size of eQTL mapping", metavar="integer"),
  make_option(c("--chunk"), type="character", default=NULL, 
              help="chunk file to extract gene cis-eqtl mapping parquet file for jaxQTL output only", metavar="character"),
  make_option(c("--susieL"), type="integer", default=10, 
              help="L used in susie", metavar="number"),
  make_option(c("--seed"), type="integer", default=2024, 
              help="Seed", metavar="number"),
  make_option(c("-o", "--out"), type="character", default=".",
              help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

setwd(opt$out)

# example: not converge
# celltype="allcells";chr=13;chunk="chunk_1";model="nb";N=982;gene="ENSG00000198033";outdir=".";L=10;seed=12;method="wald"; software="jaxqtl"

celltype <- opt$celltype # "CD4_ET"
chr <- opt$chr
chunk <- opt$chunk # "chunk_1"
method <- opt$method 
model <- opt$model  # "nb", "lm"
software <- opt$software # jaxqtl, tqtl, or jaxqtl_lm_score
N <- opt$sampleN
gene <- opt$gene  # "ENSG00000187608"

outdir <- paste0(opt$out, "/", software)

L <- opt$susieL
seed <- opt$seed

set.seed(seed)

# specify eqtl director
if (software == "jaxqtl" & model == "nb"){
  print("process for jaxqtl")
  
  res_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/"
  eqtl_dir <- paste0(res_dir, celltype, "/chr", chr, "/")
  
  R_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/nb/", celltype, "/")
  
  zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")
  
  # read in correlation (r), cell-type specific because projecting out diff expression PCs
  ld_suffix <- "ld_wt"

}else {
  # linear model
  print("process for tqtl")
  
  res_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl_new/"
  eqtl_dir <- paste0(res_dir, celltype, "/")
  
  # R_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/ld/")
  R_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/lm/", celltype, "/")  
  
  zscore_file <- paste0(eqtl_dir, "chr", chr, ".cis_qtl_pairs.", chr, ".parquet")
  
  # read in correlation (r), cell-type specific because projecting out diff expression PCs
  ld_suffix <- "ld"
}

# read in correlation (r), cell-type specific because projecting out diff expression PCs
R <- data.table::fread(paste0(R_dir, gene, ".", ld_suffix, ".tsv.gz"), header = F)
R <- as.matrix(R)
print(paste("read R: ", ld_suffix))

# use Wald first
Z <- read_parquet(zscore_file)
print(paste("read parquet: ", zscore_file))

Z <- Z %>% filter(phenotype_id == gene)

# check orders of SNP for tensorqtl
if (software != "jaxqtl"){
  # format tqtl output
  Z <- Z %>% mutate(snp = variant_id,
                    variant_id = gsub("^chr|_b37$", "", variant_id)) %>% 
    separate(variant_id, into=c("chrom", "pos", "ref", "alt"), sep="_") %>% 
    mutate(pos = as.integer(pos))
  
  # plink bim files are sorted
  check_unsorted <- is.unsorted(Z$pos)
  if (check_unsorted){
    Z <- Z %>% arrange(pos)
  }
}

if (software == "jaxqtl" && sum(Z$converged == FALSE) > 0){
  # switch to score statistics
  method <- "score"
  zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")
  Z <- read_parquet(zscore_file)
  print(paste("read parquet: ", zscore_file))
  
  Z <- Z %>% filter(phenotype_id == gene)
  if (sum(Z$converged == FALSE) > 0){
    print("score test stat has non-converge result")
    q()
  }
}

# pull vector of sumstats # check if there is missing value and convergence
z_scores <- Z$slope/Z$slope_se

# check dimension between Z and R have the same number of cis-SNPs
if (nrow(R) != length(z_scores)){
  stop("cis-SNPs in Z score not match LD mat dimension")
}

p <- length(z_scores)

# diagnostic

png(paste0(outdir, "/plot/", gene, ".", celltype, ".", model, ".", method, ".Z.png"))
susie_plot(z_scores, y = "z")
dev.off()

# by default, assume L=10; min_abs_corr=0.5
# estimate_residual_variance=FALSE (prefered for rss) will fix resid_var to be 1
est_var <- FALSE

run_susie_rss <- function(x, y) {
  fitted_rss <- NULL
  tryCatch(
    {
      fitted_rss <- susie_rss(z = z_scores, R = R, n = N, L = L, 
                              estimate_residual_variance=est_var, 
                              refine=TRUE, max_iter = 10000, 
                              check_prior = TRUE)
    },
    error = function(e) {
      cat("Variance estimate too large; turn off prior check \n")
      fitted_rss <- susie_rss(z = z_scores, R = R, n = N, L = L, 
                              estimate_residual_variance=est_var, 
                              refine=TRUE, max_iter = 10000, 
                              check_prior = FALSE)
    }
  )
  return(fitted_rss)
}

fitted_rss <- run_susie_rss()

# output marginal PIP
# post_mu <- colSums(fitted_rss$alpha * fitted_rss$mu)

# if not converged, change seed and do one more time
if (fitted_rss$converged == FALSE){
  set.seed(seed+1)
  fitted_rss <- run_susie_rss()
}


if (fitted_rss$converged == TRUE){
  cs_df <- summary(fitted_rss)$cs
  
  # get marginal PIP
  out <- Z %>% 
    select(snp, chrom, pos) %>% 
    mutate(pip = fitted_rss$pip)
  
  total_cs <- nrow(cs_df)
  
  if (!is.null(cs_df) && total_cs > 0){
    for (i in 1:total_cs){
      idx <- as.integer(strsplit(as.vector(cs_df$variable[i]), ",")[[1]])
      tmp_vec <- rep(0, length(z_scores))
      tmp_vec[idx] <- 1
      out[paste0("cs", i)] <- tmp_vec
    }
    cs_df %>% 
      write_tsv(paste0(outdir, "/", gene, ".", celltype, ".", model, ".L", L, ".estvar", as.character(est_var), ".", method, ".cs.tsv"))
  }else{
    # if no CS, then put 0 for all SNPs
    out <- out %>% mutate(cs1 = 0) 
  }
  
  out %>% 
    write_tsv(paste0(outdir, "/", gene, ".", celltype, ".", model, ".L", L, ".estvar", as.character(est_var), ".", method, ".tsv.gz"))
  
}else{
  print("Not converged")
}

# create susie plot

png(paste0(outdir, "/plot/", gene, ".", celltype, ".", model, ".L", L, ".", method, ".pip.png"))
susie_plot(fitted_rss, y="PIP", 
           main = paste0(celltype, ": ", gene))
dev.off()


# diagnostic
# estiamte lambda
# note: A larger s means there is a strong inconsistency between z scores and LD matrix.
lambda <- estimate_s_rss(z_scores, R, n=N)

tibble(lambda = lambda) %>% 
  write_tsv(paste0(outdir, "/", gene, ".", celltype, ".", model, ".L", L, ".", method, ".lambda.tsv"))

# check observed Z and expected Z
# label red if log LR > 2 and abs(z) > 2
condz_in <- kriging_rss(z_scores, R, n=N)
png(paste0(outdir, "/plot/", gene, ".", celltype, ".", model, ".", method, ".obsZ_expZ.png"))
condz_in$plot
dev.off()
