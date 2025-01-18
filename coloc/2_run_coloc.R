# run coloc susie

# fine-map for each genotype
# use r4sc env

# coloc assumptions
# H0: neither trait has a genetic association in the region
# H1: only trait 1 has a genetic association in the region
# H2: only trait 2 has a genetic association in the region
# H3: both traits are associated, but with different causal variants
# H4: both traits are associated and share a single causal variant

library(optparse)
library(tidyverse)
library(data.table)
library(coloc)
library(susieR)
library(arrow)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/coloc")

option_list <- list(
  make_option(c("--celltype"), type="character", default=NULL, 
              help="cell type", metavar="character"),
  make_option(c("--gwas"), type="character", default=NULL, 
              help="gwas file", metavar="character"),
  make_option(c("--gene"), type="character", default=NULL, 
              help="gene phenotype id", metavar="character"),
  make_option(c("--chr"), type="character", default=NULL, 
              help="chromosome of gene", metavar="character"),
  make_option(c("--model"), type="character", default=NULL, 
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

celltype <- opt$celltype # "CD4_ET"
chr <- opt$chr
chunk <- opt$chunk # "chunk_1"
method <- opt$method # "wald", "score"
model <- opt$model  # "nb", "lm"
N <- opt$sampleN
gene <- opt$gene  # "ENSG00000187608"
gwas_file <- opt$gwas 
outdir <- opt$out

L <- opt$susieL
seed <- opt$seed

set.seed(seed)

#### read summary stats ####
# use Wald first
# celltype <- "CD4_NC"
# chr <- 5
# model <- "nb"
# method <- "wald"
# gene <- "ENSG00000134352"
# chunk <- "chunk_6"
# outdir <- "RA"
# gwas_file <- "RA_Ishigaki2022.tsv.gz"

eqtl_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/"
eqtl_dir <- paste0(eqtl_dir, celltype, "/chr", chr, "/")

eqtl_zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")

print(paste("read parquet: ", eqtl_zscore_file))
eqtl_Z <- read_parquet(eqtl_zscore_file)
eqtl_Z <- eqtl_Z %>% filter(phenotype_id == gene)

# switch to score statistics if any non-converged wald stats
if (sum(eqtl_Z$converged == FALSE) > 0){
  # switch to score statistics
  method <- "score"
  eqtl_zscore_file <- paste0(eqtl_dir, chunk, ".", model, ".cis_qtl_pairs.", chr, ".", method, ".parquet")
  eqtl_Z <- read_parquet(eqtl_zscore_file)
  print(paste("read parquet: ", eqtl_zscore_file))
  
  eqtl_Z <- eqtl_Z %>% filter(phenotype_id == gene)
  if (sum(eqtl_Z$converged == FALSE) > 0){
    print("score test stat has non-converge result")
    q()
  }
}
nsnps_eqtl <- nrow(eqtl_Z)

# read gwas result
gwas_Z <- fread(paste0("data/", gwas_file), header=T)

# inner join with eqtl Z
Z_merge <- eqtl_Z %>% 
  mutate(id = row_number(),
         eqtl_maf = ifelse(af > 0.5, 1-af, af)) %>% 
  inner_join(gwas_Z %>% select(-pos), by="snp") %>% 
  arrange(id)

nsnps_merge <- nrow(Z_merge)

#### Read LD ####
# Caution: these are LD for ALT alleles, which is the effect allele

R_wt_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/nb/", celltype, "/")
R_dir <- paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_cis_ld/ld/", celltype, "/")  

# read in correlation (r), cell-type specific because projecting out diff expression PCs
R_wt <- data.table::fread(paste0(R_wt_dir, gene, ".ld_wt.tsv.gz"), header = F)
R_wt <- as.matrix(R_wt)

R <- data.table::fread(paste0(R_dir, gene, ".ld.raw.tsv.gz"), header = F)
R <- as.matrix(R)
print(paste("read correlation matrices"))

# check dimension between Z and R have the same number of cis-SNPs
if (nrow(R) != nrow(eqtl_Z)){
  stop("cis-SNPs in Z score not match LD mat dimension")
}

if (nrow(R_wt) != nrow(eqtl_Z)){
  stop("cis-SNPs in Z score not match wt LD mat dimension")
}

# subset to merged snps
R <- R[Z_merge$id, Z_merge$id]
R_wt <- R_wt[Z_merge$id, Z_merge$id]

if (nrow(R) != nrow(Z_merge)){
  stop("cis-SNPs in merged Z score not match LD mat dimension")
}

if (nrow(R_wt) != nrow(Z_merge)){
  stop("cis-SNPs in merged Z score not match wt LD mat dimension")
}

#### convert to coloc data structure and diagnotic ####
# convert df to coloc data structure

# eqtl: beta, varbeta, N, type="quant", LD, snp, position
# get cell type sample size
if (celltype == "allcells"){
  celltype_N <- 982
}else{
  celltype_N <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/celltype_samplesize.tsv") %>% 
    rename(cell = celltype) %>% 
    filter(cell == celltype) %>% 
    pull(N) 
}

eqtl_D <- Z_merge %>% 
  mutate(slope_se2 = slope_se^2) %>% 
  select(beta = slope, 
         varbeta = slope_se2,
         snp,
         MAF = eqtl_maf,
         position = pos)

colnames(R_wt) <- eqtl_D$snp
rownames(R_wt) <- eqtl_D$snp

eqtl_D <- as.list(eqtl_D)
eqtl_D$type <- "quant"
eqtl_D$N <- celltype_N
eqtl_D$LD <- R_wt

check_dataset(eqtl_D, warn.minp=1e-10)

# gwas data
# beta, varbeta, N, type="cc", LD, snp, position
gwas_D <- Z_merge %>% 
  mutate(se2 = se^2) %>% 
  select(beta, 
         varbeta = se2,
         snp,
         MAF = eqtl_maf,
         position = pos)


gwas_D <- as.list(gwas_D)
gwas_D$type <- "cc"
gwas_D$N <- max(Z_merge$N)
gwas_D$s <- max(Z_merge$N_cases) / max(Z_merge$N)

colnames(R) <- eqtl_D$snp
rownames(R) <- eqtl_D$snp
gwas_D$LD <- R

check_dataset(gwas_D, warn.minp=1e-10)

# check alignment
png(paste0(outdir, "/plot/", gene, ".", celltype, ".check_align_ld_wt.png"))
prop_pos_eqtl <- check_alignment(eqtl_D)
dev.off()
png(paste0(outdir, "/plot/", gene, ".", celltype, ".check_align_ld.png"))
prop_pos_gwas <- check_alignment(gwas_D)
dev.off()

# plot coloc
png(paste0(outdir, "/plot/", gene, ".", celltype, ".data.png"))
par(mfrow=c(2,1))
plot_dataset(gwas_D, main="GWAS")
plot_dataset(eqtl_D, main="eQTL")
dev.off()

# check Z
condz_in <- kriging_rss(eqtl_D$beta/sqrt(eqtl_D$varbeta), R_wt, n=celltype_N)
png(paste0(outdir, "/plot/", gene, ".", celltype, ".eqtl.obsZ_expZ.png"))
condz_in$plot
dev.off()

condz_in <- kriging_rss(gwas_D$beta/sqrt(gwas_D$varbeta), R, n=max(Z_merge$N))
png(paste0(outdir, "/plot/", gene, ".", celltype, ".gwas.obsZ_expZ.png"))
condz_in$plot
dev.off()

##### run susie and coloc ####
susie_eqtl <- runsusie(eqtl_D, maxit = 10000, L = 10, 
               estimate_residual_variance=TRUE, 
               check_prior=FALSE,
               refine=TRUE)

if (susie_eqtl$converged != TRUE){
  susie_eqtl <- runsusie(eqtl_D, maxit = 10000, L = 10, 
                         estimate_residual_variance=FALSE, 
                         check_prior=FALSE,
                         refine=TRUE)
}

susie_gwas <- runsusie(gwas_D, maxit = 10000, L = 10, 
                       estimate_residual_variance=TRUE, 
                       check_prior=FALSE,
                       refine=TRUE)

if (susie_gwas$converged != TRUE){
  susie_gwas <- runsusie(gwas_D, maxit = 10000, L = 10, 
                         estimate_residual_variance=FALSE, 
                         check_prior=FALSE,
                         refine=TRUE)
}

coloc.res <- coloc.susie(susie_eqtl, susie_gwas)
coloc.res$summary

if (!is.null(coloc.res$summary)){
  
  coloc.res$summary %>% 
    left_join(tibble(hit1 = names(susie_eqtl$pip),
                     pip_eqtl = susie_eqtl$pip),
              by="hit1") %>% 
    left_join(tibble(hit2 = names(susie_gwas$pip),
                     pip_gwas = susie_gwas$pip),
              by="hit2") %>% 
    mutate(celltype = celltype, gene = gene) %>% 
    write_tsv(paste0(outdir, "/", gene, ".", celltype, ".coloc.tsv"))
  
  png(paste0(outdir, "/plot/", gene, ".", celltype, ".senitivity.png"))
  for (i in 1:nrow(coloc.res$summary)){
    sensitivity(coloc.res,"H4 >= 0.9",row=i,dataset1=eqtl_D, dataset2=gwas_D)
    sensitivity(coloc.res,"H4 >= 0.9",row=i,dataset1=eqtl_D,dataset2=gwas_D)
  }
  dev.off()
}
