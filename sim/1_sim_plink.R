## create simulated genotype data and write in plink bed files (assume no missing)
# https://cran.rstudio.com/web/packages/genio/vignettes/genio.html

library(genio)

library(optparse)

option_list <- list(
  make_option(c("-m", "--mloci"), type="integer", default=NULL, 
              help="Number of loci", metavar="number"),
  make_option(c("-n", "--nind"), type="integer", default=NULL, 
              help="Number of individual", metavar="number"),
  make_option(c("-p", "--maf"), type="numeric", default=NULL, 
              help="minor allele frequency", metavar="number"),
  make_option(c("-s", "--seed"), type="integer", default=NULL, 
              help="Seed", metavar="number"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

set.seed(opt$seed)

# Data dimensions.
# Choose non-multiples of 4 to test edge cases of BED parsers.
# Number of loci.
m_loci <- opt$mloci
# Number of individuals
n_ind <- opt$nind
# Overall allele frequency
# (we don't do any actual genetics in this vignette,
# so just set to a reasonable value)
p <- opt$maf

# Total number of genotypes
n_data <- n_ind * m_loci
# Draw random genotypes from Binomial
X <- rbinom(n_data, 2, p)

# Turn into matrix (m x n)
X <- matrix(X, nrow = m_loci, ncol = n_ind)

# Inspect the first 10 individuals at the first 10 loci
X[1:10, 1:10]

######### bim #########
# We have to specify the number of loci
bim <- make_bim(n = m_loci) # all data in single chromosome

# Inspect the default values
bim
#> # A tibble: 10,001 × 6
#>      chr    id  posg   pos   alt   ref
#>    <dbl> <int> <dbl> <int> <dbl> <dbl>
#>  1     1     1     0     1     2     1
#>  2     1     2     0     2     2     1
#>  3     1     3     0     3     2     1
#>  4     1     4     0     4     2     1
#>  5     1     5     0     5     2     1
#>  6     1     6     0     6     2     1
#>  7     1     7     0     7     2     1
#>  8     1     8     0     8     2     1
#>  9     1     9     0     9     2     1
#> 10     1    10     0    10     2     1
#> # … with 9,991 more rows

# Let's add the "chr" prefix to the chromosome values,
# so we recognize them when we see them later.
# bim$chr <- paste0('chr', bim$chr)
# Make SNP IDs look like "rs" IDs
bim$id <- paste0('rs', bim$id)
# Make positions 1000 bigger
bim$pos <- bim$pos * 10
# Select randomly between Cs and Gs for the reference alleles
bim$ref <- sample(c('C', 'G'), m_loci, replace = TRUE)
# Select randomly between As and Ts for the alternative alleles
bim$alt <- sample(c('A', 'T'), m_loci, replace = TRUE)

# Inspect the table with our changes
bim


######### fam #########

# Specify the number of individuals
fam <- make_fam(n = n_ind)

# Inspect the default values
fam$fam <- 0
fam
#> # A tibble: 1,001 × 6
#>      fam    id   pat   mat   sex pheno
#>    <int> <int> <dbl> <dbl> <dbl> <dbl>
#>  1     1     1     0     0     0     0
#>  2     2     2     0     0     0     0
#>  3     3     3     0     0     0     0
#>  4     4     4     0     0     0     0
#>  5     5     5     0     0     0     0
#>  6     6     6     0     0     0     0
#>  7     7     7     0     0     0     0
#>  8     8     8     0     0     0     0
#>  9     9     9     0     0     0     0
#> 10    10    10     0     0     0     0
#> # … with 991 more rows

# Add prefixes to families and IDs to recognize them later
# fam$fam <- paste0('fam', fam$fam)
# fam$id <- paste0('id', fam$id)
# Sex values are usually 1 and 2
# fam$sex <- sample(1:2, n_ind, replace = TRUE)
# Let's make phenotypes continuous.
# Draw independently from Standard Normal.
fam$pheno <- -9
# Let's leave maternal and paternal IDs as missing

# Inspect again
fam

# Add column and row names from bim and fam tables we just created.
rownames(X) <- bim$id
colnames(X) <- fam$id
# Inspect again the first 10 individuals and loci
X[1:10, 1:10]


####### write plink triplet files #######
# Write genotypes, along with the BIM and FAM files we created.
# Omiting them would result in writing the original dummy version of these tables,
# before we edited them.


# write out region file
library(tidyverse)
tibble(chr = gsub("chr","", unique(bim$chr)), 
       lb = min(bim$pos)-1, 
       ub = max(bim$pos)+1) %>% 
  write_tsv(paste0(gsub("geno", "region", opt$out), "_cis_region.txt"), col_names = F)

