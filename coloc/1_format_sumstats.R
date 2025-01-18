# format sum stats

library(tidyverse)
library(data.table)
options(scipen=999)

# https://wanggroup.org/compbio_tutorial/allele_qc.html
allele.qc = function(a1,a2,ref1,ref2) {
  # a1 and a2 are the first data-set
  # ref1 and ref2 are the 2nd data-set
  # Make all the alleles into upper-case, as A,T,C,G:
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  # Strand flip, to change the allele representation in the 2nd data-set
  strand_flip = function(ref) {
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip
  }
  flip1 = strand_flip(ref1)
  flip2 = strand_flip(ref2)
  snp = list()
  # Remove strand ambiguous SNPs (scenario 3)
  snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  # snp[["keep"]] = rep(T, length(a1))
  # Remove non-ATCG coding
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  # as long as scenario 1 is involved, sign_flip will return TRUE
  snp[["sign_flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  # as long as scenario 2 is involved, strand_flip will return TRUE
  snp[["strand_flip"]] = (a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)
  # remove other cases, eg, tri-allelic, one dataset is A C, the other is A G, for example.
  exact_match = (a1 == ref1 & a2 == ref2) 
  snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] = F
  return(snp)
}

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/coloc/data")

## read all snps in Onek1k (after filter)
# 10    chr10_90127_C_T_b37     0     90127      T      C
allsnps <- fread("/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_updateid/hwe_maf0.05_snp/allchr.bim", header=F)
colnames(allsnps) <- c("chr", "snp", "cM", "pos", "alt", "ref")

#### EUR_all_auto-10-2021.txt.gz (RA) ####

# [File format (header)]
# - SNP: chr_pos_ref_alt (hg19)
# - Beta: effect size estimates (the effect allele is the alternative allele)
# - SE: S.E. of the effect size estimate
# - Pval: P value

df <- fread("EUR_all_auto-10-2021.txt.gz")
df <- df %>% separate(SNP, into=c("chr", "pos", "gwas_ref", "gwas_alt"), sep="_") %>% 
   filter(str_length(gwas_ref)<2 & str_length(gwas_alt)<2)
  
colnames(df) <- tolower(colnames(df))

df_merge <- df %>% mutate(chr = as.integer(chr), 
              pos = as.integer(pos)) %>% 
  inner_join(allsnps, by=c("chr", "pos"))

qc <- allele.qc(df_merge$gwas_ref, df_merge$gwas_alt, df_merge$ref, df_merge$alt)

# 4484327 snps
df_merge <- df_merge %>%
  mutate(keep = qc$keep,
         sign_flip = qc$sign_flip,
         strand_flip = qc$strand_flip) %>%
  filter(keep == TRUE) 

sum(df_merge$sign_flip) # no sign flip
sum(df_merge$strand_flip) # no strand flip

# RA: total 22,350 cases
# chr       pos snp          ref   alt      beta     se    pval     N N_cases
df_merge %>% select(chr, pos, snp, ref, alt, beta, se, pval) %>% 
  mutate(pval = as.numeric(pval),
         N = 97173, 
         N_cases = 22350) %>% 
  fwrite("./RA_Ishigaki2022.tsv.gz", sep="\t")

#### GCST90014023_buildGRCh38.tsv.gz (T1D) ####

# read all snps with rsid
allsnps <- fread("/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_rsid/allchr.bim", header=F)
colnames(allsnps) <- c("chr", "snp", "cM", "pos", "alt", "ref")

# caution: this is build 38
df <- fread("GCST90014023_buildGRCh38.tsv.gz")

df_merge <- df %>% 
  rename(snp=variant_id, chr=chromosome, pos=base_pair_location) %>% 
  select(-pos) %>% 
  mutate(chr=as.integer(chr)) %>% 
  inner_join(allsnps, by=c("snp", "chr"))

qc <- allele.qc(df_merge$other_allele, df_merge$effect_allele, df_merge$ref, df_merge$alt)

# 10135042 snps
df_merge <- df_merge %>%
  mutate(keep = qc$keep,
         sign_flip = qc$sign_flip,
         strand_flip = qc$strand_flip) %>%
  filter(keep == TRUE) 

sum(df_merge$sign_flip) # no sign flip
sum(df_merge$strand_flip) # no strand flip

# take max total N as sample size: 520580
# cases: 18942
# chr       pos snp          ref   alt      beta     se    pval     N N_cases
df_merge %>% 
  mutate(snp = paste0("chr", chr, "_", pos, "_", ref, "_", alt, "_b37"),
         beta = ifelse(sign_flip == TRUE, beta * -1, beta),
         maf = ifelse(effect_allele_frequency > 0.5, 1-effect_allele_frequency, effect_allele_frequency)) %>% 
  select(chr, pos, snp, 
         ref=other_allele, alt=effect_allele, beta, se=standard_error, pval=p_value, 
         sample_size, maf) %>% 
  mutate(N = max(sample_size), 
         N_cases = 18942) %>% 
  fwrite("./T1D_Chious2021.tsv.gz", sep="\t")

#### Leukemia ####
# chr     pos     snp     ref     alt     beta    se      pval    N       N_cases

# write out bed file for liftover
# df %>% mutate(chr=paste0("chr",chromosome), pos_1 = base_pair_location-1) %>% 
#   select(chr, pos_1, base_pair_location, everything()) %>% 
#   fwrite("GCST90014023_buildGRCh38.bed", col.names=F, sep="\t")

# next use /project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich/liftover.sh change from hg38 -> hg19
# /project/nmancuso_8/elezhang/projects/jaxqtl/code/enrich/liftover.sh GCST90014023_buildGRCh38.bed 38 19 
# beta, varbeta, N, snp

# gwas_in <- "/project/gazal_569/DATA/ldsc/sumstats_EUR_raw/come_PASS/"
# gwas_out <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/gwas/steven/"
# geno_dir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/geno_n982_info_08_rsid/"
# 
# setwd(gwas_out)
# 
# trait_list <- read_tsv("trait24.tsv", F) %>% pull(X1)
# 
# for (trait in trait_list){
#   # read in sumstats
#   print(trait)
#   gwas_df <- fread(paste0(gwas_in, trait, ".sumstats.gz"), header=T)
#   if (grepl("rs", colnames(gwas_df)[1])){
#     gwas_df <- fread(paste0(gwas_in, trait, ".sumstats.gz"), header=F, skip = 1)
#     colnames(gwas_df) <- c("snp", "a1", "a2", "n", "z")
#     print("header use diff tab")
#   }
#   colnames(gwas_df) <- tolower(colnames(gwas_df))
#   
#   gwas_update <- data.frame()
#   for (idx in 1:22){
#     bim <- fread(paste0(geno_dir, "chr", idx, ".bim"), header=F)
#     colnames(bim) <- c("chr", "snp", "cM", "pos", "alt", "ref")
#     tmp <- gwas_df %>% inner_join(bim, by="snp")
#     
#     gwas_update <- bind_rows(gwas_update, tmp)
#   }
#   
#   gwas_update <- gwas_update %>% 
#     rename(rsid = snp) %>% 
#     mutate(snp = paste0("chr", chr, "_", pos, "_", ref, "_", alt, "_b37"))
#   
#   # change d1 coding to d2
#   qc <- allele.qc(gwas_update$a1, gwas_update$a2, gwas_update$ref, gwas_update$alt)
#   
#   gwas_update <- gwas_update %>% 
#     mutate(keep = qc$keep,
#            sign_flip = qc$sign_flip,
#            strand_flip = qc$strand_flip) %>% 
#     filter(keep == TRUE) %>% 
#     mutate(z = ifelse(sign_flip == TRUE, z * -1, z))
#   
#   gwas_update %>% fwrite(paste0(gwas_out, trait, ".tsv.gz"), sep="\t")
# }
