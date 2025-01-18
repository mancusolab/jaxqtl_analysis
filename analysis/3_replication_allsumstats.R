# replication using all sumstats
library(tidyverse)
library(data.table)
library(glue)
library(arrow)
library(qvalue)

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
  # ! here we allows strand ambiguous SNPs
  # snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]] = rep(TRUE, length(a1))
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

study_sig <- FALSE # use study specific FDR

indir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha/all_celltype/"

##### eQTLgen ####
if (study_sig == TRUE){
  # significant (FDR<0.05) cis-eQTL results
  eqtlgen <- fread("/project/nmancuso_8/data/eQTLGEN/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz") %>% 
    filter(FDR < 0.05)
}else{
  # full cis-eQTL results:
  eqtlgen <- fread("/project/nmancuso_8/data/eQTLGEN/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
}

for (method in c("nb", "lm")){
  prefix <- glue("jaxqtl_{method}_leadsnp") # jaxqtl_lm_leadsnp
  query <- fread(glue("{indir}/{prefix}.tsv.gz"))
  if (study_sig == TRUE){
    p_cutoff <- 1
  }else{
    p_cutoff <- 0.05/nrow(query)
  }
  
  joined <- query %>%
    left_join(eqtlgen %>% 
                filter(Pvalue < p_cutoff) %>% 
                select(Gene, SNPChr, SNPPos, AssessedAllele, OtherAllele, Zscore_eqtlgen=Zscore), 
              by=c("phenotype_id"="Gene", "chr"="SNPChr", "pos"="SNPPos"))
  
  aligns <- allele.qc(joined$ref,joined$alt,joined$OtherAllele, joined$AssessedAllele)
  
  n_rep <- joined %>% mutate(keep = aligns$keep) %>% 
    filter(keep == TRUE & !is.na(Zscore_eqtlgen)) %>% nrow(); n_rep
  
  print(glue("{method} rep: {n_rep}/{nrow(query)}"))
}

# study wide FDR: 14406/18907; linear: 12765/16654
# p<#eqtl cutoff: NB: 14197/18907; linear: 12606/16654

##### GTEx ####

if (study_sig == TRUE){
  gtex <- fread("/project/nmancuso_8/data/GTEx/GTEXv8/from_web/eQTL_results/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz") %>%
    mutate(phenotype_id = gsub("\\..*", "", gene_id),
           variant_id = gsub("_b38", "", variant_id)) %>% 
    separate(variant_id, into=c("chr", "pos_38", "ref", "alt"), sep="_") %>%
    mutate(chr = as.integer(gsub("chr", "", chr)),
           pos_38 = as.integer(pos_38)) %>% 
    select(phenotype_id, chr, pos_38, pval_nominal, gtex_ref=ref, gtex_alt=alt, slope_gtex = slope) %>%
    filter(!is.na(chr))
}else{
  gtex_dir <- "/project/nmancuso_8/data/GTEx/GTEXv8/from_web/eQTL_results/allsumstats/Whole_Blood"
  
  gtex <- data.frame()
  for (i in 1:22){
    df <- read_parquet(glue("{gtex_dir}/Whole_Blood.v8.EUR.allpairs.chr{i}.parquet"))
    gtex <- bind_rows(gtex, df)
    print(i)
  }
  
  gtex <- gtex %>% mutate(phenotype_id = gsub("\\..*", "", phenotype_id),
                          variant_id = gsub("_b38", "", variant_id)) %>% 
    separate(variant_id, into=c("chr", "pos_38", "ref", "alt"), sep="_") %>%
    mutate(chr = as.integer(gsub("chr", "", chr)),
           pos_38 = as.integer(pos_38)) %>% 
    select(phenotype_id, chr, pos_38, pval_nominal, gtex_ref=ref, gtex_alt=alt, slope_gtex = slope) %>%
    filter(!is.na(chr))
}

for (method in c("nb", "lm")){
  prefix <- glue("jaxqtl_{method}_leadsnp") # jaxqtl_lm_leadsnp
  query <- fread(glue("{indir}/{prefix}.tsv.gz"))
  
  if (study_sig == TRUE){
    p_cutoff <- 1
  }else{
    p_cutoff <- 0.05/nrow(query)
  }
  
  joined <- query %>%
    inner_join(gtex %>% 
                 filter(pval_nominal < p_cutoff) %>% 
                 select(phenotype_id, chr, pos_38, gtex_ref, gtex_alt, slope_gtex), 
               by=c("phenotype_id", "chr", "pos_38"))
  
  aligns <- allele.qc(joined$ref,joined$alt,joined$gtex_ref, joined$gtex_alt)
  
  n_rep <- joined %>% mutate(keep = aligns$keep) %>% 
    filter(keep == TRUE & !is.na(slope_gtex)) %>% nrow(); n_rep
  
  print(glue("{method} rep: {n_rep}/{nrow(query)}"))
}

# study wide FDR: NB: 9493/18907; linear: 8447/16654
# p<#eqtl cutoff: NB: 8053/18907; linear: 7212/16654

##### eQTL catalogue ####


#catlog <- "/project/nmancuso_8/data/eQTL-Catalogue/allsumstats"
catlog <- "/project/nmancuso_8/data/eQTL-Catalogue/sig_eQTL"

# run local
meta <- read_tsv(glue("/project/nmancuso_8/data/eQTL-Catalogue/allsumstats/CT_metadata"))

CT_groups <- unique(meta$group)

for (prefix in c("nb", "lm")){
  prefix <- glue("jaxqtl_{prefix}_leadsnp")
  query <- fread(glue("{indir}/{prefix}.tsv.gz")) %>% 
    mutate(group = case_when(celltype %in% c("CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_ET", "CD8_NC", "CD8_S100B") ~ "T_cell",
                             celltype %in% c("B_IN", "B_Mem", "Plasma") ~ "B_cell",
                             celltype %in% c("NK", "NK_R") ~ "NK",
                             celltype %in% c("Mono_C", "Mono_NC") ~ "monocyte",
                             celltype %in% c("DC") ~ "DC"))
  # p_cutoff <- 0.05/nrow(query)
  # p_cutoff <- 1
  rep_res <- tibble(group = CT_groups, n_query = NA, n_rep = NA)
  
  for (CT_group in CT_groups){
    print(CT_group)
    ref <- data.frame()
    list <- meta %>% filter(group == CT_group) %>% pull(out_name)
    for (i in list){
      # get permuted p values
      tmp <- read_tsv(glue("{catlog}/{i}.tsv.gz"), col_types="cciiciinnnn") %>% 
        mutate(qval = qvalue(p_beta, fdr.level = 0.05)$qvalue) %>% 
        filter(chromosome %in% 1:22)
      ref <- bind_rows(ref, tmp)
      print(i)
    }
    
    ref <- ref %>% separate(variant, into = c("chr", "pos", "ref", "alt"), sep="_") %>% 
      filter(str_length(ref) == 1 & str_length(alt) == 1) %>% 
      distinct(molecular_trait_id, chromosome, position, ref, alt)
    
    joined <- query %>%
      filter(group == CT_group) %>%
      inner_join(ref %>% 
                   select(phenotype_id=molecular_trait_id, chr=chromosome, pos_38=position, rep_ref=ref, rep_alt=alt), 
                 by=c("phenotype_id", "chr", "pos_38")); joined
    
    aligns <- allele.qc(joined$ref,joined$alt,joined$rep_ref, joined$rep_alt)
    
    n_rep <- joined %>% mutate(keep = aligns$keep) %>% 
      filter(keep == TRUE) %>% nrow(); n_rep
    
    rep_res$n_query[rep_res$group == CT_group] <- query %>% filter(group == CT_group) %>% nrow()
    rep_res$n_rep[rep_res$group == CT_group] <- n_rep
  }
  
  rep_res %>% mutate(rate = n_rep / n_query) %>% write_tsv(glue("{indir}/{prefix}.CT_rep.tsv"))
}

