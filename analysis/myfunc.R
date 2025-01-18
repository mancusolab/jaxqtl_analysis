### myfuncs

library(dplyr)

compare_metric <- function(joined = cisgenes_summary, metric, alternative="t", 
                           set1="jaxqtl", set2="tqtl"){
  wilcox.test(cisgenes_summary %>% filter(jaxqtl_only == set1) %>% pull(metric),
              cisgenes_summary %>% filter(jaxqtl_only == set2) %>% pull(metric),
              alternative = alternative) %>% broom::tidy()
}

# show condition of evaluating lines of code
show_condition <- function(code) {
  tryCatch(code,
           error = function(c) "error",
           warning = function(c) "warning",
           message = function(c) "message"
  )
}

# test median or median of two vector
test_diff <- function(vec1, vec2, var_name, alt="t", names= c("NB", "tqtl"), method = "np"){
  vec1 <- vec1[is.finite(vec1)]
  vec2 <- vec2[is.finite(vec2)]
  print(summary(vec1))
  print(summary(vec2))
  if (method == "np"){
    res <- wilcox.test(vec1, vec2, alternative = alt) %>% broom::tidy()
  }else if (method == "t"){
    res <- t.test(vec1, vec2, alternative = alt) %>% broom::tidy()
  }
  boxplot(vec1,
          vec2,
          names = names,
          main = var_name)
  return(res)
}

# meta analysis across cell types
meta_annot <- function(dat, which_annot, which_method){
  meta_dat <- dat %>% 
    filter(annot == which_annot & method == which_method) %>% 
    mutate(TE = log(OR), lower = log(OR_L), upper = log(OR_U),
           seTE = NA) 
  
  m.gen_bin <- metagen(TE = TE,
                       seTE = seTE,
                       lower = lower,
                       upper = upper,
                       studlab = snplist,
                       data = meta_dat,
                       sm = "OR",
                       method.tau = "PM",
                       fixed = FALSE,
                       random = TRUE,
                       title = "Enrichment (Pre-calculated)")
  return(m.gen_bin)
}

# simple inverse variance meta analysis
meta_fix <- function(TE, seTE){
  wts <- 1/(seTE^2)
  sum_wts <- sum(wts)
  TE_fixed <- sum(TE * wts)/sum_wts
  TE_fixed_se <- sqrt(1/sum_wts)
  return(list(meta_TE = TE_fixed,
                meta_seTE = TE_fixed_se))
}
# annotation

annot_list_specific <- c("CTCF_Hoffman",
                         "DHS_Trynka", "DGF_ENCODE", "TFBS_ENCODE",
                         "H3K27ac_Hnisz", "H3K9ac_Trynka", "H3K4me1_Trynka", "H3K27ac_PGC2", "H3K4me3_Trynka",
                         "SuperEnhancer_Hnisz", "Enhancer_Andersson",
                         "Human_Enhancer_Villar",
                         "Human_Promoter_Villar")

annot_list_agnostic <- c("Coding_UCSC",
                         "Promoter_UCSC","UTR_3_UCSC", "UTR_5_UCSC",
                         "TSS_Hoffman")


# align alleles

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

## create shared jacard index matrix

create_jaccard_mat <- function(dat, which_col){
  # shared eQTLs
  share_M <- matrix(0, nrow = 14, ncol = 14)
  
  share_M_celltype <- unique(dat$celltype)
  colnames(share_M) <- share_M_celltype
  rownames(share_M) <- share_M_celltype
  
  for (i in 1:14){
    for (j in 1:14){
      row_cell <- share_M_celltype[i]
      col_cell <- share_M_celltype[j]
      if (i == j){
        share_M[i, j] <- 1
      }else{
        row_cell_cs <- dat %>% filter(celltype == row_cell) %>% distinct(!!sym(which_col)) %>% pull(!!sym(which_col))
        col_cell_cs <- dat %>% filter(celltype == col_cell) %>% distinct(!!sym(which_col)) %>% pull(!!sym(which_col))
        both_cs <- intersect(row_cell_cs, col_cell_cs)
        union_cs <- union(row_cell_cs, col_cell_cs)
        share_M[i, j] <- length(both_cs) / length(union_cs)
      }
    }
    print(i)
  }
  
  share_M
  share_M <- share_M[match(celltypes_grouped, rownames(share_M)), match(celltypes_grouped, rownames(share_M))]
  return(share_M)
}
