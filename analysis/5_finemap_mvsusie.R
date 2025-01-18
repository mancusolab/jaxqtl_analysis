library(mvsusieR)
library(corrplot)
load("../result/mvsusie/result/mashr_finemap.RData")

#### mvsusie results ####

## cross-trait PIPs and CS
pip_res <- read_tsv("../result/mvsusie/result/allres_cs.tsv.gz")
pip_res <- pip_res %>% mutate(snp=gsub("chr|_b37", "", snp)) %>% 
  separate(snp, into=c("chr", "pos", "ref", "alt"), sep='_') %>% 
  left_join(gene_lookup %>% select(phenotype_id, tss=end)) %>% 
  mutate(pos = as.integer(pos),
         tss_dist = pos - tss)
lfsr_res <- read_tsv("../result/mvsusie/result/allres_lfsr.tsv.gz")

pip_res %>% distinct(cs, phenotype_id) # total 5939 CS
pip_res %>% distinct(phenotype_id) # for 3642 egenes

# size of CS
# 1020/5939 has 1-SNP
pip_res %>% group_by(phenotype_id, cs) %>% count() %>% pull(n) %>% summary()
pip_res %>% group_by(phenotype_id, cs) %>% count() %>% pull(n) %>% table()

## trait-wise (average) lfsr, smaller means high confidence in the sign of effect
lfsr_res

## plots for 5939 CSs identified from 3642 eGenes
## 1. number of significant traits per CS
lfsr_thresh <- 0.01

# ENSG00000134352 (IL6ST), ENSG00000164512

lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  add_count(phenotype_id, cs) %>% 
  filter(phenotype_id == "ENSG00000164512")

lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  group_by(phenotype_id, cs) %>% 
  count() %>% 
  ungroup() %>% 
  count(n) %>% 
  mutate(prop = nn/sum(nn))

lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  group_by(phenotype_id, cs) %>% 
  count() %>% ungroup() %>% 
  left_join(pLI_score %>% select(phenotype_id, pLI, loeuf), by="phenotype_id") %>%
  left_join(eds_score %>% select(phenotype_id, EDS), by="phenotype_id") %>%
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol, gene_type), by = "phenotype_id") %>% 
  left_join(rvis_score, by=c("phenotype_id")) %>% 
  mutate(specific = as.factor(ifelse(n < 2, "specific", 
                                     ifelse(n == 14, "ubiquitous", "others")))) %>% 
  ggplot(aes(x = specific, y = loeuf)) + geom_boxplot()

lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  group_by(phenotype_id, cs) %>% 
  count() %>% ungroup() %>% 
  ggplot(aes(x = n)) + geom_histogram() +
  scale_x_continuous(breaks=seq(1,14, 1))+
  axis_theme

ggsave2(paste0(plots_dir, "Fig_mvsusie_sigcell_ct.png"),
        width = onecol_width, height = onecol_height, units = 'cm', dpi = 300)


## 2. upset plot
upset_tmp <- lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  mutate(cs = paste0(phenotype_id, "_", cs)) %>% 
  left_join(cellmeta, by="celltype")
lt <- list()
for (i in unique(upset_tmp$group)){
  lt[[i]] <- upset_tmp %>% filter(group == i) %>% pull(cs)
}

m <- make_comb_mat(lt)
UpSet(m)


## 3. pairwise sharing 

share_M <- matrix(0, nrow = 14, ncol = 14)
share_M_celltype <- unique(lfsr_res$celltype)
colnames(share_M) <- share_M_celltype
rownames(share_M) <- share_M_celltype

lfsr_res_sig <- lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
  mutate(cs = paste0(phenotype_id, "_", cs))

for (i in 1:14){
  for (j in 1:14){
    row_cell <- share_M_celltype[i]
    col_cell <- share_M_celltype[j]
    if (i == j){
      share_M[i, j] <- 1
    }else{
      row_cell_cs <- lfsr_res_sig %>% filter(celltype == row_cell) %>% distinct(cs) %>% pull(cs)
      col_cell_cs <- lfsr_res_sig %>% filter(celltype == col_cell) %>% distinct(cs) %>% pull(cs)
      both_cs <- intersect(row_cell_cs, col_cell_cs)
      union_cs <- union(row_cell_cs, col_cell_cs)
      share_M[i, j] <- length(both_cs) / length(union_cs)
    }
  }
  print(i)
}

share_M

library(RColorBrewer)
# library(grDevices)

# swap orders
share_M <- share_M[c(1:4, 6:8, 5, 12, 13, 9:11, 14), c(1:4, 6:8, 5, 12, 13, 9:11, 14)]

corrplot(share_M, method = c('color'), 
         type = "upper",
         is.corr = F, 
         col.lim = c(0, 1),
         COL1('YlOrRd', 200)
)

library("ggplotify")
p_share_M <- as.grob(~corrplot(share_M, method = c('color'), 
                            type = "upper",
                            is.corr = F, 
                            col.lim = c(0, 1),
                            COL1('YlOrRd', 200)
))

# pair-wise sharing
plot_grid(
  p_share_M,
  align = 'h',
  hjust = -1, # -1
  nrow = 1,
  axis="l",
  labels = c(""), label_size = 10, vjust=1
)
ggsave(paste0(plots_dir, "Fig_mvsusie_pairwise.png"),
        width = onecol_width*2, height = onecol_height*2, units = 'cm', dpi = 300)

## eQTLs in cell-type specific CSs
pip_res %>% left_join(lfsr_res %>% filter(lfsr < lfsr_thresh) %>% 
                         group_by(phenotype_id, cs) %>% 
                         count() %>% 
                         ungroup() %>% 
                         mutate(specific = ifelse(n < 2, 1, 0)) %>% 
                         mutate(cs = gsub("lfsr_","", cs)),
                       by=c("phenotype_id", "cs")) %>% 
  filter(!is.na(specific)) %>% 
  filter(pip >= 0.9) %>%
  ggplot(aes(x = as.factor(n), y=abs(tss_dist))) +
  xlab("number of significant cell types")+
  geom_boxplot(outlier.size = 0.5) + axis_theme

ggsave(paste0(plots_dir, "Fig_mvsusie_lfsr0.01_pip0.9_tss.png"),
       width = onecol_width*1.2, height = onecol_height, units = 'cm', dpi = 300)

