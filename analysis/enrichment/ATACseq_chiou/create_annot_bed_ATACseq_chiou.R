## prepare bed file for enrichment analysis

library(tidyverse)
library(data.table)

indir <- "/project/nmancuso_8/data/LDSC/baseline_annotation/bed/"
outdir <- "/project/nmancuso_8/elezhang/projects/jaxqtl/data/annotation/bed/"

allfiles <- list.files(indir)

bedfiles <- c()

for (file in allfiles){
  dat <- fread(paste0(indir, file), fill=TRUE, header=FALSE)
  print(colnames(dat))
  dat %>% mutate(V1 = gsub("chr", "", V1)) %>%
    fwrite(paste0(outdir, file), col.names = FALSE, sep = "\t")
  print(file)

  if (ncol(dat) == 3){
    bedfiles <- c(bedfiles, file)
  }
}

tibble(X = bedfiles) %>% write_tsv(paste0(outdir, "bedfile_list"), col_names = F)

# look at Chiou et.al 2021 supp table 6 for marker genes

# dat <- readxl::read_excel("/project/nmancuso_8/data/ANNOTATIONS/CHIOU_ET_AL/chiou_ATACseq_allcelltypes.xlsx", 
#                           skip=2) %>% 
#   janitor::clean_names()

dat <- readxl::read_excel("/project/nmancuso_8/data/ANNOTATIONS/CHIOU_ET_AL/41586_2021_3552_MOESM6_ESM.xlsx", 
                          skip=2) %>% 
  janitor::clean_names()


# cytotoxic_nk: NK
dat %>% 
  filter(cytotoxic_nk > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_NK.bed"), col_names = FALSE)

# adaptive_nk: NK_R
dat %>% 
  filter(adaptive_nk > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_NK_R.bed"), col_names = FALSE)

dat %>% 
  filter(cytotoxic_nk > 0 | adaptive_nk > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_NK_all.bed"), col_names = FALSE)

# activated_cd4_t, regulatory_t: CD4 family
dat %>% 
  filter(activated_cd4_t > 0 | regulatory_t > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_CD4_all.bed"), col_names = FALSE)

# memory_cd8_t, cytotoxic_cd8_t: CD8 family
dat %>% 
  filter(memory_cd8_t > 0 | cytotoxic_cd8_t > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_CD8_all.bed"), col_names = FALSE)

# naive_b:  B_IN
dat %>% 
  filter(naive_b > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_B_IN.bed"), col_names = FALSE)

# memory_b: B_Mem
dat %>% 
  filter(memory_b > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_B_Mem.bed"), col_names = FALSE)

dat %>% 
  filter(naive_b > 0 | memory_b > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_B_all.bed"), col_names = FALSE)

# classical_monocyte: Mono_C
dat %>% 
  filter(classical_monocyte > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_Mono_C.bed"), col_names = FALSE)

# non_classical_monocyte: Mono_NC
dat %>% 
  filter(non_classical_monocyte > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_Mono_NC.bed"), col_names = FALSE)

# conventional_dendritic: DC
dat %>% 
  filter(conventional_dendritic > 0) %>% 
  select(chrom, start, end) %>% 
  write_tsv(paste0(outdir, "ATACseq_DC.bed"), col_names = FALSE)

