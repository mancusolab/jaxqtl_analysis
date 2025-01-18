#!/bin/bash

# read in and gather results
fdr_val <- 0.05
threshold <- 0.01

annot_dir <- "../../jaxqtl_project/data/OneK1K/annotation/"

# individual cell type proportion
# individuals per cell type are sequenced from the same pool (individual_pool.tsv.gz on hpc)
ind_celltype_ct <- read_tsv("../../jaxqtl_project/data/OneK1K/pheno/celltype16/metadata/celltype_ind_counts.tsv.gz") %>% 
  filter(!cell_label %in% c("Platelets", "Erythrocytes")) %>% 
  mutate(celltype = gsub(" ", "_", cell_label)) %>% 
  select(-cell_label)

celltype_order <- ind_celltype_ct %>% group_by(celltype) %>% summarize(ct = sum(counts)) %>% 
  arrange(desc(ct)) %>% pull(celltype)

celltypes_grouped <- c("CD4_NC", "CD4_ET", "CD4_SOX4",
                       "CD8_NC", "CD8_ET", "CD8_S100B",
                       "NK", "NK_R", "Plasma",
                       "B_Mem", "B_IN",
                       "Mono_C", "Mono_NC",
                       "DC")

# read in gene annotation data

# gtf: keep only autosomal chr
# !!!! don't use start and end for TSS from this file because TSS is not processed
gtf_raw <- fread("../../jaxqtl_project/data/gene_info/Homo_sapiens.GRCh37.82.gtf.parse.gz") %>% 
  filter(type == "transcript") %>%
  mutate(gene_length = abs(start - end) + 1) %>% 
  select(-c(start, end))

# only keep autosomal chr
gtf <- fread("../../jaxqtl_project/data/gene_info/Homo_sapiens.GRCh37.82.bed.gz") %>% 
  left_join(gtf_raw, by = c("chr"="seqnames", "gene_id")) %>% 
  rename(phenotype_id = gene_id) %>% 
  filter(chr %in% as.character(1:22)) %>% 
  mutate(chr = as.integer(chr))

# add gene cytoband info to gtf file
gene_band <- fread("../../jaxqtl_project/data/gene_info/ensembl_gene_cytoband.tsv.gz") %>% 
  select(-external_gene_name) %>% 
  filter(chromosome_name %in% as.character(1:22)) %>% 
  mutate(chromosome_name = as.integer(chromosome_name))

gtf <- gtf %>% left_join(gene_band, by=c("phenotype_id"="ensembl_gene_id", "chr"="chromosome_name"))

# chr6 not MHC
window <- 500000
notMHC <- gtf %>% 
  filter(chr == 6) %>% 
  filter((end + window) < 25000000 | (end - window) > 34000000) %>% 
  distinct(phenotype_id)

MHC_genes <- gtf %>% 
  filter(chr == 6) %>% 
  anti_join(notMHC) %>% 
  distinct(phenotype_id) %>% pull(phenotype_id)
  # write_tsv("../data/gene_info/MHC_genes.tsv", col_names = F)

# MAPT region (large inversion on chr17)
MAPT_genes <- gtf %>% filter(chr == "17" & band == "q21.31") %>% pull(phenotype_id)
# tibble(genes=MAPT_genes) %>% write_tsv("../data/gene_info/MAPT_genes.tsv", col_names = F)

# hgnc <- fread("../../jaxqtl_project/data/gene_info/hgnc_complete_set_2024-03-06.txt") %>% 
#   separate(location, into = c("chr", "to_rm"), sep="q|p") %>% 
#   mutate(chr = as.integer(chr))
# 
# hgnc_new <- fread("../../jaxqtl_project/data/gene_info/hgnc_2024-04-09.txt") %>% 
#   janitor::clean_names()

# consistent with v82 TSS; GeneSymbol used by onek1k data
gene_lookup <- read_tsv("../../jaxqtl_project/data/gene_info/gene_autosome_lookup.tsv.gz") %>% 
  select(phenotype_id=Geneid, GeneSymbol) %>% 
  left_join(gtf, by=c("phenotype_id"))

gene_lookup %>% add_count(GeneSymbol) %>% filter(n<2) # 31350 unique GeneSymbol
gene_lookup %>% add_count(GeneSymbol) %>% filter(n>1) %>% distinct(GeneSymbol) # 172-84 = 88 duplicate genes

# 79 genes have duplicate gene symbol; 
# gene_lookup %>%
#   mutate(GeneSymbol = gsub("-","\\.", GeneSymbol)) %>%
#   add_count(chr, GeneSymbol) %>% filter(n>1) %>% distinct(chr, GeneSymbol) %>%
#   distinct(chr, phenotype_id)

# cell meta data
cellmeta <- read_tsv("../../jaxqtl_project/data/OneK1K/pheno/celltype16/metadata/celltype_14.tsv", F) %>% 
  rename(celltype = X1) %>% 
  mutate(celltype = gsub(" ", "_",celltype),
         lineage = c("CD4", "Sub-Lymphoid", "CD4", "CD8", "CD8", "B cells", "CD8", "B cells",
                     "Sub-Lymphoid", "Monocyte", "Monocyte", "Myeloid", "B cells", "CD4"),
         group = c("T cells", "NK", rep("T cells", 3), "B cells", "T cells", "B cells",
                   "NK", rep("Monocyte", 2), "Myeloid", "B cells", "T cells"))

celltypes <- cellmeta$celltype

# mean_counts, max_min, proportion
cell_counts <- read_tsv("../../jaxqtl_project/data/OneK1K/pheno/celltype16_new/metadata/celltype_prop.tsv.gz")
celltype_N <- read_tsv("../../jaxqtl_project/data/OneK1K/pheno/celltype16_new/metadata/celltype_samplesize.tsv")
cell_counts <- cell_counts %>% 
  filter(!cell_label %in% c("Erythrocytes", "Platelets")) %>% 
  mutate(prop = counts/total_counts, cell_label = gsub(" ", "_", cell_label)) %>% 
  rename(celltype = cell_label) %>% 
  group_by(celltype) %>% 
  summarize(mean_ct = mean(counts), 
            min_ct = min(counts),
            max_ct = max(counts),
            mean_prop = mean(prop),
            total_cell = sum(counts)) %>% 
  left_join(celltype_N)

# pLI score (for 19224 protein coding genes) and EDS score
# now change to LOEUF < 0.6 suggestive of more intolerance
# readme: https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/constraint/README.txt
pLI_score <- read_tsv(paste0(annot_dir, "gnomad.v4.0.constraint_metrics.tsv")) %>% 
  group_by(gene) %>% 
  summarize(loeuf = mean(lof.oe_ci.upper, na.rm = T),
            pLI = mean(lof.pLI, na.rm = T)) %>% 
  ungroup() %>% 
  inner_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by=c("gene" = "GeneSymbol"))

eds_score <- readxl::read_excel(paste0(annot_dir, "Wang_EDS.xlsx")) %>% 
  rename(phenotype_id=GeneSymbol) %>% 
  distinct(phenotype_id, .keep_all = T) %>%
  select(phenotype_id, EDS)

# percentile can be calcuated by ecdf(score)$scpre * 100
rvis_score <- readxl::read_excel(paste0(annot_dir, "RVIS_petrovsi.xlsx")) %>% 
  rename(RVIS="Residual Variation Intolerance Score",
         RVIS_percentile="Residual Variation Intolerance Score Percentile") %>% 
  left_join(gene_lookup %>% select(phenotype_id, GeneSymbol), by=c(`HGNC gene`="GeneSymbol"))

# ABC score
# abc <- fread(paste0(annot_dir, "abc_v3/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"))

# check duplicate variants with the same ma_count
# note: for this multi-allele snps, it is fine if the alt allele is diff, because the effect is the same
# for 
geno_rsid <- fread("../data/OneK1K/geno/allvar_rsid.tsv.gz") %>% 
  select(-c(cM)) %>% 
  # select(-c(cM, variant_id)) %>% 
  add_count(chr, pos, ref)

geno_rsid_single <- geno_rsid %>% filter(n<2) %>% select(-n)
geno_rsid_double <- geno_rsid %>% filter(n>1) %>% select(-n)

# re-annotate rsid (same ref but different alt)
geno_rsid_double <- geno_rsid_double %>% 
  group_by(chr, pos, ref) %>% 
  mutate(rsid = first(rsid)) %>% 
  ungroup()

geno_rsid <- bind_rows(geno_rsid_single, geno_rsid_double) %>% 
  mutate(variant_id = paste0("chr", chr, "_", pos, "_", ref, "_", alt, "_b37"))

# cell type match annotation
atacseq_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "ATACseq_chiou")
epimap_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "EpiMap")
scent_map <- readxl::read_excel(paste0(annot_dir, "celltype_match.xlsx"), sheet = "SCENT")
