### plot annotations

library(dplyr)
library(plotgardener)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(cicero)

# gwas RA
gwas <- data.table::fread("../data/gwas/Ishigaki2022/sumstats/RA_Ishigaki2022.tsv.gz") %>% 
  mutate(chr = paste0("chr", chr))
eqtl <- data.table::fread("../result/coloc/RA/jaxqtl.cis_qtl_pairs.nb.wald.tsv.gz")

assembly_name <- "hg19"
chr <- "chr5"
gene <- "IL6ST"
ensemid <- "ENSG00000134352"
celltype <- "CD8_NC"
leadsnp <- "chr5_55444683_G_A_b37"

df_eqtl <- eqtl %>% 
  filter(phenotype_id == ensemid) %>%
  # filter(phenotype_id == ensemid & celltype == .env$celltype) %>%
  filter(converged == TRUE)

# take cis-eQTLs from this gene
df_eqtl <- df_eqtl %>% 
  dplyr::select(chrom, pos, p=pval_nominal, snp, celltype) %>% 
  mutate(chrom = paste0("chr", chrom)) %>% 
  mutate(sig = ifelse(snp == leadsnp, 1, 0)) %>% 
  as_data_frame(); nrow(df_eqtl)

# take gwas summary stats for this gene (same variants as eQTLs)
df_gwas <- gwas %>% 
  dplyr::select(chrom=chr, pos, p=pval, snp) %>% 
  filter(snp %in% unique(df_eqtl$snp)) %>% 
  mutate(sig = ifelse(snp == leadsnp, 1, 0)) %>% 
  as_data_frame(); nrow(df_gwas)

# take colocalized causal eQTL
rsid <- "rs7731626" # "rs488540" # "rs2304204"
lead_pos <- 55444683 # 58280942 # 114884018 # 91002927 # 159701708 # A > G

tss_start <- 55290821
chr_start <- tss_start - 500000
chr_end <- tss_start + 500000

col1 <- "#882E72"
col2 <- "#5289C7" # "#36648B"
col3 <- "#F6C141"
col4 <- "#4EB265"
col5 <-  "#F1932D"

indir <- paste0("../result/annotation/Epigen")
# indir <- paste0("../result/annotation/Benaglio_2022/genome_browser_tracks")
plots_dir <- paste0("../result/annotation/figures/")

# df_chipseq1 <- readBigwig(paste0(indir, "/pbmc1_15.naive_cd4_t.scale_1e6.bw"),
#                           chrom=chr, chromstart=chr_start, chromend=chr_end)
# 
# df_chipseq2 <- readBigwig(paste0(indir, "/pbmc1_15.naive_cd8_t.scale_1e6.bw"),
#                           chrom=chr, chromstart=chr_start, chromend=chr_end)
# 
# df_chipseq3 <- readBigwig(paste0(indir, "/pbmc1_15.mem_b.scale_1e6.bw"),
#                           chrom=chr, chromstart=chr_start, chromend=chr_end)
# 
# df_chipseq4 <- readBigwig(paste0(indir, "/pbmc1_15.adaptive_NK.scale_1e6.bw"),
#                           chrom=chr, chromstart=chr_start, chromend=chr_end)
# 
# df_chipseq5 <- readBigwig(paste0(indir, "/pbmc1_15.cMono.scale_1e6.bw"),
#                           chrom=chr, chromstart=chr_start, chromend=chr_end)
# 

df_chipseq1 <- readBigwig(paste0(indir, "/Tcell/H3K27ac_CD4_T_ENCFF313TWH.bigWig"),
                          chrom=chr, chromstart=chr_start, chromend=chr_end)


df_chipseq2 <- readBigwig(paste0(indir, "/Tcell/H3K27ac_CD8_T_ENCFF635YOQ.bigWig"),
                          chrom=chr, chromstart=chr_start, chromend=chr_end)

df_chipseq3 <- readBigwig(paste0(indir, "/Bcell/H3K27ac_Bcell_ENCFF071MEQ.bigWig"),
                          chrom=chr, chromstart=chr_start, chromend=chr_end)

df_chipseq4 <- readBigwig(paste0(indir, "/NK/H3K27ac_NK_ENCFF814VKT.bigWig"),
                          chrom=chr, chromstart=chr_start, chromend=chr_end)

df_chipseq5 <- readBigwig(paste0(indir, "/Mono_C/H3K27ac_Mono_ENCFF468QFA.bigWig"),
                          chrom=chr, chromstart=chr_start, chromend=chr_end)


# GM12878 is a lymphoblastoid cell line

# HiC <- data.table::fread("../result/annotation/HiC/Rao_2014.GM12878.hg19.peakachu-merged.loops") %>%
#   dplyr::select(V1:V6)
# colnames(HiC) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
# HiC

# PCHiC data
# bait is promoter (eg. enhancer),oe is the other promoter interacting region (PIR)
pchic <- data.table::fread("../result/annotation/Javierre_2016/PCHiC_peak_matrix_cutoff5.tsv") %>% 
  #dplyr::filter(tB >=5) %>% 
  dplyr::select(chrom1 = baitChr, start1 = baitStart, end1 = baitEnd,
                chrom2 = oeChr, start2= oeStart, end2 = oeEnd,
                nCD4, nCD8) %>% 
  mutate(chrom1 = paste0("chr", chrom1),
         chrom2 = paste0("chr", chrom2)) %>% 
  tidyr::gather(key = celltype, value = score, c(nCD4, nCD8))

HiC <- pchic

# # co-accessability
# # position are on read as forward strand
# coac <- fread("../result/annotation/Benaglio_2022/cicero/t.filtered_coaccessible_sites4s.txt")
# coac <- coac %>% 
#   filter(Peak2 %in% c("chr5_55443409_55444796") &
#            Peak1 %in% c("chr5_55273640_55274515",
#                         "chr5_55278757_55279291")) %>% 
#   separate(Peak1, into=c("chrom1", "start1", "end1"), sep="_") %>% 
#   separate(Peak2, into=c("chrom2", "start2", "end2"), sep="_") %>% 
#   mutate(start1 = as.integer(start1), start2 = as.integer(start2),
#          end1 = as.integer(end1), end2 = as.integer(end2))
# HiC <- coac


# read GTF file for cicero
# raw: Homo_sapiens.GRCh37.82.gtf.gz
# after collapse: Homo_sapiens.GRCh37.82.genes.gtf
gene_anno <- rtracklayer::readGFF("../data/gene_info/Homo_sapiens.GRCh37.82.genes.gtf")

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

params_c <- pgParams(chrom = chr, chromstart = chr_start,
                     chromend = chr_end, 
                     assembly = assembly_name,
                     x = 3.5, width = 1.5, default.units = "inches")

ctcf_range1 <- pgParams(range = c(min(df_chipseq1$score), max(df_chipseq1$score)), 
                        assembly = assembly_name)

ctcf_range2 <- pgParams(range = c(min(df_chipseq2$score), max(df_chipseq2$score)),
                        assembly = assembly_name)

ctcf_range3 <- pgParams(range = c(min(df_chipseq3$score), max(df_chipseq3$score)),
                        assembly = assembly_name)

ctcf_range4 <- pgParams(range = c(min(df_chipseq4$score), max(df_chipseq4$score)),
                        assembly = assembly_name)

ctcf_range5 <- pgParams(range = c(min(df_chipseq5$score), max(df_chipseq5$score)),
                        assembly = assembly_name)


