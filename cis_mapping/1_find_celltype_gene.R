## find cell type and gene list
library(tidyverse)
library(data.table)

yazar_res <- readxl::read_excel("../data/OneK1K/yazar_result/science.abf3041_tables_s6_to_s19.xlsx", sheet = 5,
                                skip = 2) %>% 
  janitor::clean_names()

yazar_donor_cell_ct <- readxl::read_excel("../data/OneK1K/yazar_result/science.abf3041_tables_s6_to_s19.xlsx", sheet = 2,
                                       skip = 2) %>% 
  janitor::clean_names()


# gtf file
gtf <- read_tsv("../../jaxqtl/example/data/Homo_sapiens.GRCh37.87.bed.gz") %>% 
  rename(phenotype_id = gene_id)

yazar_res %>% distinct(gene_ensembl_id) %>% 
  inner_join(gtf, by = c("gene_ensembl_id" = "phenotype_id"))

yazar_res %>% filter(cell_type == "NK" & gene_ensembl_id == "ENSG00000178607") %>% View
yazar_res %>% filter(cell_type == "NK" & gene_ensembl_id == "ENSG00000132359")

yazar_res %>% filter(cell_type == "NK") %>% 
  group_by(chromosome) %>% 
  summarize(n=n()) %>% View

# 160 cis-genes on chr17 from yazar results
yazar_res %>% filter(cell_type == "NK") %>% 
  distinct(gene_ensembl_id, .keep_all = T) %>% 
  filter(chromosome==17) %>% 
  select(gene_ensembl_id) %>% 
  write_tsv("./cis_mapping/NK_chr17_genelist.tsv", col_names = F)

# 2340 cis-genes in total from yazar results
yazar_res %>% filter(cell_type == "NK") %>% 
  distinct(gene_ensembl_id, .keep_all = T) %>% 
  select(chromosome, gene_ensembl_id, position) %>% 
  left_join(gtf, by = c("gene_ensembl_id" = "phenotype_id")) %>% 
  mutate(dist = abs(position - start)) %>% filter(dist < 500000) %>% 
  write_tsv("./cis_mapping/NK_cis_genelist.t2sv", col_names = F)

yazar_res %>% filter(cell_type == "NK") %>% 
  distinct(gene_ensembl_id, .keep_all = T) %>% 
  select(gene_ensembl_id, position) %>% 
  inner_join(jaxqtl, by=c("gene_ensembl_id" = "phenotype_id")) %>% View
  filter(is.na(pval_beta))

yazar_donor_cell_ct %>% filter(cell_type == "NK") %>% 
  pull(number_of_cells) %>% 
  sum()


n981 <- read_tsv("../data/OneK1K/yazar_result/donor_n981_celltype_count.tsv") %>% 
  janitor::clean_names()

table(res$V1)

# count number of eSNPs after 5 rounds of conditional analysis
yazar_res %>% group_by(cell_type) %>% summarise(n = n()) %>% View
n981 %>% distinct(cell_type) %>% View

CD4nc <- c("central memory CD4-positive, alpha-beta T cell", "naive thymus-derived CD4-positive, alpha-beta T cell")
CD4et <- c("effector memory CD4-positive, alpha-beta T cell")
Bmem <- c("memory B cell")
NK <- c("natural killer cell")
plasma <- c("plasmablast")
dc <- c("dendritic cell", "plasmacytoid dendritic cell", "conventional dendritic cell")

n981 %>% filter(cell_type %in% NK) %>% pull(n) %>% sum()

