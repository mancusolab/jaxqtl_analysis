# add band info for genes
library(biomaRt)
library(glue)
library(RSQLite)
library(data.table)

# https://lab-notes.hakyimlab.org/post/2021-07-28-how-to-get-the-cytogenetic-band-of-a-gene/
# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl",
                      # mirror = "useast",
                      GRCh = 37)
## get the gene annotation with cytoband  from biomart
anno_gene <- getBM(attributes =c("ensembl_gene_id","external_gene_name","chromosome_name","band"),mart=ensembl,
                   verbose = T)
head(anno_gene)

anno_gene %>% write_tsv("../../jaxqtl_project/data/gene_info/ensembl_gene_cytoband.tsv.gz")

gtf <- fread("../../jaxqtl_project/data/gene_info/Homo_sapiens.GRCh37.82.bed.gz") %>% 
  rename(phenotype_id = gene_id)



