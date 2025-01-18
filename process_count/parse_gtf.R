library(rtracklayer)
library(data.table)
my_obj <- import("../data/gene_info/Homo_sapiens.GRCh37.87.gtf.gz")

as.data.frame(my_obj) %>% fwrite("../data/gene_info/Homo_sapiens.GRCh37.87.gtf.parse.gz", sep="\t")
