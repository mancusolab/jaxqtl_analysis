# extract iHS score for given population

library(tidyverse)

pops="ESN|GWD|LWK|MSL|YRI|ACB|ASW|CLM|MXL|PEL|PUR|CDX|CHB|CHS|JPT|KHV|CEU|FIN|GBR|IBS|TSI|BEB|GIH|ITU|PJL|STU"
pops=str_split_1(pops, "\\|")
ceu_idx=which(pops == "CEU")

df=read_tsv("chr10.1kg.p3.allPops.iHS.txt", skip=3)

# get Derived allele frequency
df %>% separate(DAF, into=paste0("V",1:length(pops)), sep="\\|")
