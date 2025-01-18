# calculate expression PCs

library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/")

celltypes <- read_tsv("../../pheno_meta/celltype_14.tsv", F) %>% pull(X1)
celltypes <- gsub(" ", "_", celltypes)
celltypes <- c(celltypes, "allcells")

out <- data.frame(celltype = celltypes, PC1 = NA, PC2 = NA, PC3 = NA, PC4 = NA, PC5 = NA)

for (i in seq_along(out$celltype)){
  cell <- out$celltype[i]
  df <- fread(paste0(cell, ".bed.gz"), header=T)
  df <- t(df[, 5:ncol(df)])
  df <- df[,colSums(df) > 0]
  res <- prcomp(df, center = TRUE, scale. = TRUE)
  summary(res)
  eigs <- res$sdev^2
  var_pcs <- eigs / sum(eigs)
  
  out[i, 2:6] <- var_pcs[1:5]
  print(i)
}

out %>% write_tsv("./metadata/count_5PCs_var.tsv")
