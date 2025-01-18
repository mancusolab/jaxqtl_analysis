library(tidyverse)
library(data.table)

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl")

celltype <- "NK"
params <- read_tsv(paste0("./data/pheno/celltype16/metadata/genelist/", celltype, "/params"), col_names =  F)
colnames(params) <- c("chr", "chunk")

outdir <- paste0(paste0("./result/cis/celltype16/", celltype, "_fitnull/"))

Y_df <- fread(paste0("./data/pheno/celltype16/", celltype, ".bed.gz"), header = TRUE)

genelist <- read_tsv("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16/NK/notconverged_181.genelist", 
                     col_names =  F) %>% 
  pull(X1)

offset <- colSums(Y_df[, 5:986])

# all results
# !!! don't overwrite python output
for (i in  1:nrow(params)){
  chr <- params$chr[[i]]
  chunk <- params$chunk[[i]]

  onedat <- fread(paste0(outdir, "chr", chr, "/", chunk, ".fit.mu.tsv.gz"), header=TRUE)
  Y_sub <- Y_df %>% filter(Geneid %in% colnames(onedat))
  Gene_order <- Y_sub$Geneid

  onedat <- onedat %>% select(all_of(Gene_order))

  Ymat <- Y_sub %>% select(-c("#Chr", "start", "end","Geneid")) %>% t()

  if (nrow(Ymat) == nrow(onedat)){
    resid <- (Ymat - onedat)
    resid %>% fwrite(paste0(outdir, "chr", chr, "/", chunk, ".fit.resid.man.tsv.gz"), sep="\t")
    print(paste0("Finish ", i))
  }else{
    print("number of individuals not the same")
  }
}


# not converged
allres <- tibble(id = colnames(Y_df)[5:ncol(Y_df)])
for (i in  1:nrow(params)){
  chr <- params$chr[[i]]
  chunk <- params$chunk[[i]]
  
  onedat <- fread(paste0(outdir, "chr", chr, "/", chunk, ".fit.resid.tsv.gz"), header=TRUE)
  
  onedat <- onedat %>% select(genelist[genelist %in% colnames(onedat)])
  allres <- bind_cols(allres, onedat)
  if (ncol(onedat) > 0){
    print("Found")
  }
  print(paste0("Finish ", i))
}

allres %>% fwrite(paste0(outdir, "notconverged_181.resid.tsv.gz"), sep="\t")



# calculate variance of residuals and mean of mu
allres <- tibble()
for (i in 1:nrow(params)){
  chr <- params$chr[[i]]
  chunk <- params$chunk[[i]]
  
  dat_mu <- fread(paste0(outdir, "chr", chr, "/", chunk, ".pois.fit.mu.tsv.gz"), header=TRUE)
  dat_resid <- fread(paste0(outdir, "chr", chr, "/", chunk, ".pois.fit.resid.tsv.gz"), header=TRUE)
  
  dat_mu <- dat_mu %>% select(colnames(dat_resid))
  if (all.equal(colnames(dat_mu), colnames(dat_resid))){
    allres <- bind_rows(allres, 
                        tibble(fitmu_mean = apply(dat_mu, 2, mean), 
                               resid_var = apply(dat_resid, 2, var),
                               fitmu_mean_rm_offset = apply(dat_mu / offset, 2, mean),
                               resid_var_rm_offset = apply(dat_resid / offset, 2, var),
                               phenotype_id = colnames(dat_resid))) 
    print(paste0("Finish ", i))
  }
}

allres %>% write_tsv("./")
