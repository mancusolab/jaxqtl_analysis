### count data
library(tidyverse)
library(data.table)
library(moments)

wkdir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new"; 
setwd(wkdir)

celltype_path="./metadata/celltype_14.tsv"

# args <- commandArgs(trailingOnly=TRUE)
# wkdir <- args[1] # wkdir="/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16"
# celltype_path <- args[2] # celltype_path="./metadata/celltype_14.tsv"
# percent_threshold <- as.numeric(args[3]) # 0.0
# chunk <- args[3] # 50

percent_threshold <- 0
chunk <- 50

# genelength
genelen <- fread("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/Homo_sapiens.GRCh37.82.collapse.genelength.gz") %>%
  rename(Geneid = gene)

# read cell type file
all_celltype <- read_tsv(celltype_path, F) %>% 
  rename(celltype = X1) %>% 
  mutate(celltype = str_replace(celltype, " ", "_"))

allcelltypes <- c(all_celltype$celltype, "allcells")

for (celltype in allcelltypes){
  dat <- fread(paste0(celltype, ".bed.gz"), header = TRUE)
  # align gene order with those in dat
  genelen <- dat %>% select(Geneid) %>% left_join(genelen, by = c("Geneid"))
  
  count <- dat[,5:ncol(dat)] # gene x iid
  
  libsize <- colSums(count)
  rate <- t(t(count)/libsize) # adjust library size
  rate_mean <- rowMeans(rate)
  rate_var <- apply(rate, 1, var)
  rate_ratio <- rate_mean/rate_var
  
  count_sum <- rowSums(count) # read counts across people
  express_percent <- rowMeans(count > 0)
  express_percent_6 <- rowMeans(count >= 6)
  count_mean <- rowMeans(count)
  count_var <- apply(count, 1, var)
  ratio <- count_mean / count_var
  
  skew <- apply(count, 1, skewness)
  kurt <- apply(count, 1, kurtosis)
  
  # calculate TPM (use mean gene length here)
  # see ref: https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w
  count_divide_len <- count/genelen$mean # divide each row by a value
  sum_gene <- colSums(count_divide_len, na.rm = TRUE) # sum over genes per individual
  TPM <- t(t(count_divide_len) * 1000000/sum_gene) # colSums(TPM) should all be equal to 1e6
  tpm_0.1 <- rowMeans(TPM >= 0.1)
  
  gene_summary <- tibble(chr = dat$`#Chr`,
                         Geneid = dat$Geneid,
                         count_sum = count_sum,
                         express_percent = express_percent,
                         express_percent_6 = express_percent_6,
                         count_mean = count_mean,
                         count_var = count_var,
                         mean_var_ratio = ratio,
                         tpm_0.1 = tpm_0.1,
                         rate_mean = rate_mean,
                         rate_var = rate_var,
                         rate_mean_var_ratio = rate_ratio,
                         skew = skew,
                         kurt = kurt)
  
  gene_summary %>% write_tsv(paste0("./metadata/", celltype, "_gene_summary.tsv.gz"))
  
  print(paste0("Finish: ", celltype))
}


for (celltype in allcelltypes){
  gene_summary <- read_tsv(paste0("./metadata/", celltype, "_gene_summary.tsv.gz"))
  gene_summary <- gene_summary %>% filter(express_percent > percent_threshold)
  
  if (!(dir.exists(paste0("./metadata/genelist/", celltype)))){
    dir.create(paste0("./metadata/genelist/", celltype))
  }
  
  for (chr_idx in unique(gene_summary$chr)){
    gene_summary_chr <- gene_summary %>% filter(chr == chr_idx)
    n <- nrow(gene_summary_chr)
    if (n < chunk){
      r <- rep(1, n)
    }else{
      r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
    }
    d <- split(gene_summary_chr,r)
    
    if (!(dir.exists(paste0("./metadata/genelist/", celltype, "/chr", chr_idx)))){
      dir.create(paste0("./metadata/genelist/", celltype, "/chr", chr_idx))
    }
    
    for (i in 1:length(d)){
      d[[i]] %>% 
        select(Geneid) %>% 
        write_tsv(paste0("./metadata/genelist/", celltype, "/chr", chr_idx, "/chunk_", i), col_names = F)
    }
    print(chr_idx)
  }
}

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/genelist")

params <- data.frame()
for (celltype in allcelltypes){
  allchr <- list.files(celltype)
  for (chr in allchr){
    chr_idx <- gsub("chr", "", chr)
    allfiles <- list.files(paste0(celltype, "/",chr))
    tmp <- tibble(cell_type = celltype, chr_col = chr_idx, files = allfiles)
    params <- bind_rows(params, tmp)
    print(chr)
  }
  print(celltype)
}

params %>% write_tsv("params_all_celltype", col_names = F)

# cell type
for (celltype in allcelltypes){
  params %>% filter(cell_type == celltype) %>% 
    write_tsv(paste0(celltype, "/params"), col_names=F)
  print(celltype)
}

# make result directory

setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new")

for (celltype in allcelltypes){
  if (!(dir.exists(paste0("./", celltype)))){
    dir.create(paste0("./", celltype))
  }
  
  for (chr_idx in 1:22){
    if (!(dir.exists(paste0("./", celltype, "/chr", chr_idx)))){
      dir.create(paste0("./", celltype, "/chr", chr_idx))
    } 
  }
}

# make directory for tensorqtl
setwd("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_tensorqtl_new")
for (celltype in allcelltypes){
  if (!(dir.exists(paste0("./", celltype)))){
    dir.create(paste0("./", celltype))
  }
}

