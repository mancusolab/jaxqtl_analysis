library(tidyverse)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("-m", "--method"), type="character", default=NULL, 
              help="jaxqtl or tqtl", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=NULL, 
              help="filter PIP >= cutoff", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

method=opt$method
pip_cutoff=opt$cutoff

if (method == "jaxqtl"){
  prefix="nb"
  model="nb"
}else{
  prefix=method # "tqtl", jaxqtl_lm_score
  model="lm" # for tqtl or jaxqtl_lm_score
}

params=read_tsv(paste0("../../finemap/code/", prefix ,"_allegenes.tsv"), F) %>% select(X1:X3)
colnames(params)=c("phenotype_id", "chr", "celltype")  #, "N", "chunk")

allfiles=list.files(method)

celltypes=unique(params$celltype)

# make bed file for SNPs with PIP >= cutoff
if (pip_cutoff > 0){
  for (cell in celltypes){
    
    params_cell=params %>% filter(celltype==cell)
    out=data.frame()
    
    for (idx in 1:nrow(params_cell)){
      
      gene=params_cell$phenotype_id[idx]
      prefix=paste0(gene, ".", cell)
      
      file=allfiles[grepl(glob2rx(paste0(prefix, ".", model, ".*tsv.gz")), allfiles)]
      
      if (length(file) > 0){
        
        df=fread(paste0(method,"/", file), header=T)
        
        df_pass=df %>% filter(pip >= pip_cutoff)
        
        if (nrow(df_pass)>0){
          out=bind_rows(out,df_pass %>% mutate(pos_1 = pos-1) %>%
                          select(chrom, pos_1, pos, snp, pip) %>%
                          mutate(gene=gene))
          print("found PIP > threshold")
        }
      }
    }
    print(cell)
    out %>% write_tsv(paste0(method, "/bed/", cell, ".PIP", pip_cutoff, ".bed"), col_names=F)
  }
}else{
  # make bed file for SNPs in C
  for (cell in celltypes){
    
    params_cell=params %>% filter(celltype==cell)
    out=data.frame()
    
    for (idx in 1:nrow(params_cell)){
      print(idx)
      gene=params_cell$phenotype_id[idx]
      prefix=paste0(gene, ".", cell)
      
      file=allfiles[grepl(glob2rx(paste0(prefix, ".", model, "*.tsv.gz")), allfiles)]
      cs_file=allfiles[grepl(glob2rx(paste0(prefix, ".", model, "*.cs.tsv")), allfiles)]
      
      if (length(file) > 0 && length(cs_file) > 0){
        
        df=fread(paste0(method,"/", file), header=T)
        
        cs_union=rowSums(df %>% select(starts_with("cs")))>0
        
        out=bind_rows(out, df[cs_union,] %>% mutate(pos_1 = pos-1) %>%
                        select(chrom, pos_1, pos, snp, pip) %>%
                        mutate(gene=gene))
      }
    }
    print(cell)
    out %>% write_tsv(paste0(method, "/bed/", cell, ".cs.bed"), col_names=F)
  }
}

# concatenate file in output directory using combine.sh