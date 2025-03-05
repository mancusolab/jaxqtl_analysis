## clean up sldsc results
library(tidyverse)
library(data.table)

method <- "jaxqtl"

setwd(paste0("/project/nmancuso_8/elezhang/projects/jaxqtl/result/finemap/sldsc_wald_label/", method, "/sldsc"))

# res <- read_tsv(paste0(method, "_sldsc_traits107_union.label.tsv.gz"))
res <- read_tsv(paste0(method, "_sldsc_traits107_union.label.newselection.tsv.gz"))


colnames(res)[ncol(res)] <- "filename"
res <- res %>% filter(Category != "Category")

ld <- res %>% filter(grepl("baselineLD", filename))
nold <- res %>% filter(!grepl("baselineLD", filename))

nrow(ld) + nrow(nold) == nrow(res)

ld <- ld %>% separate(filename, into=c("trait", "suffix"), sep=".baselineLD.", remove = F) %>% 
  mutate(suffix = gsub(".results$", "", suffix)) %>% 
  separate(suffix, into=c("annot", "celltype"), sep="\\.")

nold <- nold %>% separate(filename, into=c("trait", "suffix"), sep=".baseline.", remove = F) %>% 
  mutate(suffix = gsub(".results$", "", suffix)) %>% 
  separate(suffix, into=c("annot", "celltype"), sep="\\.")

bind_rows(ld, nold) %>% 
  #write_tsv(paste0(method,"_sldsc_traits107_union.label.tsv.gz"))
  write_tsv(paste0(method,"_sldsc_traits107_union.label.newselection.tsv.gz"))


# gather h2 from sldsc output .log
allfiles <- list.files()
allfiles <- allfiles[grepl(".log$", allfiles)]
allfiles <- allfiles[!grepl("baselineLD", allfiles)]

outdf <- tibble(filename = allfiles,
                h2g = NA,
                h2g_se = NA)

# TODO: calculate h2 and h2Z (look for traits with h2Z > 6)
for (idx in 1:nrow(outdf)){
  file <- outdf$filename[[idx]]
  log <- read.table(file,h=F,fill=T)
  if (length(as.character(log[which(log$V4=="h2:"),5])) > 0){
    outdf$h2g[[idx]] <- as.numeric(as.character(log[which(log$V4=="h2:"),5]))
    h2g_se <- as.character(log[which(log$V4=="h2:"),6])
    h2g_se <- as.numeric(gsub("\\(|\\)", "", h2g_se))
    outdf$h2g_se[[idx]] <- h2g_se
  }
  print(idx)
}

# outdf %>% write_tsv(paste0(method, "_sldsc_traits107_h2_baseline.tsv.gz"))
outdf %>% write_tsv(paste0(method, "_sldsc_traits107_h2_baseline.newselection.tsv.gz"))

