###################################
#     process raw count data      #
###################################

library(tidyverse)
library(Seurat)

dat <- readRDS("../data/OneK1K/Count.rds")
dim(dat) # 36571 genes x 1248980 cells 

dat[['RNA']]@data # raw count data
dat[['RNA']]@meta.features # meta feature for genes: 36571 genes x 9 features

# switch slots, put rawdata into "counts" slot 
count.data <- GetAssayData(object = dat[["RNA"]], slot = "data")
SetAssayData(object = dat, slot = "counts", new.data = count.data)

# normalize by relative count
dat <- NormalizeData(dat, normalization.method = "RC", scale.factor = 1e6)
levels(x = dat)


###################################
#     qc statistics      # 
# already done
# 1. find singlets: 1,420,567 -> 1,295,408
# 2. remove outlier cells -> 1,272,518
###################################

# 
# average read depth: should be ~34000
dat[['SCT']]@
