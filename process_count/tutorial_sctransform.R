# tutorail for sctransform
library(Seurat)
library(SeuratData)
library(tidyverse)
library(sctransform)

InstallData("ifnb")

LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

ctrl <- ifnb.list[["CTRL"]]
stim <- ifnb.list[["STIM"]]


### my data

pbmc[['percent.mt']]
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE, conserve.memory = FALSE)

ctrl <- SCTransform(ctrl, verbose = FALSE)
ctrl[['SCT']]
ctrl
