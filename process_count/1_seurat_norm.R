# Script information ------------------------------------------------------

# title: Normalize OneK1k data
# author: Jose Alquicira Hernandez
# date: 2019/07/10
# description: Data normalization is performed using the new SCTransform method
# from Seurat.
# See:
# - https://satijalab.org/seurat/v3.0/sctransform_vignette.html

# Import libraries --------------------------------------------------------

# Primary
library("tidyverse")

# Secondary
library("Seurat")


# Set output --------------------------------------------------------------

output <- set_output("2019-08-23", "norm_sctransform")

# Read data ---------------------------------------------------------------

inicio("Reading gene expression data")
input <- "../data/OneK1K/Count.rds"
data <- readRDS(input)
data[['RNA']]
data[['percent.mt']]
data[['pool']] <- as.factor(data@meta.data$pool_number)
head(data[['RNA']]@misc)
head(data@meta.data$percent.mt)

data[['RNA']]

# Normalize data ----------------------------------------------------------

print("Normalizing data using SCTransform")
options(future.globals.maxSize = 1024**3 * 2500)
data <- SCTransform(object = data, vars.to.regress = c("percent.mt", "pool"), conserve.memory = TRUE)
fin()

# Save output -------------------------------------------------------------

print("Saving normalized data")
saveRDS(data, file = output, "norm.RDS")
fin()

# Session info ------------------------------------------------------------
print_session(output)