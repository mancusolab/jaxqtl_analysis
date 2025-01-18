# need --mem-per-cpu 60Gb for this script
import os
import jax
from jaxqtl.io.pheno import SingleCellFilter, H5AD
import numpy as np

jax.config.update("jax_enable_x64", True)
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/data/")

raw_count_path = "./pheno/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz"  # n=981
celltype_path = "./pheno_meta/celltype_16.tsv"
out_dir = "./pheno/celltype16_new/"

# update column names
SingleCellFilter.id_col = "individual"
SingleCellFilter.celltype_col = "cell_label"
SingleCellFilter.mt_col = "percent.mt"
SingleCellFilter.geneid_col = "Geneid"

# update filtering metric so that every cell is kept
SingleCellFilter.min_cells = -np.Inf
SingleCellFilter.min_genes = -np.Inf
SingleCellFilter.max_genes = np.Inf
SingleCellFilter.percent_mt = 100.0
SingleCellFilter.bulk_method = "sum"


# Prepare input #
# For given cell type, create bed files from h5ad file
pheno_reader = H5AD()
dat = pheno_reader(raw_count_path)
print("Finish read in raw data.")

# relabel part of 870_871 cells as 966_967
dat.obs['individual'] = dat.obs.individual.astype(str)
mask = (dat.obs.individual == "870_871") & (dat.obs.latent == 1)
dat.obs.loc[mask, 'individual'] = "966_967"
dat.obs['individual'] = dat.obs.individual.astype('category')

# norm factor 1e5 for linear model
dat_df = pheno_reader.process(dat, SingleCellFilter,
                                divide_size_factor=False)
print("Finish processing raw data.")

# 28196 genes
# now use collapsed model v82
# old use Homo_sapiens.GRCh37.87.bed.gz
print("Writing bed files")
pheno_reader.write_bed(
    dat_df,
    gtf_bed_path="./pheno_meta/Homo_sapiens.GRCh37.82.bed.gz",
    out_dir=out_dir,
    celltype_path=celltype_path,
    autosomal_only=True
)

# import hdf5plugin
# dat_df.reset_index().to_csv("../data/pheno/yazar2022.pseudo.fix1e5.mean.gz", sep="\t", index=False)
# 
# # make file for after filtered
# pheno_reader.write_bed(
#     count_df,
#     gtf_bed_path="../data/pheno_meta/Homo_sapiens.GRCh37.87.bed.gz",
#     out_dir="../data/pheno/bed/NK/",
#     celltype_path="../data/pheno_meta/NK.tsv",
#     autosomal_only=True,
# )
