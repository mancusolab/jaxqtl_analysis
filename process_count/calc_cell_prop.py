# calculate cell type proportion per person

# need --mem-per-cpu 60Gb for this script
import os
from jaxqtl.io.pheno import SingleCellFilter, H5AD
import pandas as pd
import numpy as np

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

df = dat.obs
df = df[~df['cell_label'].isin(['Erythrocytes','Platelets'])]

group_counts = df.groupby(['individual','cell_label']).size().reset_index(name='counts')
total_counts = df.groupby(['individual']).size().reset_index(name='total_counts')
out = pd.merge(group_counts, total_counts, how="left", on="individual")

out.to_csv(f"{out_dir}/celltype_prop.tsv.gz", sep="\t", index=False)
