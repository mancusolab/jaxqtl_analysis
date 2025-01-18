import os
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/code")

import pandas as pd
import jax.numpy as jnp

from jaxqtl.families.distribution import Poisson, Gaussian
from jaxqtl.io.covar import covar_reader
from jaxqtl.io.geno import PlinkReader
from jaxqtl.io.pheno import PheBedReader
from jaxqtl.io.readfile import create_readydata

from jaxqtl.log import get_log
from jaxqtl.map import map_cis, map_cis_nominal

from jax.config import config
from numpy import random
import numpy as np

platform = "cpu"
config.update("jax_enable_x64", True)
config.update("jax_platform_name", platform)

geno_path = "../data/geno_n981/chr"  # prefix
covar_path = "../data/features/donor_features.all.6PC.tsv"
pheno_dir = "../data/pheno/bed/"

family = Poisson()

log = get_log()

# raw genotype data and impute for genotype data
geno_reader = PlinkReader()
pheno_reader = PheBedReader()

covar = covar_reader(covar_path)

celltype = "natural_killer_cell"
chr='17'

# read pheno, geno
pheno_path = pheno_dir + celltype + ".bed.gz"
pheno = pheno_reader(pheno_path)

geno_path = geno_path + chr
geno, bim, sample_info = geno_reader(geno_path)

# read gene list
genelist = pd.read_csv(f"../data/pheno_meta/NK_chr{chr}_genelist.tsv", header=None).iloc[:, 0].to_list()
# genelist = ['ENSG00000178607']  # ENR1

# for chr in [str(i) for i in range(1, 23)]:
# TODO: we only need read pheno and covar data once
dat = create_readydata(
    geno,
    bim,
    pheno,
    covar,
    autosomal_only=True
)

# remove variants with no variation
maf_threshold = 0.0
dat.filter_geno(maf_threshold)

# filter genes with no expressions at all
dat.filter_gene(geneexpr_percent_cutoff=0.0)

# before filter gene list, calculate library size and set offset
total_libsize = jnp.array(dat.pheno.count.sum(axis=1))[:, jnp.newaxis]
offset_eta = jnp.log(total_libsize)

# add expression PCs to covar, genotype PC should appended to covar outside jaxqtl
dat.add_covar_pheno_PC(k=2)

dat.filter_gene(gene_list=genelist)

if isinstance(family, Gaussian):
    dat.transform_y(y0=1.0, log_y=True)

# permute outcome values only and do cis scan
number_perm = 1
true_df = dat.pheno.count.copy()
random.seed(123)

genelist = dat.pheno_meta.gene_map.phenotype_id.to_list()

for idx in range(0, number_perm):
    for gene in genelist:
        dat.pheno.count[gene] = dat.pheno.count[gene].iloc[random.permutation(len(true_df))].values

    out_prefix = f"../result/cis/perm/{celltype}.chr{chr}.pois.offset.perm.nominal"

    map_cis_nominal(dat, family=family, standardize=False, log=log,
                    window=500000, out_path=out_prefix,
                    offset_eta=offset_eta, robust_se=True)
