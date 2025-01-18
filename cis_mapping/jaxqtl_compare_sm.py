import os
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/code")

import pandas as pd
import statsmodels.api as sm
import numpy as np

from jaxqtl.families.distribution import Poisson, Gaussian
from jaxqtl.infer.solve import QRSolve, CholeskySolve
from jaxqtl.io.covar import covar_reader
from jaxqtl.io.geno import PlinkReader
from jaxqtl.io.pheno import PheBedReader
from jaxqtl.io.readfile import create_readydata

from jaxqtl.log import get_log
from jaxqtl.map import map_cis, map_cis_nominal
from jaxqtl.infer.glm import GLM

from jax.config import config
import jax.numpy as jnp

platform = "cpu"
config.update("jax_enable_x64", True)
config.update("jax_platform_name", platform)

geno_path = "../data/geno_n981/chr"  # prefix
covar_path = "../data/features/donor_features.all.6PC.tsv"
pheno_dir = "../data/pheno/bed/NK_meansum/"
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
# genelist = pd.read_csv(f"../data/pheno_meta/NK_chr{chr}_genelist.nomissing.tsv", header=None).iloc[:, 0].to_list()
# genelist = ['ENSG00000178607']

# for chr in [str(i) for i in range(1, 23)]:
# TODO: we only need read pheno and covar data once
dat = create_readydata(
    geno,
    bim,
    pheno,
    covar,
    autosomal_only=True
)

library_size = dat.pheno.count.sum(axis=1)
library_size = jnp.log(jnp.array(library_size))[:, jnp.newaxis]

# maf_threshold = 0.0
# dat.filter_geno(maf_threshold, chr)

if isinstance(family, Gaussian):
    dat.transform_y(y0=1.0, log_y=True)

# add expression PCs to covar, genotype PC should appended to covar outside jaxqtl
dat.add_covar_pheno_PC(k=2)

# ERN1
M = jnp.hstack((dat.covar, dat.geno[:, 213943][:, jnp.newaxis]))
M = jnp.hstack((jnp.ones((M.shape[0], 1)), M))
y = jnp.array(dat.pheno.count['ENSG00000178607'])[:, jnp.newaxis]

M = jnp.hstack((dat.covar, dat.geno[:, 33653][:, jnp.newaxis]))
M = jnp.hstack((jnp.ones((M.shape[0], 1)), M))
y = jnp.array(dat.pheno.count['ENSG00000132507'])[:, jnp.newaxis]


# full model fit to compare Wald p to Score p
glmstate = GLM(X=M,y=y,family=family,append=False,maxiter=100, solver=CholeskySolve()).fit(offset_eta=library_size)

sm_mod = sm.GLM(np.array(y), np.array(M), family=sm.families.Poisson(), offset=np.array(library_size).reshape((len(library_size),))).fit()
