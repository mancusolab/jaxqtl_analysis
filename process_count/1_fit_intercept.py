import os
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/code")

import pandas as pd

import jax.numpy as jnp

from jaxqtl.families.distribution import Poisson, Gaussian
from jaxqtl.io.pheno import PheBedReader
from jaxqtl.infer.glm import GLM

from jax.config import config

platform = "cpu"
config.update("jax_enable_x64", True)
config.update("jax_platform_name", platform)

pheno_dir = "../data/pheno/bed/"
celltype_list_path = "../data/pheno_meta/celltype.tsv"

# get a list of files
celltype_list = pd.read_csv(celltype_list_path, header=None, sep="\t").iloc[:, 0].to_list()

# read raw count
pheno_reader = PheBedReader()

celltype = "natural_killer_cell"

pheno_path = pheno_dir + celltype + ".bed.gz"
pheno = pheno_reader(pheno_path)

pheno.drop(["chr", "start", "end"], axis=1, inplace=True)
mean_df = pd.DataFrame({'phenotype_id': pheno.index})

# before filter gene list, calculate library size and set offset
pheno = jnp.array(pheno.T)  # transpose to sample x genes
total_libsize = pheno.sum(axis=1)[:, jnp.newaxis]
offset_eta = jnp.log(total_libsize)

# intercept only
M = jnp.ones((len(pheno), 1))

glm = GLM(
    family=Poisson(),
    maxiter=100,
)

mean_df['exp_intercept'] = jnp.zeros((len(pheno.T),))

mu_df = pd.DataFrame(jnp.zeros(pheno.shape))

for i in range(0, len(pheno.T)):
    if pheno[:, i].sum() == 0:
        continue
    glmstate = glm.fit(M, pheno[:, i][:, jnp.newaxis], offset_eta=offset_eta)
    mean_df.loc[i, 'exp_intercept'] = jnp.exp(glmstate.beta[0])
    mu_df.loc[:, i] = glmstate.mu
    if i % 100 == 0:
        print(f"finish total {i}")

mean_df.to_csv(f"../result/mean_expr/{celltype}.exp_intercept.allgenes.gz", index=False, header=True, sep='\t')
mu_df.to_csv(f"../result/mean_expr/{celltype}.mu.allgenes.gz", index=False, header=True, sep='\t')