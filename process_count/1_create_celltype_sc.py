# need --mem-per-cpu 60Gb for this script
# use this script to create single cell data from OneK1K for GLMM

import argparse as ap
import os
import sys
import logging

from jaxqtl.io.pheno import SingleCellFilter, H5AD
import numpy as np
import pandas as pd
import scipy as scp
import scanpy as sc


def get_logger(name, path=None):
    """get logger for factorgo progress"""
    logger = logging.getLogger(name)
    if not logger.handlers:
        # Prevent logging from propagating to the root logger
        logger.propagate = 0
        console = logging.StreamHandler()
        logger.addHandler(console)

        log_format = "[%(asctime)s - %(levelname)s] %(message)s"
        date_format = "%Y-%m-%d %H:%M:%S"
        formatter = logging.Formatter(fmt=log_format, datefmt=date_format)
        console.setFormatter(formatter)

        if path is not None:
            disk_log_stream = open("{}.log".format(path), "w")
            disk_handler = logging.StreamHandler(disk_log_stream)
            logger.addHandler(disk_handler)
            disk_handler.setFormatter(formatter)

    return logger


def _main(args):
    argp = ap.ArgumentParser(description="")  # create an instance
    argp.add_argument(
        "--celltype",
        type=str,
        help="select which cell type, eg. NK",
    )
    argp.add_argument(
        "--outpath",
        type=str,
        help="output directory",
    )
    argp.add_argument(
        "-k",
        type=int,
        help="Number of expression PCs to compute",
    )
    argp.add_argument(
        "--genelist",
        type=str,
        help="path to gene list (tsv)",
    )
    
    args = argp.parse_args(args)

    os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/data/")

    raw_count_path = "./pheno/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz"  # n=982

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
    
    dat.obs['cell_label'] = dat.obs['cell_label'].str.replace(' ', '_')
    # subset data by cell type and gene list
    celltype = args.celltype # NK
    outpath = args.outpath # "../saigeqtl/input/pheno"
    
    dat = dat[dat.obs.cell_label == celltype]
    # compute PC across all genes
    sc.tl.pca(dat, n_comps=args.k, zero_center=True)  # PCA on covariance matrix (no standardization)
    np.savetxt(f"{outpath}/{celltype}.sc.2PC.gz", dat.obsm['X_pca'], delimiter="\t")

    # calculate individual-cell offset
    dat.obs['ind_cell_offset'] = dat.X.sum(axis=1)
    
    # pull up
    genelist = pd.read_csv(args.genelist,sep="\t", header=None, names=["phenotype_id", "chr", "tss"])
    dat = dat[:, genelist.phenotype_id.values] # change col names here
    
    # remove genes not expressed at all
    # dat = dat[dat.X.sum(axis=1) > 0]

    scp.io.mmwrite(f"{outpath}/{celltype}.sc.mtx", dat.X)
    dat.var.reset_index().to_csv(f"{outpath}/{celltype}.sc.genelist.gz", sep="\t", index=False)
    dat.obs.reset_index().to_csv(f"{outpath}/{celltype}.sc.cellmeta.tsv.gz", sep="\t", index=False)

    # dat.X = dat.X.todense() # convert to dense matrix; caution: takes too much memory


def run_cli():
    return _main(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
