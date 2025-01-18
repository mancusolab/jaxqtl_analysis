#!/usr/bin/env python3
import argparse as ap
import logging
import sys
from typing import List

import numpy as np
import pandas as pd
import statsmodels.api as sm

import jax
import jax.numpy as jnp

from jaxtyping import ArrayLike

from jaxqtl.io.geno import PlinkReader
from jaxqtl.io.pheno import PheBedReader
from jaxqtl.io.readfile import create_readydata, ReadyDataState
from jaxqtl.log import get_log, get_logger
from jaxqtl.map.utils import _get_geno_info, _setup_G_y


def cov_scan_sm(G: ArrayLike, y: ArrayLike):
    """
    run GLM across variants in a flanking window of given gene
    cis-widow: plus and minus W base pairs, total length 2*cis_window
    """
    beta_vec = []
    bse_vec = []
    p_vec = []

    # print(sm_res.summary())
    for snp in G.T:
        X = np.hstack((jnp.ones_like(y), snp[:, np.newaxis]))
        sm_glm = sm.Logit(
            np.array(y),
            np.array(X),
        )
        sm_res = sm_glm.fit(maxiter=100)
        beta_vec.append(sm_res.params[-1])
        bse_vec.append(sm_res.bse[-1])
        p_vec.append(sm_res.pvalues[-1])

    return beta_vec, bse_vec, p_vec


def _cis_window_cutter(geno: ArrayLike, bim: pd.DataFrame, chrom: str, start: int, end: int) -> Tuple[Array, pd.DataFrame]:
    """
    return variant list in cis for given gene
    Map is a pandas data frame

    plink bim file is 1-based
    the map file is hg19,
    emsemble use 1-based
    vcf file is one-based

    gene_name = 'ENSG00000250479', start: 24110630
    GenomicRanges example: https://biocpy.github.io/GenomicRanges/

    Returns:
        Genotype matrix for cis-variants
    """
    var_info = bim

    cis_var_info = var_info.loc[
        (var_info["chrom"] == str(chrom)) & (var_info["pos"] >= start) & (var_info["pos"] <= end)
    ]

    # subset G to cis variants (nxp)
    G_tocheck = jnp.take(geno, jnp.array(cis_var_info.i), axis=1)

    # check monomorphic: G.T[:, [0]] find first occurrence on all genotype, var x 1
    mono_var = (G_tocheck.T == G_tocheck.T[:, [0]]).all(1)  # bool (var, ), show whether given variant is monomorphic
    not_mono_var = jnp.invert(mono_var)  # reverse False and True (same as "~" operator)
    G = G_tocheck[:, not_mono_var]  # take genotype that are NOT monomorphic
    cis_var_info = cis_var_info.loc[not_mono_var.tolist()]

    # note: if no variants taken, then G has shape (n,0), cis_var_info has shape (0, 7); both 2-dim
    return G, cis_var_info


def map_annot_sm(
    geno: ReadyDataState,
    annot: pd.DataFrame,

    celltype: str,
    pip_prefix: List,
    pip_suffix: str,
    log=None,
    window: int = 500000,
    verbose: bool = True,
):
    """eQTL Mapping for all cis-SNP gene pairs

    append_intercept: add a column of ones in front for intercepts in design matrix
    standardize: on covariates

    Returns:
        score test statistics and p value (no effect estimates)
        write out parquet file by chrom for efficient data storage and retrieval
    """

    if log is None:
        log = get_log()

    pos_df = pheno[["chr", "start", "end"]]
    gene_info = GeneMetaData(pos_df)

    out_columns = ["phenotype_id", "chrom", "annot", "slope", "slope_se", "pval_nominal", "model_converged"]

    phenotype_id = []
    chrom_list = []
    slope = []
    slope_se = []
    nominal_p = []
    converged = []

    for gene in gene_info:
        gene_name, chrom, start_min, end_max = gene
        lstart = max(0, start_min - window)
        rend = end_max + window

        # read geno on one chromosome and pull cis G (nxM) and y for this gene
        _, var_df = _cis_window_cutter(geno, bim, str(chrom), lstart, rend)

        # TODO: use var_df to subset cis variant in annot file
        cov = jnp.array(annot.iloc[var_df.i.values])

        # skip if no cis SNPs found
        if cov.shape[1] == 0:
            if verbose:
                log.info(
                    "No covariate found for %s over region %s:%s-%s. Skipping.",
                    gene_name,
                    str(chrom),
                    str(lstart),
                    str(rend),
                )
            continue

        if verbose:
            log.info(
                "Performing cis-cov scan (statsmodel) for %s over region %s:%s-%s",
                gene_name,
                str(chrom),
                str(lstart),
                str(rend),
            )

        # read in pip_y
        beta, bse, pval = cov_scan_sm(G, y)

        if verbose:
            log.info(
                "Finished cis-qtl scan (statsmodel) for %s over region %s:%s-%s",
                gene_name,
                str(chrom),
                str(lstart),
                str(rend),
            )

        var_df["phenotype_id"] = gene_name
        var_df["tss"] = start_min
        var_df_all = pd.concat([var_df_all, var_df], ignore_index=True)
        gene_mapped_list.loc[len(gene_mapped_list)] = [gene_name, chrom, start_min]

        # combine results
        af.append(g_info.af)
        ma_count.append(g_info.ma_count)

        nominal_p.append(sm_p)
        Z.append(chi2)
        num_var_cis.append(var_df.shape[0])

        # combine results
        phenotype_id.append(gene_name)
        slope.append(result.beta.item())
        slope_se.append(result.se.item())
        nominal_p.append(result.p.item())
        converged.append(result.converged.item())  # whether full model converged
        alpha.append(result.alpha.item())
        chrom_list.append(chrom)

    # write result
    result_out = [phenotype_id, chrom_list, slope, slope_se, nominal_p, converged, alpha]

    result_df = pd.DataFrame(result_out, index=out_columns).T

    return result_df


def main(args):
    argp = ap.ArgumentParser(description="")  # create an instance
    argp.add_argument("-covar", type=str, help="Covariate path")
    argp.add_argument("-pheno", type=str, help="Pheno path")
    argp.add_argument("-genelist", type=str, help="Path to gene list (no header)")
    argp.add_argument(
        "--platform",
        "-p",
        type=str,
        choices=["cpu", "gpu", "tpu"],
        help="platform, cpu, gpu or tpu",
    )
    argp.add_argument("-window", type=int, default=500000)
    argp.add_argument("-max-iter", type=int, default=1000)
    argp.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Verbose for logger",
    )
    argp.add_argument("-out", type=str, help="out file prefix")

    args = argp.parse_args(args)  # a list a strings

    jax.config.update("jax_enable_x64", True)
    jax.config.update("jax_platform_name", args.platform)

    log = get_logger(__name__, args.out)
    if args.verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)

    # raw genotype data and impute for genotype data
    geno_reader = PlinkReader()
    geno, bim, sample_info = geno_reader(args.geno)

    annot = pd.read_csv(args.annot, header=True, sep="\t")

    pos_df = pd.read_csv(args.genelist, header=None, sep="\t").iloc[:, 0].to_list()

    gene_info = GeneMetaData(pos_df)

    out_df = map_nominal_covar(
        dat,
        family=family,
        test=WaldTest(),
        standardize=args.standardize,
        robust_se=args.robust,
        log=log,
        max_iter=args.max_iter,
        prop_cutoff=args.prop_cutoff,
    )
    out_df.to_csv(args.out + ".below" + str(args.prop_cutoff) + ".cis_wald.tsv.gz", sep="\t", index=False)

    return 0


def run_cli():
    return main(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
