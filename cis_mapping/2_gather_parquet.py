# gather parquet results and pull cis genes results

import pandas as pd
import glob
import os

cell_meta = pd.read_csv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno_meta/celltype_14.tsv",
                        header=None)
cell_meta = cell_meta.iloc[:,0].values
cell_meta = [x.replace(" ", "_") for x in cell_meta]

indir = "/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis"
gtex_1per = "1per"
methods = ['celltype16', 'celltype16_tensorqtl']

for celltype in cell_meta:
  allgenes = []
  jqtl_genes = pd.read_csv(f"{indir}/celltype16/{celltype}/jaxqtl_allres_cis_score.{celltype}.cisgenes.tsv", sep='\t', header=None,
  names=['gene']).gene.values
  tqtl_genes = pd.read_csv(f"{indir}/celltype16_tensorqtl/{celltype}/tqtl_allres.{celltype}.{gtex_1per}.cisgenes.tsv", sep='\t', header=None,
  names=['gene']).gene.values
  allgenes = list(set(jqtl_genes) | set(tqtl_genes))
  
  jqtl_df = pd.DataFrame()
  tqtl_df = pd.DataFrame()
  for chr in range(1,23):
    jqtl_pheno_path = f"{indir}/celltype16/{celltype}/chr{chr}"
    os.chdir(jqtl_pheno_path)
    jqtl_allfiles = glob.glob(f"*.wald.parquet")
    print(chr)
    for file in jqtl_allfiles:
      df = pd.read_parquet(file)
      jqtl_df = pd.concat([jqtl_df, df.loc[df['phenotype_id'].isin(allgenes)]])
    
    tqtl_pheno_path = f"{indir}/celltype16_tensorqtl/{celltype}"
    os.chdir(tqtl_pheno_path)
    tqtl_allfiles = glob.glob(f"*.parquet")
    print(chr)
    for file in tqtl_allfiles:
      df = pd.read_parquet(file)
      tqtl_df = pd.concat([tqtl_df, df.loc[df['phenotype_id'].isin(allgenes)]])
  
  jqtl_df.to_csv(f"{indir}/celltype16/{celltype}/jaxqtl_bothcisgenes_cis_score.slope.tsv.gz", sep='\t', index=False)
  tqtl_df.to_csv(f"{indir}/celltype16_tensorqtl/{celltype}/tqtl_bothcisgenes_cis_score.slope.tsv.gz", sep='\t', index=False)
