import pandas as pd
import glob
import os

model = "pois"
method = "wald"

# os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha")
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new")

# genelist = pd.read_csv("./all_celltype/allcisgenes.tsv.gz", sep="\t")
# genelist = pd.read_csv("./all_celltype/lm_score_allcisgenes.tsv.gz", sep="\t")
genelist = pd.read_csv("./all_celltype/pois_score_allcisgenes.tsv.gz", sep="\t")

genelist.columns = ["phenotype_id", "celltype", "snp"]

cell_meta = genelist.celltype.unique()

alldat = pd.DataFrame()
for celltype in cell_meta:
  genelist_cell = genelist[genelist['celltype'] == celltype]
  for idx in range(1,23):
      allfiles=glob.glob(f"{celltype}/chr{idx}/chunk_*.{model}.*.{method}.parquet")
      for file in allfiles:
         onedat=pd.read_parquet(file)
         alldat=pd.concat([alldat, onedat.merge(genelist_cell, on=['phenotype_id','snp'])])
         print(file)
  
alldat.to_csv(f"all_celltype/jaxqtl.cis_qtl_pairs.{model}.{method}.egenes.tsv.gz", sep='\t', index=False)


# tensorqtl files
genelist.columns = ["phenotype_id", "celltype", "variant_id"]

alldat = pd.DataFrame()
for celltype in cell_meta:
  genelist_cell = genelist[genelist['celltype'] == celltype]
  for idx in range(1,23):
      allfiles=glob.glob(f"{celltype}/chr{idx}.cis_qtl_pairs.{idx}.parquet")
      for file in allfiles:
         onedat=pd.read_parquet(file)
         alldat=pd.concat([alldat, onedat.merge(genelist_cell, on=['phenotype_id','variant_id'])])
         print(file)

alldat.to_csv(f"all_celltype/tqtl.cis_qtl_pairs.egenes.tsv.gz", sep='\t', index=False)


# poisson
model = "pois"
method = "score"
idx=22
celltype="CD4_NC"

for celltype in ["CD4_NC","B_IN", "Plasma"]:
  alldat=pd.DataFrame()
  allfiles=glob.glob(f"{celltype}/chr{idx}/chunk_*.{model}.cis_qtl_*.{method}.parquet")
  for file in allfiles:
     onedat=pd.read_parquet(file)
     alldat=pd.concat([alldat, onedat])
     print(file)
  alldat.to_csv(f"all_celltype/jaxqtl.cis_qtl_pairs.{celltype}.{model}.{method}.tsv.gz", sep='\t', index=False)

# gather slopes and slope_se from parquet files

model = "nb"
method = "wald"

os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha")

genelist = pd.read_csv("./all_celltype/finemap_eqtl_pip0.5.tsv", sep="\t")
n_rows = genelist.shape[0]

alldat = pd.DataFrame()
for idx in range(n_rows):
  celltype=genelist['celltype'][idx]
  chunk=genelist['chunk'][idx]
  chr=genelist['chr'][idx]
  file=glob.glob(f"{celltype}/chr{chr}/{chunk}.{model}.*.{method}.parquet")
  onedat=pd.read_parquet(file)
  onedat['celltype']=celltype
  alldat=pd.concat([alldat, onedat.merge(genelist, on=['phenotype_id','snp', 'celltype'])])
  print(file)

alldat
alldat.drop_duplicates().to_csv(f"all_celltype/jaxqtl.finemap_eqtl_pip0.5.{model}.{method}.tsv.gz", sep='\t', index=False)

# get slope estimates for each fine-mapped eQTL in each cell type

model = "nb"
method = "wald"

os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha")

genelist = pd.read_csv("./all_celltype/finemap_eqtl_pip0.5.tsv", sep="\t")
gene_chunk = pd.read_csv("./all_celltype/finemap_eqtl_pip0.5_genes_celltype.tsv", sep="\t")

n_rows = gene_chunk.shape[0]
cell_meta = genelist.celltype.unique()

varlist = genelist.drop_duplicates(subset=['phenotype_id', 'snp'])
varlist = varlist[['phenotype_id', 'snp']]

alldat = pd.DataFrame()
for idx in range(n_rows):
  celltype=gene_chunk['celltype'][idx]
  chunk=gene_chunk['chunk'][idx]
  chr=gene_chunk['chr'][idx]
  file=glob.glob(f"{celltype}/chr{chr}/{chunk}.{model}.*.{method}.parquet")
  onedat=pd.read_parquet(file)
  onedat['celltype']=celltype
  alldat=pd.concat([alldat, onedat.merge(varlist, on=['phenotype_id','snp'])])
  print(f"{idx}: {file}")
  
alldat.to_csv(f"./all_celltype/jaxqtl.finemap_eqtl_pip0.5_allest.{model}.{method}.tsv.gz", sep='\t', index=False)

# get results for genes that found only by NB

model = "nb"
method = "wald"

os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha")

genelist = pd.read_csv("./all_celltype/eqtl_nb_only.tsv", sep="\t")
varlist = genelist.drop_duplicates(subset=['phenotype_id'])
varlist = varlist[['phenotype_id']]

gene_dict = pd.read_csv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/genelist/chunk_gene_dict.tsv.gz", 
                        sep="\t")
gene_chunk = gene_dict.merge(varlist, on=['phenotype_id'])

n_rows = gene_chunk.shape[0]

alldat = pd.DataFrame()
for idx in range(n_rows):
  celltype=gene_chunk['celltype'][idx]
  chunk=gene_chunk['chunk'][idx]
  chr=gene_chunk['chr'][idx]
  file=glob.glob(f"{celltype}/chr{chr}/{chunk}.{model}.*.{method}.parquet")
  onedat=pd.read_parquet(file)
  onedat['celltype']=celltype
  alldat=pd.concat([alldat, onedat.merge(varlist, on=['phenotype_id'])])
  print(f"{idx}: {file}")
  
alldat.to_csv(f"./all_celltype/eqtl_nb_only.all_pairs.tsv.gz", sep='\t', index=False)



# grab poisson results
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new")

# genelist = pd.read_csv("./all_celltype/allcisgenes.tsv.gz", sep="\t")
# genelist = pd.read_csv("./all_celltype/lm_score_allcisgenes.tsv.gz", sep="\t")
genelist = pd.read_csv("./all_celltype/pois_score_allcisgenes.tsv.gz", sep="\t")

genelist.columns = ["phenotype_id", "celltype", "snp"]

model = "pois"
method = "wald"

gene_dict = pd.read_csv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/genelist/chunk_gene_dict.tsv.gz", 
                        sep="\t")
varlist = gene_dict.merge(genelist, on=['phenotype_id', 'celltype'])
gene_chunk = varlist.drop_duplicates(subset=['chunk', 'celltype', 'chr'])
gene_chunk = gene_chunk.reset_index(drop=True)

n_rows = gene_chunk.shape[0]

alldat = pd.DataFrame()
for idx in range(n_rows):
  celltype=gene_chunk['celltype'][idx]
  chunk=gene_chunk['chunk'][idx]
  chr=gene_chunk['chr'][idx]
  file=glob.glob(f"{celltype}/chr{chr}/{chunk}.*.{model}.*.{method}.parquet")
  onedat=pd.read_parquet(file)
  onedat['celltype']=celltype
  alldat=pd.concat([alldat, onedat.merge(varlist, on=['phenotype_id', 'celltype', 'snp'])])
  print(f"{idx}: {file}")
  
alldat.to_csv(f"all_celltype/jaxqtl.cis_qtl_pairs.{model}.{method}.egenes.tsv.gz", sep='\t', index=False)


# grab nb bulk (allcells) results
os.chdir("/project/nmancuso_8/elezhang/projects/jaxqtl/result/cis/celltype16_new_fixalpha")

# genelist = pd.read_csv("./all_celltype/allcisgenes.tsv.gz", sep="\t")
# genelist = pd.read_csv("./all_celltype/lm_score_allcisgenes.tsv.gz", sep="\t")
genelist = pd.read_csv("./all_celltype/nb_allcells_score_allcisgenes.tsv.gz", sep="\t")

genelist.columns = ["phenotype_id", "celltype", "snp"]

model = "nb"
method = "wald"

gene_dict = pd.read_csv("/project/nmancuso_8/elezhang/projects/jaxqtl/data/pheno/celltype16_new/metadata/genelist/chunk_gene_dict.tsv.gz", 
                        sep="\t")
varlist = gene_dict.merge(genelist, on=['phenotype_id', 'celltype'])
gene_chunk = varlist.drop_duplicates(subset=['chunk', 'celltype', 'chr'])
gene_chunk = gene_chunk.reset_index(drop=True)

n_rows = gene_chunk.shape[0]

alldat = pd.DataFrame()
for idx in range(n_rows):
  celltype=gene_chunk['celltype'][idx]
  chunk=gene_chunk['chunk'][idx]
  chr=gene_chunk['chr'][idx]
  file=glob.glob(f"{celltype}/chr{chr}/{chunk}.{model}.*.{method}.parquet")
  onedat=pd.read_parquet(file)
  onedat['celltype']=celltype
  alldat=pd.concat([alldat, onedat.merge(varlist, on=['phenotype_id', 'celltype', 'snp'])])
  print(f"{idx}: {file}")
  
alldat.to_csv(f"all_celltype/jaxqtl.allcells.cis_qtl_pairs.{model}.{method}.egenes.tsv.gz", sep='\t', index=False)
