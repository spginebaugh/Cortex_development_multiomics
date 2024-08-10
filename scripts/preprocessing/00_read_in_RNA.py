# %%
import os

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy

import pybiomart
import scipy.sparse
# %%
meta = pd.read_csv("data/raw_data/GSE162170_rna_cell_metadata.txt.gz", compression='gzip', sep='\t')
meta.index = meta['Cell.ID']
meta.index.name = None

# %%
adata = sc.read_csv("data/raw_data/GSE162170_rna_counts.tsv.gz", delimiter='\t', dtype="int32")

# %%
adata = adata.transpose()

# %% 
sum(adata.obs_names==meta.index)
adata.obs = meta
# %%
bmt_queries = sc.queries.biomart_annotations(
    "hsapiens", 
    ["ensembl_gene_id", "start_position", "end_position", "chromosome_name","gene_biotype","hgnc_symbol"])

# %%
bmt_queries = bmt_queries.drop_duplicates(subset="ensembl_gene_id")

# %%
adata.var["ensembl_gene_id"] = adata.var_names

# %%
var_mat = adata.var

# %%
var_mat = var_mat.merge(bmt_queries, how = "left", on = "ensembl_gene_id")
# %%
var_mat["gene_name"] = var_mat["hgnc_symbol"]
var_mat.loc[var_mat["hgnc_symbol"].isna(), "gene_name"] = var_mat.loc[var_mat["hgnc_symbol"].isna(), "ensembl_gene_id"]

# %%
var_mat.index = var_mat["gene_name"]
var_mat.index.name = None

adata.var = var_mat
# %%
adata.var_names_make_unique()

# %%
adata.X = scipy.sparse.csr_matrix(adata.X.copy(), dtype=np.int32)
adata.layers["counts"] = adata.X.copy()

# %%
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# %%
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# %%
sc.pl.violin(adata, ['n_genes_by_counts', "total_counts", "pct_counts_mt"],
             jitter=0.4, groupby = 'Sample.ID')
# %%
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# %% filtering
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# %%
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

# %%
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="Sample.ID")

# %%
sc.tl.pca(adata)
# %%
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# %%
sc.pl.umap(adata, color = ["Tissue.ID","Assay","seurat_clusters"])