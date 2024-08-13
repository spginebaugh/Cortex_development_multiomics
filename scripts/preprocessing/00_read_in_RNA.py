# Import Libraries
import os

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy

import pybiomart
import scipy.sparse

# set file paths
meta_path = "data/raw_data/GSE162170_rna_cell_metadata.txt.gz"
counts_path = "data/raw_data/GSE162170_rna_counts.tsv.gz"
write_h5ad_path = "data/processed_data/filtered_rna.h5ad"

# Import metadata
meta = pd.read_csv(meta_path, compression='gzip', sep='\t')
meta.index = meta['Cell.ID']
meta.index.name = None

# Import counts
adata = sc.read_csv(counts_path, delimiter='\t', dtype="int32")
adata = adata.transpose()

# Add metadata to anndata object
adata.obs = meta

# Change gene names from ensembl id to gene symbols
bmt_queries = sc.queries.biomart_annotations(
    "hsapiens", 
    ["ensembl_gene_id", "start_position", "end_position", "chromosome_name","gene_biotype","hgnc_symbol"])

bmt_queries = bmt_queries.drop_duplicates(subset="ensembl_gene_id")
adata.var["ensembl_gene_id"] = adata.var_names

var_mat = adata.var
var_mat = var_mat.merge(bmt_queries, how = "left", on = "ensembl_gene_id")
var_mat["gene_name"] = var_mat["hgnc_symbol"]
var_mat.loc[var_mat["hgnc_symbol"].isna(), "gene_name"] = var_mat.loc[var_mat["hgnc_symbol"].isna(), "ensembl_gene_id"]

var_mat.index = var_mat["gene_name"]
var_mat.index.name = None

# Add gene info to anndata object
adata.var = var_mat
adata.var_names_make_unique()

# Add QC and normalize counts
adata.X = scipy.sparse.csr_matrix(adata.X.copy(), dtype=np.int32)
adata.layers["counts"] = adata.X.copy()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# filtering
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# normalize
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# save to h5ad
adata.write_h5ad(write_h5ad_path)
