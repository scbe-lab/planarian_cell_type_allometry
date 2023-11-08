import scanpy as sc
import numpy as np
import anndata
import gzip
from scipy.io import mmwrite
adata = sc.read_h5ad('/mnt/sda/alberto/colabos/sizes/data/smed_size_analysis_202306.h5ad')
adata.obs.to_csv('/mnt/sda/alberto/colabos/sizes/data/sizes_Idents.csv')
adata.var.to_csv('/mnt/sda/alberto/colabos/sizes/data/sizes_genes.csv') # var
adata.raw.var.to_csv('/mnt/sda/alberto/colabos/sizes/data/sizes_genes_all.csv') # var, all
