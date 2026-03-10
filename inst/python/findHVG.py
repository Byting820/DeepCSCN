import scanpy
import numpy as np
import pandas as pd

def findHVG(expr, num_HVG):
    """
    Find highly variable genes (HVGs).
    :param expr: (pd.DataFrame) Cell-by-gene expression matrix.
    :param num_HVG: (int) The number of highly variable genes.
    :return: (pd.DataFrame) Cell-by-HVG expression matrix.
    """
    adata = scanpy.AnnData(expr.copy())
    scanpy.pp.log1p(adata)
    hvg_idx = scanpy.pp.highly_variable_genes(adata, n_top_genes=num_HVG, inplace=False)
    if len(np.where(hvg_idx.highly_variable.values)[0]) > num_HVG:
        idx_by_mean = np.argsort(hvg_idx[hvg_idx.highly_variable == True].means.values)[::-1]
        hvg_idx = idx_by_mean[:num_HVG]
    else:
        hvg_idx = np.where(hvg_idx.highly_variable.values)[0]
    return expr.iloc[:, hvg_idx]
