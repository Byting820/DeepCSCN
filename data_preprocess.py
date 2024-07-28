'''
Description:
    Pre-process experimental data, including removing non-expressed cells/genes and find highly variable genes (HVGs).
'''

import scanpy
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix

# ===============================
#          LOAD RAW DATA
# ===============================
def _formatStatsDict(stats_dict):
    return pd.DataFrame(data=stats_dict, index=["# Cells", "# Genes", "Sparsity", "Max Cnt"]).T

def loadMouseCortex():
    '''
    Load mouse cortex raw data.
    :return:
    (dict) A dictionary of experimental data. Keys are data names (protocols) and values are
           cell-by-gene count expression matrices.
    (pd.DataFrame) A table of dataset statistics.
    '''
    # Load cell names and gene names
    cell_names = pd.read_csv("data/mouse_cortex/cell.names.new.txt", header=None).values.squeeze()
    gene_names = pd.read_csv("data/mouse_cortex/genes.count.txt", header=None).values.squeeze()
    num_genes = len(gene_names)
    num_cells = len(cell_names)
    # Load UMI
    expr_mat = pd.read_csv("data/mouse_cortex/count.umis.txt.gz", compression='gzip', sep=' ', skiprows=1).values
    expr_mat = coo_matrix((expr_mat[:, 2], (expr_mat[:, 0]-1, expr_mat[:, 1]-1)), shape=(num_genes, num_cells)).toarray()
    # Split by experiment and protocols
    exp_types = ["Cortex1", "Cortex2"]
    protocol_types = ["Smart_seq2", "10xChromium"]
    expr_mat_dict = {}    
    expr_stats_dict = {}
    for e_t in exp_types:
        for p_t in protocol_types:
            tmp_idx = np.where([True if (e_t in each) and (p_t in each) else False for each in cell_names])[0]
            tmp_mat = pd.DataFrame(data=expr_mat[:, tmp_idx].T, index=cell_names[tmp_idx], columns=gene_names)
            expr_mat_dict["{}-{}".format(e_t, p_t)] = tmp_mat
            expr_stats_dict["{}-{}".format(e_t, p_t)] = [
                tmp_mat.shape[0],
                tmp_mat.shape[1],
                np.count_nonzero(tmp_mat.values)/np.prod(tmp_mat.shape),
                np.max(tmp_mat.values)
            ]
    return expr_mat_dict, _formatStatsDict(expr_stats_dict)


def loadPBMC():
    '''
    Load human PBMC raw data.
    :return:
    (dict) A dictionary of experimental data. Keys are data names (protocols) and values are
            cell-by-gene count expression matrices.
    (pd.DataFrame) A table of dataset statistics.
    '''
    # Load cell names and gene names
    cell_names = pd.read_csv("data/PBMC/cells.umi.new.txt", header=None).values.squeeze()
    gene_names = pd.read_csv("data/PBMC/genes.umi.txt", header=None).values.squeeze()
    num_genes = len(gene_names)
    num_cells = len(cell_names)
    # Load UMI
    expr_mat = pd.read_csv("data/PBMC/counts.umi.txt.gz", compression='gzip', sep=' ', skiprows=1).values
    expr_mat = coo_matrix((expr_mat[:, 2], (expr_mat[:, 0] - 1, expr_mat[:, 1] - 1)), shape=(num_genes, num_cells)).toarray()
    # Split by experiment and protocols
    exp_types = ["pbmc1", "pbmc2"]
    protocol_types = ["Drop", "inDrops"]
    expr_mat_dict = {}
    expr_stats_dict = {}
    for e_t in exp_types:
        for p_t in protocol_types:
            tmp_idx = np.where([True if (e_t in each) and (p_t in each) else False for each in cell_names])[0]
            if len(tmp_idx) == 0:
                print("No samples for {}-{}".format(e_t, p_t))
                continue
            tmp_mat = pd.DataFrame(data=expr_mat[:, tmp_idx].T, index=cell_names[tmp_idx], columns=gene_names)
            expr_mat_dict["{}-{}".format(e_t, p_t)] = tmp_mat
            expr_stats_dict["{}-{}".format(e_t, p_t)] = [
                tmp_mat.shape[0],
                tmp_mat.shape[1],
                np.count_nonzero(tmp_mat.values) / np.prod(tmp_mat.shape),
                np.max(tmp_mat.values)
            ]
    return expr_mat_dict, _formatStatsDict(expr_stats_dict)


# ===============================
#          PRE-PROCESSING
# ===============================
def qcFilter(expr):
    '''
    Remove non-expressed cells and genes with quality control.
    :param expr: (pd.DataFrame): Cell-by-gene expression matrix.
    :return: (pd.DataFrame): Cell-by-gene expression matrix after quality control.
    '''
    adata = scanpy.AnnData(expr)
    scanpy.pp.filter_cells(adata, min_counts=1)
    scanpy.pp.filter_genes(adata, min_counts=1)
    return adata.to_df()

def findHVG(expr, num_HVG):
    '''
    Find highly variable genes (HVGs).
    :param expr: (pd.DataFrame): Cell-by-gene expression matrix.
    :param num_HVG: (int) The number of highly variable genes.
    :return: (pd.DataFrame): Cell-by-HVG expression matrix.
    '''
    adata = scanpy.AnnData(expr.copy())
    scanpy.pp.log1p(adata)
    hvg_idx = scanpy.pp.highly_variable_genes(adata, n_top_genes=num_HVG, inplace=False)
    if len(np.where(hvg_idx.highly_variable.values)[0]) > num_HVG:
        idx_by_mean = np.argsort(hvg_idx[hvg_idx.highly_variable == True].means.values)[::-1]
        hvg_idx = idx_by_mean[:num_HVG]
    else:
        hvg_idx = np.where(hvg_idx.highly_variable.values)[0]
    return expr.iloc[:, hvg_idx]


# ===============================
#         MAIN FUNCTION
# ===============================
def preprocessMouseCortex():
    '''
    Pre-process mouse cortex data.
    '''
    print("=" * 70)
    print("Loading data for mouse cortex...")
    expr_mat_dict, expr_stats_dict = loadMouseCortex()
    print(expr_stats_dict)
    # -----
    print("=" * 70)
    print("Quality control...")
    qc_expr_mat_dict = {each: qcFilter(expr_mat_dict[each]) for each in expr_mat_dict} 
    qc_expr_stats_dict = {each: [
        qc_expr_mat_dict[each].shape[0],
        qc_expr_mat_dict[each].shape[1],
        np.count_nonzero(qc_expr_mat_dict[each].values) / np.prod(qc_expr_mat_dict[each].shape),
        np.max(qc_expr_mat_dict[each].values)
    ] for each in qc_expr_mat_dict}
    print(_formatStatsDict(qc_expr_stats_dict))
    # -----
    print("=" * 70)
    print("Find highly variable genes...")
    hvg100_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=100) for each in qc_expr_mat_dict}
    hvg500_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=500) for each in qc_expr_mat_dict}
    hvg1000_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=1000) for each in qc_expr_mat_dict}
    print("[ 100 HVGs ]")
    hvg100_expr_stats_dict = {each: [
        hvg100_expr_mat_dict[each].shape[0],
        hvg100_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg100_expr_mat_dict[each].values) / np.prod(hvg100_expr_mat_dict[each].shape),
        np.max(hvg100_expr_mat_dict[each].values)
    ] for each in hvg100_expr_mat_dict}
    print(_formatStatsDict(hvg100_expr_stats_dict))
    print("[ 500 HVGs ]")
    hvg500_expr_stats_dict = {each: [
        hvg500_expr_mat_dict[each].shape[0],
        hvg500_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg500_expr_mat_dict[each].values) / np.prod(hvg500_expr_mat_dict[each].shape),
        np.max(hvg500_expr_mat_dict[each].values)
    ] for each in hvg500_expr_mat_dict}
    print(_formatStatsDict(hvg500_expr_stats_dict))
    print("[ 1000 HVGs ]")
    hvg1000_expr_stats_dict = {each: [
        hvg1000_expr_mat_dict[each].shape[0],
        hvg1000_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg1000_expr_mat_dict[each].values) / np.prod(hvg1000_expr_mat_dict[each].shape),
        np.max(hvg1000_expr_mat_dict[each].values)
    ] for each in hvg1000_expr_mat_dict}
    print(_formatStatsDict(hvg1000_expr_stats_dict))
    # -----
    print("=" * 70)
    print("Saving data...")
    save_path = "data/mouse_cortex/processed/expr"
    for each in qc_expr_mat_dict:
        hvg100_expr_mat_dict[each].to_csv("{}/{}-100hvg.csv".format(save_path, each))
        hvg500_expr_mat_dict[each].to_csv("{}/{}-500hvg.csv".format(save_path, each))
        hvg1000_expr_mat_dict[each].to_csv("{}/{}-1000hvg.csv".format(save_path, each))


def preprocessPBMC():
    '''
    Pre-process human PBMC data.
    '''
    print("=" * 70)
    print("Loading data for PBMC...")
    expr_mat_dict, expr_stats_dict = loadPBMC()
    print(expr_stats_dict)
    # -----
    print("=" * 70)
    print("Quality control...")
    qc_expr_mat_dict = {each: qcFilter(expr_mat_dict[each]) for each in expr_mat_dict}
    qc_expr_stats_dict = {each: [
        qc_expr_mat_dict[each].shape[0],
        qc_expr_mat_dict[each].shape[1],
        np.count_nonzero(qc_expr_mat_dict[each].values) / np.prod(qc_expr_mat_dict[each].shape),
        np.max(qc_expr_mat_dict[each].values)
    ] for each in qc_expr_mat_dict}
    print(_formatStatsDict(qc_expr_stats_dict))
    # for each in qc_expr_mat_dict:
    #     qc_expr_mat_dict[each].to_csv("{}/{}.csv".format(save_path, each))

    # -----
    print("=" * 70)
    print("Find highly variable genes...")
    hvg100_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=100) for each in qc_expr_mat_dict}
    hvg500_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=500) for each in qc_expr_mat_dict}
    hvg1000_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=1000) for each in qc_expr_mat_dict}
    # hvg2000_expr_mat_dict = {each: findHVG(qc_expr_mat_dict[each], num_HVG=2000) for each in qc_expr_mat_dict}
    print("[ 100 HVGs ]")
    hvg100_expr_stats_dict = {each: [
        hvg100_expr_mat_dict[each].shape[0],
        hvg100_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg100_expr_mat_dict[each].values) / np.prod(hvg100_expr_mat_dict[each].shape),
        np.max(hvg100_expr_mat_dict[each].values)
    ] for each in hvg100_expr_mat_dict}
    print(_formatStatsDict(hvg100_expr_stats_dict))
    print("[ 500 HVGs ]")
    hvg500_expr_stats_dict = {each: [
        hvg500_expr_mat_dict[each].shape[0],
        hvg500_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg500_expr_mat_dict[each].values) / np.prod(hvg500_expr_mat_dict[each].shape),
        np.max(hvg500_expr_mat_dict[each].values)
    ] for each in hvg500_expr_mat_dict}
    print(_formatStatsDict(hvg500_expr_stats_dict))
    print("[ 1000 HVGs ]")
    hvg1000_expr_stats_dict = {each: [
        hvg1000_expr_mat_dict[each].shape[0],
        hvg1000_expr_mat_dict[each].shape[1],
        np.count_nonzero(hvg1000_expr_mat_dict[each].values) / np.prod(hvg1000_expr_mat_dict[each].shape),
        np.max(hvg1000_expr_mat_dict[each].values)
    ] for each in hvg1000_expr_mat_dict}
    print(_formatStatsDict(hvg1000_expr_stats_dict))

    # print("[ 2000 HVGs ]")
    # hvg2000_expr_stats_dict = {each: [
    #     hvg2000_expr_mat_dict[each].shape[0],
    #     hvg2000_expr_mat_dict[each].shape[1],
    #     np.count_nonzero(hvg2000_expr_mat_dict[each].values) / np.prod(hvg2000_expr_mat_dict[each].shape),
    #     np.max(hvg2000_expr_mat_dict[each].values)
    # ] for each in hvg2000_expr_mat_dict}
    # print(_formatStatsDict(hvg2000_expr_stats_dict))

    # -----
    print("=" * 70)
    print("Saving data...")
    save_path = "data/PBMC/processed/expr/"
    for each in qc_expr_mat_dict:
        hvg100_expr_mat_dict[each].to_csv("{}/{}-100hvg.csv".format(save_path, each))
        hvg500_expr_mat_dict[each].to_csv("{}/{}-500hvg.csv".format(save_path, each))
        hvg1000_expr_mat_dict[each].to_csv("{}/{}-1000hvg.csv".format(save_path, each))
        # hvg2000_expr_mat_dict[each].to_csv("{}/{}-2000hvg.csv".format(save_path, each))


# ===============================

if __name__ == '__main__':
    # preprocessMouseCortex()
    preprocessPBMC()
