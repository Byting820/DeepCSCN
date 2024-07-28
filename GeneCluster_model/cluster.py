import os
import pickle
import time
import faiss
import numpy as np

import torch
from torch.utils.data import Dataset


__all__ = ['Kmeans']

def preprocess_features(npdata, pca=20):
    """ PCA
    Preprocess an array of features.
    Args:
        npdata (np.array N * ndim): features to preprocess
        pca (int): dim of output
    Returns:
        np.array of dim N * pca: data PCA-reduced, whitened and L2-normalized
    """
    _, ndim = npdata.shape
    npdata =  npdata.astype('float32')

    # Apply PCA-whitening with Faiss
    mat = faiss.PCAMatrix(ndim, pca, eigen_power=-0.5)  
    mat.train(npdata)  
    assert mat.is_trained
    npdata = mat.apply_py(npdata)

    # L2 normalization
    row_sums = np.linalg.norm(npdata, axis=1)  
    npdata = npdata / row_sums[:, np.newaxis]

    return npdata   


def run_kmeans(x, nmb_clusters, verbose=False):
    """Runs kmeans on 1 GPU.
    Args:
        x: data
        nmb_clusters (int): number of clusters
    Returns:
        list: ids of data in each cluster
    """
    n_data, d = x.shape

    # faiss implementation of k-means
    clus = faiss.Clustering(d, nmb_clusters)

    # Change faiss seed at each k-means so that the randomly picked
    # initialization centroids do not correspond to the same feature ids
    # from an epoch to another.
    # clus.seed = np.random.randint(1234)    

    clus.niter = 20   

    clus.max_points_per_centroid = 1000   
     
    res = faiss.StandardGpuResources()
    flat_config = faiss.GpuIndexFlatConfig()
    flat_config.useFloat16 = False
    flat_config.device = 0
    index = faiss.GpuIndexFlatL2(res, d, flat_config)

    # perform the training
    clus.train(x, index)  
    _, I = index.search(x, 1)  
    stats = clus.iteration_stats
    losses = np.array([stats.at(i).obj for i in range(stats.size())])
    if verbose:
        print('k-means loss evolution: {0}'.format(losses))

    return [int(n[0]) for n in I], losses[-1]


class Kmeans(object):
    def __init__(self, nmb_clusters):
        self.nmb_clusters = nmb_clusters

    def run(self, data, verbose=False):
        """Performs k-means clustering.
            Args:
                x_data (np.array N * dim): data to cluster
        """
        end = time.time()
        
        # PCA-reducing, whitening and L2-normalization
        self.xb = preprocess_features(data)

        # cluster the data
        I, loss = run_kmeans(self.xb, self.nmb_clusters, verbose) 
        self.I = I
       
        self.data_lists = [[] for i in range(self.nmb_clusters)]
        for i in range(len(data)):
            self.data_lists[I[i]].append(i)

        if verbose:
            print('k-means time: {0:.0f} s'.format(time.time() - end))

        return loss
    
    
class Logger(object):
    """ Class to update every epoch to keep trace of the results
    Methods:
        - log() log and save
    """

    def __init__(self, path):
        self.path = path
        self.data = []

    def log(self, train_point):
        self.data.append(train_point)
        with open(os.path.join(self.path), 'wb') as fp:
            pickle.dump(self.data, fp, -1)   


def arrange_clustering(data_lists):
    pseudolabels = []
    data_indexes = []
    for cluster, datas in enumerate(data_lists):
        data_indexes.extend(datas)
        pseudolabels.extend([cluster] * len(datas))
    indexes = np.argsort(data_indexes)  
    return np.asarray(pseudolabels)[indexes]    
    