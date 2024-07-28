import os
import pickle

import pandas as pd
import numpy as np

import torch
import models
from torch.utils.data import Dataset
from torch.utils.data.sampler import Sampler


def load_model(path):
    '''
    load model from checkpoint path
    path: checkpoint file path
    '''
    if os.path.isfile(path):
        print("=> loading checkpoint '{}'".format(path))
        checkpoint = torch.load(path)

        # size of the top layer
        N = checkpoint['state_dict']['top_layer.bias'].size()

        # build skeleton of the model
        model = models.__dict__[checkpoint['arch']](out=int(N[0]))

        # deal with a dataparallel table
        def rename_key(key):
            if not 'module' in key:
                return key
            return ''.join(key.split('.module'))

        checkpoint['state_dict'] = {rename_key(key): val
                                    for key, val
                                    in checkpoint['state_dict'].items()}

        # load weights
        model.load_state_dict(checkpoint['state_dict'])
        print("Loaded")
    else:
        model = None
        print("=> no checkpoint found at '{}'".format(path))
    return model


class ExpressionDataset(Dataset):
    '''
    custom Dataset class for classifier training.
    __len__: returns the number of samples in our dataset.
    __getitem__:  loads and returns a sample from the dataset at the given index idx.
    '''
    def __init__(self, data_path):
        '''
        label_file: label
        '''
        self.datas = pd.read_csv(data_path, sep='\t', index_col=0)    ##load datas
        self.datas = np.expand_dims(np.array(self.datas), axis=1)
        self.labels = np.zeros((len(self.datas)), dtype='float32')
    
    def __len__(self):
        return len(self.datas)

    def __getitem__(self, idx):
        data = self.datas[idx]
        label = self.labels[idx]
        data = torch.FloatTensor(data)
        label = torch.LongTensor(np.array(label))
        return data, label
    
    
class ReassignedDataset(Dataset):

    def __init__(self, data_indexes, pseudolabels, dataset):
        self.datas, self.labels = self.make_dataset(data_indexes, pseudolabels, dataset)
    
    def make_dataset(self, data_indexes, pseudolabels, dataset):
        label_to_idx = {label: idx for idx, label in enumerate(set(pseudolabels))}
        datas = np.zeros(dataset.shape, dtype='float32')
        labels = np.zeros((len(pseudolabels)), dtype='float32')
        for j, idx in enumerate(data_indexes):  
            data = dataset[idx]  
            pseudolabel = label_to_idx[pseudolabels[j]]
            datas[idx] = data    
            labels[idx] = pseudolabel
        return datas, labels     
    
    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        data = self.datas[idx]
        label = self.labels[idx]
        data = torch.FloatTensor(data)
        label = torch.LongTensor(np.array(label))
        return data, label
    
    
def pseudolabels_assign(datas_lists, dataset):
    """Creates a dataset from clustering, with clusters as labels.
    Args:
        data_lists (list of list): for each cluster, the list of data indexes
                                    belonging to this cluster
        dataset (list): initial dataset
    Returns:
        ReassignedDataset(torch.utils.data.Dataset): a dataset with clusters as
                                                     labels
    """
    assert datas_lists is not None
    pseudolabels = []
    data_indexes = []
    for cluster, datas in enumerate(datas_lists):
        data_indexes.extend(datas)
        pseudolabels.extend([cluster] * len(datas))
    return ReassignedDataset(data_indexes, pseudolabels, dataset)


class UnifLabelSampler(Sampler):
    """Samples elements uniformely accross pseudolabels.
        Args:
            N (int): size of returned iterator.
            images_lists: dict of key (target), value (list of data with this target)
    """

    def __init__(self, N, data_lists):
        self.N = N
        self.data_lists = data_lists
        self.indexes = self.generate_indexes_epoch()

    def generate_indexes_epoch(self):
        nmb_non_empty_clusters = 0
        for i in range(len(self.data_lists)):
            if len(self.data_lists[i]) != 0:
                nmb_non_empty_clusters += 1

        size_per_pseudolabel = int(self.N / nmb_non_empty_clusters) + 1
        res = np.array([])

        for i in range(len(self.data_lists)):
            # skip empty clusters
            if len(self.data_lists[i]) == 0:
                continue
            indexes = np.random.choice(
                self.data_lists[i],
                size_per_pseudolabel,
                replace=(len(self.data_lists[i]) <= size_per_pseudolabel)
            )
            res = np.concatenate((res, indexes))

        np.random.shuffle(res)
        res = list(res.astype('int'))
        if len(res) >= self.N:
            return res[:self.N]
        res += res[: (self.N - len(res))]
        return res

    def __iter__(self):
        return iter(self.indexes)

    def __len__(self):
        return len(self.indexes)
    
    
    
class AverageMeter(object):
    """Computes and stores the average and current value"""
    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count
        
    
