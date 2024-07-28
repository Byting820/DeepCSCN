import os
import sys
import time
import argparse
from tkinter import Variable
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import torch.nn.functional as F
from torch.autograd import Variable
from sklearn.cluster import KMeans
from torch.optim import lr_scheduler

import torch
from torch import nn
import torch.nn.parallel
import torch.backends.cudnn as cudnn
from sklearn.metrics.cluster import normalized_mutual_info_score,homogeneity_score

from torch.utils.data import dataset
from sklearn.metrics import silhouette_score,calinski_harabasz_score,adjusted_rand_score

import models
import cluster
from utils import ExpressionDataset, pseudolabels_assign, UnifLabelSampler, AverageMeter

def parse_args():
    parser = argparse.ArgumentParser()
    
    # working_dir = os.path.dirname(os.path.abspath('main_copy2.py')) 
    parser.add_argument('--nmb_cluster', type=int, default=10,
                        help='number of cluster for k-means (default: 100)')
    parser.add_argument('--lr', type=float, default=0.1,
                        help='learning rate (default: 0.05)')
    parser.add_argument('--momentum', type=float, default=0.9,
                        help='momentum (default: 0.9)')
    parser.add_argument('--wd', default=-5, type=float,
                        help='weight decay pow (default: -5)')
    parser.add_argument('--data_path', metavar='PATH', help='path to dataset',
                        default='count.csv')

    parser.add_argument('--batch', type=int, default=32,
                        help='batch size')
    parser.add_argument('--workers', type=int, default=20,
                        help='number of data loading workers')
    parser.add_argument('--epochs', type=int, default=200,
                        help='epoch number')
    parser.add_argument('--clustering', type=str, choices=['Kmeans'],
                        default='Kmeans', help='clustering algorithm (default: Kmeans)')
    parser.add_argument('--arch', type=str, choices=['sppCNN'],
                        default='sppCNN', help='feature embedding model architecture')
    parser.add_argument('--reassign', type=float, default=1.,
                        help="""how many epochs of training between two consecutive
                        reassignments of clusters (default: 1)""")
    parser.add_argument('--ckpt_path', type=str, default='train_res', help='')
    parser.add_argument('--checkpoints', type=int, default=1000,
                        help='how many iterations between two checkpoints (default: 1000)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        default=True, help='verbose mode')

    return parser.parse_args()


def compute_features(dataloader, model, N):
    end = time.time()
    model.eval()  # close train model
    for batch, (input_tensor, _) in enumerate(dataloader):  
        #input_tensor = input_tensor.to()
        with torch.no_grad():
            aux = model(input_tensor).data.cpu().numpy() 
        
        if batch == 0:
            features = np.zeros((N, aux.shape[1]), dtype='float32')  

        aux = aux.astype('float32')  
        if batch < len(dataloader) - 1:  
            features[batch * args.batch: (batch + 1) * args.batch] = aux
        else:

            features[batch * args.batch:] = aux

        time_cost = time.time() - end
        end = time.time()

        if (batch % 200) == 0:
            print('{0} / {1}\t'
                  'Time cost: {time:.3f} s'
                  .format(batch, len(dataloader), time=time_cost))
    return features     
    

def Visualization(features, data_lists=None): 
    """ Cluster visualization """
   
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(features)

    if data_lists is None:
        kmeans = KMeans(n_clusters=args.nmb_cluster,random_state=10)     
        kmeans.fit(embedding)
        gene_label = kmeans.predict(embedding)
        
    else:
        gene_label = np.zeros(len(features))      
        for i in range(len(data_lists)):
            gene_label[data_lists[i]] = [i for m in range(len(data_lists[i]))]   

    
    picture = plt.scatter(embedding[:, 0], embedding[:, 1], s=5, c=gene_label)             
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar(boundaries=np.arange(11)-0.5).set_ticks(np.arange(10))
    plt.title('UMAP projection of the dataset')
    # plt.gca().legend().remove()
    return picture


class FocalLoss(nn.Module):
    def __init__(self, alpha=None, gamma=2, reduction='mean'):
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction
        
    def forward(self, output, targets):
        ce_loss = F.cross_entropy(output, targets, reduction='none')
        pt = torch.exp(-ce_loss)
        if self.alpha is not None:
            alpha_t = self.alpha[targets]
            focal_loss = alpha_t * (1 - pt) ** self.gamma * ce_loss
        else:
            focal_loss = (1 - pt) ** self.gamma * ce_loss
        if self.reduction == 'mean':
            return torch.mean(focal_loss)
        elif self.reduction == 'sum':
            return torch.sum(focal_loss)
        else:
            return focal_loss
        

def main(args):
    # model
    # os.environ["CUDA_VISIBLE_DEVICES"] = "0"
    model = models.__dict__[args.arch](out=args.nmb_cluster)
    fea_dim = int(model.top_layer.weight.size()[1])  # torch.Size([100, 200])
    model.top_layer = None
    model.features = torch.nn.DataParallel(model.features)  
    #if torch.cuda.is_available():
     #   device = torch.device("cuda")
    #else:
    #    device = torch.device("cpu")
    # device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    #model.features.to(device)
    
    model.cuda()
    cudnn.benchmark = True
    
    # load data
    end = time.time()
    dataset = ExpressionDataset(args.data_path)

    if args.verbose:
        print('Load dataset: {0:.2f} s'.format(time.time() - end))

    dataloader = torch.utils.data.DataLoader(dataset,
                                            batch_size=args.batch,
                                            num_workers=args.workers,
                                            pin_memory=True)


    deepcluster = cluster.__dict__[args.clustering](args.nmb_cluster) 
    if args.verbose:
            print('Cluster the features with {}'.format(args.clustering))
    
    
    # create optimizer
    optimizer = torch.optim.SGD(
        filter(lambda x: x.requires_grad, model.parameters()),
        lr=args.lr,
        momentum=args.momentum,
        weight_decay=10**args.wd,
    )

    scheduler = lr_scheduler.ExponentialLR(optimizer, gamma=0.9)

    # define loss function
    # loss_fn = nn.CrossEntropyLoss().cuda()
    # loss_fn = FocalLoss(args.nmb_cluster).cuda()
    loss_fn = FocalLoss(gamma=2, reduction='mean').cuda()
    
    # creating checkpoint repo
    ckpt_path = os.path.join(args.ckpt_path, 'checkpoints')
    if not os.path.isdir(ckpt_path):
        os.makedirs(ckpt_path)
    
    cluster_log = cluster.Logger(os.path.join(args.ckpt_path,'clusters'))

    last_update = 0
    patience = 20     
    for epoch in range(args.epochs):
        end = time.time()
        
        model.top_layer = None
        model.classifier = nn.Sequential(*list(model.classifier.children())[:-1]) 
        features = compute_features(dataloader, model, len(dataset))
        print(features.shape)
        # print(np.isnan(features))
        np.save(args.ckpt_path+'/features.npy',features)

        clustering_loss = deepcluster.run(features, verbose=args.verbose) 
        
        if epoch % 20 == 0:
            visual = Visualization(features, deepcluster.data_lists)

            plt.savefig(os.path.join(args.ckpt_path, 'Epoch{}_cluster_picture.png'.format(epoch)))
            plt.close()


        if args.verbose:
            print('Assign pseudo labels')            
        train_dataset = pseudolabels_assign(deepcluster.data_lists,
                                            dataset.datas)

        sampler = UnifLabelSampler(int(args.reassign * len(train_dataset)),   
                                   deepcluster.data_lists)
        
        train_dataloader = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=args.batch,
            num_workers=args.workers,
            sampler=sampler,
            pin_memory=True)
        
        # set last fully connected layer
        mlp = list(model.classifier.children())
        mlp.append(nn.ReLU(inplace=True).cuda())  
        model.classifier = nn.Sequential(*mlp)
        model.top_layer = nn.Linear(fea_dim, len(deepcluster.data_lists))
        model.top_layer.weight.data.normal_(0, 0.01)
        model.top_layer.bias.data.zero_()
        model.top_layer.cuda()

         
        loss = train(train_dataloader, model, loss_fn, optimizer, epoch)
        scheduler.step()
        
        gene_list = pd.DataFrame({'clus_id':train_dataset.labels})
        gene_list['clus_id'] = gene_list['clus_id'].astype(int)
        gene_list.index = pd.read_csv(args.data_path,sep = '\t',index_col = 0).index
        gene_list.to_csv(args.ckpt_path+'/genecluster_res.txt')

        print(train_dataset.datas.shape)
        print(train_dataset.labels.shape)
        
        silhouette_avg = silhouette_score(features,train_dataset.labels)
        ch_score = calinski_harabasz_score(features,train_dataset.labels)
        print('silhouette_avg: {:.3f}\n'
              'ch_score: {:.3f}\n'
              .format(silhouette_avg,ch_score))


        # print log
        if args.verbose:
            print('###### Epoch [{0}] ###### \n'
                  'Time: {1:.3f} s\n'
                  'Clustering loss: {2:.3f} \n'
                  'ConvNet loss: {3:.3f}\n'
                  .format(epoch, time.time() - end, clustering_loss, loss))
            
            try:
                NMI = normalized_mutual_info_score(
                        cluster.arrange_clustering(deepcluster.data_lists),
                        cluster.arrange_clustering(cluster_log.data[-1])
                    )
                print('NMI against previous assignment: {0:.3f}'.format(NMI))
            except IndexError:
                pass

            cluster_log.log(deepcluster.data_lists)
            #print(cluster_log.data)
            print('####################### \n')
            
        # save running checkpoint
        torch.save({'epoch': epoch + 1,
		            'arch': args.arch,
                    'clustering': args.clustering,
                    'state_dict': model.state_dict(),
                    'optimizer' : optimizer.state_dict()},
                   os.path.join(args.ckpt_path, 'checkpoint.pth.tar'))

        # # early stopping
        if epoch == 0:
             best_loss = float('inf')
        else:
             if loss < best_loss:
                 best_loss = loss
                 last_update = epoch
             else:
                 if epoch - last_update > patience:
                     print("Early stopping!")
                     break        


    
def train(dataloader, model, loss_fn, optimizer, epoch):

    losses = AverageMeter()
    batch_time = AverageMeter()
    
    # switch to train mode
    model.train()

    # create an optimizer for the last fc layer
    optimizer_tl = torch.optim.SGD(
        model.top_layer.parameters(),
        lr=args.lr,
        weight_decay=10**args.wd,
    )

    for batch, (input_tensor, target) in enumerate(dataloader):
        end = time.time()

        # save checkpoint
        n = len(dataloader) * epoch + batch
        if n % args.checkpoints == 0:
            path = os.path.join(
                args.ckpt_path,
                'checkpoints',
                'checkpoint_' + str(n / args.checkpoints) + '.pth.tar',
            )
            if args.verbose:
                print('Save checkpoint at: {0}'.format(path))
            torch.save({
                'epoch': epoch + 1,
		        'arch': args.arch,
                'clustering': args.clustering,
                'state_dict': model.state_dict(),
                'optimizer' : optimizer.state_dict()
            }, path)

        target = target.cuda()
        input_var = torch.autograd.Variable(input_tensor.cuda())
        target_var = torch.autograd.Variable(target)

        output = model(input_var)
        # probs = F.softmax(output,dim=1)
        # # pred = np.argmax(probs)
        # print(probs)

        # loss = loss_fn(output, target_var)
        loss = loss_fn(output, target_var)

        # record loss
        losses.update(loss.item(), input_tensor.size(0))

        # compute gradient and do SGD step
        optimizer.zero_grad()
        optimizer_tl.zero_grad()
        loss.backward()    
        optimizer.step()
        optimizer_tl.step()

        # measure elapsed time
        batch_time.update(time.time() - end)

        if args.verbose and (batch % 200) == 0:
            print('Epoch: [{0}][{1}/{2}]\t'
                  'Time: {batch_time.val:.3f} ({batch_time.avg:.3f})\t'
                  'Loss: {loss.val:.4f} ({loss.avg:.4f})'
                  .format(epoch, batch, len(dataloader),
                          batch_time=batch_time,
                          loss=losses))

    return losses.avg
    
if __name__ == '__main__':
    start = time.time()
    args = parse_args()
    sys.stdout = open(args.ckpt_path + '/result.log', mode = 'w')
    main(args)
    total_time = (time.time()-start)/60
    print('model_run_time' + str(total_time) + 'min')
