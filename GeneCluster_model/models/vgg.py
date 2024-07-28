import torch
import torch.nn as nn
import math
from random import random as rd
import torch.nn.functional as F
import warnings
warnings.filterwarnings("ignore")

# __all__ = [ 'VGG', 'vgg16']
# out_put_size = [1,2,3]
out_put_size = [1,3,6,8,10,16]

class VGG(nn.Module):

    def __init__(self, features, num_classes):
        super(VGG, self).__init__()
        self.features = features
        self.classifier = nn.Sequential(
            nn.Linear(512 * sum(out_put_size), 2048),
            nn.ReLU(True),
            nn.Dropout(0.5),
            nn.Linear(2048, 2048),
            nn.ReLU(True)
        )
        self.top_layer = nn.Linear(2048, num_classes)   
        self._initialize_weights()
        

    def forward(self, x):
        x = self.features(x)
        x = self.classifier(x)
        if self.top_layer:
            x = self.top_layer(x)
        return x

    def _initialize_weights(self):
        """ Weight initialization"""       
        for y,m in enumerate(self.modules()): 
            if isinstance(m, nn.Conv1d): 
                #print(y)
                n = m.kernel_size[0] * m.kernel_size[0] * m.out_channels
                for i in range(m.out_channels):
                    m.weight.data[i].normal_(0, math.sqrt(2. / n))
                if m.bias is not None: 
                    m.bias.data.zero_()
            elif isinstance(m, nn.BatchNorm1d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()
            elif isinstance(m, nn.Linear):
                m.weight.data.normal_(0, 0.01)
                m.bias.data.zero_()


def make_layers(input_dim, batch_norm):
    """ Different convolutional layers are defined according to different configuration files """
    layers = []   
    in_channels = input_dim
    cfg = [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 'M', 512, 512, 512, 'M', 512, 512, 512, 'SPP']  
    for v in cfg:
        if v == 'M':    
            layers += [nn.MaxPool1d(kernel_size=2, stride=2)]
        elif v == 'SPP':
            layers += [SPP_1d(out_put_size)]
        else:           
            conv1d = nn.Conv1d(in_channels, v, kernel_size=3, padding=1)   
            if batch_norm:   
                layers += [conv1d, nn.BatchNorm1d(v), nn.ReLU(inplace=True)]   
            else:
                layers += [conv1d, nn.ReLU(inplace=True)]   
            in_channels = v
    return nn.Sequential(*layers)   

class SPP_1d(nn.Module):
    '''
    1D Spatial pyramid pool layer
    '''
    def __init__(self, out_pool_size): 
        '''
        out_pool_size:The output size of each layer of the pool pyramid
        '''
        super(SPP_1d, self).__init__()
        self.out_pool_size = out_pool_size
        
    def forward(self, x):
        N, C, L = x.size()
        for i in range(len(self.out_pool_size)):
            wid = int(math.ceil(L / self.out_pool_size[i]))   
            pad = wid * self.out_pool_size[i] - L             
            maxpool = nn.MaxPool1d(wid, stride=wid)
            x_pad = F.pad(x, (0, pad), "constant", 0)
            output = maxpool(x_pad)
            if (i == 0):
                spp = output.view(N, -1)
            else:
                spp = torch.cat((spp, output.view(N, -1)), 1)
        return spp

def vgg16(bn=True, out=100):
    """instantiation"""
    dim = 1
    model = VGG(make_layers(dim, bn), out)
    return model

