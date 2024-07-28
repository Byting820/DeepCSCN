import math
from typing import Tuple

import torch
from torch import nn, Tensor
import torch.nn.functional as F
from torch.nn import TransformerEncoder, TransformerEncoderLayer
from torch.utils.data import dataset

config = [(256, 3, 1, 1), (512, 3, 1, 1), (1024, 3, 1, 1),(1024, 3, 1, 1), (512, 3, 1, 1), (256, 3, 1, 1), "SPP"]
out_put_size = [1,3,6,8,10,16]

class SPPNet(nn.Module):
    def __init__(self, features, num_classes):
        super(SPPNet, self).__init__()
        self.features = features
        self.classifier = nn.Sequential(nn.Dropout(0.2),
                            nn.Linear(256 * sum(out_put_size), 2048), 
                            nn.ReLU(inplace=False),
                            nn.Dropout(0.2),
                            nn.Linear(2048, 2048),
                            nn.ReLU(inplace=False))
                            
        self.top_layer = nn.Linear(2048, num_classes)
        self._initialize_weights()
        
    def forward(self, x):
        x = self.features(x)
        x = self.classifier(x)
        if self.top_layer:
            x = self.top_layer(x)
        return x

    def _initialize_weights(self):
        for y, m in enumerate(self.modules()):
            if isinstance(m, nn.Conv1d):            
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

          

def make_layers_features(config, input_dim, bn):
    layers = []
    in_channels = input_dim
    for i in config:
        if i == 'MP':
            layers += [nn.MaxPool1d(kernel_size=3, stride=2)]
        elif i == 'SPP':
            layers += [SPP_1d(out_put_size)]
        else:
            conv1d = nn.Conv1d(in_channels, i[0], kernel_size=i[1],
                               stride=i[2], padding= i[3])
            if bn: 
                layers += [conv1d, nn.BatchNorm1d(i[0]), nn.ReLU(inplace=False)]  
            else:
                layers += [conv1d, nn.ReLU(inplace=False)]
            in_channels = i[0]
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
            

             
def sppCNN(bn=True, out=100):
    dim = 1
    model = SPPNet(make_layers_features(config, dim, bn=bn), out)
    return model

