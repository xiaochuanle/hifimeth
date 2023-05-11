import torch
from torch.nn import Linear
import torch.nn.functional as F 
from torch_geometric.nn import GraphConv, GatedGraphConv
from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool
from torch_geometric.nn import BatchNorm
from torch_geometric.data import Data, Dataset 
from torch_geometric.loader import DataLoader

import graph_conv_backbone
from graph_conv_backbone import GNN_graph_conv_backbone

class KmerModel(torch.nn.Module):
    def __init__(self, num_node_features, hidden_channels = 128, num_classes = 2):
        super(KmerModel, self).__init__()

        self.num_node_features = num_node_features
        self.num_classes = num_classes

        self.net1 = GNN_graph_conv_backbone(num_node_features, hidden_channels, num_classes)
        self.lin1 = Linear(hidden_channels, num_classes)
        #self.net2 = GNN_graph_conv_backbone(num_node_features, hidden_channels, num_classes)

        #self.lin = Linear(2 * hidden_channels, num_classes)   

        n_params = sum(p.numel() for p in self.parameters())
        print('Total params: %d' % (n_params))
        n_trainable_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print('Total trainable params: %d' % (n_trainable_params))

    def forward(self, x, edge_index, batch):
        '''
        num_node_features = self.num_node_features
        x1, x2 = torch.split(x, num_node_features, dim=-1)
        x1 = self.net1(x1, edge_index, batch)
        x2 = self.net2(x2, edge_index, batch)    

        x = torch.cat([x1, x2], dim=-1)    

        x = self.lin(x)
        
        return x 
        '''

        x = self.net1(x, edge_index, batch)
        x = self.lin1(x)
        return x
