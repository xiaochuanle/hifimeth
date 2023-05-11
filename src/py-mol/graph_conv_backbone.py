import torch
from torch.nn import Linear
import torch.nn.functional as F 
from torch_geometric.nn import GraphConv, GatedGraphConv
from torch_geometric.nn import global_mean_pool, global_max_pool, global_add_pool
from torch_geometric.nn import BatchNorm
from torch_geometric.data import Data, Dataset 
from torch_geometric.loader import DataLoader

class GNN_graph_conv_backbone(torch.nn.Module):
    def __init__(self, num_node_features, hidden_channels, num_classes):
        super(GNN_graph_conv_backbone, self).__init__()

        aggr = 'mean'
        self.conv1 = GraphConv(num_node_features, hidden_channels, aggr=aggr)
        self.conv2 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv3 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv4 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv5 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv6 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv7 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv8 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv9 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv10 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv11 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        '''
        self.conv12 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv13 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv14 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        self.conv15 = GraphConv(hidden_channels, hidden_channels, aggr=aggr)
        '''
        self.conv = GraphConv(hidden_channels, hidden_channels, aggr=aggr)

        self.bn1 = BatchNorm(num_node_features)
        self.bn2 = BatchNorm(hidden_channels)
        self.bn3 = BatchNorm(hidden_channels)
        self.bn4 = BatchNorm(hidden_channels)
        self.bn5 = BatchNorm(hidden_channels)
        self.bn6 = BatchNorm(hidden_channels)
        self.bn7 = BatchNorm(hidden_channels)
        self.bn8 = BatchNorm(hidden_channels)
        self.bn9 = BatchNorm(hidden_channels)
        self.bn10 = BatchNorm(hidden_channels)
        self.bn11 = BatchNorm(hidden_channels)
        '''
        self.bn12 = BatchNorm(hidden_channels)
        self.bn13 = BatchNorm(hidden_channels)
        self.bn14 = BatchNorm(hidden_channels)
        self.bn15 = BatchNorm(hidden_channels)
        '''
        self.bn = BatchNorm(hidden_channels)

    def forward(self, x, edge_index, batch):
        x = self.bn1(x)
        x = self.conv1(x, edge_index)
        x = x.relu()

        x = self.bn2(x)
        x = self.conv2(x, edge_index)
        x = x.relu()

        x = self.bn3(x)
        x = self.conv3(x, edge_index)
        x = x.relu()

        x = self.bn4(x)
        x = self.conv4(x, edge_index)
        x = x.relu()

        x = self.bn5(x)
        x = self.conv5(x, edge_index)
        x = x.relu()

        x = self.bn6(x)
        x = self.conv6(x, edge_index)
        x = x.relu()

        x = self.bn7(x)
        x = self.conv7(x, edge_index)
        x = x.relu()

        x = self.bn8(x)
        x = self.conv8(x, edge_index)
        x = x.relu()

        x = self.bn9(x)
        x = self.conv9(x, edge_index)
        x = x.relu()

        x = self.bn10(x)
        x = self.conv10(x, edge_index)
        x = x.relu()

        x = self.bn11(x)
        x = self.conv11(x, edge_index)
        x = x.relu()

        '''
        x = self.bn12(x)
        x = self.conv12(x, edge_index)
        x = x.relu()

        x = self.bn13(x)
        x = self.conv13(x, edge_index)
        x = x.relu()

        x = self.bn14(x)
        x = self.conv14(x, edge_index)
        x = x.relu()

        x = self.bn15(x)
        x = self.conv15(x, edge_index)
        x = x.relu()
        '''

        x = self.bn(x)
        x = self.conv(x, edge_index)

        x = global_mean_pool(x, batch)    

        return x
