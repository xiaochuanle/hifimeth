import torch
import torch.nn as nn
import numpy as np
import torch.nn.functional as F 

from sample_dataset import KMER_BASE_FEATURES

class DNAModNet(torch.nn.Module):
    def __init__(self, kmer):
        super(DNAModNet, self).__init__()

        hd0 = KMER_BASE_FEATURES
        hd1 = 128
        hd2 = 128
        hd3 = 128
        hd4 = 96
        hd5 = 96
        hd6 = 96
        hd7 = 64
        hd8 = 64

        ks1 = 13
        ks2 = 3
        ks3 = 3
        ks4 = 3
        ks5 = 3
        ks6 = 3
        ks7 = 3
        ks8 = 3

        self.bn0 = nn.BatchNorm1d(hd0)
        self.convs = nn.Sequential(
            nn.Conv1d(hd0, hd1, ks1, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd1),
            nn.ReLU(),

            nn.Conv1d(hd1, hd2, ks2, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd2),
            nn.ReLU(),

            nn.Conv1d(hd2, hd3, ks3, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd3),
            nn.ReLU(),

            nn.Conv1d(hd3, hd4, ks4, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd4),
            nn.ReLU(),

            nn.Conv1d(hd4, hd5, ks5, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd5),
            nn.ReLU(),

            nn.Conv1d(hd5, hd6, ks6, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd6),
            nn.ReLU(),

            nn.Conv1d(hd6, hd7, ks7, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd7),
            nn.ReLU(),

            nn.Conv1d(hd7, hd8, ks8, stride = 2, padding = 1, bias = False),
            nn.BatchNorm1d(hd8),
            nn.ReLU(),
        )
        X = torch.zeros(1, KMER_BASE_FEATURES, kmer)
        X = self.convs(X)
        conv_output_dim = X.view(-1).shape[0]
        self.fc1 = nn.Linear(conv_output_dim, 256)
        self.fc2 = nn.Linear(256, 2)

        num_params = sum(p.numel() for p in self.parameters())
        print(f"Total params: {num_params}")
        num_trainable_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print(f"Total trainable params: {num_trainable_params}")

    def forward(self, X):
        X = X.permute(0, 2, 1)
        X = self.bn0(X)
        X = self.convs(X)

        X = torch.flatten(X, 1)
        X = self.fc1(X)
        X = F.relu(X)
        X = self.fc2(X)

        return X