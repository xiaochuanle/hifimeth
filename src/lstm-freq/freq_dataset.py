# -*- coding: utf-8 -*-
import os
import random
import re
import sys
import time
import math

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

mod_site_features = 23
mod_site_input_features = 21

def assemble_one_graph_data(sample, num_sites, has_label):

    num_nodes = num_sites
    nodes = np.zeros([num_nodes, mod_site_features], dtype=np.float32)
    idx = 0
    for i in range(num_sites):
        nodes[i][0:mod_site_input_features] = sample[idx:idx+mod_site_input_features]
        nodes[i][mod_site_input_features] = 0.0
        nodes[i][mod_site_input_features + 1] = 1.0
        idx += mod_site_input_features
    target_pos = num_sites // 2
    nodes[target_pos][mod_site_input_features] = 1.0
    nodes[target_pos][mod_site_input_features+1] = 0.0

    t_nodes = torch.from_numpy(nodes)

    if has_label:
        y = np.array([sample[-1]], dtype=np.float32)
        t_y = torch.from_numpy(y)
        return t_nodes, t_y
    else:
        return t_nodes

class FreqFeaturesDataset(Dataset):
    def __init__(self, features_path, num_sites, has_label):
        self.features_path = features_path
        self.num_sites = num_sites
        self.has_label = has_label
        super(FreqFeaturesDataset, self).__init__()

        self.samples = np.loadtxt(features_path)
        self.n_samples = np.size(self.samples, 0)
        self.sample_elems = np.size(self.samples, 1)
        self.node_features = mod_site_features
        print('has_label: %d,sample elems: %d, num_sites %d' % (has_label, self.sample_elems, num_sites))
        print('Load %d samples' % (self.n_samples))

        n_pos, n_neg = 0, 0
        for i in range(self.n_samples):
            if self.samples[i][-1] > 0.5:
                n_pos += 1
            else:
                n_neg += 1
        print('Pos samples: %d neg samples: %d' % (n_pos, n_neg))

    def __len__(self):
        return self.n_samples

    def __getitem__(self, idx):
        return assemble_one_graph_data(self.samples[idx], self.num_sites, self.has_label)

    def print_sample(self, idx):
        sample = self.samples[idx]
        print('dist: %f, label: %f' % (sample[0], sample[-1]))
        print('mod at site 0')
        print(sample[1:1+mod_site_features])

def dump_usage(program):
    print('')
    print('USAGE:')
    print('%s features-path' % (program))

def main(argv):
    n_param = len(argv)
    if n_param < 2:
        dump_usage(argv[0])
        sys.exit(1)

    features_path = argv[1]
    batch_size = 256
    num_sites = 11
    has_label = 1

    dataset = FreqFeaturesDataset(features_path, num_sites, has_label) 
    dataloader = DataLoader(dataset, batch_size=512, shuffle=True, num_workers=32)

    n_samples = dataset.n_samples
    for i in range(10):
        dataset.print_sample(i)
    for i in range(n_samples - 10, n_samples):
        dataset.print_sample(i)

    batch_i = 0
    for x, y in tqdm(dataloader):
        pass

    '''
    for bi, batch in enumerate(dataloader):
        print('batch %d, samples: %d' % (bi, len(batch)))
        g0 = batch[0]
        print(g0.x)
    '''

if __name__ == "__main__":
    main(sys.argv)
