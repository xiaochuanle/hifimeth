# -*- coding: utf-8 -*-
import os
import random
import re
import sys
import time
import math

import numpy as np
import torch
from torch_geometric.data import Data, Dataset 
from torch_geometric.loader import DataLoader
from tqdm import tqdm

kmer_base_features = 15

'''
The lossy encoding for IPD and pulsewidth values into the available 256 codepoints is as follows (codec v1):
Frames	            Encoding
0 .. 63	0,          1, .. 63
64, 66, .. 190	    64, 65, .. 127
192, 196 .. 444	    128, 129 .. 191
448, 456, .. 952	192, 193 .. 255
'''
def fill_codev1_to_frame_table():

    CodeV1ToFrameTableSize = 256
    codev1_to_frame_table = [0] * CodeV1ToFrameTableSize
    idx = 0

    for i in range(64):
        codev1_to_frame_table[idx] = i 
        idx += 1

    for i in range(64, 128):
        codev1_to_frame_table[idx] = (i - 64) * 2 + 64
        idx += 1

    for i in range(128, 192):
        codev1_to_frame_table[idx] = (i - 128) * 4 + 192
        idx += 1

    for i in range(192, 256):
        codev1_to_frame_table[idx] = (i - 192) * 8 + 448
        idx += 1

    return codev1_to_frame_table

def phred_quality_score_to_prob(phred_quality_scores):
    # QUAL = (-10logP) + 33
    phred_quality_scores -= 33.0
    phred_quality_scores /= -10.0
    return 1.0 - np.power(10.0, phred_quality_scores)

dna_encode_table = [ [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0] ]

def assemble_one_graph_data(sample, num_positions, num_node_features, codev1_to_frame_table, t_edge_index, has_label):
    idx = 0

    qual = np.array(sample[idx:idx+num_positions], dtype=np.float32)
    qual = phred_quality_score_to_prob(qual)
    idx += num_positions

    fseq = np.array(sample[idx:idx+num_positions], dtype=np.uint8)
    idx += num_positions
    
    fipd = np.array([codev1_to_frame_table[x] for x in sample[idx:idx+num_positions]], dtype=np.float32)
    fipd /= 952.0
    idx += num_positions

    fpw = np.array([codev1_to_frame_table[x] for x in sample[idx:idx+num_positions]], dtype=np.float32)
    fpw /= 952.0
    idx += num_positions   

    label = None 
    if has_label:
        label = np.array([sample[idx]], dtype=np.uint8)
        idx += 1

    ripd = np.array([codev1_to_frame_table[x] for x in sample[idx:idx+num_positions]], dtype = np.float32)
    ripd /= 952.0
    idx += num_positions

    rpw = np.array([codev1_to_frame_table[x] for x in sample[idx:idx+num_positions]], dtype = np.float32)
    rpw /= 952.0

    num_nodes = num_positions
    nodes = np.zeros([num_nodes, num_node_features], dtype=np.float32)
    motif_pos = num_positions // 2
    motif_len = 2
    for i in range(num_positions):
        nodes[i][0:4] = dna_encode_table[fseq[i]]
        nodes[i][4:8] = dna_encode_table[3-fseq[i]]
        nodes[i][8] = qual[i]
        nodes[i][9] = fipd[i]
        nodes[i][10] = fpw[i]
        nodes[i][11] = ripd[num_positions - 1 - i]
        nodes[i][12] = rpw[num_positions - 1 - i]
        nodes[i][13] = 1.0
        nodes[i][14] = 0.0

    for i in range(motif_pos, motif_pos + motif_len):
        nodes[i][kmer_base_features-2] = 0.0
        nodes[i][kmer_base_features-1] = 1.0

    t_nodes = torch.from_numpy(nodes)

    if has_label:
        y = np.array([label[0]], dtype=np.int64)
        t_y = torch.from_numpy(y)
        return Data(x=t_nodes, edge_index=t_edge_index, y=t_y, num_nodes=num_nodes)
    else:
        return Data(x=t_nodes, edge_index=t_edge_index, num_nodes=num_nodes)

class KmerFeaturesDataset(Dataset):
    def __init__(self, features_path, num_positions, has_label):
        self.data_root = os.path.dirname(features_path)
        self.features_path = features_path
        self.num_positions = num_positions
        self.sample_elems = num_positions * 6
        if has_label:
            self.sample_elems += 1
        print('has_label: %d,sample elems: %d, num_positions %d' % (has_label, self.sample_elems, num_positions))
        self.has_label = has_label
        super(KmerFeaturesDataset, self).__init__(self.data_root)

        self.samples = np.load(features_path)

        self.n_samples = np.size(self.samples, 0)
        self.node_features = kmer_base_features
        print('Load %d samples' % (self.n_samples))

        num_edges = num_positions - 1
        edge_index = np.zeros([2, num_edges], dtype=np.int64)
        ei = 0
        for i in range(1, num_positions):
            edge_index[0][ei] = i 
            edge_index[1][ei] = i - 1
            ei += 1
        assert ei == num_edges
        self.t_edge_index = torch.from_numpy(edge_index)

        self.codev1_to_frame_table = fill_codev1_to_frame_table()

    @property
    def raw_file_names(self):
        return [self.features_path]

    @property
    def processed_file_names(self):
        return [self.features_path]

    def download(self):
        pass 

    def len(self):
        return self.n_samples

    def get(self, idx):
        return assemble_one_graph_data(self.samples[idx], self.num_positions, self.node_features, self.codev1_to_frame_table, self.t_edge_index, self.has_label)

    def print_sample(self, idx):
        sample = self.samples[idx]
        vocab = "ACGT"
        kmer = [vocab[x] for x in sample[self.num_positions:self.num_positions+20]]
        print('sample %d\t' % (idx))
        print(kmer)


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
    num_positions = 400
    label_size = 2
    batch_size = 256
    has_label = 1

    dataset = KmerFeaturesDataset(features_path, num_positions, has_label) 
    dataloader = DataLoader(dataset, batch_size=512, shuffle=True, num_workers=32)

    n_samples = dataset.n_samples
    for i in range(10):
        dataset.print_sample(i)
    for i in range(n_samples - 10, n_samples):
        dataset.print_sample(i)

    batch_i = 0
    for batch in tqdm(dataloader):
        pass

    '''
    for bi, batch in enumerate(dataloader):
        print('batch %d, samples: %d' % (bi, len(batch)))
        g0 = batch[0]
        print(g0.x)
    '''

if __name__ == "__main__":
    main(sys.argv)
