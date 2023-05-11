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

import kmer_features_dataset
from kmer_features_dataset import kmer_base_features
from kmer_features_dataset import fill_codev1_to_frame_table
from kmer_features_dataset import dna_encode_table

input_base_features = 6

def phred_quality_score_to_prob(phred_quality_scores):
    # QUAL = (-10logP) + 33
    phred_quality_scores = float(phred_quality_scores)
    phred_quality_scores -= 33.0
    phred_quality_scores /= -10.0
    return 1.0 - math.pow(10.0, phred_quality_scores)

def build_base_feature(xlist, codev1_to_frame_table):
    assert len(xlist) == input_base_features
    f = [0.0] * kmer_base_features
    f[0:4] = dna_encode_table[xlist[0]]
    f[4:8] = dna_encode_table[3-xlist[0]]
    f[8] = phred_quality_score_to_prob(xlist[1])
    f[9] = codev1_to_frame_table[xlist[2]]/952.0
    f[10] = codev1_to_frame_table[xlist[3]]/952.0
    f[11] = codev1_to_frame_table[xlist[4]]/952.0
    f[12] = codev1_to_frame_table[xlist[5]]/952.0
    f[13] = 1.0
    f[14] = 0.0
    return f

class ReadGraphDataset(Dataset):
    def __init__(self, features_path, sample_offset_path, kmer_size, motif, has_label):
        self.data_root = os.path.dirname(features_path)
        self.features_path = features_path
        self.sample_offset_path = sample_offset_path
        self.kmer_size = kmer_size 
        self.motif = motif 
        self.has_label = has_label

        super(ReadGraphDataset, self).__init__(self.data_root)

        f = open(features_path, "rb")
        strb = f.read()
        f.close()
        self.features = np.frombuffer(strb, dtype=np.uint8)
        self.features = self.features.reshape([-1, input_base_features])
        del strb 

        self.samples = np.loadtxt(sample_offset_path)
        self.n_samples = len(self.samples)

        print('Load %d features including %d samples' % (self.features.shape[0], self.n_samples))

        num_edges = kmer_size - 1
        edge_index = np.zeros([2, num_edges], dtype=np.int64)
        ei = 0
        for i in range(1, kmer_size):
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
        half_kmer_size = int(self.kmer_size / 2)
        feature_idx = int(self.samples[idx][1])
        fidx_from = int(feature_idx - half_kmer_size)
        fidx_to = int(feature_idx + half_kmer_size)
        nodes = [ build_base_feature(self.features[p], self.codev1_to_frame_table) for p in range(fidx_from, fidx_to) ]
        motif_size = len(self.motif)
        nodes = np.array(nodes, dtype=np.float32)
        for i in range(half_kmer_size, half_kmer_size + motif_size):
            nodes[i][kmer_base_features-2] = 0.0
            nodes[i][kmer_base_features-1] = 1.0
        t_nodes = torch.from_numpy(nodes)
        return Data(x=t_nodes, edge_index=self.t_edge_index, num_nodes=self.kmer_size)

    def print_sample(self, idx):
        vocab = "ACGT"
        half_kmer_size = self.kmer_size / 2
        read_idx = int(self.samples[idx][0])
        feature_idx = int(self.samples[idx][1])
        fidx_from = int(feature_idx - half_kmer_size)
        fidx_to = int(feature_idx + half_kmer_size)

        vocab = "ACGT"
        kmer = []
        for i in range(fidx_from, fidx_to):
            c = vocab[int(self.features[i][0])]
            kmer.append(c)
        print('sample %d, read_id = %d, kmer = %s' % (idx, read_idx, kmer))

def dump_usage(program):
    print('')
    print('USAGE:')
    print('%s features-path sample_offset_path kmer-size' % (program))

def main(argv):
    n_param = len(argv)
    if n_param < 4:
        dump_usage(argv[0])
        sys.exit(1)

    features_path = argv[1]
    sample_offset_path = argv[2]
    kmer_size = int(argv[3])
    motif = 'CG'
    has_label = 0

    dataset = ReadGraphDataset(features_path, sample_offset_path, kmer_size, motif, kmer_size)
    dataloader = DataLoader(dataset, batch_size = 512, shuffle = False, num_workers = 8)
    for i in range(10):
        dataset.print_sample(i)
    for i in range(dataset.n_samples-10, dataset.n_samples):
        dataset.print_sample(i)

    batch_i = 0
    for batch in tqdm(dataloader):
        g0 = batch[0]
        if batch_i < 10:
            print(g0.x)
        batch_i += 1

    '''
    for bi, batch in enumerate(dataloader):
        print('batch %d, samples: %d' % (bi, len(batch)))
        g0 = batch[0]
        print(g0.x)
    '''

if __name__ == "__main__":
    main(sys.argv)