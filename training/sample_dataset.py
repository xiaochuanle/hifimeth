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

INPUT_BASE_FEATURES = 5
KMER_BASE_FEATURES = 8

kMaxKineticsValue = 952.0

'''
The lossy encoding for IPD and pulsewidth values into the available 256 codepoints is as follows (codec v1):
Frames              Encoding
0 .. 63 0,          1, .. 63
64, 66, .. 190      64, 65, .. 127
192, 196 .. 444     128, 129 .. 191
448, 456, .. 952        192, 193 .. 255
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

    codev1_to_frame_table = [ 1.0*v/kMaxKineticsValue for v in codev1_to_frame_table ]

    return codev1_to_frame_table

codev1_to_frame_table = fill_codev1_to_frame_table()

def load_features(features_path):
    f = open(features_path, 'rb')
    strb = f.read()
    f.close()
    features = np.frombuffer(strb, dtype=np.uint8)
    del strb 
    nf = features.shape[0]
    print(f"features-shape: {features.shape}")
    assert (nf % INPUT_BASE_FEATURES) == 0
    return features

def load_samples(samples_path):
    samples = np.loadtxt(samples_path, dtype = np.uint64)
    print(f"Sample data shape: {samples.shape}")
    return samples

def load_offsets(offsets_path):
    dtype = [ ('offset', np.int64), ('id', np.int32), ('size', np.int32), ('fn', np.int32), ('rn', np.int32) ]
    offsets = np.loadtxt(offsets_path, dtype=dtype, delimiter='\t')
    print(f"offsets data shape: {offsets.shape}")
    return offsets

dna_one_hot_encoding = [
    [ 1.0, 0.0, 0.0, 0.0 ],
    [ 0.0, 1.0, 0.0, 0.0 ],
    [ 0.0, 0.0, 1.0, 0.0 ],
    [ 0.0, 0.0, 0.0, 1.0 ]
]

def assemble_one_sample_features(features, samples, offsets, kmer, idx):
    qid = samples[idx][0]
    qoff = samples[idx][1]
    label = int(samples[idx][2])

    base_offset = np.uint64(offsets[qid][0]) * INPUT_BASE_FEATURES
    qsize = offsets[qid][2]

    fipd_offset = np.uint64(base_offset + qsize)
    fpw_offset = np.uint64(fipd_offset + qsize)
    ripd_offset = np.uint64(fpw_offset + qsize)
    rpw_offset = np.uint64(ripd_offset + qsize)

    p1 = np.uint64(base_offset + qoff)
    c1 = features[p1]
    assert c1 == 1 or c1 == 2, f"c = {c1}"

    hk = kmer // 2
    ri = np.uint64(0) if qoff >= hk else np.uint64(hk - qoff)
    qfrom = np.uint64(qoff-hk) if qoff >= hk else np.uint64(0)
    qto = np.uint64(qoff + 1 + hk) if qoff + 1 + hk <= qsize else np.uint64(qsize)
    rj = np.uint64(ri + qto - qfrom)
    F = np.zeros([kmer, KMER_BASE_FEATURES], dtype = np.float32)

    bf, bt = np.uint64(base_offset+qfrom), np.uint64(base_offset+qto)

    a1f, a1t = np.uint64(fipd_offset+qfrom), np.uint64(fipd_offset+qto)
    a1 = [ codev1_to_frame_table[x] for x in features[a1f:a1t] ]

    a2f, a2t = np.uint64(fpw_offset+qfrom), np.uint64(fpw_offset+qto)
    a2 = [ codev1_to_frame_table[x] for x in features[a2f:a2t] ]

    a3f, a3t = np.uint64(ripd_offset+qfrom), np.uint64(ripd_offset+qto)
    a3 = [ codev1_to_frame_table[x] for x in features[a3f:a3t] ]

    a4f, a4t = np.uint64(rpw_offset+qfrom), np.uint64(rpw_offset+qto)
    a4 = [ codev1_to_frame_table[x] for x in features[a4f:a4t] ]

    if c1 == 1:
        F[ri:rj, 0:4] = [ dna_one_hot_encoding[b] for b in features[bf:bt] ]
        F[ri:rj, 4]   = a1 
        F[ri:rj, 5]   = a2 
        F[ri:rj, 6]   = a3
        F[ri:rj, 7]   = a4
    else:
        F[ri:rj, 0:4] = [ dna_one_hot_encoding[int(3-b)] for b in features[bf:bt] ]
        F[ri:rj, 4]   = a3
        F[ri:rj, 5]   = a4
        F[ri:rj, 6]   = a1 
        F[ri:rj, 7]   = a2       
        F = np.flip(F, axis=0).copy()

    assert F[hk][1] == 1.0

    label = 1 if label > 0 else 0
    return torch.from_numpy(F), label

class KmerFeaturesDataset(Dataset):
    def __init__(self, features, samples_path, offsets_path, kmer_size):
        super(KmerFeaturesDataset, self).__init__()

        self.features = features
        self.samples = load_samples(samples_path)
        self.offsets = load_offsets(offsets_path)
        self.kmer_size = kmer_size

    def __len__(self):
        return self.samples.shape[0]

    def __getitem__(self, idx):
        return assemble_one_sample_features(self.features, self.samples, self.offsets, self.kmer_size, idx)
