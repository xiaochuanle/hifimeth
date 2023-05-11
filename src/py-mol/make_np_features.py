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

def dump_usage(program):
    print('')
    print('USAGE:')
    print('%s num-positions has-label raw-features-path [saved-samples]' % (program))

def bin_to_np_array(kmer_size, has_label, raw_features_path, saved_samples, np_features_path):
    f = open(raw_features_path, "rb")
    strb = f.read()
    f.close()

    sample_elems = kmer_size * 6
    if has_label:
        sample_elems += 1
    print('sample elements: %d' % (sample_elems))

    file_size = len(strb)
    if file_size % sample_elems != 0:
        print('File size does not match sample elements (%d)' % (sample_elems))
        sys.exit(1)        

    print('Add elements to np array');
    array = np.frombuffer(strb, dtype=np.uint8)
    array = array.reshape([-1, sample_elems])
    if saved_samples > 0:
        array = array[0:saved_samples]
    print('Done')
    print('save %d samples' % (np.size(array, 0)))
    np.save(np_features_path, array);

def main(argv):
    n_param = len(argv)
    if n_param < 4:
        dump_usage(argv[0])
        sys.exit(1)

    kmer_size = (int)(argv[1])
    has_label = (int)(argv[2])
    raw_features_path = argv[3]
    np_features_path = raw_features_path + '.npy'
    saved_samples = 0
    if n_param == 5:
        saved_samples = (int)(argv[4])

    bin_to_np_array(kmer_size, has_label, raw_features_path, saved_samples, np_features_path)


if __name__ == "__main__":
    main(sys.argv)