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

import kmer_model
from kmer_model import KmerModel

import kmer_util
from kmer_util import eval_model

import kmer_features_dataset
from kmer_features_dataset import KmerFeaturesDataset
from kmer_features_dataset import kmer_base_features

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

def dump_usage(program):
    print('')
    print('USAGE:')
    print('%s model-path kmer-size features' % (program))

def main(argv):
    n_param = len(argv)
    if n_param < 4:
        dump_usage(argv[0])
        sys.exit(1)

    model_path = argv[1]
    kmer_size = int(argv[2])
    features_path = argv[3]

    model = KmerModel(num_node_features = kmer_base_features)
    model.load_state_dict(torch.load(model_path, map_location=device))
    print(model)
    model.to(device) 
    model.eval()
    criterion = torch.nn.CrossEntropyLoss()

    data = KmerFeaturesDataset(features_path, kmer_size, 1)
    data_loader = DataLoader(data, batch_size = 512, shuffle = False, num_workers = 16)

    eval_model(model, criterion, data_loader, -1, device)

if __name__ == "__main__":
    main(sys.argv)
