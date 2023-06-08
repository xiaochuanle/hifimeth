# -*- coding: utf-8 -*-

import argparse
import os
import sys
import numpy as np
import random
import time
import torch
from tqdm import tqdm
from torch_geometric.data import Data, Dataset
from torch_geometric.loader import DataLoader

import kmer_features_dataset
from kmer_features_dataset import kmer_base_features

import read_graph
from read_graph import ReadGraphDataset

import kmer_model
from kmer_model import KmerModel

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
#device = torch.device("cpu")

torch.multiprocessing.set_sharing_strategy('file_system')

def get_time_string():
    localtime = time.asctime(time.localtime(time.time()))
    return localtime

def load_num_chunks(features_dir):
    path = features_dir + '/num_chunks'
    num_chunks = -1
    with open(path, "r") as f:
        line = f.readline();
        num_chunks = int(line)
    return num_chunks

def chunk_is_processed(features_dir, chunk_id):
    path = features_dir + '/chunk_' + str(chunk_id) + '.bin.done'
    return os.path.exists(path)

def chunk_make_processed(features_dir, chunk_id):
    path = features_dir + '/chunk_' + str(chunk_id) + '.bin.done'
    file = open(path, "w")
    file.close()

def call_5mc(args):
    model_path = args.model 
    kmer_size = args.kmer_size
    batch_size = args.batch_size
    dl_workers = args.dl_workers
    features_dir = args.features
    motif = 'CG'

    print('model: %s' % (model_path))
    print('kmer size: %d' % (kmer_size))
    print('batch size: %d' % (batch_size))
    print('dl_workers: %d' % (dl_workers))
    print('features directory: %s' % (features_dir))

    model = KmerModel(num_node_features = kmer_base_features)
    model.load_state_dict(torch.load(model_path, map_location=device))
    #print(model)
    model.to(device) 
    model.eval()

    num_chunks = load_num_chunks(features_dir)
    print('number of chunks: %d' % (num_chunks))

    for x in range(num_chunks):
        if chunk_is_processed(features_dir, x):
            print('chunk %d is proceseed, skip it.' % (x))
            continue

        features_path = features_dir + '/chunk_' + str(x) + '.bin'
        sample_offset_path = features_dir + '/chunk_' + str(x) + '.samples'
        print(features_path)
        print(sample_offset_path)
        now = get_time_string()
        print('[%s] process %s' % (now, features_path))

        data = ReadGraphDataset(features_path, sample_offset_path, kmer_size, motif, 0)
        data_loader = DataLoader(data, batch_size = batch_size, shuffle = False, num_workers = dl_workers)

        n_samples = len(data)
        prob_list = np.zeros([n_samples], dtype=np.float32)
        prob_idx = 0

        with torch.no_grad():
            for data in tqdm(data_loader):
                data.to(device)
                out = model(data.x, data.edge_index, data.batch)
                out = torch.softmax(out, 1)
                out.to('cpu')
                n = len(out)
                for i in range(n):
                    prob_list[prob_idx+i] = out[i][1]
                prob_idx += n 

        prob_path = features_dir + '/chunk_' + str(x) + '.prob';
        now = get_time_string()
        print('[%s] probability path: %s' % (now, prob_path))
        np.savetxt(prob_path, prob_list, fmt="%.4f", delimiter="\n")

        del data_loader
        del data 

        chunk_make_processed(features_dir, x)

def main():
    parser = argparse.ArgumentParser("call 5mc")

    parser.add_argument('--model', type=str, required=True,
                        help="Path to the model")
    parser.add_argument('--kmer_size', type=int, default=400, required=False,
                        help="len of kmer. Default: 400")
    parser.add_argument('--batch_size', type=int, default=512, required=False,
                        help="numbe of examples in one batch. Default: 512")
    parser.add_argument('--dl_workers', type=int, default=16, required=False,
                        help="number of CPU threads")
    parser.add_argument('--features', type=str, required=True,
                        help="directory containint example features")                    

    args = parser.parse_args()

    call_5mc(args)

if __name__ == "__main__":
    if device == torch.device("cpu"):
        print('CUDA is not available')
    else:
        print('CUDA is available')

    main()
